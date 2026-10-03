#!/usr/bin/env python
"""
orca_summary.py
===============
ORCA の .out から必要な項目だけを抜き出して、小さな JSON / CSV にまとめる。
(巨大な .out をチャットに貼らずに済ませ、トークンを節約するためのスクリプト)

依存パッケージなし (標準ライブラリのみ)。OPI や ORCA 本体も不要。
ORCA 6.1.1 の出力で動作確認。CRLF 改行にも対応。

使い方:
  python orca_summary.py calc.out                       # JSON を標準出力へ
  python orca_summary.py calc.out other.out -o s.json   # 複数ファイル -> JSON ファイル
  python orca_summary.py DIR -r --csv summary.csv       # DIR 以下の *.out を再帰的に -> CSV
  python orca_summary.py calc.out --brief               # 1 行の要約だけ表示

出力の単位:
  エネルギーは Eh、双極子は *_au (原子単位) と *_debye を別キーで保持、
  振動数は cm^-1、軌道エネルギーは eV、時間は秒。
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

FLOAT = r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?"


def _last(pattern, text, flags=0):
    """最後にマッチした最初のグループを返す (なければ None)。"""
    m = None
    for m in re.finditer(pattern, text, flags):
        pass
    return m.group(1) if m else None


def _f(x):
    return float(x) if x is not None else None


def parse_orbitals(text):
    """最後の ORBITAL ENERGIES ブロックから HOMO/LUMO (closed shell 想定) を返す。"""
    idx = text.rfind("ORBITAL ENERGIES")
    if idx < 0:
        return None
    rows = []
    for line in text[idx:].splitlines()[4:]:
        m = re.match(rf"\s*(\d+)\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$", line)
        if not m:
            if rows:
                break
            continue
        rows.append((float(m.group(2)), float(m.group(3)), float(m.group(4))))
    occ = [r for r in rows if r[0] > 0.0]
    virt = [r for r in rows if r[0] == 0.0]
    if not occ or not virt:
        return None
    homo, lumo = occ[-1], virt[0]
    return {
        "homo_ev": homo[2],
        "lumo_ev": lumo[2],
        "gap_ev": round(lumo[2] - homo[2], 4),
        "open_shell_or_fractional": any(0.0 < r[0] < 2.0 and r[0] != 1.0 for r in rows) or None,
    }


def parse_frequencies(text):
    idx = text.rfind("VIBRATIONAL FREQUENCIES")
    if idx < 0:
        return None
    freqs = []
    for line in text[idx:].splitlines():
        m = re.match(rf"\s*\d+:\s+({FLOAT})\s+cm\*\*-1", line)
        if m:
            freqs.append(float(m.group(1)))
        elif freqs and line.strip() == "":
            break
    nonzero = [x for x in freqs if abs(x) > 1e-3]
    imag = [x for x in nonzero if x < 0]
    out = {
        "n_modes_nonzero": len(nonzero),
        "n_imaginary": len(imag),
        "lowest_cm1": min(nonzero) if nonzero else None,
        "imaginary_cm1": imag[:5] or None,
    }
    for key, pat in (
        ("temperature_K", rf"Temperature\s+\.\.\.\s+({FLOAT}) K"),
        ("electronic_energy_eh", rf"Electronic energy\s+\.\.\.\s+({FLOAT}) Eh"),
        ("zpe_eh", rf"Zero point energy\s+\.\.\.\s+({FLOAT}) Eh"),
        ("enthalpy_eh", rf"Total Enthalpy\s+\.\.\.\s+({FLOAT}) Eh"),
        ("entropy_term_eh", rf"Final entropy term\s+\.\.\.\s+({FLOAT}) Eh"),
        ("gibbs_eh", rf"Final Gibbs free energy\s+\.\.\.\s+({FLOAT}) Eh"),
        ("g_minus_eel_eh", rf"G-E\(el\)\s+\.\.\.\s+({FLOAT}) Eh"),
    ):
        out[key] = _f(_last(pat, text))
    return out


def summarise(path):
    data = Path(path).read_bytes()
    if data.startswith(b"\xff\xfe"):
        raw = data.decode("utf-16-le", errors="replace")
    elif data.startswith(b"\xfe\xff"):
        raw = data.decode("utf-16-be", errors="replace")
    elif data.startswith(b"\xef\xbb\xbf"):
        raw = data.decode("utf-8-sig", errors="replace")
    else:
        raw = data.decode("utf-8", errors="replace")
    text = raw.replace("\r\r\n", "\n").replace("\r\n", "\n").replace("\r", "\n")

    s = {"file": str(path)}
    s["orca_version"] = _last(r"Program Version\s+(\S+)", text)

    # 入力の "! ..." 行 (ORCA が出力冒頭にエコーしたもの)
    kw = re.findall(r"^\|\s*\d+>\s*!\s*(.+?)\s*$", text, re.M)
    s["keywords"] = " ".join(kw) if kw else None

    s["n_atoms"] = int(_last(r"Number of atoms\s+\.+\s+(\d+)", text) or 0) or None
    if s["n_atoms"] is None:  # XTB2 など: プロパティ出力 "# Number of atoms / # / 118" の形式
        na = _last(r"#\s*Number of atoms\s*\n#\s*\n\s*(\d+)", text)
        s["n_atoms"] = int(na) if na else None
    ch = _last(r"Total Charge\s+Charge\s+\.+\s+(-?\d+)", text)
    mu = _last(r"Multiplicity\s+Mult\s+\.+\s+(\d+)", text)
    s["charge"] = int(ch) if ch is not None else None
    s["multiplicity"] = int(mu) if mu is not None else None

    s["terminated_normally"] = "ORCA TERMINATED NORMALLY" in text
    m = re.findall(
        r"TOTAL RUN TIME:\s+(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", text
    )
    if m:
        d, h, mi, se, ms = map(int, m[-1])
        s["run_time_s"] = round(d * 86400 + h * 3600 + mi * 60 + se + ms / 1000, 1)
    else:
        s["run_time_s"] = None

    n_conv = len(re.findall(r"SCF CONVERGED AFTER", text))
    n_fail = len(re.findall(r"SCF NOT CONVERGED", text))
    s["scf_converged"] = (n_conv > 0 and n_fail == 0) if (n_conv or n_fail) else None
    s["scf_cycles_last"] = _last(r"SCF CONVERGED AFTER\s+(\d+)\s+CYCLES", text)
    if s["scf_cycles_last"] is not None:
        s["scf_cycles_last"] = int(s["scf_cycles_last"])

    s["final_energy_eh"] = _f(_last(rf"FINAL SINGLE POINT ENERGY\s+({FLOAT})", text))

    # 構造最適化
    cyc = _last(r"GEOMETRY OPTIMIZATION CYCLE\s+(\d+)", text)
    if cyc is not None:
        s["opt"] = {
            "cycles": int(cyc),
            "converged": "THE OPTIMIZATION HAS CONVERGED" in text,
        }
    else:
        s["opt"] = None

    s["freq"] = parse_frequencies(text)

    # 双極子 (a.u. と Debye を別キーで保持)
    dip = _last(rf"Total Dipole Moment\s+:\s+({FLOAT}\s+{FLOAT}\s+{FLOAT})", text)
    dau = _last(rf"Magnitude \(a\.u\.\)\s+:\s+({FLOAT})", text)
    dde = _last(rf"Magnitude \(Debye\)\s+:\s+({FLOAT})", text)
    if dau or dde:
        s["dipole"] = {
            "total_au": [float(x) for x in dip.split()] if dip else None,
            "magnitude_au": _f(dau),
            "magnitude_debye": _f(dde),
        }
    else:
        s["dipole"] = None

    s["orbitals"] = parse_orbitals(text)

    # 警告・エラー
    warns = re.findall(r"^\s*WARNING:\s*(.+?)\s*$", text, re.M)
    uniq = list(dict.fromkeys(w[:110] for w in warns))
    s["warnings"] = {"count": len(warns), "first": uniq[:3]} if warns else None
    err = re.findall(r"^.*(?:ORCA finished by error termination|aborting the run).*$", text, re.M)
    errmsg = re.findall(r"^\s*ERROR\s*!*\s*$\n\s*(.+)$", text, re.M)
    s["error"] = ([e.strip()[:150] for e in err[:1]] + [e.strip()[:150] for e in errmsg[:2]]) or None
    return s


def brief(s):
    parts = [Path(s["file"]).name]
    parts.append("OK" if s["terminated_normally"] else "NOT-TERMINATED")
    if s["final_energy_eh"] is not None:
        parts.append(f"E={s['final_energy_eh']:.8f}")
    if s["scf_converged"] is False:
        parts.append("SCF-FAIL")
    if s["opt"]:
        parts.append("opt=" + ("conv" if s["opt"]["converged"] else "NOT-conv") + f"({s['opt']['cycles']})")
    if s["freq"]:
        parts.append(f"imag={s['freq']['n_imaginary']}")
        if s["freq"]["gibbs_eh"] is not None:
            parts.append(f"G={s['freq']['gibbs_eh']:.8f}")
    if s["run_time_s"] is not None:
        parts.append(f"{s['run_time_s'] / 60:.1f}min")
    return "  ".join(parts)


def flatten(s):
    row = {}
    for k, v in s.items():
        if v is None:
            continue
        if isinstance(v, dict):
            for k2, v2 in v.items():
                row[f"{k}.{k2}"] = json.dumps(v2, ensure_ascii=False) if isinstance(v2, (list, dict)) else v2
        elif isinstance(v, list):
            row[k] = json.dumps(v, ensure_ascii=False)
        else:
            row[k] = v
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help=".out ファイル、またはフォルダ")
    ap.add_argument("-r", "--recursive", action="store_true", help="フォルダ内の *.out を再帰的に探す")
    ap.add_argument("-o", "--output", help="JSON の出力先 (省略時は標準出力)")
    ap.add_argument("--csv", help="CSV の出力先 (平坦化した表)")
    ap.add_argument("--brief", action="store_true", help="1 ファイル 1 行の短い要約だけ表示")
    a = ap.parse_args()

    files = []
    for p in map(Path, a.paths):
        if p.is_dir():
            files += sorted(p.rglob("*.out") if a.recursive else p.glob("*.out"))
        else:
            files.append(p)
    if not files:
        sys.exit("no .out files found")

    results = []
    for f in files:
        try:
            results.append(summarise(f))
        except Exception as e:  # 1 ファイルの失敗で全体を止めない
            results.append({"file": str(f), "parse_error": f"{type(e).__name__}: {e}"})

    if a.brief:
        for s in results:
            print(brief(s) if "parse_error" not in s else f"{s['file']}  PARSE-ERROR {s['parse_error']}")
    else:
        # None のキーを落として、さらに小さくする
        def prune(o):
            if isinstance(o, dict):
                return {k: prune(v) for k, v in o.items() if v is not None}
            return o

        payload = [prune(s) for s in results]
        payload = payload[0] if len(payload) == 1 else payload
        txt = json.dumps(payload, ensure_ascii=False, indent=1)
        if a.output:
            Path(a.output).write_text(txt, encoding="utf-8")
            print(f"wrote {a.output} ({len(results)} files)")
        else:
            print(txt)

    if a.csv:
        rows = [flatten(s) for s in results]
        cols = list(dict.fromkeys(k for r in rows for k in r))
        with open(a.csv, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {a.csv}")


if __name__ == "__main__":
    main()
