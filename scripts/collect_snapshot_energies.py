#!/usr/bin/env python3
"""extract_snapshots.py で作って実行したORCA出力(.out)から ΔE を集計する。

  ΔE_A-S   = E(A+nS) - E(A) - E(nS)
  ΔE_shell = E(nS) - n * E(S)         (--es と --nsolv を渡したとき)
  ΔE_total = ΔE_A-S + ΔE_shell

E(S) は溶媒1分子のエネルギー(Eh)。どのジオメトリのE(S)を使うか(緩和済みの気相最適化構造など)は
研究側で決めて --es に渡す。値はこのスクリプトでは仮定しない。

使い方:
  python collect_snapshot_energies.py snap_water --out water_dE.csv
  python collect_snapshot_energies.py snap_water --es -76.3123 --nsolv 10 --out water_dE.csv
"""
import argparse
import csv
import glob
import os
import re

HARTREE_TO_KCAL = 627.5094740631
PAT = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")


def read_text(path):
    """.out をバイト列で読み、エンコーディングを判別して文字列にする。

    Windows PowerShell 5.1 の `>` リダイレクトは UTF-16 LE (BOMあり) で保存するため、
    open() の既定エンコーディングでは行が読めない。BOM/NUL を見て UTF-16 を判別する。
    """
    with open(path, "rb") as f:
        data = f.read()
    if data[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in data[:400]:
        return data.decode("utf-16", errors="replace")
    return data.decode("utf-8-sig", errors="replace")


def read_energy(path):
    if not os.path.exists(path):
        return None
    last = None
    for ln in read_text(path).splitlines():
        m = PAT.search(ln)
        if m:
            last = float(m.group(1))
    return last


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--es", type=float, default=None, help="溶媒1分子のエネルギー (Eh)")
    ap.add_argument("--nsolv", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    rows, missing = [], 0
    for d in sorted(glob.glob(os.path.join(a.outdir, "f*"))):
        if not os.path.isdir(d):
            continue
        tag = os.path.basename(d)
        e = {t: read_energy(os.path.join(d, f"{tag}_{t}.out")) for t in ("AnS", "A", "nS")}
        if any(v is None for v in e.values()):
            missing += 1
            continue
        row = {"frame": int(tag[1:]), "E_AnS_Eh": e["AnS"], "E_A_Eh": e["A"], "E_nS_Eh": e["nS"],
               "dE_AS_kcal": round((e["AnS"] - e["A"] - e["nS"]) * HARTREE_TO_KCAL, 3)}
        if a.es is not None:
            row["dE_shell_kcal"] = round((e["nS"] - a.nsolv * a.es) * HARTREE_TO_KCAL, 3)
            row["dE_total_kcal"] = round(row["dE_AS_kcal"] + row["dE_shell_kcal"], 3)
        rows.append(row)

    print(f"集計 {len(rows)} フレーム, 未完了/欠損 {missing}")
    if rows:
        v = [r["dE_AS_kcal"] for r in rows]
        print(f"dE_AS (kcal/mol): 平均 {sum(v) / len(v):.2f}, 最小 {min(v):.2f}, 最大 {max(v):.2f}")
    if a.out and rows:
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"書き出し: {a.out}")


if __name__ == "__main__":
    main()
