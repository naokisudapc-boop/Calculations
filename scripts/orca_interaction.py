#!/usr/bin/env python
"""
orca_interaction.py
===================
結晶断片ダイマーの会合エネルギーを .out から自動計算する (人が桁を書き写さない)。

    dE_int = E(complex) - E(monoA) - E(monoB)

make_crystal_fragments.py が作った命名規則のファイルを探す:
    <tag>_complex_<level>_SP.out, <tag>_monoA_<level>_SP.out, <tag>_monoB_<level>_SP.out
    <tag>_pull_+X.XA_<level>_SP.out   (あれば引き離し曲線も計算)
<level> は XTB2, r2SCAN3c など。レベルごとに別の表を出す。

同じフォルダの orca_summary.py (文字コード判定つき) を使うので、
scripts フォルダに両方を置いておくこと。標準ライブラリのみ。

使い方 (CrystalFragment フォルダで):
  python ..\\..\\scripts\\orca_interaction.py .                 # 表を表示
  python ..\\..\\scripts\\orca_interaction.py . --level XTB2    # レベルを限定
  python ..\\..\\scripts\\orca_interaction.py . --json dEint.json --csv dEint.csv

チェック内容 (問題があれば dE は出さず、理由を表示):
  * 3 ファイルが揃っている / 正常終了 / SCF 収束 / エネルギーがある
  * 原子数・電荷が complex = monoA + monoB になっている
  注意 (dE は出すが表示):
  * monoA と monoB のエネルギー差が 1e-5 Eh を超える
    (結晶中の対称等価な分子どうしなら、ほぼ同じエネルギーになるはず)
  * |dE| が 60 kcal/mol を超える (中性ダイマーとしては異常。フラグメント取り違えの疑い)

結晶中の 1 分子あたりのペア和 (--summary を指定したときだけ):
  make_crystal_fragments.py が出力した summary.csv の multiplicity (同種ダイマーの隣接分子数) を読み、
      E_pair = 1/2 * Σ (multiplicity_i * dE_i)
  を計算する。次の場合は合計を出さず ERROR を表示して終了コード 1 を返す:
    * summary.csv に載っているダイマーの .out が無い / 異常終了 / SCF 未収束
    * multiplicity が欠けている・不正
    * summary.csv に無いダイマーの結果が存在する
    * 警告つきのダイマーがある (--allow-warnings で許可)
    * 単一点以外 (最適化・振動解析) の出力が混ざっている
  python orca_interaction.py . --level XTB2 --summary summary.csv [--expect -54.47 --tol 0.01]
  (--expect: 期待値 kcal/mol との照合。許容差 --tol を超えたら終了コード 1)
  注意: 合計は近接ダイマー (summary.csv の cutoff 内) だけの単純なペア和。多体効果・分子の変形・遠距離寄与は含まない。

換算: 1 Eh = 627.509474 kcal/mol = 2625.499638 kJ/mol
"""
import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_summary import summarise  # noqa: E402

EH_TO_KCAL = 627.509474
EH_TO_KJ = 2625.499638
AB_TOL_EH = 1e-5
BIG_KCAL = 60.0   # 中性ダイマーの会合エネルギーとしては異常に大きい値

NAME = re.compile(r"^(?P<tag>.+?)_(?P<role>complex|monoA|monoB)_(?P<level>.+?)_SP\.out$")
PULL = re.compile(r"^(?P<tag>.+?)_pull_(?P<d>[+-]?\d+(?:\.\d+)?)A_(?P<level>.+?)_SP\.out$")
ROLES = ("complex", "monoA", "monoB")


def collect(root):
    groups = {}
    for p in sorted(Path(root).rglob("*.out")):
        m = NAME.match(p.name)
        if m:
            groups.setdefault((m["level"], m["tag"]), {"pull": {}})[m["role"]] = p
            continue
        m = PULL.match(p.name)
        if m:
            groups.setdefault((m["level"], m["tag"]), {"pull": {}})["pull"][float(m["d"])] = p
    return groups


def check(name, s):
    """1 ファイルの致命的な問題を文字列のリストで返す。"""
    probs = []
    if "parse_error" in s:
        return [f"{name}: parse error"]
    if not s.get("terminated_normally"):
        probs.append(f"{name}: 正常終了していない")
    if s.get("scf_converged") is False:
        probs.append(f"{name}: SCF 未収束")
    if s.get("final_energy_eh") is None:
        probs.append(f"{name}: エネルギーが見つからない")
    return probs


def analyse(level, tag, g):
    res = {"dimer": tag, "level": level, "status": "ok", "problems": [], "notes": []}
    missing = [r for r in ROLES if r not in g]
    if missing:
        res["status"] = "incomplete"
        res["problems"].append("ファイル無し: " + ", ".join(missing))
        return res

    S = {}
    for r in ROLES:
        try:
            S[r] = summarise(g[r])
        except Exception as e:  # noqa: BLE001
            S[r] = {"parse_error": f"{type(e).__name__}: {e}"}
        res["problems"] += check(r, S[r])
    if res["problems"]:
        res["status"] = "failed"
        return res

    n = {r: S[r].get("n_atoms") for r in ROLES}
    if None not in n.values() and n["complex"] != n["monoA"] + n["monoB"]:
        res["problems"].append(f"原子数が合わない: complex {n['complex']} != {n['monoA']} + {n['monoB']}")
    q = {r: S[r].get("charge") for r in ROLES}
    if None not in q.values() and q["complex"] != q["monoA"] + q["monoB"]:
        res["problems"].append(f"電荷が合わない: complex {q['complex']} != {q['monoA']} + {q['monoB']}")
    if res["problems"]:
        res["status"] = "failed"
        return res

    Ec, Ea, Eb = (S[r]["final_energy_eh"] for r in ROLES)
    dE = Ec - Ea - Eb
    res["single_point"] = all(S[r].get("opt") is None and S[r].get("freq") is None for r in ROLES)
    res["keywords"] = S["complex"].get("keywords")
    res.update(
        E_complex_eh=Ec, E_monoA_eh=Ea, E_monoB_eh=Eb,
        dE_eh=dE, dE_kcal_mol=dE * EH_TO_KCAL, dE_kJ_mol=dE * EH_TO_KJ,
    )
    if abs(dE) * EH_TO_KCAL > BIG_KCAL:
        res["status"] = "warning"
        res["notes"].append(f"|dE| が {BIG_KCAL:g} kcal/mol を超えている (フラグメントの取り違えを疑う)")
    if abs(Ea - Eb) > AB_TOL_EH:
        res["status"] = "warning"
        res["notes"].append(f"monoA と monoB のエネルギー差 {abs(Ea - Eb):.2e} Eh (> {AB_TOL_EH:g})")

    # 引き離し曲線 (あれば)
    if g["pull"]:
        curve = [{"d_A": 0.0, "dE_kcal_mol": dE * EH_TO_KCAL, "vs_crystal_kcal_mol": 0.0}]
        for d, p in sorted(g["pull"].items()):
            try:
                s = summarise(p)
            except Exception:  # noqa: BLE001
                continue
            if check("pull", s):
                res["notes"].append(f"pull {d:+.1f} Å は未完了/失敗のため除外")
                continue
            e = s["final_energy_eh"]
            curve.append({
                "d_A": d,
                "dE_kcal_mol": (e - Ea - Eb) * EH_TO_KCAL,
                "vs_crystal_kcal_mol": (e - Ec) * EH_TO_KCAL,
            })
        res["curve"] = curve
    return res


def read_summary(path):
    """summary.csv -> (ダイマー名の順序リスト, {ダイマー名: multiplicity 文字列})"""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        if "dimer" not in cols:
            raise ValueError(f"{path}: 'dimer' 列がありません")
        if "multiplicity" not in cols:
            raise ValueError(
                f"{path}: 'multiplicity' 列がありません "
                "(新しい make_crystal_fragments.py で summary.csv を作り直してください)"
            )
        order, mult = [], {}
        for row in rd:
            d = (row["dimer"] or "").strip()
            if d:
                order.append(d)
                mult[d] = (row.get("multiplicity") or "").strip()
    return order, mult


def pair_sum(level, rows, order, mult_raw, allow_warnings=False):
    """1/2 * Σ(multiplicity * dE)。条件を満たさなければ computed=False とエラー一覧を返す。"""
    errs = []
    byd = {r["dimer"]: r for r in rows}
    mult = {}
    for d in order:
        try:
            m = int(mult_raw.get(d, ""))
            if m < 1:
                raise ValueError
            mult[d] = m
        except ValueError:
            errs.append(f"multiplicity is missing or invalid for {d}.")
    for d in order:
        r = byd.get(d)
        if r is None:
            errs.append(f"{d}: no {level} output files found.")
        elif "dE_eh" not in r:
            errs.append(f"{d}: output missing or not usable ({r['status']}: {'; '.join(r['problems'])}).")
        else:
            if r["notes"] and not allow_warnings:
                errs.append(f"{d}: has warnings ({'; '.join(r['notes'])}); use --allow-warnings to include.")
            if not r.get("single_point", True):
                errs.append(f"{d}: output is not a plain single point (opt/freq found).")
    for d in byd:
        if d not in order:
            errs.append(f"{d}: has {level} results but is not listed in the summary (multiplicity unknown).")
    if errs:
        return {"computed": False, "level": level, "errors": errs}
    total = sum(mult[d] * byd[d]["dE_eh"] for d in order) / 2
    kws = sorted({byd[d].get("keywords") for d in order if byd[d].get("keywords")})
    return {
        "computed": True, "level": level, "dimer_types": len(order),
        "neighbor_molecules": sum(mult.values()), "keywords": kws,
        "E_pair_eh": total, "kcal_mol": total * EH_TO_KCAL, "kJ_mol": total * EH_TO_KJ,
        "multiplicity": {d: mult[d] for d in order},
    }


def fmt_pair_sum(ps):
    if not ps["computed"]:
        out = [f"[{ps['level']}] Pairwise crystal interaction estimate"]
        out += [f"ERROR: {e}" for e in ps["errors"]]
        out.append("Pairwise sum was NOT calculated.")
        return "\n".join(out)
    kw = f" (output keywords: {'; '.join(ps['keywords'])})" if ps["keywords"] else ""
    return "\n".join([
        "Pairwise crystal interaction estimate",
        "--------------------------------------",
        f"Dimer types       : {ps['dimer_types']}",
        f"Neighbor molecules: {ps['neighbor_molecules']}",
        f"Method            : {ps['level']}{kw}",
        "Geometry          : crystal_fixed (assumed: inputs from make_crystal_fragments.py; not verifiable from .out)",
        "Calculation       : single_point (checked: no opt/freq in outputs)",
        "",
        "1/2 Σ(multiplicity × ΔEint)",
        f"= {ps['kcal_mol']:.2f} kcal/mol",
        f"= {ps['kJ_mol']:.1f} kJ/mol",
        "Scope: pair sum over listed near neighbours only (no many-body, deformation or long-range terms).",
    ])


def fmt_table(level, rows):
    lines = [f"[{level}]  dE_int = E(complex) - E(monoA) - E(monoB)"]
    hdr = f"{'dimer':<6}{'E_complex (Eh)':>20}{'E_monoA (Eh)':>20}{'E_monoB (Eh)':>20}{'dE (Eh)':>17}{'kcal/mol':>10}{'kJ/mol':>10}  note"
    lines.append(hdr)
    ok = sorted([r for r in rows if "dE_eh" in r], key=lambda r: r["dE_eh"])
    bad = [r for r in rows if "dE_eh" not in r]
    for r in ok:
        note = "; ".join(r["notes"]) if r["notes"] else ""
        lines.append(
            f"{r['dimer']:<6}{r['E_complex_eh']:>20.11f}{r['E_monoA_eh']:>20.11f}{r['E_monoB_eh']:>20.11f}"
            f"{r['dE_eh']:>17.11f}{r['dE_kcal_mol']:>10.2f}{r['dE_kJ_mol']:>10.2f}  {note}"
        )
    for r in bad:
        lines.append(f"{r['dimer']:<6}  -- {r['status']}: " + " / ".join(r["problems"]))
    for r in ok:
        if "curve" in r:
            pts = "  ".join(f"{c['d_A']:+.1f}Å:{c['vs_crystal_kcal_mol']:+.2f}" for c in r["curve"])
            lines.append(f"  {r['dimer']} 引き離し (結晶配置との差, kcal/mol): {pts}")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", nargs="?", default=".", help="検索するフォルダ (再帰的に探す)")
    ap.add_argument("--level", help="計算レベルを限定 (例: XTB2, r2SCAN3c)")
    ap.add_argument("--json", help="小さな JSON の出力先")
    ap.add_argument("--csv", help="CSV の出力先")
    ap.add_argument("--summary", help="make_crystal_fragments.py の summary.csv (multiplicity を読んでペア和を計算)")
    ap.add_argument("--allow-warnings", action="store_true", help="警告つきのダイマーもペア和に含める")
    ap.add_argument("--expect", type=float, help="ペア和の期待値 (kcal/mol)。--summary と併用")
    ap.add_argument("--tol", type=float, default=0.01, help="--expect の許容差 (kcal/mol, 既定 0.01)")
    a = ap.parse_args()

    groups = collect(a.root)
    if a.level:
        groups = {k: v for k, v in groups.items() if k[0] == a.level}
    if not groups:
        sys.exit("対象の .out が見つかりません (命名規則: <tag>_complex_<level>_SP.out など)")

    by_level = {}
    for (level, tag), g in sorted(groups.items()):
        by_level.setdefault(level, []).append(analyse(level, tag, g))

    for level, rows in by_level.items():
        print(fmt_table(level, rows))
        print()

    pair = {}
    exit_code = 0
    if a.summary:
        try:
            order, mult_raw = read_summary(a.summary)
        except (OSError, ValueError) as e:
            sys.exit(f"ERROR: {e}\nPairwise sum was NOT calculated.")
        for level, rows in by_level.items():
            ps = pair[level] = pair_sum(level, rows, order, mult_raw, a.allow_warnings)
            print(fmt_pair_sum(ps))
            if not ps["computed"]:
                exit_code = 1
            elif a.expect is not None:
                ok = abs(ps["kcal_mol"] - a.expect) <= a.tol
                ps["expect_kcal_mol"], ps["expect_tol"], ps["expect_pass"] = a.expect, a.tol, ok
                print(f"Check vs expected {a.expect:.2f} ± {a.tol:g} kcal/mol: " + ("PASS" if ok else "FAIL"))
                if not ok:
                    exit_code = 1
            print()

    if a.json:
        out = {}
        for level, rows in by_level.items():
            d = {}
            for r in rows:
                if "dE_eh" in r:
                    e = {
                        "dE_kcal_mol": round(r["dE_kcal_mol"], 2),
                        "dE_kJ_mol": round(r["dE_kJ_mol"], 2),
                        "dE_eh": round(r["dE_eh"], 11),
                    }
                    if r["notes"]:
                        e["notes"] = r["notes"]
                    if "curve" in r:
                        e["detach_kcal_mol"] = {f"{c['d_A']:+.1f}A": round(c["vs_crystal_kcal_mol"], 2) for c in r["curve"]}
                else:
                    e = {"status": r["status"], "problems": r["problems"]}
                d[r["dimer"]] = e
            out[level] = {"single_point": True, "dimers": d}
            if level in pair:
                ps = pair[level]
                out[level]["pair_sum"] = (
                    {"computed": False, "errors": ps["errors"]} if not ps["computed"] else {
                        "computed": True,
                        "formula": "1/2 * sum(multiplicity * dE_int)",
                        "dimer_types": ps["dimer_types"],
                        "neighbor_molecules": ps["neighbor_molecules"],
                        "kcal_mol": round(ps["kcal_mol"], 2),
                        "kJ_mol": round(ps["kJ_mol"], 1),
                        "geometry": "crystal_fixed (assumed)",
                        "multiplicity": ps["multiplicity"],
                        **({"expect_pass": ps["expect_pass"]} if "expect_pass" in ps else {}),
                    }
                )
        Path(a.json).write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
        print(f"wrote {a.json}")

    if a.csv:
        cols = ["level", "dimer", "status", "E_complex_eh", "E_monoA_eh", "E_monoB_eh",
                "dE_eh", "dE_kcal_mol", "dE_kJ_mol", "notes", "problems"]
        with open(a.csv, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            for rows in by_level.values():
                for r in rows:
                    row = {k: r.get(k) for k in cols}
                    row["notes"] = "; ".join(r["notes"])
                    row["problems"] = "; ".join(r["problems"])
                    w.writerow(row)
        print(f"wrote {a.csv}")

    if exit_code:
        sys.exit(exit_code)


if __name__ == "__main__":
    main()
