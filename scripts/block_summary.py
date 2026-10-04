#!/usr/bin/env python3
"""ΔE_A-S (collect_snapshot_energies.py の出力) と traj CSV (analyze_md_trajectory.py の出力) を
フレーム番号で結合し、ブロックごとの平均・傾向・接触数との関係を表示する。

使い方:
  python block_summary.py etoh_dE.csv etoh_traj.csv --block 20
  python block_summary.py water_dE.csv water_traj.csv --block 20 --tail 50

出力:
  - 結合数 (ΔE行数 / traj行数 / 結合数)
  - ブロック表: frames, n, dE 平均±SD (kcal/mol), n4 平均, n_far 平均, dE/n4 (平均の比)
  - 全体と末尾 --tail フレームでの r(dE, n4), 傾き (kcal/mol per contact), dE の時間傾向 (kcal/mol per ps)

注意: 1フレーム = 0.1 ps を仮定 (--ps-per-frame で変更)。フレーム間は自己相関があるため、
SD は独立標本の誤差ではない。時間傾向がゼロに近いときだけ「定常」とみなせる。
"""
import argparse
import csv
import math
import statistics as st


def load(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def linfit(x, y):
    n = len(x)
    if n < 3:
        return float("nan"), float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx == 0 or syy == 0:
        return float("nan"), float("nan")
    return sxy / sxx, sxy / math.sqrt(sxx * syy)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dE_csv")
    ap.add_argument("traj_csv")
    ap.add_argument("--block", type=int, default=20, help="1ブロックのフレーム数")
    ap.add_argument("--tail", type=int, default=50, help="末尾の解析フレーム数")
    ap.add_argument("--ps-per-frame", type=float, default=0.1)
    a = ap.parse_args()

    e_rows, t_rows = load(a.dE_csv), load(a.traj_csv)
    tm = {int(r["frame"]): r for r in t_rows if r.get("n_atoms_ok", "1") == "1"}
    rows = []
    for r in e_rows:
        k = int(r["frame"])
        if k in tm:
            rows.append({"frame": k, "de": float(r["dE_AS_kcal"]),
                         "n4": float(tm[k]["n_in_shell"]), "far": float(tm[k]["n_far"])})
    rows.sort(key=lambda r: r["frame"])
    print(f"ΔE行数 {len(e_rows)} / traj行数 {len(t_rows)} / 結合 {len(rows)}")
    if not rows:
        return

    print(f"\n{'frames':>11} {'n':>3} {'dE mean±SD':>14} {'n4':>6} {'n_far':>6} {'dE/n4':>8}")
    f0 = rows[0]["frame"]
    blocks = {}
    for r in rows:
        blocks.setdefault((r["frame"] - f0) // a.block, []).append(r)
    for b in sorted(blocks):
        g = blocks[b]
        de = [r["de"] for r in g]
        n4 = st.mean(r["n4"] for r in g)
        sd = st.pstdev(de) if len(de) > 1 else 0.0
        ratio = f"{st.mean(de) / n4:8.2f}" if n4 > 0 else "     n/a"
        print(f"{g[0]['frame']:>5}-{g[-1]['frame']:<5} {len(g):>3} "
              f"{st.mean(de):7.2f}±{sd:5.2f} {n4:6.2f} {st.mean(r['far'] for r in g):6.2f} {ratio}")

    def report(label, sub):
        de = [r["de"] for r in sub]
        n4 = [r["n4"] for r in sub]
        fr = [r["frame"] * a.ps_per_frame for r in sub]
        s_n, r_n = linfit(n4, de)
        s_t, _ = linfit(fr, de)
        print(f"{label:>10}: n={len(sub):3d}  mean dE={st.mean(de):7.2f}  "
              f"r(dE,n4)={r_n:6.3f}  slope(dE vs n4)={s_n:6.2f} kcal/mol/contact  "
              f"dE時間傾向={s_t:6.2f} kcal/mol/ps")

    print()
    report("全体", rows)
    report(f"末尾{a.tail}", rows[-a.tail:])


if __name__ == "__main__":
    main()
