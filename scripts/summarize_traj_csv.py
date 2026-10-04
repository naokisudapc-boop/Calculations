#!/usr/bin/env python3
"""analyze_md_trajectory.py (v2) が出した *_traj.csv を溶媒ごとに並べて要約する。

使い方 (名前=パス を複数):
  python summarize_traj_csv.py Water=Water\\MD\\water_traj.csv MeOH=MeOH\\MD\\10\\meoh_traj.csv

表示 (有効フレームのみ):
  n_frames, n3.5 / n4 / n5 (平均±標準偏差), n_far(平均), max_mol_min_dist(平均 / 最大),
  n4 の前半・後半フレーム平均 (傾向の確認)
注意: フレームは自己相関があるため、標準偏差は独立標本の誤差ではない。
"""
import csv
import statistics as st
import sys


def col(rows, name):
    return [float(r[name]) for r in rows if r.get(name, "") != ""]


def ms(v):
    return f"{st.mean(v):.2f}±{st.pstdev(v):.2f}"


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    out = []
    for arg in sys.argv[1:]:
        name, _, path = arg.partition("=")
        with open(path, newline="") as f:
            rows = [r for r in csv.DictReader(f) if r.get("n_atoms_ok") == "1"]
        need = ["n_in_shell", "n_in_shell_3.5", "n_in_shell_5", "n_far", "max_mol_min_dist"]
        miss = [c for c in need if c not in rows[0]]
        if miss:
            sys.exit(f"{path}: 列が足りません {miss} (v2スクリプトで再実行してください)")
        n4 = col(rows, "n_in_shell")
        h = len(n4) // 2
        out.append((
            name, len(rows), ms(col(rows, "n_in_shell_3.5")), ms(n4), ms(col(rows, "n_in_shell_5")),
            f"{st.mean(col(rows, 'n_far')):.2f}",
            f"{st.mean(col(rows, 'max_mol_min_dist')):.1f}/{max(col(rows, 'max_mol_min_dist')):.1f}",
            f"{st.mean(n4[:h]):.2f}/{st.mean(n4[h:]):.2f}",
        ))
    head = ("solvent", "frames", "n3.5", "n4", "n5", "n_far", "maxd mean/max", "n4 前半/後半")
    widths = [max(len(str(r[i])) for r in out + [head]) for i in range(len(head))]
    for r in [head] + out:
        print("  ".join(str(c).ljust(w) for c, w in zip(r, widths)))


if __name__ == "__main__":
    main()
