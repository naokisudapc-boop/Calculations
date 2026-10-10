#!/usr/bin/env python3
"""PMF window 入力の生成 (新規ファイルのみ。ORCA は実行しない。既存は上書きしない)。

使い方 (HalfShell で):
  python ..\\..\\..\\..\\scripts\\pmf_make_windows.py 5 10 20 --run 300        # 予備試験
  python ..\\..\\..\\..\\scripts\\pmf_make_windows.py 0 0.5 ... --run 1500 --out PMF_win_prod

基準構造: PMF_relax_f014\\PMF_f0014_relaxed_r0.xyz (緩和済み, v->+x 回転済み, 殻COM=原点)
各窓: A+DMSO (0..98) を +x へ r だけ剛体並進。殻(99..629)は不動。
Target(r) = CV0 + r  (CV0 は基準構造で再計算)
熱浴: CSVR 298.15 K Timecon 50 fs Region active, Timestep 1 fs, Spring k (kJ/mol/A^2)
窓ごとに別フォルダ・別名・別乱数シード。launch スクリプト(run_windows.ps1)も生成する。
"""
import argparse
import os
import sys

import numpy as np

BASE = r"PMF_relax_f014\PMF_f0014_relaxed_r0.xyz"
NA, NS, NTOT = 59, 40, 630
SH0 = NA + NS
MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "S": 32.065}


def read_xyz(path):
    L = open(path, encoding="utf-8").read().splitlines()
    n = int(L[0].split()[0])
    sym, X = [], []
    for l in L[2:2 + n]:
        p = l.split()
        sym.append(p[0]); X.append([float(v) for v in p[1:4]])
    return sym, np.array(X)


def tag(r):
    return f"{r:05.1f}".replace(".", "p")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("r", type=float, nargs="+")
    ap.add_argument("--run", type=int, default=300)
    ap.add_argument("--k", type=float, default=50.0)
    ap.add_argument("--out", default="PMF_win_pilot")
    a = ap.parse_args()
    if not os.path.exists(BASE):
        sys.exit(f"見つからない: {BASE}")
    sym, X0 = read_xyz(BASE)
    if len(sym) != NTOT:
        sys.exit("原子数が630でない")
    m = np.array([MASS[s] for s in sym])
    com = lambda X, sl: (m[sl, None] * X[sl]).sum(0) / m[sl].sum()
    CV0 = float(com(X0, slice(0, NA))[0] - com(X0, slice(SH0, NTOT))[0])
    if os.path.exists(a.out):
        sys.exit(f"{a.out} が既にある。上書きしないので中止")
    os.makedirs(a.out)
    base = os.path.abspath(a.out)
    ps = ["$ErrorActionPreference = 'Continue'", "$procs = @()"]
    print(f"CV0 = {CV0:.8f} A  (基準 {BASE})")
    print("r[A]   Target[A]    A/DMSO-殻 最短距離[A]   フォルダ")
    for i, r in enumerate(a.r):
        t, d = tag(r), os.path.join(base, f"r{tag(r)}")
        os.makedirs(d)
        X = X0.copy()
        X[:SH0, 0] += r
        tgt = CV0 + r
        dmin = float(np.linalg.norm(X[:SH0, None] - X[None, SH0:], axis=-1).min())
        name = f"win_r{t}"
        md = ["%md", "  Timestep 1.0_fs", "  Manage_Region active Define 0..98",
              f"  Randomize {101 + i}", "  Initvel 298.15_K Region active",
              "  Thermostat CSVR 298.15_K Timecon 50.0_fs Region active",
              "  Manage_Colvar Define 1 Distance X Group 0..58 Group 99..629",
              f"  Restraint Add Colvar 1 Harmonic Spring {a.k} Target {tgt:.8f}_A",
              f'  Dump Position Stride 10 Filename "{name}_traj.xyz"', f"  Run {a.run}", "end", "", "* xyz 0 1"]
        with open(os.path.join(d, name + ".inp"), "w", newline="\n") as f:
            f.write(f"# PMF window r={r} A: Target=CV0+r={tgt:.8f} (CV0={CV0:.8f}); A+DMSO translated +x by r; shell frozen\n")
            f.write("\n".join(["! XTB2 MD", "", "%pal", "  nprocs 1", "end", "", "%maxcore 2000", ""] + md) + "\n")
            for s, (x, y, z) in zip(sym, X):
                f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
            f.write("*\n")
        ps.append(f'$procs += Start-Process -FilePath "C:\\ORCA_6.1.1\\orca.exe" -ArgumentList "{name}.inp" '
                  f'-WorkingDirectory "{d}" -RedirectStandardOutput "{d}\\{name}.out" '
                  f'-RedirectStandardError "{d}\\{name}.err" -NoNewWindow -PassThru')
        print(f"{r:5.1f}  {tgt:10.6f}   {dmin:8.3f}               {d}")
    ps += ['$procs | Wait-Process', 'Get-ChildItem -Recurse "' + base + '" -Filter *.out | ForEach-Object { $_.FullName; Select-String -Path $_.FullName -Pattern "ORCA TERMINATED NORMALLY" | Select-Object -First 1 }']
    open(os.path.join(base, "run_windows.ps1"), "w", encoding="utf-8").write("\n".join(ps) + "\n")
    print(f"\n起動スクリプト: {a.out}\\run_windows.ps1 (まだ実行していない)")


if __name__ == "__main__":
    main()
