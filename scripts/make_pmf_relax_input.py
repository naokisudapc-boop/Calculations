#!/usr/bin/env python3
"""[1/2] H+DMSO のみ緩和する ORCA 入力を新規作成する。ORCA は実行しない。既存ファイルは変更しない。

固定: 殻(99..629) と A の重原子。 active: A の H (元素記号から自動抽出) + DMSO 59..98。
緩和中は Colvar / Restraint を使わない。 Minimize は LBFGS, StepLimit 0.1 A, MaxGrad/RMSGrad 明示。

HalfShell フォルダで実行:  python ..\\..\\..\\..\\scripts\\make_pmf_relax_input.py
出力: PMF_DMSO_f014\\relax\\ (PMF_f0014_relax.inp, before.xyz)
"""
import os
import shutil
import sys

SRC = r"PMF_DMSO_f014\PMF_f0014_rotated_r0.xyz"
OUT = r"PMF_DMSO_f014\relax"
NA, NS, NTOT = 59, 40, 630
STEPS, MAXGRAD, RMSGRAD, STEPLIMIT = 300, 5.0, 1.0, 0.1   # MaxGrad/RMSGrad: kJ/mol/A (ORCA既定と同じ), StepLimit: A


def read_xyz(p):
    L = open(p, encoding="utf-8").read().splitlines()
    n = int(L[0])
    at = [l.split() for l in L[2:2 + n]]
    return [a[0] for a in at], [(float(a[1]), float(a[2]), float(a[3])) for a in at]


def main():
    if not os.path.exists(SRC):
        sys.exit(f"見つからない: {SRC} (HalfShell フォルダで実行)")
    if os.path.exists(OUT):
        sys.exit(f"{OUT} が既にある。上書きしないので中止")
    sym, xyz = read_xyz(SRC)
    assert len(sym) == NTOT, len(sym)
    hA = [i for i in range(NA) if sym[i] == "H"]
    assert len(hA) == 25, len(hA)
    assert all(sym[NA + 10 * m + 1] == "S" for m in range(4))
    active = ", ".join(str(i) for i in hA) + f", {NA}..{NA + NS - 1}"
    os.makedirs(OUT)
    shutil.copyfile(SRC, os.path.join(OUT, "before.xyz"))
    inp = [
        "# relax: A-H + 4 DMSO active; shell and A heavy atoms frozen; no Colvar/Restraint",
        "",
        "! XTB2 MD",
        "",
        "%pal",
        "  nprocs 1",
        "end",
        "",
        "%maxcore 2000",
        "",
        "%md",
        f"  Manage_Region active Define {active}",
        f"  Minimize LBFGS Steps {STEPS} MaxGrad {MAXGRAD} RMSGrad {RMSGRAD} StepLimit {STEPLIMIT}",
        '  Dump Position Stride 0 Filename "PMF_f0014_relaxed_final.xyz" Replace',
        "end",
        "",
        "* xyz 0 1",
    ]
    with open(os.path.join(OUT, "PMF_f0014_relax.inp"), "w", newline="\n") as f:
        f.write("\n".join(inp) + "\n")
        for s, (x, y, z) in zip(sym, xyz):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        f.write("*\n")
    print(f"active = A-H {len(hA)} 原子 + DMSO {NS} 原子 = {len(hA) + NS} 原子 (固定 {NTOT - len(hA) - NS})")
    print(f"A-H 番号: {hA}")
    print(f"出力: {OUT}\\PMF_f0014_relax.inp  (ORCA 未実行)")
    print("Minimize 後の最終構造は Dump Position Stride 0 で PMF_f0014_relaxed_final.xyz に出る想定 (未検証)")


if __name__ == "__main__":
    main()
