#!/usr/bin/env python3
"""H + DMSO のみの事前緩和 (r=0, frame14)。ORCA は実行しない。既存ファイルは変更・上書きしない。

  python ..\\..\\..\\..\\scripts\\pmf_relax_r0.py make      # 緩和入力を PMF_relax_f014 に生成
  (ORCA を手動で実行: orca PMF_relax_f014\\relax.inp > relax.out)
  python ..\\..\\..\\..\\scripts\\pmf_relax_r0.py analyze   # 緩和後の検証 + r=0 50step MD 試験入力を生成

HalfShell で実行する。入力元は PMF_DMSO_f014\\PMF_f0014_rotated_r0.xyz (v -> +x 回転済み, 殻COM=原点)。
固定: 殻(99..629) と A の重原子。活性: A の H + DMSO(59..98)。Colvar/Restraint なし。
"""
import os
import re
import sys

import numpy as np

SRC = r"PMF_DMSO_f014\PMF_f0014_rotated_r0.xyz"
DIR = "PMF_relax_f014"
NA, NS, NTOT = 59, 40, 630
SH0 = NA + NS
MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "S": 32.065}
MAXGRAD, RMSGRAD, STEPLIM, MAXSTEPS = 5.0, 1.0, 0.1, 500   # kJ/mol/A, kJ/mol/A, A (ORCA既定と同水準を明示)
K_SPRING, RUN_STEPS, TEMP, TIMECON = 50.0, 50, "298.15_K", "100.0_fs"
HEADER = ["! XTB2 MD", "", "%pal", "  nprocs 1", "end", "", "%maxcore 2000", ""]


def read_xyz_last(path):
    L = open(path, encoding="utf-8").read().splitlines()
    n = int(L[0].split()[0])
    nfr = len(L) // (n + 2)
    blk = L[(nfr - 1) * (n + 2) + 2: nfr * (n + 2)]
    sym, X = [], []
    for l in blk:
        p = l.split()
        sym.append(p[0]); X.append([float(v) for v in p[1:4]])
    return sym, np.array(X), nfr


def write_coords(f, sym, X):
    for s, (x, y, z) in zip(sym, X):
        f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def check_layout(sym):
    ok = len(sym) == NTOT
    ok &= all(sym[NA + 10 * m: NA + 10 * m + 10] == ["C", "S", "C", "O"] + ["H"] * 6 for m in range(4))
    return ok


def active_list(sym):
    return [i for i in range(NA) if sym[i] == "H"]


def mk():
    if not os.path.exists(SRC):
        sys.exit(f"見つからない: {SRC}")
    if os.path.exists(DIR):
        sys.exit(f"{DIR} が既にある。上書きしないので中止")
    sym, X, _ = read_xyz_last(SRC)
    if not check_layout(sym):
        sys.exit("原子数/原子順が想定と違う。中止")
    hA = active_list(sym)
    os.makedirs(DIR)
    md = [
        "%md",
        f"  Manage_Region active Define {', '.join(map(str, hA))}, {NA}..{SH0 - 1}",
        '  Dump Position Stride 0 Format XYZ Filename "relax_start.xyz" Replace',
        f"  Minimize LBFGS Steps {MAXSTEPS} MaxGrad {MAXGRAD} RMSGrad {RMSGRAD} StepLimit {STEPLIM}_A",
        '  Dump Position Stride 0 Format XYZ Filename "relax_end.xyz" Replace',
        "end", "", "* xyz 0 1"]
    with open(os.path.join(DIR, "relax.inp"), "w", newline="\n") as f:
        f.write("# pre-relaxation r=0 frame14: shell(99..629)+A heavy FIXED; A-H + DMSO(59..98) active; no restraint\n")
        f.write("\n".join(HEADER + md) + "\n")
        write_coords(f, sym, X)
        f.write("*\n")
    print(f"生成: {DIR}\\relax.inp  (active = A-H {len(hA)} 個 + DMSO 40 = {len(hA) + 40} 原子)")
    print(f"A-H 原子番号: {hA}")
    print("次: ORCA を手動実行 (例) orca .\\PMF_relax_f014\\relax.inp > .\\PMF_relax_f014\\relax.out")


def com(X, m, sl):
    return (m[sl, None] * X[sl]).sum(0) / m[sl].sum()


def hbonds(sym, X):
    out = {}
    heavy = [i for i in range(NA) if sym[i] != "H"]
    for h in active_list(sym):
        d = np.linalg.norm(X[heavy] - X[h], axis=1)
        j = heavy[int(d.argmin())]
        out[h] = (sym[j], j, float(d.min()))
    return out


def min_pair(X, act):
    """活性原子と全原子の最短距離 (結合込み, ORCAの表示と同定義)"""
    best = (9e9, -1, -1)
    for i in act:
        d = np.linalg.norm(X - X[i], axis=1)
        d[i] = 9e9
        j = int(d.argmin())
        if d[j] < best[0]:
            best = (float(d[j]), min(i, j), max(i, j))
    return best


def analyze():
    p0, p1 = os.path.join(DIR, "relax_start.xyz"), os.path.join(DIR, "relax_end.xyz")
    out = os.path.join(DIR, "relax.out")
    for p in (SRC, p0, p1):
        if not os.path.exists(p):
            sys.exit(f"見つからない: {p} (ORCA 実行後に実行する)")
    tag = "PMF_f0014_r000_relaxed_test"
    if os.path.exists(os.path.join(DIR, tag + ".inp")):
        sys.exit("試験入力が既にある。上書きしないので中止")
    sym, X0, _ = read_xyz_last(SRC)
    s1, X1, nfr = read_xyz_last(p1)
    if s1 != sym:
        sys.exit("緩和後 xyz の元素順が元と違う。中止")
    m = np.array([MASS[s] for s in sym])
    hA = set(active_list(sym))
    idx = {
        "A 重原子 (固定)": [i for i in range(NA) if i not in hA],
        "殻 99..629 (固定)": list(range(SH0, NTOT)),
        "A-H (活性)": sorted(hA),
        "DMSO 59..98 (活性)": list(range(NA, SH0)),
    }
    L = ["== 緩和の収束状況 (relax.out の該当行) =="]
    if os.path.exists(out):
        T = open(out, encoding="utf-8", errors="replace").read().splitlines()
        hit = [t for t in T if re.search(r"onverg|MaxGrad|RMSGrad|Minimi", t)]
        L += hit[-12:] if hit else ["(該当行なし。末尾15行)"] + T[-15:]
    else:
        L.append("relax.out なし")
    L += ["", f"relax_end.xyz フレーム数: {nfr} (最終フレームを使用)", "", "== 最大変位 [A] =="]
    for k, v in idx.items():
        d = np.linalg.norm(X1[v] - X0[v], axis=1)
        L.append(f"{k:22s} max {d.max():.6f}  mean {d.mean():.6f}  (n={len(v)})")
    fixed_ok = all(np.linalg.norm(X1[idx[k]] - X0[idx[k]], axis=1).max() < 1e-5
                   for k in ("A 重原子 (固定)", "殻 99..629 (固定)"))
    L.append(f"固定原子が動いていない (<1e-5 A): {'OK' if fixed_ok else 'NG'}")

    b0, b1 = hbonds(sym, X0), hbonds(sym, X1)
    L += ["", "== A の H-重原子距離 [A] (前 -> 後) =="]
    for el in ("N", "O", "C"):
        hs = [h for h in b0 if b0[h][0] == el]
        if not hs:
            continue
        if el in ("N", "O"):
            for h in hs:
                L.append(f"H{h:<2d}-{el}{b0[h][1]:<2d}  {b0[h][2]:.4f} -> {b1[h][2]:.4f}")
        else:
            a = np.array([b0[h][2] for h in hs]); b = np.array([b1[h][2] for h in hs])
            L.append(f"C-H (n={len(hs)})  前 min/mean/max {a.min():.4f}/{a.mean():.4f}/{a.max():.4f} "
                     f"-> 後 {b.min():.4f}/{b.mean():.4f}/{b.max():.4f}")
    act = sorted(hA) + list(range(NA, SH0))
    q0, q1 = min_pair(X0, act), min_pair(X1, act)
    L += ["", "== 活性原子の最短原子間距離 (結合を含む) ==",
          f"前 {q0[0]:.4f} A (atoms {q0[1]}-{q0[2]}, {sym[q0[1]]}-{sym[q0[2]]})",
          f"後 {q1[0]:.4f} A (atoms {q1[1]}-{q1[2]}, {sym[q1[1]]}-{sym[q1[2]]})"]

    sh = slice(SH0, NTOT)
    cv_old = float(com(X0, m, slice(0, NA))[0] - com(X0, m, sh)[0])
    cv_new = float(com(X1, m, slice(0, NA))[0] - com(X1, m, sh)[0])
    L += ["", f"CV0 前 {cv_old:.8f} A", f"CV0 後 (新) {cv_new:.8f} A", f"Target(r=0) = {cv_new:.8f} A  (window: CV0 + r)"]

    with open(os.path.join(DIR, "PMF_f0014_relaxed_r0.xyz"), "w", newline="\n") as f:
        f.write(f"{NTOT}\nrelaxed r=0 frame14 (A-H + DMSO relaxed)\n")
        write_coords(f, sym, X1)
    md = ["%md", "  Timestep 1.0_fs", f"  Initvel {TEMP}", f"  Thermostat NHC {TEMP} Timecon {TIMECON}",
          f"  Manage_Region active Define {', '.join(map(str, sorted(hA)))}, {NA}..{SH0 - 1}",
          "  Manage_Colvar Define 1 Distance X Group 0..58 Group 99..629",
          f"  Restraint Add Colvar 1 Harmonic Spring {K_SPRING} Target {cv_new:.8f}_A",
          f'  Dump Position Stride 1 Filename "{tag}_traj.xyz"', f"  Run {RUN_STEPS}", "end", "", "* xyz 0 1"]
    L += ["", "== r=0 試験入力 (MD 50 step) =="]
    with open(os.path.join(DIR, tag + ".inp"), "w", newline="\n") as f:
        f.write(f"# relaxed r=0 test: active = A-H + DMSO (shell, A heavy frozen); CV0={cv_new:.8f}; Target=CV0+r\n")
        f.write("\n".join(HEADER + md) + "\n")
        write_coords(f, sym, X1)
        f.write("*\n")
    L += md[:-3] + ["", f"出力: {DIR}\\{tag}.inp, PMF_f0014_relaxed_r0.xyz, report.txt", "ORCA(MD)は実行していない。"]
    open(os.path.join(DIR, "report.txt"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    {"make": mk, "analyze": analyze}.get(sys.argv[1] if len(sys.argv) > 1 else "", lambda: sys.exit(__doc__))()
