#!/usr/bin/env python3
"""PMF試験入力 (frame14, r=0, 50 steps) を新規作成する。ORCA は実行しない。既存ファイルは変更・上書きしない。

やること:
  1. SHS_DMSO_f014/SHS_f0014_r+0.0A_XTB2_SP.inp を読む (読み取り専用)
  2. 原子数 630 と原子順 A=0..58 / DMSO=59..98 / 殻=99..629 を検証
  3. v (selection.json) を +x へ回す最小回転 R を作り、殻の質量重心(COM)を軸に全原子を回転し、
     殻 COM を原点へ平行移動
  4. 回転後の CV0 = COM_x(A) - COM_x(殻) を質量加重で計算 (Distance X の初期値)
  5. 新規フォルダ PMF_DMSO_f014 に 試験 inp / 回転後 xyz / 変換情報 json / report.txt を書く

使い方 (HalfShell フォルダで):  python ..\\..\\..\\..\\scripts\\make_pmf_test_r0.py
Target(window) = CV0 + r  (r そのものを Target にしない)
"""
import json
import os
import sys

import numpy as np

SRC = r"SHS_DMSO_f014\SHS_f0014_r+0.0A_XTB2_SP.inp"
SEL = r"SHS_DMSO_f014\selection.json"
OUT = "PMF_DMSO_f014"
NA, NS, NTOT = 59, 40, 630
MASS = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "S": 32.065}  # 近似質量 (COMの重みのみ)
K_SPRING = 50.0   # kJ/mol/A^2 (ORCA既定値)
RUN_STEPS = 50
TEMP = "298.15_K"
TIMECON = "100.0_fs"


def read_inp(path):
    raw = open(path, "rb").read()
    txt = raw.decode("utf-16") if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else raw.decode("utf-8-sig")
    lines = txt.splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("* xyz"))
    at, j = [], i0 + 1
    while j < len(lines) and lines[j].strip() != "*":
        p = lines[j].split()
        if len(p) == 4:
            at.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        j += 1
    return at


def main():
    for p in (SRC, SEL):
        if not os.path.exists(p):
            sys.exit(f"見つからない: {p} (HalfShell フォルダで実行する)")
    if os.path.exists(OUT):
        sys.exit(f"{OUT} が既にある。上書きしないので中止")

    at = read_inp(SRC)
    sym = [a[0] for a in at]
    X = np.array([[a[1], a[2], a[3]] for a in at])
    chk = []
    chk.append(("原子数 == 630", len(at) == NTOT, len(at)))
    ok_d = all(sym[NA + 10 * m: NA + 10 * m + 10] == ["C", "S", "C", "O"] + ["H"] * 6 for m in range(4))
    chk.append(("DMSO x4 が 59..98 に C,S,C,O,H*6 で並ぶ", ok_d, ""))
    from collections import Counter
    cA = Counter(sym[:NA])
    chk.append(("A(0..58) が C25H25N5O4", dict(cA) == {"C": 25, "H": 25, "N": 5, "O": 4}, dict(cA)))
    chk.append(("殻(99..) に S が無い/原子数531", len(at) - NA - NS == 531 and "S" not in sym[NA + NS:], len(at) - NA - NS))

    v = np.array(json.load(open(SEL, encoding="utf-8"))["v"], float)
    chk.append(("|v| == 1", abs(np.linalg.norm(v) - 1) < 1e-6, np.linalg.norm(v)))
    v = v / np.linalg.norm(v)

    # v -> +x の最小回転 (Rodrigues)
    ex = np.array([1.0, 0.0, 0.0])
    k = np.cross(v, ex)
    c = float(v @ ex)
    s2 = float(k @ k)
    if s2 < 1e-12:
        sys.exit("v が x 軸と平行/反平行。別処理が必要")
    Kx = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    R = np.eye(3) + Kx + Kx @ Kx * ((1 - c) / s2)
    chk.append(("R 直交 & det=+1", np.allclose(R @ R.T, np.eye(3), atol=1e-12) and abs(np.linalg.det(R) - 1) < 1e-12,
                f"det={np.linalg.det(R):.12f}"))
    chk.append(("R v == +x", np.allclose(R @ v, ex, atol=1e-12), (R @ v).tolist()))

    m = np.array([MASS[s] for s in sym])
    shell = slice(NA + NS, NTOT)
    pivot = (m[shell, None] * X[shell]).sum(0) / m[shell].sum()
    Y = (X - pivot) @ R.T          # 殻COMが原点、v が +x
    # 内部距離の保存確認
    d0 = np.linalg.norm(X[:NA, None] - X[None, :NA], axis=-1)
    d1 = np.linalg.norm(Y[:NA, None] - Y[None, :NA], axis=-1)
    chk.append(("A内部距離が回転で不変 (max差<1e-9)", np.abs(d0 - d1).max() < 1e-9, np.abs(d0 - d1).max()))

    comA = (m[:NA, None] * Y[:NA]).sum(0) / m[:NA].sum()
    comS = (m[shell, None] * Y[shell]).sum(0) / m[shell].sum()
    CV0 = float(comA[0] - comS[0])
    chk.append(("CV0 > 0 (A が殻の +x 側)", CV0 > 0, CV0))

    os.makedirs(OUT)
    tag = "PMF_f0014_r000_test"
    with open(os.path.join(OUT, "PMF_f0014_rotated_r0.xyz"), "w", newline="\n") as f:
        f.write(f"{NTOT}\nframe14 r=0 rotated so v->+x; shell COM at origin\n")
        for s, (x, y, z) in zip(sym, Y):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

    inp = [
        "# PMF test: frame 14, r=0 window, rotated (v -> +x), 4 DMSO + A active, shell frozen, XTB2 MD",
        f"# CV0 = {CV0:.8f} A ; Target = CV0 + r ; Spring {K_SPRING} kJ/mol/A^2",
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
        "  Timestep 1.0_fs",
        f"  Initvel {TEMP}",
        f"  Thermostat NHC {TEMP} Timecon {TIMECON}",
        "  Manage_Region active Define 0..98",
        "  Manage_Colvar Define 1 Distance X Group 0..58 Group 99..629",
        f"  Restraint Add Colvar 1 Harmonic Spring {K_SPRING} Target {CV0:.8f}",
        f'  Dump Position Stride 1 Filename "{tag}_traj.xyz"',
        f"  Run {RUN_STEPS}",
        "end",
        "",
        "* xyz 0 1",
    ]
    with open(os.path.join(OUT, tag + ".inp"), "w", newline="\n") as f:
        f.write("\n".join(inp) + "\n")
        for s, (x, y, z) in zip(sym, Y):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        f.write("*\n")

    json.dump(dict(source=SRC, v=v.tolist(), R=R.tolist(), pivot_shell_COM_original=pivot.tolist(),
                   CV0=CV0, comA_rot=comA.tolist(), comShell_rot=comS.tolist(), spring=K_SPRING,
                   masses="approximate average masses, COM weights only (ORCA internal masses may differ ~1e-4)",
                   target_rule="Target = CV0 + r"),
              open(os.path.join(OUT, "pmf_transform.json"), "w", encoding="utf-8"), indent=1)

    lines = ["== チェック =="] + [f"[{'OK' if ok else 'NG'}] {n}  {val}" for n, ok, val in chk]
    lines += ["", f"v = {v.tolist()}", "R (v -> +x) =", *[str(r) for r in R.tolist()],
              f"元座標での殻COM(回転中心) = {pivot.tolist()}",
              f"回転後 A-COM = {comA.tolist()}", f"回転後 殻-COM = {comS.tolist()}",
              f"CV0 = {CV0:.8f} A", f"Target(r=0) = {CV0:.8f} A", f"Target(r) = {CV0:.8f} + r",
              "", f"出力: {OUT}\\{tag}.inp / PMF_f0014_rotated_r0.xyz / pmf_transform.json / report.txt",
              "ORCA は実行していない。"]
    open(os.path.join(OUT, "report.txt"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if not all(ok for _, ok, _ in chk):
        print("\n!! NG あり。試験入力は使わず確認すること")


if __name__ == "__main__":
    main()
