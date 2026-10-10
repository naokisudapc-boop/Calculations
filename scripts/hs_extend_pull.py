#!/usr/bin/env python3
"""HalfShell の引き離し曲線を延長する (XTB2, nprocs 1)。

既存の +0.0 Å / +5.0 Å の .inp から A (先頭 N_A 原子) の実際の座標変位を読み、
v を丸めずに +8 / +12 / +20 Å の入力を作る。計算後は E_bind(r) を表にする。

  python hs_extend_pull.py make    [--distances 8 12 20]
  python hs_extend_pull.py collect [--ea Eh] [--eshell Eh]

make が行う検査 (1つでも外れたら中止):
  - 2つの .inp の原子数・元素記号の並びが同じ
  - 先頭 N_A(=59) 原子の変位が全て同一ベクトル (許容 1e-6 Å)
  - 59原子目以降の変位が 0 (許容 1e-6 Å)
  - A が先頭59原子であること (元素組成が A_mono_XTB2_SP.inp と一致、あれば)
出力: HS_cluster_+8.0A_XTB2_SP.inp 等と run_extend.ps1 (すでに正常終了した .out はスキップ)
"""
import argparse
import glob
import os
import re
import sys

N_A = 59
HARTREE_TO_KCAL = 627.5094740631
PAT = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
TOL = 1e-6


def read_text(path):
    with open(path, "rb") as f:
        d = f.read()
    if d[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in d[:400]:
        return d.decode("utf-16", errors="replace")
    return d.decode("utf-8-sig", errors="replace")


def read_inp(path):
    """-> (header_lines, atoms[(sym,x,y,z)], footer_lines)。 '* xyz' から '*' までが座標。"""
    lines = read_text(path).splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("* xyz"))
    atoms, j = [], i0 + 1
    while j < len(lines) and lines[j].strip() != "*":
        p = lines[j].split()
        if len(p) == 4:
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        j += 1
    return lines[: i0 + 1], atoms, lines[j:]


def energy(path):
    if not os.path.exists(path):
        return None
    last = None
    for ln in read_text(path).splitlines():
        m = PAT.search(ln)
        if m:
            last = float(m.group(1))
    return last


def fmt(r):
    return f"{r:.1f}"


def make(a):
    p0, p5 = a.ref0, a.ref5
    h0, at0, f0 = read_inp(p0)
    h5, at5, f5 = read_inp(p5)
    if len(at0) != len(at5):
        sys.exit(f"原子数が違う: {len(at0)} vs {len(at5)}")
    if [x[0] for x in at0] != [x[0] for x in at5]:
        sys.exit("元素記号の並びが違う")
    n = len(at0)
    d = [(b[1] - c[1], b[2] - c[2], b[3] - c[3]) for c, b in zip(at0, at5)]
    ref = d[0]
    bad_a = [i for i in range(N_A) if max(abs(d[i][k] - ref[k]) for k in range(3)) > TOL]
    bad_s = [i for i in range(N_A, n) if max(abs(x) for x in d[i]) > TOL]
    print(f"原子数 {n} / 先頭{N_A}原子の変位 (+5.0 - +0.0) = ({ref[0]:.8f}, {ref[1]:.8f}, {ref[2]:.8f})")
    print(f"  |変位| = {sum(x * x for x in ref) ** 0.5:.8f} Å (期待 5.0)")
    if bad_a or bad_s:
        sys.exit(f"検査失敗: Aの変位が不揃い {len(bad_a)} 原子 (例 {bad_a[:5]}), "
                 f"殻が動いている {len(bad_s)} 原子 (例 {bad_s[:5]})。作成を中止。")
    print("  検査OK: 先頭59原子が同一ベクトルで平行移動、残り %d 原子は不動" % (n - N_A))
    if os.path.exists(a.mono):
        _, am, _ = read_inp(a.mono)
        comp = lambda L: sorted(x[0] for x in L)
        ok = comp(am) == comp(at0[:N_A])
        print(f"  A_mono の元素組成と先頭{N_A}原子: {'一致' if ok else '不一致!'}")
        if not ok:
            sys.exit("A_mono と先頭59原子の組成が違う。中止。")
    v = tuple(x / 5.0 for x in ref)  # 1 Å あたりの変位 (実測値 / 5)
    print(f"  v (実測/5) = ({v[0]:.10f}, {v[1]:.10f}, {v[2]:.10f})")

    # 入力のヘッダは +0.0 Å (np1) を流用し、1行目のコメントだけ差し替え
    outs = []
    for r in a.distances:
        name = f"HS_cluster_+{fmt(r)}A_XTB2_SP"
        head = list(h0)
        head[0] = f"# half shell: A displaced +{fmt(r)} A along v (measured from +0.0/+5.0 inps), 9 fixed neighbours, XTB2 SP, nprocs 1"
        lines = head[:]
        for i, (s, x, y, z) in enumerate(at0):
            if i < N_A:
                x, y, z = x + v[0] * r, y + v[1] * r, z + v[2] * r
            lines.append(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}")
        lines.append("*")
        with open(name + ".inp", "w", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        outs.append(name)
        print(f"  書き出し: {name}.inp")
    ps = ["# 延長した HalfShell 計算を順に実行 (正常終了済みはスキップ)",
          'param([string]$Orca = "C:\\ORCA_6.1.1\\orca.exe")',
          "Set-Location $PSScriptRoot"]
    for name in outs:
        ps.append(f'$o = "{name}.out"')
        ps.append('if (-not ((Test-Path $o) -and (Select-String -Path $o -Pattern "ORCA TERMINATED NORMALLY" -Quiet))) {')
        ps.append(f'  & $Orca "{name}.inp" 2>&1 | Out-File -FilePath $o -Encoding utf8')
        ps.append("}")
    with open("run_extend.ps1", "w", newline="\r\n") as f:
        f.write("\n".join(ps) + "\n")
    print("  書き出し: run_extend.ps1  (実行: .\\run_extend.ps1)")


def collect(a):
    ea = a.ea if a.ea is not None else energy(a.mono_out)
    es = a.eshell if a.eshell is not None else energy(a.shell_out)
    if ea is None or es is None:
        sys.exit("E_A / E_shell が取れない。--ea / --eshell で指定してください。")
    print(f"E_A = {ea:.12f} Eh ({'引数' if a.ea is not None else a.mono_out})")
    print(f"E_shell = {es:.12f} Eh ({'引数' if a.eshell is not None else a.shell_out})")
    rows = []
    for p in glob.glob("HS_cluster_+*A_*.out"):
        m = re.match(r"HS_cluster_\+(\d+\.\d)A_", os.path.basename(p))
        if not m or "1core" in p:
            continue
        r = float(m.group(1))
        txt = read_text(p)
        e = energy(p)
        if e is None or "ORCA TERMINATED NORMALLY" not in txt:
            print(f"  (スキップ: {p} は未完了)")
            continue
        rows.append((r, e, os.path.basename(p)))
    # 同じ r が複数あれば np1 を優先 (表示のみ。重複は両方出す)
    rows.sort()
    if not rows:
        sys.exit("対象の .out が見つからない")
    e0 = next((e for r, e, _ in rows if r == 0.0), None)
    print(f"\n{'r(Å)':>6} {'E_cluster (Eh)':>20} {'ΔE_pull':>9} {'E_bind(r)':>10}  file")
    for r, e, fn in rows:
        eb = (e - es - ea) * HARTREE_TO_KCAL
        dp = (e - e0) * HARTREE_TO_KCAL if e0 is not None else float("nan")
        print(f"{r:6.1f} {e:20.12f} {dp:9.2f} {eb:10.2f}  {fn}")
    print("\nE_bind(r) = E_cluster(r) - E_shell - E_A (kcal/mol)。剛体・凍結・XTB2・気相。")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("make")
    m.add_argument("--ref0", default="HS_cluster_+0.0A_SP_np1.inp")
    m.add_argument("--ref5", default="HS_cluster_+5.0A_XTB2_SP.inp")
    m.add_argument("--mono", default="A_mono_XTB2_SP.inp")
    m.add_argument("--distances", type=float, nargs="+", default=[8.0, 12.0, 20.0])
    c = sub.add_parser("collect")
    c.add_argument("--ea", type=float, default=None)
    c.add_argument("--eshell", type=float, default=None)
    c.add_argument("--mono-out", dest="mono_out", default="A_mono_XTB2_SP.out")
    c.add_argument("--shell-out", dest="shell_out", default="HS_shell_XTB2_SP.out")
    a = ap.parse_args()
    make(a) if a.cmd == "make" else collect(a)


if __name__ == "__main__":
    main()
