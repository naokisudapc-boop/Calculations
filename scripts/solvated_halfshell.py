#!/usr/bin/env python3
"""溶媒付きHalfShell (Model A: Aの周辺だけ溶媒) の構造を、既存MDスナップショットから作る。

設計 (ユーザー承認 2026-10-06):
  MDフレーム(Apixaban 59原子 + 溶媒 n 分子)の A を、HalfShell(+0.0 Å)の A (先頭59原子) へ
  Kabsch で重ね合わせ、その回転・並進を溶媒にも適用する。溶媒を人工配置しない。
  各溶媒分子は次の理由で除外し、理由を必ず記録する:
    far   : A との最小原子間距離 > --max-as
    behind: A重心からの COM を v へ射影した値 < --min-proj (結晶側/露出側でない)
    clash : pull 経路上 (r = distances の全点) で殻9分子との最小原子間距離 < --clash
  残った溶媒は A と一緒に剛体移動 (r Å × v)。殻は動かさない。
  v は HS_cluster_+0.0A_SP_np1.inp と HS_cluster_+5.0A_XTB2_SP.inp の A 変位 /5 (hs_extend_pull.py と同じ)。

使い方 (HalfShell フォルダで実行。numpy 必要):
  python ..\\..\\..\\..\\..\\scripts\\solvated_halfshell.py scan  TRAJ.xyz
  python ...\\solvated_halfshell.py make  TRAJ.xyz --frame 150 --outdir SHS_DMSO_f0150
  python ...\\solvated_halfshell.py collect --outdir SHS_DMSO_f0150
  python ...\\solvated_halfshell.py proxy TRAJ.xyz --frame 10 --keep 8 9 10 --outdir SHS_DMSO_f10   # S_all / A+S_keep / S_keep (既存ファイルは上書きしない)

原子数が 700 を超えたら XTB の stack overflow 対策で XTBEXE=C:\\ORCA_6.1.1\\xTB_bigstack\\xtb.exe を設定する
(生成する run.ps1 に含めてある)。

注意:
  * これは rigid / 凍結 / 気相 / GFN2-xTB の単一点。自由エネルギーではない。
  * E_bind_solv(r) = E[A(r)+S+Shell] - E[Shell] - E[A+S]。溶媒和エネルギー自体は含まず、
    「溶媒を伴う A と結晶殻の相互作用」。真空 E_bind(r) との差 ΔΔE(r) も同じ意味で読むこと。
  * MD の A と結晶の A は立体配座が異なりうる。重ね合わせ RMSD を必ず確認 (既定 >1.0 Å で中止)。
"""
import argparse
import csv
import json
import os
import re
import sys

import numpy as np

N_A = 59
N_SOLV_ATOMS = 10  # DMSO
HARTREE_TO_KCAL = 627.5094740631
PAT = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")


def read_text(path):
    with open(path, "rb") as f:
        d = f.read()
    if d[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in d[:400]:
        return d.decode("utf-16", errors="replace")
    return d.decode("utf-8-sig", errors="replace")


def read_inp(path):
    lines = read_text(path).splitlines()
    i0 = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("* xyz"))
    atoms, j = [], i0 + 1
    while j < len(lines) and lines[j].strip() != "*":
        p = lines[j].split()
        if len(p) == 4:
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        j += 1
    return lines[: i0 + 1], atoms


def read_xyz_frames(path):
    lines = read_text(path).splitlines()
    i, frames = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        atoms = []
        for ln in lines[i + 2: i + 2 + n]:
            p = ln.split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        if len(atoms) != n:
            break
        frames.append(atoms)
        i += 2 + n
    return frames


def arr(atoms):
    return np.array([[a[1], a[2], a[3]] for a in atoms], float)


def kabsch(P, Q):
    """P を Q へ重ねる回転 R と並進 (x' = (x - Pc) @ R + Qc) と RMSD。"""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = (Vt.T @ D @ U.T).T  # row-vector 形式
    rmsd_all = np.sqrt((((P - Pc) @ R + Qc - Q) ** 2).sum(1).mean())
    return R, Pc, Qc, rmsd_all


def mindist(X, Y):
    return float(np.sqrt(((X[:, None, :] - Y[None, :, :]) ** 2).sum(-1)).min())


COV = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}


def heavy_graph(X, syms):
    idx = [i for i, s in enumerate(syms) if s != "H"]
    g = {i: set() for i in idx}
    for p in range(len(idx)):
        for q in range(p + 1, len(idx)):
            i, j = idx[p], idx[q]
            if np.linalg.norm(X[i] - X[j]) < 1.25 * (COV[syms[i]] + COV[syms[j]]):
                g[i].add(j)
                g[j].add(i)
    return idx, g


def wl_colors(idx, g, syms, rounds=5):
    import hashlib
    col = {i: syms[i] for i in idx}
    for _ in range(rounds):
        col = {i: hashlib.md5((col[i] + "|" + ",".join(sorted(col[j] for j in g[i]))).encode()).hexdigest()
               for i in idx}
    return col


def heavy_isomorphisms(symA, XA, symB, XB, limit=64):
    """結晶A(重原子) -> MDのA(重原子) の、元素と結合を保つ全対応 (対称性で複数) を最大 limit 個返す。"""
    ia, ga = heavy_graph(XA, symA)
    ib, gb = heavy_graph(XB, symB)
    if len(ia) != len(ib):
        return []
    ca, cb = wl_colors(ia, ga, symA), wl_colors(ib, gb, symB)
    if sorted(ca.values()) != sorted(cb.values()):
        return []
    order, seen = [], set()
    start = min(ia, key=lambda i: sum(1 for j in ia if ca[j] == ca[i]))
    queue = [start]
    seen.add(start)
    while queue:
        u = queue.pop(0)
        order.append(u)
        for w in sorted(ga[u]):
            if w not in seen:
                seen.add(w)
                queue.append(w)
    sols = []

    def bt(k, m, used):
        if len(sols) >= limit:
            return
        if k == len(order):
            sols.append(dict(m))
            return
        a = order[k]
        for b in ib:
            if b in used or cb[b] != ca[a]:
                continue
            if all((n in ga[a]) == (m[n] in gb[b]) for n in m):
                m[a] = b
                used.add(b)
                bt(k + 1, m, used)
                del m[a]
                used.discard(b)

    bt(0, {}, set())
    return sols


class Ctx:
    pass


def load_context(a):
    c = Ctx()
    c.head, at0 = read_inp(a.ref0)
    _, at5 = read_inp(a.ref5)
    if len(at0) != len(at5):
        sys.exit("ref0/ref5 の原子数が違う")
    X0, X5 = arr(at0), arr(at5)
    d = X5[:N_A] - X0[:N_A]
    if np.abs(d - d[0]).max() > 1e-6 or np.abs(X5[N_A:] - X0[N_A:]).max() > 1e-6:
        sys.exit("ref0→ref5 で A が一様に動いていない/殻が動いている。中止")
    c.v = d[0] / 5.0
    c.sym0 = [x[0] for x in at0]
    c.XA = X0[:N_A]
    c.XS = X0[N_A:]
    c.n_shell = len(c.XS)
    c.isos = None
    c.heavy_idx = [i for i in range(N_A) if c.sym0[i] != "H"]
    return c


def analyse_frame(c, atoms, a, nsolv):
    exp = N_A + N_SOLV_ATOMS * nsolv
    if len(atoms) != exp:
        return None, f"原子数 {len(atoms)} != {exp}"
    sym = [x[0] for x in atoms]
    if sorted(sym[:N_A]) != sorted(c.sym0[:N_A]):
        return None, "MDのAと結晶Aで元素組成が違う"
    X = arr(atoms)
    if c.isos is None:
        c.isos = heavy_isomorphisms(c.sym0[:N_A], c.XA, sym[:N_A], X[:N_A])
        print(f"重原子の同型写像: {len(c.isos)} 個 (対称性による。各フレームで RMSD 最小のものを採用)")
        if not c.isos:
            return None, "重原子の結合グラフが結晶AとMDのAで一致しない (結合判定または分子の違い)"
    best = None
    for iso in c.isos:
        P = X[[iso[i] for i in c.heavy_idx]]
        Q = c.XA[c.heavy_idx]
        R_, Pc_, Qc_, rm = kabsch(P, Q)
        if best is None or rm < best[0]:
            best = (rm, R_, Pc_, Qc_)
    rmsd_h, R, Pc, Qc = best
    rmsd = rmsd_h
    Y = (X - Pc) @ R + Qc
    A = Y[:N_A]
    Acen = c.XA.mean(0)
    rows = []
    for m in range(nsolv):
        idx = slice(N_A + m * N_SOLV_ATOMS, N_A + (m + 1) * N_SOLV_ATOMS)
        S = Y[idx]
        com = S.mean(0)
        proj = float((com - Acen) @ c.v / np.linalg.norm(c.v))
        mAS = mindist(S, A)
        path_min = min(mindist(S + c.v * r, c.XS) for r in a.distances)
        r0 = mindist(S, c.XS)
        if mAS > a.max_as:
            st = "far"
        elif proj < a.min_proj:
            st = "behind"
        elif path_min < a.clash:
            st = "clash"
        else:
            st = "keep"
        rows.append(dict(mol=m + 1, min_A_S=round(mAS, 3), min_shell_S_r0=round(r0, 3),
                         min_shell_S_path=round(path_min, 3), proj_v=round(proj, 3), status=st))
    return dict(rmsd_all=float(rmsd), rmsd_heavy=rmsd_h, rows=rows, Y=Y), None


def cmd_scan(a, c):
    frames = read_xyz_frames(a.xyz)
    print(f"frames={len(frames)}  v=({c.v[0]:.6f},{c.v[1]:.6f},{c.v[2]:.6f})")
    print("frame  rmsd_heavy  keep  far  behind  clash  min_A_S(keep)")
    for k in range(a.start, len(frames), a.step):
        res, err = analyse_frame(c, frames[k], a, a.nsolv)
        if res is None:
            print(f"{k:5d}  skip: {err}")
            break
        st = [r["status"] for r in res["rows"]]
        kp = [r["min_A_S"] for r in res["rows"] if r["status"] == "keep"]
        print(f"{k:5d}  {res['rmsd_heavy']:9.3f}  {st.count('keep'):4d} {st.count('far'):4d} {st.count('behind'):6d} {st.count('clash'):6d}  "
              f"{(min(kp) if kp else float('nan')):.2f}")


def write_inp(path, head, comment, syms, X):
    h = list(head)
    h[0] = comment
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(h) + "\n")
        for s, (x, y, z) in zip(syms, X):
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")
        f.write("*\n")


def cmd_make(a, c):
    frames = read_xyz_frames(a.xyz)
    k = len(frames) - 1 if a.frame is None else a.frame
    atoms = frames[k]
    res, err = analyse_frame(c, atoms, a, a.nsolv)
    if res is None:
        sys.exit(err)
    print(f"frame {k}: RMSD(all)={res['rmsd_all']:.3f} Å, RMSD(heavy)={res['rmsd_heavy']:.3f} Å")
    if res["rmsd_heavy"] > a.max_rmsd and not a.force:
        sys.exit(f"RMSD(heavy) が {a.max_rmsd} Å を超える。A の立体配座が結晶と違いすぎる。--force で続行可 (非推奨)")
    keep = [r["mol"] for r in res["rows"] if r["status"] == "keep"]
    print("mol  min_A_S  min_shell_S(r=0)  min_shell_S(path)  proj_v   status")
    for r in res["rows"]:
        print(f"{r['mol']:3d}  {r['min_A_S']:7.3f}  {r['min_shell_S_r0']:15.3f}  {r['min_shell_S_path']:16.3f}  {r['proj_v']:7.3f}  {r['status']}")
    if not keep:
        sys.exit("残る溶媒分子が 0。--min-proj / --clash / --max-as を見直すか別フレームを使う")
    os.makedirs(a.outdir, exist_ok=True)
    Y = res["Y"]
    sym = [x[0] for x in atoms]
    solv_idx = []
    for m in keep:
        solv_idx += list(range(N_A + (m - 1) * N_SOLV_ATOMS, N_A + m * N_SOLV_ATOMS))
    sym_S = [sym[i] for i in solv_idx]
    S0 = Y[solv_idx]
    # 構成: [A(59, 結晶座標) | 溶媒(MD→結晶座標) | 殻(531)]。A は結晶の A を使う (RMSD 上記)
    symA = c.sym0[:N_A]
    symShell = c.sym0[N_A:]
    names = []
    for r in a.distances:
        nm = f"SHS_f{k:04d}_r+{r:.1f}A_XTB2_SP"
        X = np.vstack([c.XA + c.v * r, S0 + c.v * r, c.XS])
        write_inp(os.path.join(a.outdir, nm + ".inp"), c.head,
                  f"# solvated half shell (Model A): frame {k}, {len(keep)} DMSO rigid with A, A displaced +{r:.1f} A along v",
                  symA + sym_S + symShell, X)
        names.append(nm)
    nm_as = f"SHS_f{k:04d}_AS_XTB2_SP"
    write_inp(os.path.join(a.outdir, nm_as + ".inp"), c.head, f"# reference E[A+S]: frame {k}, same geometry as r=0",
              symA + sym_S, np.vstack([c.XA, S0]))
    names.append(nm_as)
    n_tot = N_A + len(sym_S) + c.n_shell
    with open(os.path.join(a.outdir, "run.ps1"), "w", newline="\r\n") as f:
        f.write('param([string]$Orca = "C:\\ORCA_6.1.1\\orca.exe")\nSet-Location $PSScriptRoot\n')
        if n_tot > 700:
            f.write('$env:XTBEXE = "C:\\ORCA_6.1.1\\xTB_bigstack\\xtb.exe"  # 700原子超の stack overflow 対策\n')
        for nm in names:
            f.write(f'$o = "{nm}.out"\n'
                    f'if (-not ((Test-Path $o) -and (Select-String -Path $o -Pattern "ORCA TERMINATED NORMALLY" -Quiet))) {{\n'
                    f'  & $Orca "{nm}.inp" 2>&1 | Out-File -FilePath $o -Encoding utf8\n}}\n')
    with open(os.path.join(a.outdir, "selection.json"), "w", encoding="utf-8") as f:
        json.dump(dict(traj=a.xyz, frame=k, kept=keep, rmsd_heavy=res["rmsd_heavy"], rmsd_all=res["rmsd_all"],
                       v=list(map(float, c.v)), params=dict(max_as=a.max_as, min_proj=a.min_proj, clash=a.clash,
                       distances=a.distances), rows=res["rows"], n_atoms_cluster=n_tot,
                       e_shell_file=a.shell_out), f, ensure_ascii=False, indent=1)
    with open(os.path.join(a.outdir, "selection.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(res["rows"][0].keys()))
        w.writeheader()
        w.writerows(res["rows"])
    print(f"\n残した溶媒 {len(keep)}/{a.nsolv} 分子、クラスター {n_tot} 原子、出力 {len(names)} 入力 -> {a.outdir}")
    print("実行: " + os.path.join(a.outdir, "run.ps1") + "  (ORCA は各自 PowerShell で実行)")


def cmd_proxy(a, c):
    """新 solvation proxy 用の3入力 (既存ファイルは一切上書きしない):
      SHS_f{k}_Sall_XTB2_SP   : 溶媒 全 nsolv 分子のみ            (S_all)
      SHS_f{k}_ASkeep_XTB2_SP : A(結晶座標) + keep 溶媒            (A+S_keep)
      SHS_f{k}_Skeep_XTB2_SP  : keep 溶媒のみ                      (S_keep)
    座標・原子順は make の AS 入力と同一 ([A(59, 結晶A) | 溶媒(MD→結晶座標)])。
    溶媒は MD の元の分子順 (--keep は昇順に並べる)。凍結・XTB2 SP・nprocs は ref0 のヘッダーを流用。
    """
    frames = read_xyz_frames(a.xyz)
    k = len(frames) - 1 if a.frame is None else a.frame
    atoms = frames[k]
    res, err = analyse_frame(c, atoms, a, a.nsolv)
    if res is None:
        sys.exit(err)
    print(f"frame {k}: RMSD(heavy)={res['rmsd_heavy']:.3f} Å")
    if res["rmsd_heavy"] > a.max_rmsd and not a.force:
        sys.exit(f"RMSD(heavy) が {a.max_rmsd} Å を超える。--force で続行可 (非推奨)")
    keep_auto = [r["mol"] for r in res["rows"] if r["status"] == "keep"]
    keep = sorted(a.keep) if a.keep else keep_auto
    if any(m < 1 or m > a.nsolv for m in keep) or not keep:
        sys.exit(f"--keep が不正: {keep}")
    print(f"keep(使用)={keep}  keep(自動選別)={keep_auto}" + ("" if keep == keep_auto else "  <-- 不一致"))
    Y = res["Y"]
    sym = [x[0] for x in atoms]

    def sidx(mols):
        out = []
        for m in mols:
            out += list(range(N_A + (m - 1) * N_SOLV_ATOMS, N_A + m * N_SOLV_ATOMS))
        return out

    i_all, i_keep = sidx(range(1, a.nsolv + 1)), sidx(keep)
    symA = c.sym0[:N_A]
    jobs = [
        ("Sall", [sym[i] for i in i_all], Y[i_all], f"S_all ({a.nsolv} solvent, frame {k})"),
        ("ASkeep", symA + [sym[i] for i in i_keep], np.vstack([c.XA, Y[i_keep]]),
         f"A+S_keep (keep {keep}, frame {k})"),
        ("Skeep", [sym[i] for i in i_keep], Y[i_keep], f"S_keep (keep {keep}, frame {k})"),
    ]
    os.makedirs(a.outdir, exist_ok=True)
    made = {}
    for tag, s_, X_, label in jobs:
        nm = f"SHS_f{k:04d}_{tag}_XTB2_SP"
        pth = os.path.join(a.outdir, nm + ".inp")
        if os.path.exists(pth) and not a.force:
            sys.exit(f"既に存在 (上書きしない): {pth}")
        write_inp(pth, c.head, f"# solvation proxy input: {label}; same coordinates as the AS input", s_, X_)
        made[tag] = (nm + ".inp", s_, X_)
        print(f"  {nm}.inp  atoms={len(s_)}")

    # 整合性: 既存 AS 入力 / A_mono 入力 / 部分集合の座標一致
    def cmp(label, path, s_, X_):
        if not os.path.exists(path):
            print(f"  (比較対象なし) {label}: {path}")
            return
        _, at = read_inp(path)
        same_sym = [x[0] for x in at] == list(s_)
        dev = float(np.abs(arr(at) - X_).max()) if len(at) == len(s_) else float("nan")
        print(f"  {label}: n={len(at)} vs {len(s_)}  元素順一致={same_sym}  最大座標差={dev:.2e} Å")

    cmp("ASkeep vs 既存AS", os.path.join(a.outdir, f"SHS_f{k:04d}_AS_XTB2_SP.inp"), made["ASkeep"][1], made["ASkeep"][2])
    cmp("ASkeep[:59] vs A_mono", a.amono_inp, made["ASkeep"][1][:N_A], made["ASkeep"][2][:N_A])
    ok1 = made["Skeep"][1] == made["ASkeep"][1][N_A:] and bool(np.array_equal(made["Skeep"][2], made["ASkeep"][2][N_A:]))
    sel = [i_all.index(i) for i in i_keep]
    ok2 = [made["Sall"][1][i] for i in sel] == made["Skeep"][1] and bool(np.array_equal(made["Sall"][2][sel], made["Skeep"][2]))
    print(f"  Skeep == ASkeep[59:]: {ok1}   Skeep == Sall[keep]: {ok2}")
    print("ORCA は未実行。")


def energy(path):
    if not os.path.exists(path):
        return None
    last = None
    for ln in read_text(path).splitlines():
        m = PAT.search(ln)
        if m:
            last = float(m.group(1))
    return last


def cmd_collect(a, c):
    es = energy(a.shell_out)
    if es is None:
        sys.exit("E_shell が取れない: " + a.shell_out)
    e_as = None
    rows = []
    for fn in sorted(os.listdir(a.outdir)):
        if not fn.endswith(".out"):
            continue
        if "ORCA TERMINATED NORMALLY" not in read_text(os.path.join(a.outdir, fn)):
            print(f"(未完了/異常: {fn})")
            continue
        e = energy(os.path.join(a.outdir, fn))
        if "_AS_" in fn:
            e_as = e
        else:
            m = re.search(r"_r\+(\d+\.\d)A_", fn)
            if m:
                rows.append((float(m.group(1)), e, fn))
    if e_as is None:
        sys.exit("E[A+S] (…_AS_XTB2_SP.out) が未完了")
    vac = {}
    if os.path.exists(a.vac_csv):
        for r in csv.DictReader(open(a.vac_csv, encoding="utf-8-sig")):
            if r["r_A"].replace(".", "").replace("-", "").isdigit():
                vac[round(float(r["r_A"]), 1)] = float(r["E_bind_kcal_mol"])
    print(f"E_shell={es:.12f}  E[A+S]={e_as:.12f}")
    print(" r(Å)  E_bind_solv  E_bind_vac  ΔΔE   (kcal/mol)")
    for r, e, fn in sorted(rows):
        eb = (e - es - e_as) * HARTREE_TO_KCAL
        v = vac.get(round(r, 1))
        print(f"{r:5.1f} {eb:11.2f} " + (f"{v:10.2f} {eb - v:7.2f}" if v is not None else ""))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["scan", "make", "collect", "proxy"])
    ap.add_argument("xyz", nargs="?")
    ap.add_argument("--ref0", default="HS_cluster_+0.0A_SP_np1.inp")
    ap.add_argument("--ref5", default="HS_cluster_+5.0A_XTB2_SP.inp")
    ap.add_argument("--shell-out", dest="shell_out", default="HS_shell_XTB2_SP.out")
    ap.add_argument("--vac-csv", dest="vac_csv", default="HalfShell_pull_summary.csv")
    ap.add_argument("--nsolv", type=int, default=10)
    ap.add_argument("--frame", type=int, default=None, help="省略で最終フレーム")
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--step", type=int, default=10)
    ap.add_argument("--distances", type=float, nargs="+", default=[0.0, 2.0, 5.0, 8.0, 12.0, 20.0])
    ap.add_argument("--max-as", dest="max_as", type=float, default=6.0)
    ap.add_argument("--min-proj", dest="min_proj", type=float, default=0.0)
    ap.add_argument("--clash", type=float, default=2.2)
    ap.add_argument("--max-rmsd", dest="max_rmsd", type=float, default=1.0)
    ap.add_argument("--keep", type=int, nargs="+", default=None, help="proxy: keep する溶媒分子番号 (1始まり)。省略で自動選別")
    ap.add_argument("--amono-inp", dest="amono_inp", default="A_mono_XTB2_SP.inp", help="proxy: A 単独入力 (整合性確認用)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--outdir", default="SHS_out")
    a = ap.parse_args()
    if a.cmd != "collect" and not a.xyz:
        sys.exit("trajectory xyz を指定")
    c = load_context(a)
    {"scan": cmd_scan, "make": cmd_make, "collect": cmd_collect, "proxy": cmd_proxy}[a.cmd](a, c)


if __name__ == "__main__":
    main()
