#!/usr/bin/env python3
r"""
任意の軸分子 F についての3体項解析 (GFN2-xTB を直接実行)

  dE3(A, F, Dj) = E(A,F,Dj) - E(A,F) - E(A,Dj) - E(F,Dj) + E(A) + E(F) + E(Dj)

を、shell 内の全 Dj について、全エネルギーと成分(SCC, isoES, anisoES,
anisoXC, dispersion, repulsion, rest)ごとに計算し、その合計を
「F を shell から1個抜いたときの A 中心の非加算性の変化」(leave-one-out) と比べる。

  Δnonadd(F) = [E_bind(shell) - E_bind(shell \ F)] - dE2(A,F)
  E_bind(S)  = E(A+S) - E(S) - E(A)

F が最後に足した分子(D9)なら、これは Stage08->09 の増分と同じ量になる。
F が D3 のように途中の分子でも、同じ形式で定義できる。

使い方 (PowerShell):
  python manybody_trimers.py prepare --cluster <Stage09のcluster XYZ> --focus D9  --loo --outdir d9_trimers
  python manybody_trimers.py run     --focus D9  --loo --outdir d9_trimers --jobs 2
  python manybody_trimers.py analyze --focus D9  --loo --outdir d9_trimers

  python manybody_trimers.py prepare --cluster <Stage09のcluster XYZ> --focus D3  --loo --outdir d9_trimers
  python manybody_trimers.py run     --focus D3  --loo --outdir d9_trimers --jobs 2
  python manybody_trimers.py analyze --focus D3  --loo --outdir d9_trimers

  python manybody_trimers.py compare --focus D9 --vs D3 --outdir d9_trimers

  --focus には labels の名前を書く。D2, D3 は D2a, D3a に自動で読み替える。
  同じ --outdir を使うと、モノマーと A+Dj ダイマー(軸に依存しない計算)は再利用される。

  --loo を付けると、Stage09 全体(cluster/shell)と、F を抜いた cluster/shell の
  4ジョブ(885, 826, 826, 767原子相当)を追加で計算する。付けないと、
  F=D9 のときだけ、提示済みの Stage08->09 の増分と比較する。

前提: クラスターXYZは「1分子 = 連続59原子」で、ブロック0が A。
      ブロックの並びは --labels で与える(既定は追加順 A, D2a, D3a, D1a, ...)。
"""
import argparse
import csv
import math
import os
import re
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HARTREE_TO_KCAL = 627.5095
NAT = 59
DEFAULT_LABELS = "A,D2a,D3a,D1a,D1b,D5a,D5b,D7a,D7b,D4,D6a,D6b,D8a,D8b,D9"
DEFAULT_XTB = r"C:\ORCA_6.1.1\xTB_bigstack\xtb.exe"

# --- 既存のORCA経由の結果(ユーザー提示値)。ブロック対応と環境の検証に使う ---
REF_MONO_A = -96.82920663126
REF_PAIR = {  # E(A + Dj) [Eh]
    "D2a": -193.68574044537, "D3a": -193.68510754425,
    "D1a": -193.68183227145, "D1b": -193.68183223959,
    "D5a": -193.67100990812, "D5b": -193.67100988946,
    "D7a": -193.66802373609, "D7b": -193.66802370222,
    "D4": -193.66692068363,
    "D6a": -193.66426271929, "D6b": -193.66426261695,
    "D8a": -193.66117673282, "D8b": -193.66117672052,
    "D9": -193.66100195157,
}
REF_STAGE09 = {"c_full": -1453.01955644620, "s_full": -1356.05581635188}  # [Eh]
REF_STAGE08 = {"c_loo": -1356.12753790806, "s_loo": -1259.16255964158}    # [Eh] (F=D9のときのみ)
REF_STAGE02_NONADD_EH = -0.00063466341   # = dE3(A,D2a,D3a) の恒等式
# 以前の整理による重心距離 [Å](参考表示のみ。判定には使わない)
EXPECTED_CENTROID = {
    "D2a": 5.70, "D3a": 6.29, "D1a": 10.25, "D1b": 10.25,
    "D5a": 10.33, "D5b": 10.33, "D7a": 10.73, "D7b": 10.73,
    "D4": 11.74, "D6a": 13.84, "D6b": 13.84,
    "D8a": 13.90, "D8b": 13.90, "D9": 12.03,
}
# Stage08->Stage09 の実測増分 [kcal/mol](ユーザー提示値。--loo なしで F=D9 のときの比較用)
REF_STEP_D9 = {"tot": 2.4014, "scc": 2.4014, "iso": 0.3153, "aes": 0.2255,
               "axc": 0.0008, "disp": 1.8945, "rep": 0.0}

PATTERNS = {
    "tot":  r"^\s*::\s+total energy\s+(-?\d+\.\d+)\s+Eh",
    "scc":  r"^\s*::\s+SCC energy\s+(-?\d+\.\d+)\s+Eh",
    "iso":  r"^\s*::\s+->\s+isotropic ES\s+(-?\d+\.\d+)\s+Eh",
    "aes":  r"^\s*::\s+->\s+anisotropic ES\s+(-?\d+\.\d+)\s+Eh",
    "axc":  r"^\s*::\s+->\s+anisotropic XC\s+(-?\d+\.\d+)\s+Eh",
    "disp": r"^\s*::\s+->\s+dispersion\s+(-?\d+\.\d+)\s+Eh",
    "rep":  r"^\s*::\s+repulsion energy\s+(-?\d+\.\d+)\s+Eh",
}
KEYS = ["tot", "scc", "iso", "aes", "axc", "disp", "rep", "rest"]
HEAD = {"tot": "Total", "scc": "SCC", "iso": "IsoES", "aes": "AnisoES",
        "axc": "AnisoXC", "disp": "Disp", "rep": "Repul", "rest": "rest"}
OK_TEXT = "normal termination of xtb"


# ------------------------------------------------------------------ helpers
def resolve_focus(labels, focus):
    if focus == "A":
        sys.exit("--focus に A は指定できません")
    if focus in labels:
        return focus
    if focus + "a" in labels:
        print(f"(注) --focus {focus} を {focus}a に読み替えました")
        return focus + "a"
    sys.exit(f"--focus {focus} が labels にありません: {','.join(labels)}")


def read_blocks(path, labels, nat):
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    n = int(lines[0].split()[0])
    atoms = [l.rstrip() for l in lines[2:2 + n]]
    if len(atoms) != n:
        sys.exit(f"XYZの原子数が合いません: 宣言 {n}, 実際 {len(atoms)}")
    if n != nat * len(labels):
        sys.exit(f"原子数 {n} が 分子数({len(labels)}) x {nat} = {nat*len(labels)} と一致しません。"
                 f" --labels / --natoms を確認してください")
    blocks = [atoms[i * nat:(i + 1) * nat] for i in range(len(labels))]
    # 元素組成だけを比べる。原子の並び順は問わず、XYZの座標・原子順はそのまま使う
    el0 = Counter(l.split()[0] for l in blocks[0])
    for i, b in enumerate(blocks):
        if Counter(l.split()[0] for l in b) != el0:
            sys.exit(f"ブロック{i}({labels[i]})の元素組成がブロック0と違います。"
                     f" 分子の境界(59原子ずつ)が合っていない可能性があります")
    return blocks


def centroid(block):
    xs = ys = zs = 0.0
    for l in block:
        p = l.split()
        xs += float(p[1]); ys += float(p[2]); zs += float(p[3])
    n = len(block)
    return (xs / n, ys / n, zs / n)


def build_jobs(labels, focus, loo):
    idx = {l: i for i, l in enumerate(labels)}
    A, F = idx["A"], idx[focus]
    partners = [l for l in labels if l not in ("A", focus)]
    jobs = {"m_A": [A], f"m_{focus}": [F], f"p_A_{focus}": [A, F]}
    for l in partners:
        j = idx[l]
        jobs[f"m_{l}"] = [j]
        jobs[f"p_A_{l}"] = [A, j]
        jobs[f"d_{focus}_{l}"] = [F, j]
        jobs[f"t_A_{focus}_{l}"] = [A, F, j]
    if loo:  # 大きいジョブは最後
        everything = list(range(len(labels)))
        shell = [i for i in everything if i != A]
        jobs["c_full"] = everything
        jobs["s_full"] = shell
        jobs[f"c_loo_{focus}"] = [i for i in everything if i != F]
        jobs[f"s_loo_{focus}"] = [i for i in shell if i != F]
    return jobs, partners


def job_dir(outdir, name):
    return Path(outdir) / "jobs" / name


# ------------------------------------------------------------------ prepare
def cmd_prepare(a):
    labels = a.labels.split(",")
    focus = resolve_focus(labels, a.focus)
    blocks = read_blocks(a.cluster, labels, a.natoms)
    jobs, partners = build_jobs(labels, focus, a.loo)
    c0 = centroid(blocks[0])
    print("ブロック対応表 (ブロック0 = A からの重心距離)")
    print(f"{'#':>3} {'label':<5} {'dist(Å)':>8} {'前回の目安':>10} {'差':>7}")
    for i, (lab, b) in enumerate(zip(labels, blocks)):
        d = math.dist(centroid(b), c0)
        exp = EXPECTED_CENTROID.get(lab)
        mark = "  <- 軸" if lab == focus else ""
        if exp is None:
            print(f"{i+1:>3} {lab:<5} {d:8.2f}{mark}")
        else:
            print(f"{i+1:>3} {lab:<5} {d:8.2f} {exp:10.2f} {d-exp:+7.2f}{mark}")
    print("\n注意: 距離が目安と大きく違う行は、並び(--labels)が違う可能性があります。")
    print("      最終的な対応の判定は analyze の「ペアエネルギー検証」で行います。")
    for name, members in jobs.items():
        atoms = [l for m in members for l in blocks[m]]
        d = job_dir(a.outdir, name)
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{name}.xyz").write_text(
            f"{len(atoms)}\n{name}\n" + "\n".join(atoms) + "\n", encoding="utf-8")
    print(f"\n軸 = {focus}: {len(jobs)} 個のジョブ入力を {Path(a.outdir)/'jobs'} に作成しました。")
    print(f"  トリマー {len(partners)}, {focus}+Dj ダイマー {len(partners)}, "
          f"A+Dj ダイマー {len(partners)+1}, モノマー {len(partners)+2}"
          + (", 大きいクラスター 4 (leave-one-out用)" if a.loo else ""))


# ------------------------------------------------------------------ run
def run_job(a, name):
    d = job_dir(a.outdir, name)
    xyz = d / f"{name}.xyz"
    out = d / f"{name}.out"
    if not xyz.exists():
        return name, "NO_INPUT (先に prepare を実行)"
    if out.exists() and OK_TEXT in out.read_text(errors="replace"):
        return name, "skip(完了済み)"
    xtb = str(Path(a.xtb).resolve()) if Path(a.xtb).exists() else a.xtb
    cmd = [xtb, xyz.name, "-c", "0", "-u", "0", "-P", str(a.threads),
           "--acc", "1.0", "--namespace", name]
    if a.input_tmp:
        cmd += ["--input", str(Path(a.input_tmp).resolve())]
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
    try:
        with open(out, "wb") as f:
            rc = subprocess.run(cmd, cwd=d, stdout=f, stderr=subprocess.STDOUT, env=env).returncode
    except FileNotFoundError:
        return name, f"FAILED(xtbが見つかりません: {a.xtb})"
    ok = OK_TEXT in out.read_text(errors="replace")
    return name, "ok" if ok else f"FAILED(rc={rc})"


def cmd_run(a):
    labels = a.labels.split(",")
    focus = resolve_focus(labels, a.focus)
    jobs, _ = build_jobs(labels, focus, a.loo)
    names = list(jobs)
    print(f"軸 = {focus}: {len(names)} ジョブを実行します (同時実行 {a.jobs}, xtbスレッド {a.threads})")
    failed = []
    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for k, (name, st) in enumerate(ex.map(lambda n: run_job(a, n), names), 1):
            print(f"[{k:>2}/{len(names)}] {name:<28} {st}")
            if st.startswith(("FAILED", "NO_INPUT")):
                failed.append(name)
    if failed:
        print("\n失敗したジョブ:", ", ".join(failed))
        print("各 jobs\\<name>\\<name>.out の末尾を確認してください。")
        sys.exit(1)
    print("\n全ジョブ完了。次: analyze")


# ------------------------------------------------------------------ analyze
def parse_job(outdir, name):
    p = job_dir(outdir, name) / f"{name}.out"
    if not p.exists():
        sys.exit(f"出力がありません: {p}\n(--loo を付けたなら prepare と run も --loo 付きで実行してください)")
    t = p.read_text(errors="replace")
    if OK_TEXT not in t:
        sys.exit(f"正常終了していません: {p}")
    e = {}
    for k, rx in PATTERNS.items():
        m = re.findall(rx, t, flags=re.M)
        if not m:
            sys.exit(f"{p} に '{k}' の行が見つかりません。サマリー部分を数行確認してください")
        e[k] = float(m[-1])
    e["rest"] = e["scc"] - e["iso"] - e["aes"] - e["axc"] - e["disp"]
    return e


def delta3(E, focus, l):
    t, pf, pl = E[f"t_A_{focus}_{l}"], E[f"p_A_{focus}"], E[f"p_A_{l}"]
    d, mA, mF, mL = E[f"d_{focus}_{l}"], E["m_A"], E[f"m_{focus}"], E[f"m_{l}"]
    return {k: (t[k] - pf[k] - pl[k] - d[k] + mA[k] + mF[k] + mL[k]) * HARTREE_TO_KCAL
            for k in KEYS}


def marginal_nonadd(E, focus):
    """F を shell から抜いたときの A 中心の非加算性の変化 [kcal/mol]"""
    out = {}
    for k in KEYS:
        dbind = E["c_full"][k] - E["s_full"][k] - E[f"c_loo_{focus}"][k] + E[f"s_loo_{focus}"][k]
        pair2 = E[f"p_A_{focus}"][k] - E[f"m_{focus}"][k] - E["m_A"][k]
        out[k] = (dbind - pair2) * HARTREE_TO_KCAL
    return out


def line(name, d, width=10):
    return f"{name:<8}" + "".join(f"{d[k]:+{width}.4f}" for k in KEYS)


def cmd_analyze(a):
    labels = a.labels.split(",")
    focus = resolve_focus(labels, a.focus)
    jobs, partners = build_jobs(labels, focus, a.loo)
    E = {n: parse_job(a.outdir, n) for n in jobs}
    tol = a.tol

    # ---- 検証: 既存のORCA経由の値との一致
    print(f"== 軸 = {focus} ==\n")
    print("== 検証: ORCA経由の既存値との比較 (Eh) ==")
    bad = 0

    def chk(label, diff):
        nonlocal bad
        flag = "" if abs(diff) <= tol else "  <-- 不一致"
        bad += abs(diff) > tol
        print(f"  {label:<18} 差 {diff:+.2e}{flag}")

    chk("mono A", E["m_A"]["tot"] - REF_MONO_A)
    for l in partners + [focus]:
        if l in REF_PAIR:
            chk(f"pair A+{l}", E[f"p_A_{l}"]["tot"] - REF_PAIR[l])
    if a.loo:
        for n, ref in REF_STAGE09.items():
            chk(f"{n} (Stage09)", E[n]["tot"] - ref)
        if focus == "D9":
            chk("Stage08 cluster", E["c_loo_D9"]["tot"] - REF_STAGE08["c_loo"])
            chk("Stage08 shell", E["s_loo_D9"]["tot"] - REF_STAGE08["s_loo"])
    if bad:
        print(f"\n!! {bad} 件が許容差 {tol:g} Eh を超えています。")
        print("   原因候補: (1) --labels の並びが実際のブロックと違う, (2) xtb の設定が")
        print("   ORCA経由と違う (--input-tmp で ORCA の .input.tmp を渡す), (4) 別のクラスター。")
        print("   このまま以降の結果を解釈しないでください。\n")
    else:
        print("  -> 全て許容差内。ブロック対応と計算条件は既存の結果と整合しています。\n")

    # ---- 各 Dj の3体項
    rows = [(l, delta3(E, focus, l)) for l in partners]
    tot = {k: sum(r[k] for _, r in rows) for k in KEYS}

    print(f"== 3体項 dE3(A,{focus},Dj)  [kcal/mol] ==")
    print(f"{'Dj':<8}" + "".join(f"{HEAD[k]:>10}" for k in KEYS))
    for l, r in rows:
        print(line(l, r))
    print("-" * (8 + 10 * len(KEYS)))
    print(line("合計", tot))

    # ---- 恒等式チェック: dE3(A,D2a,D3a) = Stage02 の E_nonadd
    if focus in ("D2a", "D3a") and "D2a" in labels and "D3a" in labels:
        other = "D3a" if focus == "D2a" else "D2a"
        r = dict(rows)[other]
        d_eh = r["tot"] / HARTREE_TO_KCAL - REF_STAGE02_NONADD_EH
        flag = "" if abs(d_eh) <= tol else "  <-- 不一致"
        print(f"\n恒等式チェック: dE3(A,D2a,D3a) = {r['tot']:+.4f} kcal/mol"
              f"  (Stage02 の E_nonadd -0.3983 kcal/mol との差 {d_eh:+.2e} Eh){flag}")

    # ---- 比較対象(実測増分)
    ref, src = None, None
    if a.loo:
        ref = marginal_nonadd(E, focus)
        src = f"leave-one-out: {focus}を除いた{len(labels)-1}分子クラスターとの差"
    elif focus == "D9":
        ref = dict(REF_STEP_D9)
        ref["rest"] = ref["scc"] - ref["iso"] - ref["aes"] - ref["axc"] - ref["disp"]
        src = "Stage08->09 の提示値"
    print()
    if ref is None:
        print("(比較対象の増分がありません。--loo を付けて prepare / run / analyze をやり直してください)")
    else:
        print(f"== 比較: 3体項の合計 vs 実測増分 [kcal/mol]  ({src}) ==")
        print(f"{'':<8}" + "".join(f"{HEAD[k]:>10}" for k in KEYS))
        print(line("3体合計", tot))
        print(line("実測", ref))
        print(line("残差", {k: ref[k] - tot[k] for k in KEYS}) + "   <- 実測 - 3体合計 (4体以上・その他)")
        print(f"{'3体/実測':<8}" + "".join(
            (f"{tot[k]/ref[k]:10.2f}" if abs(ref[k]) > 1e-3 else f"{'--':>10}") for k in KEYS))
        if a.loo and focus == "D9":
            print("\nStage08->09 提示値との差 (leave-one-out の計算が正しいことの確認):")
            diff = {k: ref[k] - {**REF_STEP_D9, "rest": REF_STEP_D9["scc"] - REF_STEP_D9["iso"]
                                 - REF_STEP_D9["aes"] - REF_STEP_D9["axc"] - REF_STEP_D9["disp"]}[k]
                    for k in KEYS}
            print(line("差", diff))
    print("\nrest = SCC - (IsoES+AnisoES+AnisoXC+Disp)。xtbのサマリーに個別表示されない項の合計。")

    # ---- 寄与の大きい Dj
    for key in ("tot", "disp"):
        top = sorted(rows, key=lambda x: -abs(x[1][key]))[:3]
        print(f"{HEAD[key]} への寄与が大きい上位3: " + ", ".join(f"{l} ({r[key]:+.3f})" for l, r in top))

    # ---- CSV
    out = Path(a.outdir) / f"{focus}_trimer_3body_components.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Dj"] + [f"{HEAD[k]}_kcal" for k in KEYS])
        for l, r in rows:
            w.writerow([l] + [f"{r[k]:.6f}" for k in KEYS])
        w.writerow(["SUM"] + [f"{tot[k]:.6f}" for k in KEYS])
    sm = Path(a.outdir) / f"{focus}_summary.csv"
    with open(sm, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["key", "sum3", "ref", "resid", "ratio", "ref_source"])
        for k in KEYS:
            if ref is None:
                w.writerow([k, f"{tot[k]:.6f}", "", "", "", ""])
            else:
                ratio = f"{tot[k]/ref[k]:.4f}" if abs(ref[k]) > 1e-3 else ""
                w.writerow([k, f"{tot[k]:.6f}", f"{ref[k]:.6f}", f"{ref[k]-tot[k]:.6f}", ratio, src])
    print(f"CSV: {out}\n     {sm}")


# ------------------------------------------------------------------ compare
def read_summary(outdir, focus):
    p = Path(outdir) / f"{focus}_summary.csv"
    if not p.exists():
        sys.exit(f"{p} がありません。先に {focus} で analyze を実行してください")
    with open(p, newline="", encoding="utf-8-sig") as f:
        return {r["key"]: r for r in csv.DictReader(f)}


def cmd_compare(a):
    labels = a.labels.split(",")
    f1 = resolve_focus(labels, a.focus)
    if not a.vs:
        sys.exit("compare には --vs が必要です")
    f2 = resolve_focus(labels, a.vs)
    s1, s2 = read_summary(a.outdir, f1), read_summary(a.outdir, f2)

    def val(s, key, col):
        v = s[key][col]
        return float(v) if v not in ("", None) else None

    def cell(v, fmt):
        return f"{'--':>10}" if v is None else f"{v:{fmt}}".rjust(10)

    print(f"== 軸の比較 [kcal/mol]: {f1} vs {f2} ==")
    print(f"{'成分':<9}{'項目':<10}{f1:>10}{f2:>10}")
    for key in ("tot", "disp", "iso", "aes", "rest"):
        for label, col, fmt in (("3体合計", "sum3", "+.4f"), ("実測", "ref", "+.4f"),
                                ("残差", "resid", "+.4f"), ("3体/実測", "ratio", ".2f")):
            nm = HEAD[key] if label == "3体合計" else ""
            print(f"{nm:<9}{label:<10}" + cell(val(s1, key, col), fmt) + cell(val(s2, key, col), fmt))
        print()
    print(f"参照: {f1}: {s1['tot']['ref_source']}\n      {f2}: {s2['tot']['ref_source']}")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["prepare", "run", "analyze", "compare"])
    ap.add_argument("--cluster", help="prepare: Stage09 cluster のXYZ")
    ap.add_argument("--outdir", default="trimers")
    ap.add_argument("--labels", default=DEFAULT_LABELS, help="ブロック順のラベル(カンマ区切り, 先頭は A)")
    ap.add_argument("--focus", default="D9", help="軸にする分子のラベル (D2,D3 は D2a,D3a に読み替え)")
    ap.add_argument("--vs", default=None, help="compare: 比較する2つ目の軸")
    ap.add_argument("--loo", action="store_true",
                    help="leave-one-out の実測増分を計算する (大きいクラスター4ジョブを追加)")
    ap.add_argument("--natoms", type=int, default=NAT)
    ap.add_argument("--xtb", default=DEFAULT_XTB)
    ap.add_argument("--input-tmp", default=None, help="任意: ORCAが作った xtb の .input.tmp")
    ap.add_argument("--threads", type=int, default=1, help="xtb 1ジョブあたりのスレッド数")
    ap.add_argument("--jobs", type=int, default=1, help="同時に走らせるジョブ数")
    ap.add_argument("--tol", type=float, default=5e-6, help="検証の許容差 [Eh]")
    a = ap.parse_args()
    if a.command == "prepare":
        if not a.cluster:
            sys.exit("prepare には --cluster が必要です")
        cmd_prepare(a)
    elif a.command == "run":
        cmd_run(a)
    elif a.command == "analyze":
        cmd_analyze(a)
    else:
        cmd_compare(a)


if __name__ == "__main__":
    main()
