#!/usr/bin/env python3
r"""
D9軸の3体項解析 (GFN2-xTB を直接実行)

  dE3(A, F, Dj) = E(A,F,Dj) - E(A,F) - E(A,Dj) - E(F,Dj) + E(A) + E(F) + E(Dj)

を、F(=D9) と shell 内の全 Dj について、全エネルギーと成分(SCC, isoES,
anisoES, anisoXC, dispersion, repulsion)ごとに計算する。

使い方 (PowerShell):
  python d9_trimers.py prepare --cluster <Stage09のcluster XYZ> --outdir d9_trimers
  python d9_trimers.py run     --outdir d9_trimers
  python d9_trimers.py analyze --outdir d9_trimers

  1) prepare : XYZを切り出して jobs\<name>\<name>.xyz を作り、ブロック対応表を表示する
               (対応表を目で確認してから run に進むこと)
  2) run     : 各ジョブを xtb で実行する。完走済みのジョブは自動でスキップする
  3) analyze : 成分別の3体項、合計、Stage08->09の実測増分との比較、検証結果を表示し
               CSV を書き出す

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
# 以前の整理による重心距離 [Å](参考表示のみ。判定には使わない)
EXPECTED_CENTROID = {
    "D2a": 5.70, "D3a": 6.29, "D1a": 10.25, "D1b": 10.25,
    "D5a": 10.33, "D5b": 10.33, "D7a": 10.73, "D7b": 10.73,
    "D4": 11.74, "D6a": 13.84, "D6b": 13.84,
    "D8a": 13.90, "D8b": 13.90, "D9": 12.03,
}
# Stage08->Stage09 の実測増分 [kcal/mol](ユーザー提示値。focus が D9 のときのみ比較)
REF_STEP = {"tot": 2.4014, "scc": 2.4014, "iso": 0.3153, "aes": 0.2255,
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


# ------------------------------------------------------------------ XYZ
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
    from collections import Counter
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


# ------------------------------------------------------------------ jobs
def build_jobs(labels, focus):
    if "A" not in labels or focus not in labels:
        sys.exit("labels に 'A' と focus が必要です")
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
    return jobs, partners


def job_dir(outdir, name):
    return Path(outdir) / "jobs" / name


# ------------------------------------------------------------------ prepare
def cmd_prepare(a):
    labels = a.labels.split(",")
    blocks = read_blocks(a.cluster, labels, a.natoms)
    jobs, partners = build_jobs(labels, a.focus)
    c0 = centroid(blocks[0])
    print("ブロック対応表 (ブロック0 = A からの重心距離)")
    print(f"{'#':>3} {'label':<5} {'dist(Å)':>8} {'前回の目安':>10} {'差':>7}")
    for i, (lab, b) in enumerate(zip(labels, blocks)):
        c = centroid(b)
        d = math.dist(c, c0)
        exp = EXPECTED_CENTROID.get(lab)
        if exp is None:
            print(f"{i+1:>3} {lab:<5} {d:8.2f}")
        else:
            print(f"{i+1:>3} {lab:<5} {d:8.2f} {exp:10.2f} {d-exp:+7.2f}")
    print("\n注意: 距離が目安と大きく違う行は、並び(--labels)が違う可能性があります。")
    print("      最終的な対応の判定は analyze の「ペアエネルギー検証」で行います。")
    for name, members in jobs.items():
        atoms = [l for m in members for l in blocks[m]]
        d = job_dir(a.outdir, name)
        d.mkdir(parents=True, exist_ok=True)
        (d / f"{name}.xyz").write_text(
            f"{len(atoms)}\n{name}\n" + "\n".join(atoms) + "\n", encoding="utf-8")
    print(f"\n{len(jobs)} 個のジョブ入力を {Path(a.outdir)/'jobs'} に作成しました。")
    print(f"  モノマー {len(partners)+2}, A+Dj ダイマー {len(partners)+1}, "
          f"{a.focus}+Dj ダイマー {len(partners)}, トリマー {len(partners)}")


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
    jobs, _ = build_jobs(labels, a.focus)
    names = list(jobs)
    print(f"{len(names)} ジョブを実行します (同時実行 {a.jobs}, xtbスレッド {a.threads})")
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
        sys.exit(f"出力がありません: {p}")
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


def cmd_analyze(a):
    labels = a.labels.split(",")
    jobs, partners = build_jobs(labels, a.focus)
    E = {n: parse_job(a.outdir, n) for n in jobs}

    # ---- 検証: 既存のORCA経由の値との一致
    print("== 検証: ORCA経由の既存値との比較 (Eh) ==")
    tol = a.tol
    bad = 0
    dA = E["m_A"]["tot"] - REF_MONO_A
    print(f"  mono A          差 {dA:+.2e}")
    bad += abs(dA) > tol
    for l in partners + [a.focus]:
        if l in REF_PAIR:
            d = E[f"p_A_{l}"]["tot"] - REF_PAIR[l]
            flag = "" if abs(d) <= tol else "  <-- 不一致"
            bad += abs(d) > tol
            print(f"  pair A+{l:<4}    差 {d:+.2e}{flag}")
    if bad:
        print(f"\n!! {bad} 件が許容差 {tol:g} Eh を超えています。")
        print("   原因候補: (1) --labels の並びが実際のブロックと違う, (2) xtb の設定が")
        print("   ORCA経由と違う (--input-tmp で ORCA の .input.tmp を渡す), (3) 別のクラスター。")
        print("   このまま以降の結果を解釈しないでください。\n")
    else:
        print("  -> 全て許容差内。ブロック対応と計算条件は既存の結果と整合しています。\n")

    # ---- 各 Dj の3体項
    rows = []
    for l in partners:
        r = delta3(E, a.focus, l)
        rows.append((l, r))
    tot = {k: sum(r[k] for _, r in rows) for k in KEYS}

    print(f"== 3体項 dE3(A,{a.focus},Dj)  [kcal/mol] ==")
    print(f"{'Dj':<6}" + "".join(f"{HEAD[k]:>10}" for k in KEYS))
    for l, r in rows:
        print(f"{l:<6}" + "".join(f"{r[k]:+10.4f}" for k in KEYS))
    print("-" * (6 + 10 * len(KEYS)))
    print(f"{'合計':<6}" + "".join(f"{tot[k]:+10.4f}" for k in KEYS))

    if a.focus == "D9":
        ref = dict(REF_STEP)
        ref["rest"] = ref["scc"] - ref["iso"] - ref["aes"] - ref["axc"] - ref["disp"]
        print(f"{'実測':<6}" + "".join(f"{ref[k]:+10.4f}" for k in KEYS) + "   <- Stage08->09 の増分")
        print(f"{'残り':<6}" + "".join(f"{ref[k]-tot[k]:+10.4f}" for k in KEYS) + "   <- 4体以上・その他")
        print(f"{'3体/実測':<6}" + "".join(
            (f"{tot[k]/ref[k]:10.2f}" if abs(ref[k]) > 1e-3 else f"{'--':>10}") for k in KEYS))
    print("\nrest = SCC - (IsoES+AnisoES+AnisoXC+Disp)。xtbのサマリーに個別表示されない項の合計。")

    out = Path(a.outdir) / f"{a.focus}_trimer_3body_components.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(["Dj"] + [f"{HEAD[k]}_kcal" for k in KEYS])
        for l, r in rows:
            w.writerow([l] + [f"{r[k]:.6f}" for k in KEYS])
        w.writerow(["SUM"] + [f"{tot[k]:.6f}" for k in KEYS])
    print(f"CSV: {out}")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["prepare", "run", "analyze"])
    ap.add_argument("--cluster", help="prepare: Stage09 cluster のXYZ")
    ap.add_argument("--outdir", default="d9_trimers")
    ap.add_argument("--labels", default=DEFAULT_LABELS, help="ブロック順のラベル(カンマ区切り, 先頭は A)")
    ap.add_argument("--focus", default="D9", help="軸にする分子のラベル")
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
    else:
        cmd_analyze(a)


if __name__ == "__main__":
    main()