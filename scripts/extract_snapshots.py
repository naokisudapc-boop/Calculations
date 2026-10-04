#!/usr/bin/env python3
"""MDトラジェクトリからスナップショットを取り出し、ORCA入力(A+nS / A / nS)を作る (Step 2 準備)。

前提は analyze_md_trajectory.py と同じ: 各フレームは Apixaban 59原子 -> 溶媒 n 分子(連続)。
同じ scripts フォルダに analyze_md_trajectory.py が必要。

使い方:
  python extract_snapshots.py traj.xyz --solvent water --outdir snap_water --step 2
  python extract_snapshots.py traj.xyz --solvent etoh --outdir snap_etoh --start 100 --end 200 --step 2

出力 (outdir/f0000/): f0000_AnS.xyz/.inp, f0000_A.xyz/.inp, f0000_nS.xyz/.inp と outdir/run_all.ps1
全て同一ジオメトリ(MD構造そのまま、最適化なし)の単一点。電荷0・一重項を仮定(要確認)。
"""
import argparse
import os

from analyze_md_trajectory import N_APIXABAN, SOLVENT_ATOMS, read_xyz_frames


def write_xyz(path, atoms, comment):
    with open(path, "w", newline="\n") as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for s, x, y, z in atoms:
            f.write(f"{s:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def write_inp(path, xyzname, method, charge, mult, nprocs):
    with open(path, "w", newline="\n") as f:
        f.write(f"! {method}\n")
        if nprocs > 1:
            f.write(f"%pal nprocs {nprocs} end\n")
        f.write(f"* xyzfile {charge} {mult} {xyzname}\n")


PS1 = r'''# Run all ORCA inputs under this folder; skip jobs that already terminated normally.
param([string]$Orca = "{orca}")
Get-ChildItem -Path $PSScriptRoot -Recurse -Filter *.inp | Sort-Object FullName | ForEach-Object {{
    $out = [IO.Path]::ChangeExtension($_.FullName, ".out")
    $done = (Test-Path $out) -and (Select-String -Path $out -Pattern "ORCA TERMINATED NORMALLY" -Quiet)
    if (-not $done) {{
        Push-Location $_.DirectoryName
        & $Orca $_.Name 2>&1 | Out-File -FilePath $out -Encoding utf8
        Pop-Location
    }}
}}
'''


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--solvent", required=True, choices=SOLVENT_ATOMS)
    ap.add_argument("--nsolv", type=int, default=10)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--end", type=int, default=None, help="最後のフレーム(含む)。省略で最終フレーム")
    ap.add_argument("--step", type=int, default=1)
    ap.add_argument("--method", default="XTB2", help="ORCAのシンプル入力キーワード")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--mult", type=int, default=1)
    ap.add_argument("--nprocs", type=int, default=1)
    ap.add_argument("--orca", default=r"C:\ORCA_6.1.1\orca.exe", help="run_all.ps1 に埋め込む orca.exe のパス")
    a = ap.parse_args()

    ns = SOLVENT_ATOMS[a.solvent]
    expected = N_APIXABAN + ns * a.nsolv
    frames = read_xyz_frames(a.xyz)
    end = len(frames) - 1 if a.end is None else min(a.end, len(frames) - 1)
    sel = range(a.start, end + 1, a.step)
    os.makedirs(a.outdir, exist_ok=True)

    n = 0
    for k in sel:
        atoms = frames[k][1]
        if len(atoms) != expected:
            print(f"frame {k}: 原子数 {len(atoms)} != {expected} -> スキップ")
            continue
        d = os.path.join(a.outdir, f"f{k:04d}")
        os.makedirs(d, exist_ok=True)
        parts = {"AnS": atoms, "A": atoms[:N_APIXABAN], "nS": atoms[N_APIXABAN:]}
        for tag, at in parts.items():
            base = f"f{k:04d}_{tag}"
            write_xyz(os.path.join(d, base + ".xyz"), at, f"{a.xyz} frame {k} {tag}")
            write_inp(os.path.join(d, base + ".inp"), base + ".xyz", a.method, a.charge, a.mult, a.nprocs)
        n += 1

    with open(os.path.join(a.outdir, "run_all.ps1"), "w", newline="\r\n") as f:
        f.write(PS1.format(orca=a.orca))
    print(f"{n} フレーム分 ({n * 3} ジョブ) を {a.outdir} に出力。実行: .\\run_all.ps1")


if __name__ == "__main__":
    main()
