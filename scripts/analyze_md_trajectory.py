#!/usr/bin/env python3
"""MD trajectory(髫阪・辟夂ｹ晁ｼ釆樒ｹ晢ｽｼ郢晢｣ｰxyz)邵ｺ・ｮ闔蜿･辯暮囓・｣隴ｫ繝ｻ(Step 1) v2邵ｲ繝ｻ
陷題ざ鄂ｲ(髫補悪・｢・ｺ髫ｱ繝ｻ:
  - 陷ｷ繝ｻ繝ｵ郢晢ｽｬ郢晢ｽｼ郢晢｣ｰ邵ｺ・ｮ陷ｴ貅ｷ・ｭ螳｣・ｰ繝ｻ繝ｻ Apixaban(59陷ｴ貅ｷ・ｭ繝ｻ -> 雋・ｽｶ陝占ｲ槭・陝・・・・n 陋溘・(鬨ｾ・｣驍ｯ螢ｹ繝ｻ陷ｷ蠕｡・ｸﾂ鬯・・
  - 雋・ｽｶ陝舌・陋ｻ繝ｻ・ｭ闊後・陷ｴ貅ｷ・ｭ蜈育・邵ｺ・ｯ --solvent 邵ｺ・ｧ雎趣ｽｺ邵ｺ・ｾ郢ｧ繝ｻ(water=3, meoh=6, mecn=6, etoh=9, dmso=10)

闖ｴ・ｿ邵ｺ繝ｻ蟀ｿ:
  python analyze_md_trajectory.py traj.xyz --solvent water --out water_traj.csv
  python analyze_md_trajectory.py traj.xyz --solvent water --center 0 0 0 --out water_traj.csv

陷・ｽｺ陷我ｽ朶V(郢晁ｼ釆樒ｹ晢ｽｼ郢晢｣ｰ邵ｺ譁絶・):
  frame, n_atoms_ok, min_A_S (隴崢髴醍ｬｬ逎・恪譎槫ｱｬ), molN_min (陷ｷ繝ｻ・ｺ・ｶ陝占ｲ槭・陝・・竊但邵ｺ・ｮ隴崢陝・ｸ樊ｬ｡陝・ｮ｣菫｣髴肴辨螻ｬ),
  n_in_shell (--cutoff 闔会ｽ･陷繝ｻ繝ｻ陋ｻ繝ｻ・ｭ蜈育・), n_in_shell_<c> (--extra-cutoffs 邵ｺ譁絶・),
  n_far (--far ・・・郢ｧ蛹ｻ・企ｩ包｣ｰ邵ｺ繝ｻ繝ｻ陝・・辟・, max_mol_min_dist,
  A_centroid_to_cluster_centroid (A邵ｺ・ｮ陟趣ｽｧ隶灘雀・ｹ・ｳ陜ｮ繝ｻ- 陷茨ｽｨ陷ｴ貅ｷ・ｭ闊後・陟趣ｽｧ隶灘雀・ｹ・ｳ陜ｮ繝ｻ 隴鯉ｽｧ A_com_to_center),
  A_centroid_to_sphere_center (--center 隰悶・・ｮ螢ｽ蜃ｾ邵ｺ・ｮ邵ｺ・ｿ: A邵ｺ・ｮ陟趣ｽｧ隶灘雀・ｹ・ｳ陜ｮ繝ｻ- 隲｡菫ｶ謫夐・・・ｸ・ｭ陟｢繝ｻ

雎包ｽｨ隲｢繝ｻ 鬩･讎奇ｽｿ繝ｻ髮会ｽｪ鬩･荳槫・鬩･繝ｻ邵ｺ・ｧ邵ｺ・ｯ邵ｺ・ｪ邵ｺ荳橸ｽｺ・ｧ隶灘生繝ｻ陷雁｡・ｴ豕鯉ｽｹ・ｳ陜ｮ繝ｻ(centroid)邵ｲ繝ｻ"""
import argparse
import csv
import math
import sys

N_APIXABAN = 59
SOLVENT_ATOMS = {"water": 3, "meoh": 6, "mecn": 6, "etoh": 9, "dmso": 10, "dmf": 12, "acetone": 10, "chloroform": 5}


def read_xyz_frames(path):
    with open(path) as f:
        lines = f.read().splitlines()
    i, frames = 0, []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1] if i + 1 < len(lines) else ""
        atoms = []
        for ln in lines[i + 2 : i + 2 + n]:
            p = ln.split()
            atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
        if len(atoms) != n:
            print(f"Frame truncated: expected {n} atoms, got {len(atoms)}", file=sys.stderr)
            break
        frames.append((comment, atoms))
        i += 2 + n
    return frames


def dist(a, b):
    return math.sqrt((a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2)


def centroid(atoms):
    n = len(atoms)
    return (sum(x[1] for x in atoms) / n, sum(x[2] for x in atoms) / n, sum(x[3] for x in atoms) / n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz")
    ap.add_argument("--solvent", required=True, choices=SOLVENT_ATOMS)
    ap.add_argument("--nsolv", type=int, default=10)
    ap.add_argument("--cutoff", type=float, default=4.0, help="cutoff for n_in_shell")
    ap.add_argument("--extra-cutoffs", type=float, nargs="*", default=[3.5, 5.0],
                    help="additional cutoff values for n_in_shell")
    ap.add_argument("--far", type=float, default=8.0,
                    help="distance threshold for n_far")
    ap.add_argument("--center", type=float, nargs=3, default=None,
                    metavar=("X", "Y", "Z"),
                    help="sphere center coordinates X Y Z")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    ns = SOLVENT_ATOMS[a.solvent]
    expected = N_APIXABAN + ns * a.nsolv
    frames = read_xyz_frames(a.xyz)
    print(f"{a.xyz}: {len(frames)} frames, expected atoms={expected}")

    rows = []
    for k, (_, atoms) in enumerate(frames):
        if len(atoms) != expected:
            rows.append({"frame": k, "n_atoms_ok": 0})
            continue
        apix = atoms[:N_APIXABAN]
        cc = centroid(atoms)
        ac = centroid(apix)
        row = {"frame": k, "n_atoms_ok": 1,
               "A_centroid_to_cluster_centroid": round(math.dist(ac, cc), 3)}
        if a.center is not None:
            row["A_centroid_to_sphere_center"] = round(math.dist(ac, a.center), 3)
        mins = []
        for m in range(a.nsolv):
            s = atoms[N_APIXABAN + m * ns : N_APIXABAN + (m + 1) * ns]
            d = min(dist(p, q) for p in apix for q in s)
            mins.append(d)
            row[f"mol{m + 1}_min"] = round(d, 3)
        row["min_A_S"] = round(min(mins), 3)
        row["max_mol_min_dist"] = round(max(mins), 3)
        row["n_in_shell"] = sum(d <= a.cutoff for d in mins)
        for c in a.extra_cutoffs:
            row[f"n_in_shell_{c:g}"] = sum(d <= c for d in mins)
        row["n_far"] = sum(d > a.far for d in mins)
        rows.append(row)

    good = [r for r in rows if r.get("n_atoms_ok")]
    if good:
        print(f"Valid frames: {len(good)}/{len(rows)}")
        print(f"n_in_shell({a.cutoff:g}): min={min(r["n_in_shell"] for r in good)}, "
              f"mean={sum(r["n_in_shell"] for r in good) / len(good):.2f}, "
              f"max={max(r["n_in_shell"] for r in good)}")
        print(f"n_far(>{a.far:g}): mean={sum(r["n_far"] for r in good) / len(good):.2f}")
        print(f"max_mol_min_dist: max={max(r["max_mol_min_dist"] for r in good):.3f}")
        if a.center is not None:
            print(f"A_centroid_to_sphere_center: max={max(r["A_centroid_to_sphere_center"] for r in good):.3f}")
    if a.out and rows:
        first = ("frame", "n_atoms_ok")
        keys = sorted({k for r in rows for k in r}, key=lambda s: (s not in first, s))
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"Wrote: {a.out}")



if __name__ == "__main__":
    main()




