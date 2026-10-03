#!/usr/bin/env python
"""
make_crystal_fragments.py
=========================
CIF (Apixaban crystal) -> crystal-fragment dimer ORCA inputs (plan B).

For the central molecule in a NxNxN supercell (default 5), find neighbouring molecules,
rank them by contact strength, and for the top-N dimers write ORCA inputs:

  1. interaction energy   dE_int = E(AB) - E(A) - E(B)   (A, B in dimer geometry)
       XTB2 and r2SCAN-3c single points
  2. rigid pull-apart     B is displaced along the A->B centroid vector by
       the given distances (XTB2 single points) -> crude detachment curve

Usage:
  python make_crystal_fragments.py Apixaban_FormN1.cif --out CrystalFragment --top 3

Requires: pip install ase scipy numpy

Notes:
  * Only geometry handling is done here; all energies come from ORCA.
  * X-ray H positions are unreliable. Consider relaxing H (heavy atoms fixed)
    before trusting absolute energies.
  * The molecule size is auto-detected (largest connected component);
    give --natoms 59 to force it (Apixaban C25H25N5O4 = 59 atoms).
"""
import argparse
import os
from collections import Counter

import numpy as np
from ase.io import read
from ase.neighborlist import natural_cutoffs, neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist


def find_molecules(sc, natoms=None):
    i, j = neighbor_list("ij", sc, natural_cutoffs(sc, mult=1.15))
    n = len(sc)
    g = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    ncomp, labels = connected_components(g, directed=False)
    sizes = Counter(labels)
    if natoms is None:
        # 分子の大きさ = 最大の連結成分 (セル境界で切れた断片が多いので最頻値は使えない)
        natoms = max(sizes.values())
    mols = [np.where(labels == c)[0] for c in range(ncomp) if sizes[c] == natoms]
    return mols, natoms


def write_orca(path, title, header, symbols, pos, nprocs, maxcore=2000):
    with open(path, "w", newline="\n") as f:
        f.write(f"# {title}\n\n{header}\n\n%pal\n  nprocs {nprocs}\nend\n\n")
        f.write(f"%maxcore {maxcore}\n\n* xyz 0 1\n")
        for s, p in zip(symbols, pos):
            f.write(f"{s:<2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}\n")
        f.write("*\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--out", default="CrystalFragment")
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--cutoff", type=float, default=4.0,
                    help="min interatomic distance (A) to count as neighbour")
    ap.add_argument("--natoms", type=int, default=None)
    ap.add_argument("--supercell", type=int, default=5,
                    help="supercell repeat N x N x N (neighbours cut at the edge are dropped, so keep it large)")
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--skip", type=int, default=0,
                    help="do not write files for the first N dimers (numbering is kept)")
    ap.add_argument("--levels", default="XTB2,r2SCAN3c", help="comma list: XTB2, r2SCAN3c")
    ap.add_argument("--pull", default="0.5,1.0,1.5,2.0,3.0,5.0",
                    help="rigid displacements (A) of molecule B")
    a = ap.parse_args()

    atoms = read(a.cif)
    sc = atoms.repeat((a.supercell,) * 3)
    sc.pbc = False  # do not connect molecules through the supercell boundary
    mols, natoms = find_molecules(sc, a.natoms)
    print(f"molecule size = {natoms} atoms, complete molecules in supercell: {len(mols)}")

    pos = sc.get_positions()
    sym = sc.get_chemical_symbols()
    centre = sc.get_cell().sum(axis=0) / 2
    cent = np.array([pos[m].mean(axis=0) for m in mols])
    ic = int(np.argmin(np.linalg.norm(cent - centre, axis=1)))
    A = mols[ic]

    heavy = lambda idx: np.array([k for k in idx if sym[k] != "H"])
    hA = heavy(A)
    rows = []
    for k, m in enumerate(mols):
        if k == ic:
            continue
        d = cdist(pos[A], pos[m])
        dmin = d.min()
        if dmin > a.cutoff:
            continue
        hm = heavy(m)
        dh = cdist(pos[hA], pos[hm])
        n_contact = int((dh < 3.8).sum())               # heavy-heavy contacts
        n_polar = 0                                      # N/O...N/O (H-bond-like)
        for p, q in zip(*np.where(dh < 3.4)):
            if sym[hA[p]] in ("N", "O") and sym[hm[q]] in ("N", "O"):
                n_polar += 1
        cd = float(np.linalg.norm(cent[k] - cent[ic]))
        rows.append(dict(k=k, dmin=dmin, cd=cd, n_contact=n_contact, n_polar=n_polar))

    # de-duplicate symmetry-equivalent neighbours (signature)
    sig_of = lambda r: (round(r["cd"], 1), round(r["dmin"], 1), r["n_contact"], r["n_polar"])
    mult = Counter(sig_of(r) for r in rows)   # number of neighbour molecules of each dimer type
    seen, uniq = set(), []
    for r in sorted(rows, key=lambda r: (-r["n_polar"], -r["n_contact"], r["dmin"])):
        sig = sig_of(r)
        if sig in seen:
            continue
        seen.add(sig)
        r["mult"] = mult[sig]
        uniq.append(r)

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "summary.csv"), "w") as f:
        f.write("dimer,centroid_dist_A,min_dist_A,heavy_contacts_lt3.8A,NO_contacts_lt3.4A,multiplicity\n")
        for n, r in enumerate(uniq, 1):
            f.write(f"D{n},{r['cd']:.2f},{r['dmin']:.2f},{r['n_contact']},{r['n_polar']},{r['mult']}\n")
    print(f"unique neighbour dimers: {len(uniq)}, neighbour molecules within cutoff: {len(rows)} (summary.csv written)")

    pulls = [float(x) for x in a.pull.split(",") if x]
    xtb = "! XTB2"
    dft = "! r2SCAN-3c TightSCF"
    bat = ["@echo off", "rem run all ORCA jobs (orca must be on PATH)"]
    for n, r in enumerate(uniq[: a.top], 1):
        if n <= a.skip:
            continue
        B = mols[r["k"]]
        d = os.path.join(a.out, f"D{n}")
        os.makedirs(d, exist_ok=True)
        idx_ab = np.concatenate([A, B])
        sets = {"complex": idx_ab, "monoA": A, "monoB": B}
        tag = f"D{n}"
        for lvl, hdr in (("XTB2", xtb), ("r2SCAN3c", dft)):
            if lvl not in a.levels.split(","):
                continue
            for name, idx in sets.items():
                fn = f"{tag}_{name}_{lvl}_SP.inp"
                write_orca(os.path.join(d, fn),
                           f"{tag} {name}: {lvl} single point (crystal geometry, "
                           f"centroid {r['cd']:.2f} A, min contact {r['dmin']:.2f} A)",
                           hdr, [sym[k] for k in idx], pos[idx], a.nprocs)
                bat.append(f'orca "{tag}\\{fn}" > "{tag}\\{fn[:-4]}.out"')
        # rigid pull-apart of B along A->B centroid vector (XTB2)
        v = cent[r["k"]] - pos[A].mean(axis=0)
        v /= np.linalg.norm(v)
        for dz in pulls:
            pB = pos[B] + v * dz
            allp = np.vstack([pos[A], pB])
            allsym = [sym[k] for k in A] + [sym[k] for k in B]
            fn = f"{tag}_pull_{dz:+.1f}A_XTB2_SP.inp"
            write_orca(os.path.join(d, fn),
                       f"{tag} pull-apart: molecule B displaced {dz:+.1f} A along "
                       f"centroid vector (rigid), XTB2 single point",
                       xtb, allsym, allp, a.nprocs)
            bat.append(f'orca "{tag}\\{fn}" > "{tag}\\{fn[:-4]}.out"')
    with open(os.path.join(a.out, "run_all.bat"), "w", newline="\r\n") as f:
        f.write("\n".join(bat) + "\n")
    print("done ->", a.out)


if __name__ == "__main__":
    main()
