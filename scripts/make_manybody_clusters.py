#!/usr/bin/env python
"""
make_manybody_clusters.py
=========================
Many-body (non-additivity) test inputs for the Apixaban crystal-fragment study (plan B).

Question: can the binding of the central molecule A be written as a sum of pair terms
dE_int(A,k) (= the D1-D9 dimers with their multiplicities)?

Everything is XTB2 single point, crystal geometry, rigid molecules (same protocol as the D1-D9 pulls).

Test 1 (FullShell/, cumulative shells)
    stage n = A + all neighbours of the n strongest dimer types
    E_bind(A|shell) = E(A+shell) - E(shell) - E(A)         <- binding of A to a fixed shell
    E_pair_sum      = sum_k [E(A+k) - E(A) - E(k)]         <- pairs in Pairs/
    non-additivity  = E_bind - E_pair_sum                  (3-body and higher terms that involve A)
    The shell-only energy is subtracted, so neighbour-neighbour interactions cancel exactly.

Test 2 (HalfShell/, simultaneous detachment, "surface-site" model)
    A is pulled rigidly along one common direction v; only the neighbours on the supporting side
    (projection of the A->k centroid vector on v below --tau) are kept and stay fixed.
    cluster(d):  E_bind(d) = E(A(d)+shell) - E(shell) - E(A)
    pair sum(d): sum_k [E(A(d)+k) - E(A) - E(k)]   (A displaced by the SAME vector in every pair)
    The difference shows where the pair picture fails during detachment.

Usage (in the CrystalFragment folder):
    python ..\\..\\scripts\\make_manybody_clusters.py Apixaban_FormN1_1060528.cif --out ManyBody

Requires: pip install ase scipy numpy
No ORCA job is started; run_manybody.ps1 (written to --out) runs them.
"""
import argparse
import json
import os
from collections import Counter

import numpy as np
from ase.io import read
from ase.neighborlist import natural_cutoffs, neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist

# XTB2 crystal-geometry dE_int (kcal/mol) of the nine dimer types. Used ONLY for the default
# stage order and for the weights of the automatic pull direction; no energy is taken from here.
DE_TYPE = {"D1": -14.70, "D2": -17.15, "D3": -16.75, "D4": -5.34, "D5": -7.90,
           "D6": -3.67, "D7": -6.03, "D8": -1.73, "D9": -1.62}
DEFAULT_ORDER = "D2,D3,D1,D5,D7,D4,D6,D8,D9"


def find_molecules(sc, natoms=None):
    i, j = neighbor_list("ij", sc, natural_cutoffs(sc, mult=1.15))
    n = len(sc)
    g = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    ncomp, labels = connected_components(g, directed=False)
    sizes = Counter(labels)
    if natoms is None:
        natoms = max(sizes.values())
    mols = [np.where(labels == c)[0] for c in range(ncomp) if sizes[c] == natoms]
    return mols, natoms


def write_orca(path, title, symbols, pos, nprocs, maxcore):
    with open(path, "w", newline="\n") as f:
        f.write(f"# {title}\n\n! XTB2\n\n%pal\n  nprocs {nprocs}\nend\n\n")
        f.write(f"%maxcore {maxcore}\n\n* xyz 0 1\n")
        for s, p in zip(symbols, pos):
            f.write(f"{s:<2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}\n")
        f.write("*\n")


def neighbour_table(pos, sym, mols, cent, ic, cutoff):
    """Same ranking / signature logic as make_crystal_fragments.py, so the D1-D9 labels agree."""
    A = mols[ic]
    hA = np.array([k for k in A if sym[k] != "H"])
    rows = []
    for k, m in enumerate(mols):
        if k == ic:
            continue
        d = cdist(pos[A], pos[m])
        dmin = float(d.min())
        if dmin > cutoff:
            continue
        hm = np.array([q for q in m if sym[q] != "H"])
        dh = cdist(pos[hA], pos[hm])
        n_contact = int((dh < 3.8).sum())
        n_polar = 0
        for p, q in zip(*np.where(dh < 3.4)):
            if sym[hA[p]] in ("N", "O") and sym[hm[q]] in ("N", "O"):
                n_polar += 1
        cd = float(np.linalg.norm(cent[k] - cent[ic]))
        rows.append(dict(k=k, dmin=dmin, cd=cd, n_contact=n_contact, n_polar=n_polar))
    sig_of = lambda r: (round(r["cd"], 1), round(r["dmin"], 1), r["n_contact"], r["n_polar"])
    mult = Counter(sig_of(r) for r in rows)
    label, uniq = {}, []
    for r in sorted(rows, key=lambda r: (-r["n_polar"], -r["n_contact"], r["dmin"])):
        sig = sig_of(r)
        if sig in label:
            continue
        label[sig] = f"D{len(uniq) + 1}"
        r["mult"] = mult[sig]
        uniq.append(r)
    for r in rows:
        r["type"] = label[sig_of(r)]
    return rows, uniq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--out", default="ManyBody")
    ap.add_argument("--cutoff", type=float, default=4.0)
    ap.add_argument("--natoms", type=int, default=None)
    ap.add_argument("--supercell", type=int, default=5)
    ap.add_argument("--nprocs", type=int, default=4)
    ap.add_argument("--maxcore", type=int, default=2000)
    ap.add_argument("--order", default=DEFAULT_ORDER,
                    help="dimer types in the order they are added to the cumulative shell")
    ap.add_argument("--direction", default="best",
                    help="pull direction of A for the half shell: best | auto | away:Dn | x,y,z")
    ap.add_argument("--clash-min", type=float, default=2.0,
                    help="best: smallest allowed A...kept-neighbour distance (A) during the pull")
    ap.add_argument("--tau", type=float, default=0.1,
                    help="neighbour k is kept in the half shell if (A->k unit vector).v < tau")
    ap.add_argument("--pull", default="0.5,1.0,1.5,2.0,3.0,5.0")
    a = ap.parse_args()

    atoms = read(a.cif)
    sc = atoms.repeat((a.supercell,) * 3)
    sc.pbc = False
    mols, natoms = find_molecules(sc, a.natoms)
    pos = sc.get_positions()
    sym = sc.get_chemical_symbols()
    centre = sc.get_cell().sum(axis=0) / 2
    cent = np.array([pos[m].mean(axis=0) for m in mols])
    ic = int(np.argmin(np.linalg.norm(cent - centre, axis=1)))
    A = mols[ic]
    rows, uniq = neighbour_table(pos, sym, mols, cent, ic, a.cutoff)
    print(f"molecule size = {natoms} atoms, complete molecules in supercell: {len(mols)}")
    print(f"neighbour molecules within cutoff: {len(rows)}; dimer types: {len(uniq)}")
    print("type  mult  centroid_A  min_dist_A  heavy_contacts  NO_contacts")
    for r in uniq:
        print(f"{r['type']:<5s} {r['mult']:>4d}  {r['cd']:>10.2f}  {r['dmin']:>10.2f}  {r['n_contact']:>14d}  {r['n_polar']:>11d}")

    order = [t.strip() for t in a.order.split(",") if t.strip()]
    unknown = set(order) ^ {r["type"] for r in uniq}
    if unknown:
        raise SystemExit(f"--order must list exactly the dimer types {sorted({r['type'] for r in uniq})}; mismatch: {sorted(unknown)}")
    rank = {t: i for i, t in enumerate(order)}
    rows.sort(key=lambda r: (rank[r["type"]], r["k"]))
    for j, r in enumerate(rows, 1):
        r["id"] = j

    symA = [sym[i] for i in A]
    cA = pos[A].mean(axis=0)

    def geom(shift, ks):
        """A (optionally displaced) followed by the molecules in ks."""
        s = list(symA)
        p = [pos[A] + shift]
        for k in ks:
            s += [sym[i] for i in mols[k]]
            p.append(pos[mols[k]])
        return s, np.vstack(p)

    manifest = []
    nonempty = lambda ks: [r["k"] for r in rows if r["k"] in ks]

    def emit(rel, title, s, p, **meta):
        path = os.path.join(a.out, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        write_orca(path, title, s, p, a.nprocs, a.maxcore)
        manifest.append(dict(file=rel.replace("\\", "/"), natoms=len(s), **meta))

    zero = np.zeros(3)

    # ---------------- Test 1: pairs, monomers, cumulative shells ----------------
    s, p = geom(zero, [])
    emit("FullShell/monoA_XTB2_SP.inp", "monoA (central molecule), XTB2 single point", s, p, kind="monoA")
    for r in rows:
        tag = f"{r['id']:02d}_{r['type']}"
        s, p = geom(zero, [r["k"]])
        emit(f"Pairs/Pair_{tag}_XTB2_SP.inp", f"pair A+{r['type']} (neighbour #{r['id']}), crystal geometry, XTB2 SP",
             s, p, kind="pair", neighbour=r["id"], type=r["type"])
        s = [sym[i] for i in mols[r["k"]]]
        emit(f"Pairs/Mono_{tag}_XTB2_SP.inp", f"neighbour #{r['id']} ({r['type']}) alone, XTB2 SP",
             s, pos[mols[r["k"]]], kind="mono_neighbour", neighbour=r["id"], type=r["type"])

    for n in range(2, len(order) + 1):
        types = order[:n]
        ks = [r["k"] for r in rows if r["type"] in types]
        ids = [r["id"] for r in rows if r["type"] in types]
        name = "+".join(types)
        s, p = geom(zero, ks)
        emit(f"FullShell/Stage{n:02d}_{name}_cluster_XTB2_SP.inp",
             f"stage {n}: A + shell of {name} ({len(ks)} neighbours), XTB2 SP", s, p,
             kind="stage_cluster", stage=n, types=types, neighbours=ids)
        s = []
        pp = []
        for k in ks:
            s += [sym[i] for i in mols[k]]
            pp.append(pos[mols[k]])
        emit(f"FullShell/Stage{n:02d}_{name}_shell_XTB2_SP.inp",
             f"stage {n}: shell only of {name} ({len(ks)} neighbours), XTB2 SP", s, np.vstack(pp), kind="stage_shell",
             stage=n, types=types, neighbours=ids)

    # ---------------- Test 2: half shell + common pull direction ----------------
    U = {r["k"]: (cent[r["k"]] - cA) / np.linalg.norm(cent[r["k"]] - cA) for r in rows}
    pulls = [float(x) for x in a.pull.split(",") if x]
    if a.direction == "best":
        # scan directions on a sphere: keep the neighbours on the supporting side (projection < tau),
        # require that no kept neighbour comes closer than --clash-min to A at any pull distance,
        # and maximise the summed |dE_int| of the kept neighbours.
        n_dir = 4000
        i_ = np.arange(n_dir) + 0.5
        ph = np.arccos(1 - 2 * i_ / n_dir)
        th = np.pi * (1 + 5 ** 0.5) * i_
        dirs = np.c_[np.cos(th) * np.sin(ph), np.sin(th) * np.sin(ph), np.cos(ph)]
        W = {r["k"]: (pos[A][:, None, :] - pos[mols[r["k"]]][None, :, :]).reshape(-1, 3) for r in rows}
        W2 = {k: (w ** 2).sum(axis=1) for k, w in W.items()}
        best = None
        for vv in dirs:
            kept = [r for r in rows if float(U[r["k"]] @ vv) < a.tau]
            if not kept:
                continue
            md = min(float(np.sqrt((W2[r["k"]] + 2 * d * (W[r["k"]] @ vv) + d * d).min()))
                     for r in kept for d in pulls)
            if md < a.clash_min:
                continue
            score = sum(abs(DE_TYPE[r["type"]]) for r in kept)
            if best is None or (score, md) > best[:2]:
                best = (score, md, vv)
        if best is None:
            raise SystemExit("no clash-free direction found; lower --clash-min or raise --tau")
        v = best[2]
        why = (f"best: clash-free (>= {a.clash_min:.1f} A) direction keeping the largest summed |dE_int| "
               f"({best[0]:.1f} kcal/mol, closest approach {best[1]:.2f} A)")
    elif a.direction == "auto":
        v = -sum(abs(DE_TYPE[r["type"]]) * U[r["k"]] for r in rows)
        aniso = np.linalg.norm(v) / sum(abs(DE_TYPE[r["type"]]) for r in rows)
        v = v / np.linalg.norm(v)
        why = f"auto (opposite to the binding-weighted sum of neighbour directions, anisotropy {aniso:.2f})"
    elif a.direction.startswith("away:"):
        t = a.direction.split(":", 1)[1]
        k0 = next(r["k"] for r in rows if r["type"] == t)
        v = -U[k0]
        why = f"away from the first {t} neighbour"
    else:
        v = np.array([float(x) for x in a.direction.split(",")])
        v = v / np.linalg.norm(v)
        why = "user supplied"
    keep = [r for r in rows if float(U[r["k"]] @ v) < a.tau]
    drop = [r for r in rows if float(U[r["k"]] @ v) >= a.tau]
    print(f"\nhalf-shell pull direction v = ({v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f}) [{why}]")
    print("neighbour  type  proj_on_v   kept   min_dist_to_A(d=0 / 2 / 5 A)")
    for r in rows:
        proj = float(U[r["k"]] @ v)
        dd = [cdist(pos[A] + x * v, pos[mols[r["k"]]]).min() for x in (0.0, 2.0, 5.0)]
        print(f"#{r['id']:<8d} {r['type']:<5s} {proj:>9.2f}   {'yes' if r in keep else 'no ':<4s}   "
              f"{dd[0]:.2f} / {dd[1]:.2f} / {dd[2]:.2f}")
    if not keep:
        raise SystemExit("half shell is empty; raise --tau or change --direction")
    ks = [r["k"] for r in keep]
    ids = [r["id"] for r in keep]
    s = []
    pp = []
    for k in ks:
        s += [sym[i] for i in mols[k]]
        pp.append(pos[mols[k]])
    emit("HalfShell/HS_shell_XTB2_SP.inp", f"half shell only ({len(ks)} neighbours), XTB2 SP", s, np.vstack(pp),
         kind="hs_shell", neighbours=ids)
    for d in [0.0] + pulls:
        s, p = geom(d * v, ks)
        emit(f"HalfShell/HS_cluster_{d:+.1f}A_XTB2_SP.inp",
             f"half shell: A displaced {d:+.1f} A along v, {len(ks)} fixed neighbours, XTB2 SP", s, p,
             kind="hs_cluster", d=d, neighbours=ids)
    for r in keep:
        for d in pulls:
            s, p = geom(d * v, [r["k"]])
            emit(f"HalfShell/HS_pair_{r['id']:02d}_{r['type']}_{d:+.1f}A_XTB2_SP.inp",
                 f"A displaced {d:+.1f} A along v + neighbour #{r['id']} ({r['type']}), XTB2 SP", s, p,
                 kind="hs_pair", d=d, neighbour=r["id"], type=r["type"])

    meta = dict(cif=os.path.basename(a.cif), natoms_molecule=int(natoms), order=order,
                direction=[float(x) for x in v], direction_rule=why, tau=a.tau, pulls=pulls,
                neighbours=[dict(id=r["id"], type=r["type"], centroid_A=r["cd"], min_dist_A=r["dmin"])
                            for r in rows],
                half_shell_neighbours=ids, jobs=manifest)
    with open(os.path.join(a.out, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1, ensure_ascii=False)

    ps1 = r"""# run_manybody.ps1 -- run all XTB2 single points under this folder (small jobs first).
# Jobs whose .out already ends with 'ORCA TERMINATED NORMALLY' are skipped. Run from this folder.
Get-ChildItem -Recurse -Filter '*_XTB2_SP.inp' | Sort-Object Length | ForEach-Object {
    $name = $_.Name
    $base = [IO.Path]::GetFileNameWithoutExtension($name)
    $out  = Join-Path $_.DirectoryName ($base + '.out')
    if ((Test-Path $out) -and (Select-String -Path $out -Pattern 'ORCA TERMINATED NORMALLY' -Quiet)) { "skip $name"; return }
    "run  $name"
    Push-Location $_.DirectoryName
    cmd /c "orca `"$name`" > `"$base.out`""
    Pop-Location
}
"""
    with open(os.path.join(a.out, "run_manybody.ps1"), "w", encoding="ascii", newline="\r\n") as f:
        f.write(ps1)

    kinds = Counter(m["kind"] for m in manifest)
    print(f"\ninputs written to {a.out}: {len(manifest)} files  {dict(kinds)}")
    print(f"largest job: {max(m['natoms'] for m in manifest)} atoms")
    print("manifest.json and run_manybody.ps1 written")


if __name__ == "__main__":
    main()
