from pathlib import Path
import re
import csv

BASE = Path.cwd()

DISTANCES = ["0.5", "1.0", "1.5", "2.0", "3.0", "5.0"]

MULTIPLICITY = {
    "D1": 2,
    "D2": 1,
    "D3": 1,
    "D4": 1,
    "D5": 2,
    "D6": 2,
    "D7": 2,
    "D8": 2,
    "D9": 1,
}

EH_TO_KCAL = 627.5094740631
EH_TO_KJ = 2625.499638

def read_orca(path):
    raw = path.read_bytes()

    # BOM / UTF-16 detection
    if raw.startswith(b'\xff\xfe'):
        return raw.decode("utf-16-le", errors="replace")
    if raw.startswith(b'\xfe\xff'):
        return raw.decode("utf-16-be", errors="replace")

    # UTF-16 without BOM: many null bytes
    if b'\x00' in raw[:1000]:
        try:
            text = raw.decode("utf-16-le", errors="replace")
            if "FINAL SINGLE POINT ENERGY" in text or "TOTAL ENERGY" in text:
                return text
        except Exception:
            pass

    # UTF-8
    return raw.decode("utf-8", errors="replace")


def final_energy(path):
    text = read_orca(path)

    # Preferred ORCA final energy
    hits = re.findall(
        r"FINAL SINGLE POINT ENERGY\s+"
        r"([-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",
        text
    )

    if hits:
        return float(hits[-1])

    # Fallback: XTB total energy
    hits = re.findall(
        r"::\s*total energy\s+"
        r"([-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*Eh",
        text,
        flags=re.IGNORECASE
    )

    if hits:
        return float(hits[-1])

    raise ValueError("energy not found")


rows = []
errors = []

print("==========================================")
print(" Apixaban D1-D9 rigid-pull XTB2 analysis")
print("==========================================")

for d in range(1, 10):

    dn = f"D{d}"

    complex_file = BASE / dn / f"{dn}_complex_XTB2_SP.out"
    mono_a_file = BASE / dn / f"{dn}_monoA_XTB2_SP.out"
    mono_b_file = BASE / dn / f"{dn}_monoB_XTB2_SP.out"

    try:
        e_complex_crystal = final_energy(complex_file)
        e_mono_a = final_energy(mono_a_file)
        e_mono_b = final_energy(mono_b_file)

    except Exception as e:
        errors.append(f"{dn}: reference energy error: {e}")
        continue

    crystal_de_eh = e_complex_crystal - e_mono_a - e_mono_b
    crystal_de_kcal = crystal_de_eh * EH_TO_KCAL
    crystal_de_kj = crystal_de_eh * EH_TO_KJ

    print(f"\n{dn}")
    print(f"  Crystal ΔEint = {crystal_de_kcal: .6f} kcal/mol")

    for dist in DISTANCES:

        pull_file = BASE / dn / f"{dn}_pull_+{dist}A_XTB2_SP.out"

        try:
            e_pull = final_energy(pull_file)

            de_eh = e_pull - e_mono_a - e_mono_b
            de_kcal = de_eh * EH_TO_KCAL
            de_kj = de_eh * EH_TO_KJ

            rows.append({
                "Dimer": dn,
                "Distance_A": float(dist),
                "E_complex_Eh": e_pull,
                "E_monoA_Eh": e_mono_a,
                "E_monoB_Eh": e_mono_b,
                "DeltaEint_Eh": de_eh,
                "DeltaEint_kcal_mol": de_kcal,
                "DeltaEint_kJ_mol": de_kj,
                "Multiplicity": MULTIPLICITY[dn],
                "Crystal_DeltaEint_kcal_mol": crystal_de_kcal,
                "Crystal_DeltaEint_kJ_mol": crystal_de_kj,
            })

        except Exception as e:
            errors.append(
                f"{dn} +{dist}A: {e}"
            )

# -------------------------------------------------
# CSV 1: all 54 pull points
# -------------------------------------------------

csv_path = BASE / "pull_curves_D1-D9_XTB2.csv"

fields = [
    "Dimer",
    "Distance_A",
    "E_complex_Eh",
    "E_monoA_Eh",
    "E_monoB_Eh",
    "DeltaEint_Eh",
    "DeltaEint_kcal_mol",
    "DeltaEint_kJ_mol",
    "Multiplicity",
    "Crystal_DeltaEint_kcal_mol",
    "Crystal_DeltaEint_kJ_mol",
]

with csv_path.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

# -------------------------------------------------
# CSV 2: one row per dimer
# -------------------------------------------------

summary_path = BASE / "pull_curves_D1-D9_summary.csv"

summary_fields = [
    "Dimer",
    "Multiplicity",
    "Crystal_DeltaEint_kcal_mol",
    "DeltaE_0.5_kcal_mol",
    "DeltaE_1.0_kcal_mol",
    "DeltaE_1.5_kcal_mol",
    "DeltaE_2.0_kcal_mol",
    "DeltaE_3.0_kcal_mol",
    "DeltaE_5.0_kcal_mol",
]

with summary_path.open("w", newline="", encoding="utf-8-sig") as f:

    writer = csv.DictWriter(f, fieldnames=summary_fields)
    writer.writeheader()

    for d in range(1, 10):

        dn = f"D{d}"
        rr = [
            x for x in rows
            if x["Dimer"] == dn
        ]

        if len(rr) != 6:
            continue

        values = {
            str(x["Distance_A"]): x["DeltaEint_kcal_mol"]
            for x in rr
        }

        crystal = rr[0]["Crystal_DeltaEint_kcal_mol"]

        writer.writerow({
            "Dimer": dn,
            "Multiplicity": MULTIPLICITY[dn],
            "Crystal_DeltaEint_kcal_mol": crystal,
            "DeltaE_0.5_kcal_mol": values["0.5"],
            "DeltaE_1.0_kcal_mol": values["1.0"],
            "DeltaE_1.5_kcal_mol": values["1.5"],
            "DeltaE_2.0_kcal_mol": values["2.0"],
            "DeltaE_3.0_kcal_mol": values["3.0"],
            "DeltaE_5.0_kcal_mol": values["5.0"],
        })

# -------------------------------------------------
# Pairwise crystal interaction
# -------------------------------------------------

total_pair = 0.0

for d in range(1, 10):

    dn = f"D{d}"

    rr = [
        x for x in rows
        if x["Dimer"] == dn
    ]

    if not rr:
        continue

    crystal_de = rr[0]["Crystal_DeltaEint_kcal_mol"]

    total_pair += (
        0.5
        * MULTIPLICITY[dn]
        * crystal_de
    )

# -------------------------------------------------
# Console summary
# -------------------------------------------------

print("\n==========================================")
print(" SUMMARY")
print("==========================================")

print(f"Pull points parsed : {len(rows)} / 54")

print(
    f"Pairwise crystal interaction : "
    f"{total_pair:.6f} kcal/mol"
)

print(
    f"Pairwise crystal interaction : "
    f"{total_pair * 4.184:.3f} kJ/mol"
)

print(f"\nCSV:")
print(csv_path)

print("\nSummary CSV:")
print(summary_path)

if errors:

    print("\n==========================================")
    print(" ERRORS")
    print("==========================================")

    for error in errors:
        print(error)

else:

    print("\nAll 54 pull points parsed successfully.")

# -------------------------------------------------
# Plot
# -------------------------------------------------

try:

    import matplotlib.pyplot as plt

    plt.figure(figsize=(9, 6))

    for d in range(1, 10):

        dn = f"D{d}"

        rr = sorted(
            [
                x for x in rows
                if x["Dimer"] == dn
            ],
            key=lambda x: x["Distance_A"]
        )

        if rr:

            plt.plot(
                [x["Distance_A"] for x in rr],
                [x["DeltaEint_kcal_mol"] for x in rr],
                marker="o",
                label=dn
            )

    plt.axhline(
        0,
        linewidth=0.8
    )

    plt.xlabel(
        "Rigid pull distance (Å)"
    )

    plt.ylabel(
        "ΔEint (kcal/mol)"
    )

    plt.title(
        "Apixaban crystal dimer rigid-pull curves — ORCA XTB2"
    )

    plt.legend(ncol=3)

    plt.grid(
        True,
        alpha=0.25
    )

    plt.tight_layout()

    plot_path = (
        BASE /
        "pull_curves_D1-D9_XTB2.png"
    )

    plt.savefig(
        plot_path,
        dpi=180
    )

    plt.close()

    print("\nPlot:")
    print(plot_path)

except Exception as e:

    print("\nPlot skipped:")
    print(e)
