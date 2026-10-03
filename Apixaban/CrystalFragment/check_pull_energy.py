from pathlib import Path
import re

BASE = Path.cwd()

def final_energy(path):
    text = path.read_text(encoding="utf-8", errors="ignore")

    patterns = [
        r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)(?:\s|$)",
        r"FINAL SINGLE POINT ENERGY\s*([-+]?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
    ]

    for pattern in patterns:
        hits = re.findall(pattern, text)
        if hits:
            return float(hits[-1])

    raise ValueError("FINAL SINGLE POINT ENERGY が見つかりません")

for d in range(1, 10):
    dn = f"D{d}"

    print(f"\n=== {dn} ===")

    for filename in [
        f"{dn}_complex_XTB2_SP.out",
        f"{dn}_monoA_XTB2_SP.out",
        f"{dn}_monoB_XTB2_SP.out",
        f"{dn}_pull_+0.5A_XTB2_SP.out",
    ]:
        path = BASE / dn / filename

        if not path.exists():
            print(f"[MISSING] {path}")
            continue

        try:
            energy = final_energy(path)
            print(f"[OK] {filename}")
            print(f"     E = {energy:.12f} Eh")
        except Exception as e:
            print(f"[ERROR] {filename}")
            print(f"        {e}")
