import re
from pathlib import Path
import csv

fasta_dir = Path(".")  # folder with orthogroup FASTAs
output_file = "orthogroup_screen.csv"

required_taxa = [
    "Spirogyra",
    "Mougeotia",
    "Chara",
    "Coleochaete"
]

def get_species(header):
    header = header.strip().lstrip(">")
    
    species = re.split(r"[_.|]", header)[0]
    
    return species

rows = []

for fp in fasta_dir.glob("*.fa"):
    sequences = 0
    species_set = set()
    
    with open(fp) as f:
        for line in f:
            if line.startswith(">"):
                sequences += 1
                sp = get_species(line)
                species_set.add(sp)
    
    n_species = len(species_set)
    ratio = sequences / n_species if n_species > 0 else 0
    
    taxa_present = all(any(t in sp for sp in species_set) for t in required_taxa)
    
    rows.append({
        "orthogroup": fp.name,
        "sequences": sequences,
        "species": n_species,
        "ratio": round(ratio, 3),
        "has_required_taxa": taxa_present
    })

with open(output_file, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Done. Output written to {output_file}")
