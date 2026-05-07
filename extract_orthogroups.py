from pathlib import Path

# load all sequences into memory
seqs = {}

for fasta in Path(".").glob("Species*.fa"):
    species = fasta.stem.replace("Species","")
    with open(fasta) as f:
        header = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                seqs[f"{species}_{header.split('_')[-1]}"] = line

# make output folder
outdir = Path("OG_FASTA")
outdir.mkdir(exist_ok=True)

# parse orthogroups
with open("OrthoFinder/Results_Apr26/Orthogroups/Orthogroups.txt") as f:
    for i, line in enumerate(f):
        parts = line.strip().split()
        if len(parts) < 5:
            continue

        og = f"OG{i:06d}"
        with open(outdir / f"{og}.fa", "w") as out:
            for gene in parts:
                if gene in seqs:
                    out.write(f">{gene}\n{seqs[gene]}\n")
