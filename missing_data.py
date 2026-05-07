#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

fasta = Path("concatenated_dedup.fasta")
out_csv = Path("missing_data_per_taxon.csv")

records = {}
name = None
seqs = []

with fasta.open() as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seqs)
            name = line[1:].split()[0]
            seqs = []
        else:
            seqs.append(line)

if name is not None:
    records[name] = "".join(seqs)

rows = []
for taxon, seq in records.items():
    total = len(seq)
    missing = sum(1 for c in seq if c in "-?XxNn")
    present = total - missing
    rows.append({
        "Taxon": taxon,
        "Alignment_length": total,
        "Present_sites": present,
        "Missing_or_gap_sites": missing,
        "Percent_missing": round(100 * missing / total, 2),
        "Percent_present": round(100 * present / total, 2)
    })

df = pd.DataFrame(rows).sort_values("Percent_missing", ascending=False)
df.to_csv(out_csv, index=False)

print(df)
print(f"\nWrote: {out_csv}")