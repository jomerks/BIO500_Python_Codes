#!/usr/bin/env python3

from pathlib import Path
import argparse

def read_fasta(path):
    records = {}
    name = None
    seq = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    records[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            records[name] = "".join(seq)

    return records

def wrap(seq, width=80):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="input concatenated FASTA")
parser.add_argument("-o", default="concatenated.filtered.fasta")
parser.add_argument("--min_nongap", type=float, default=0.5)
args = parser.parse_args()

records = read_fasta(args.i)

kept = {}
removed = []

for taxon, seq in records.items():
    nongap = sum(1 for c in seq if c not in "-?XNxn")
    prop = nongap / len(seq)

    if prop >= args.min_nongap:
        kept[taxon] = seq
    else:
        removed.append((taxon, prop))

with open(args.o, "w") as out:
    for taxon, seq in kept.items():
        out.write(f">{taxon}\n{wrap(seq)}\n")

print(f"Kept: {len(kept)} taxa")
print(f"Removed: {len(removed)} taxa")
print(f"Output: {args.o}")

print("\nRemoved taxa:")
for taxon, prop in removed:
    print(f"{taxon}\t{prop:.3f}")
