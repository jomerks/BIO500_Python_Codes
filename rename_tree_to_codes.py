import re

treefile = "gene_alignments/gene_alignments_renamed/gene_alignments_renamed_dedup/concat_cleaned.treefile"
outfile = "gene_alignments/gene_alignments_renamed/gene_alignments_renamed_dedup/concat_cleaned_named.tree"

mapping = {}

with open("SpeciesIDs.txt") as f:
    for line in f:
        if ":" not in line:
            continue
        sid, code = line.strip().split(":", 1)
        code = code.strip().replace(".faa", "").replace("(1)", "_1")
        mapping[sid.strip()] = code

with open(treefile) as f:
    tree = f.read()

for sid, code in sorted(mapping.items(), key=lambda x: -len(x[0])):
    tree = re.sub(rf'\b{sid}_([A-Za-z0-9]+)\b', rf'{code}_\1', tree)

with open(outfile, "w") as out:
    out.write(tree)

print(f"Wrote {outfile}")
