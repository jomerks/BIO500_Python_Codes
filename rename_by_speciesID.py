#!/usr/bin/env python3
from pathlib import Path
import re

SPECIES_ID_FILE = "SpeciesIDs.txt"
INPUT_PATTERN = "*.cleaned.aln.fa"

# Read SpeciesIDs.txt
# Expected format: 10: CATH.faa
id_to_code = {}

with open(SPECIES_ID_FILE) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        m = re.match(r"^(\d+):\s+(.+?)\.faa", line)
        if m:
            num = m.group(1)
            code = m.group(2)
            id_to_code[num] = code

print(f"Loaded {len(id_to_code)} species IDs")

for fp in sorted(Path(".").glob(INPUT_PATTERN)):
    out_fp = fp.with_name(fp.stem + ".renamed.fasta")

    used = {}

    with open(fp) as inp, open(out_fp, "w") as out:
        for line in inp:
            if line.startswith(">"):
                old = line[1:].strip().split()[0]

                # old header format: 10_16842
                species_num = old.split("_")[0]

                if species_num in id_to_code:
                    code = id_to_code[species_num]

                    # handle duplicate sequences from same species in same gene
                    used[code] = used.get(code, 0) + 1
                    if used[code] == 1:
                        new = code
                    else:
                        new = f"{code}_copy{used[code]}"

                    out.write(f">{new}\n")
                else:
                    out.write(f">{old}\n")
            else:
                out.write(line)

    print(f"Wrote {out_fp}")
