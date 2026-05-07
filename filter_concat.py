#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict

# -------- SETTINGS --------
INPUT_PATTERN = "*.cleaned.aln.renamed.fasta"   # change if needed
MIN_TAXA_PER_GENE = 30                  # keep genes with at least this many taxa
MIN_GENES_PER_TAXON = 3                 # keep taxa present in at least 3 genes
MAX_SITE_GAP_FRACTION = 0.70            # remove concat columns with >70% gaps
OUT_FASTA = "concatenated.filtered.fasta"
OUT_PARTITIONS = "partitions.filtered.txt"
# --------------------------


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).replace(" ", "")
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        seqs[name] = "".join(chunks).replace(" ", "")

    return seqs


def write_fasta(seqs, outpath, width=80):
    with open(outpath, "w") as out:
        for name, seq in seqs.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")


# 1. Read alignments
files = sorted(Path(".").glob(INPUT_PATTERN))

if not files:
    raise SystemExit(f"No files matched: {INPUT_PATTERN}")

alignments = []

for fp in files:
    seqs = read_fasta(fp)

    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise SystemExit(f"ERROR: sequences in {fp} are not same length")

    aln_len = list(lengths)[0]

    if len(seqs) >= MIN_TAXA_PER_GENE:
        alignments.append((fp.name, seqs, aln_len))
    else:
        print(f"Skipping {fp.name}: only {len(seqs)} taxa")

if not alignments:
    raise SystemExit("No alignments passed MIN_TAXA_PER_GENE filter")

# 2. Count how many genes each taxon appears in
taxon_gene_count = defaultdict(int)

for name, seqs, aln_len in alignments:
    for taxon, seq in seqs.items():
        # count only taxa with some real amino acids
        nongap = sum(1 for c in seq if c not in "-?Xx")
        if nongap > 0:
            taxon_gene_count[taxon] += 1

kept_taxa = sorted(
    t for t, n in taxon_gene_count.items()
    if n >= MIN_GENES_PER_TAXON
)

print(f"Kept genes: {len(alignments)}")
print(f"Kept taxa: {len(kept_taxa)}")

if not kept_taxa:
    raise SystemExit("No taxa passed MIN_GENES_PER_TAXON filter")

# 3. Concatenate
concat = {t: [] for t in kept_taxa}
partitions_raw = []
start = 1

for gene_name, seqs, aln_len in alignments:
    end = start + aln_len - 1
    partitions_raw.append((gene_name, start, end))

    for taxon in kept_taxa:
        concat[taxon].append(seqs.get(taxon, "-" * aln_len))

    start = end + 1

concat = {t: "".join(parts) for t, parts in concat.items()}

# 4. Remove gappy columns from final concatenated matrix
seq_len = len(next(iter(concat.values())))
keep_cols = []

for i in range(seq_len):
    column = [seq[i] for seq in concat.values()]
    gaps = sum(1 for c in column if c in "-?")
    gap_fraction = gaps / len(column)

    if gap_fraction <= MAX_SITE_GAP_FRACTION:
        keep_cols.append(i)

filtered = {
    taxon: "".join(seq[i] for i in keep_cols)
    for taxon, seq in concat.items()
}

# 5. Remove taxa that became all gaps
filtered = {
    taxon: seq for taxon, seq in filtered.items()
    if any(c not in "-?Xx" for c in seq)
}

write_fasta(filtered, OUT_FASTA)

# partitions are approximate after site filtering, so do not use them for now
with open(OUT_PARTITIONS, "w") as out:
    out.write("# Partitions not reliable after site filtering\n")
    out.write("# Use unpartitioned IQ-TREE unless you implement partition remapping\n")

print(f"Wrote: {OUT_FASTA}")
print(f"Wrote: {OUT_PARTITIONS}")
print(f"Final alignment length: {len(next(iter(filtered.values())))}")
