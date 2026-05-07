from pathlib import Path
import csv

IN_DIR = Path("mafft_alignments")
OUT_DIR = Path("mafft_alignments_cleaned")
OUT_CSV = "alignment_screen_cleaned_ranked.csv"

OUT_DIR.mkdir(exist_ok=True)

# Remove sequences with less than this fraction of the alignment present
MIN_NONGAP_FRACTION = 0.33

# Optional: do not keep alignments that drop below this many taxa
MIN_TAXA_AFTER_CLEANING = 40

def read_fasta(path):
    records = {}
    name = None
    seqs = []

    with open(path) as f:
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

    return records

def write_fasta(records, path):
    with open(path, "w") as out:
        for name, seq in records.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i+80] + "\n")

def taxon_from_header(header):
    return header.split("_")[0]

def score_alignment(records):
    seqs = list(records.values())
    headers = list(records.keys())

    aln_len = len(seqs[0])
    n_seq = len(seqs)
    taxa = {taxon_from_header(h) for h in headers}
    n_taxa = len(taxa)

    total_chars = n_seq * aln_len
    gap_chars = sum(seq.count("-") for seq in seqs)
    gap_fraction = gap_chars / total_chars if total_chars else 0

    non_gap_lengths = [len(seq.replace("-", "")) for seq in seqs]
    mean_nongap = sum(non_gap_lengths) / len(non_gap_lengths)

    fragmentary = sum(1 for x in non_gap_lengths if x < MIN_NONGAP_FRACTION * aln_len)
    fragmentary_fraction = fragmentary / n_seq if n_seq else 0
    seq_taxon_ratio = n_seq / n_taxa if n_taxa else 0

    passes_loose_screen = (
        aln_len >= 60 and
        n_taxa >= 40 and
        seq_taxon_ratio <= 1.20
    )

    passes_preferred_screen = (
        aln_len >= 80 and
        n_taxa >= 50 and
        seq_taxon_ratio <= 1.15 and
        gap_fraction <= 0.75 and
        fragmentary_fraction <= 0.20
    )

    score = (
        (gap_fraction * 2.0) +
        (fragmentary_fraction * 1.5) +
        (1 / aln_len) -
        (n_taxa / 1000)
    )

    return {
        "sequences": n_seq,
        "taxa": n_taxa,
        "seq_taxon_ratio": round(seq_taxon_ratio, 3),
        "alignment_length": aln_len,
        "gap_fraction": round(gap_fraction, 3),
        "mean_nongap_length": round(mean_nongap, 1),
        "fragmentary_sequences": fragmentary,
        "fragmentary_fraction": round(fragmentary_fraction, 3),
        "passes_loose_screen": passes_loose_screen,
        "passes_preferred_screen": passes_preferred_screen,
        "score": round(score, 5)
    }

rows = []

for fp in sorted(IN_DIR.glob("*.aln.fa")):
    records = read_fasta(fp)
    if not records:
        continue

    seqs = list(records.values())
    aln_len = len(seqs[0])

    cleaned = {}
    removed = []

    for name, seq in records.items():
        nongap_len = len(seq.replace("-", ""))
        nongap_fraction = nongap_len / aln_len if aln_len else 0

        if nongap_fraction >= MIN_NONGAP_FRACTION:
            cleaned[name] = seq
        else:
            removed.append(name)

    if len(cleaned) == 0:
        continue

    n_taxa_cleaned = len({taxon_from_header(h) for h in cleaned})

    if n_taxa_cleaned < MIN_TAXA_AFTER_CLEANING:
        continue

    out_fp = OUT_DIR / fp.name.replace(".aln.fa", ".cleaned.aln.fa")
    write_fasta(cleaned, out_fp)

    stats = score_alignment(cleaned)

    row = {
        "orthogroup": fp.stem.replace(".aln", ""),
        "cleaned_alignment_file": out_fp.name,
        "removed_sequences": len(removed),
        **stats
    }

    rows.append(row)

if not rows:
    raise SystemExit("No cleaned alignments passed the minimum taxa threshold.")

rows = sorted(
    rows,
    key=lambda r: (
        not r["passes_preferred_screen"],
        not r["passes_loose_screen"],
        r["score"],
        r["gap_fraction"],
        -r["alignment_length"],
        -r["taxa"]
    )
)

with open(OUT_CSV, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote cleaned alignments to: {OUT_DIR}")
print(f"Wrote ranked summary to: {OUT_CSV}")
print("Top 10 cleaned candidates:")

for r in rows[:10]:
    print(
        r["orthogroup"],
        "score=", r["score"],
        "len=", r["alignment_length"],
        "gap=", r["gap_fraction"],
        "taxa=", r["taxa"],
        "removed=", r["removed_sequences"]
    )
