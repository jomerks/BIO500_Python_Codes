from pathlib import Path
import csv

ALN_DIR = Path("mafft_alignments")
OUT_CSV = "alignment_screen_ranked.csv"

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

def taxon_from_header(header):
    return header.split("_")[0]

patterns = ["*.aln.fa", "*.aln.fasta", "*.fa", "*.faa", "*.fasta"]

files = []
for pattern in patterns:
    files.extend(ALN_DIR.glob(pattern))

files = sorted(set(files))

if not files:
    raise SystemExit(f"No alignment files found in {ALN_DIR}")

rows = []

for fp in files:
    records = read_fasta(fp)

    if not records:
        continue

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

    fragmentary = sum(1 for x in non_gap_lengths if x < 0.33 * aln_len)
    fragmentary_fraction = fragmentary / n_seq if n_seq else 0

    seq_taxon_ratio = n_seq / n_taxa if n_taxa else 0

    # Loose screen: removes obvious bad candidates only
    passes_loose_screen = (
        aln_len >= 60 and
        n_taxa >= 40 and
        seq_taxon_ratio <= 1.20
    )

    # Preferred screen: realistic for transcriptome data
    passes_preferred_screen = (
        aln_len >= 80 and
        n_taxa >= 50 and
        seq_taxon_ratio <= 1.15 and
        gap_fraction <= 0.75 and
        fragmentary_fraction <= 0.50
    )

    # Lower score = better
    score = (
        (gap_fraction * 2.0) +
        (fragmentary_fraction * 1.5) +
        (1 / aln_len) -
        (n_taxa / 1000)
    )

    rows.append({
        "orthogroup": fp.stem.replace(".aln", ""),
        "alignment_file": fp.name,
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
    })

if not rows:
    raise SystemExit("Files were found, but no FASTA records were read.")

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

print(f"Screened {len(rows)} alignments.")
print(f"Wrote {OUT_CSV}.")
print("Top 10 candidates:")

for r in rows[:10]:
    print(
        r["orthogroup"],
        "score=", r["score"],
        "len=", r["alignment_length"],
        "gap=", r["gap_fraction"],
        "taxa=", r["taxa"],
        "frag=", r["fragmentary_sequences"]
    )
