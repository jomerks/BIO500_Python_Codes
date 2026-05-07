from pathlib import Path
import csv

ALN_DIR = Path("mafft_alignments")
OUT_CSV = "alignment_screen.csv"

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
    # Works for CODE_geneID or Species_geneID style.
    # Change this if your headers need a different rule.
    return header.split("_")[0]

rows = []

for fp in sorted(ALN_DIR.glob("*.aln.fa")):
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

    fragmentary = sum(
        1 for x in non_gap_lengths
        if x < 0.33 * aln_len
    )

    seq_taxon_ratio = n_seq / n_taxa if n_taxa else 0

    passes = (
        aln_len >= 100 and
        gap_fraction <= 0.50 and
        seq_taxon_ratio <= 1.15 and
        fragmentary <= 5 and
        n_taxa >= 50
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
        "passes_second_screen": passes
    })

rows = sorted(
    rows,
    key=lambda r: (
        not r["passes_second_screen"],
        r["gap_fraction"],
        -r["alignment_length"],
        -r["taxa"],
        r["seq_taxon_ratio"]
    )
)

with open(OUT_CSV, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print(f"Done. Wrote {OUT_CSV}")
