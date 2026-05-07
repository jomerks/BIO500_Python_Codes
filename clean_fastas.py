from pathlib import Path
import re

allowed = set("ACDEFGHIKLMNPQRSTVWY*X")

def clean_seq(seq):
    seq = re.sub(r"\s+", "", seq).upper()
    seq = re.sub(r"[^A-Z*]", "", seq)
    return "".join(c if c in allowed else "X" for c in seq)

for path in sorted(Path(".").glob("Species*.fa")):
    out = Path("CLEAN_FASTA") / path.name
    with open(path) as fh, open(out, "w") as oh:
        header = None
        seq_parts = []

        def write_record():
            if header is None:
                return
            seq = clean_seq("".join(seq_parts))
            if len(seq) > 0:
                oh.write(f">{header}\n")
                for i in range(0, len(seq), 80):
                    oh.write(seq[i:i+80] + "\n")

        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                write_record()
                header = line[1:].strip().replace(" ", "_")
                seq_parts = []
            else:
                seq_parts.append(line)

        write_record()

print("Cleaned FASTA files written to CLEAN_FASTA/")
