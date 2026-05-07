import pandas as pd
import re
from pathlib import Path

treefile = Path("gene_alignments/gene_alignments_renamed/gene_alignments_renamed_dedup/concat_cleaned_named.tree")
csvfile = Path("taxa_classification_corrected.csv")
outfile = Path("gene_alignments/gene_alignments_renamed/gene_alignments_renamed_dedup/concat_taxonomy_colors_exact.txt")

print("Tree exists:", treefile.exists(), treefile.resolve())
print("CSV exists:", csvfile.exists(), csvfile.resolve())

tree = treefile.read_text().strip()
print("Tree length:", len(tree))
print("First 200 chars:")
print(tree[:200])

tips = re.findall(r'([A-Z0-9]{4}(?:_[A-Za-z0-9]+)+)(?=:)', tree)

print("Tips found:", len(tips))
print("First 20 tips:", tips[:20])

if len(tips) == 0:
    raise SystemExit("ERROR: No tips found. You are reading the wrong tree file or it is not Newick text.")

tips = list(dict.fromkeys(tips))

df = pd.read_csv(csvfile)
df["TaxonCode"] = df["TaxonCode"].astype(str).str.strip()
df["Subgroup"] = df["Subgroup"].astype(str).str.strip()

def assign_major_group(subgroup):
    s = subgroup.strip().lower()

    if s in {"mamiellophyceae", "nephroselmidophyceae", "pyramimonadophyceae"}:
        return "Chlorophyte algae"
    if s in {"charales", "charophytes", "chlorokybales", "coleochaetales",
             "klebsormidiales", "klebsormidiophyceae", "mesostigmatales",
             "zygnematophyceae"}:
        return "Streptophyte algae"
    if s in {"hornwort", "liverwort", "moss"}:
        return "Bryophyte"
    if s == "lycophyte":
        return "Lycophyte"
    if s == "monilophyte":
        return "Monilophyte"
    if s == "gymnosperm":
        return "Gymnosperm"
    if s in {"ana grade", "eudicot", "magnoliid", "monocot"}:
        return "Angiosperm"
    return "Unknown"

df["MajorGroup"] = df["Subgroup"].apply(assign_major_group)

colors = {
    "Chlorophyte algae": "#e9d2e1",
    "Streptophyte algae": "#d9ead3",
    "Bryophyte": "#b6d7a8",
    "Lycophyte": "#cfe2f3",
    "Monilophyte": "#9fc5e8",
    "Gymnosperm": "#d9d2e9",
    "Angiosperm": "#fce5cd",
    "Unknown": "#ffffff"
}

code_to_group = dict(zip(df["TaxonCode"], df["MajorGroup"]))

with open(outfile, "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tTaxonomic group\n")
    f.write("COLOR\t#000000\n")
    f.write("LEGEND_TITLE\tTaxonomic group\n")
    f.write("LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\n")
    f.write("LEGEND_COLORS\t#e9d2e1\t#d9ead3\t#b6d7a8\t#cfe2f3\t#9fc5e8\t#d9d2e9\t#fce5cd\n")
    f.write("LEGEND_LABELS\tChlorophyte algae\tStreptophyte algae\tBryophytes\tLycophytes\tMonilophytes (ferns)\tGymnosperms\tAngiosperms\n\n")
    f.write("DATA\n")

    for tip in tips:
        code = tip[:4]
        group = code_to_group.get(code, "Unknown")
        f.write(f"{tip}\t{colors[group]}\t{group}\n")

print("Wrote:", outfile)
