import pandas as pd
import re

treefile = "4152/OG004152.codes.tree"
csvfile = "taxa_classification_corrected.csv"
outfile = "4152/itol_taxonomy_colors_OG004152.txt"

df = pd.read_csv(csvfile)
df["TaxonCode"] = df["TaxonCode"].astype(str).str.strip()
df["MajorGroup"] = df["MajorGroup"].astype(str).str.strip()

colors = {
    "Chlorophyte algae": "#c9daf8",
    "Streptophyte algae": "#d9ead3",
    "Bryophyte": "#b6d7a8",
    "Lycophyte": "#cfe2f3",
    "Monilophyte": "#9fc5e8",
    "Gymnosperm": "#d9d2e9",
    "Angiosperm": "#fce5cd"
}

code_to_group = dict(zip(df["TaxonCode"], df["MajorGroup"]))

with open(treefile) as f:
    tree = f.read()

# Extract all tip labels from Newick tree
tips = re.findall(r'(?<=[(,])([^():,;]+):', tree)

with open(outfile, "w") as f:
    f.write("DATASET_COLORSTRIP\n")
    f.write("SEPARATOR TAB\n")
    f.write("DATASET_LABEL\tTaxonomic group\n")
    f.write("COLOR\t#000000\n")
    f.write("LEGEND_TITLE\tTaxonomic group\n")
    f.write("LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\n")
    f.write("LEGEND_COLORS\t#c9daf8\t#d9ead3\t#b6d7a8\t#cfe2f3\t#9fc5e8\t#d9d2e9\t#fce5cd\n")
    f.write("LEGEND_LABELS\tChlorophyte algae\tStreptophyte algae\tBryophytes\tLycophytes\tMonilophytes (ferns)\tGymnosperms\tAngiosperms\n\n")
    f.write("DATA\n")

    for tip in tips:
        code = tip.split("_")[0]
        group = code_to_group.get(code)

        if group is None:
            print(f"WARNING: no group found for {tip}")
            continue

        color = colors.get(group, "#ffffff")
        f.write(f"{tip}\t{color}\t{group}\n")

print(f"Wrote {outfile}")