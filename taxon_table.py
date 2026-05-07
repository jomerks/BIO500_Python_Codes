import pandas as pd
import matplotlib.pyplot as plt

# Load data
taxa = pd.read_csv("taxa_classification.csv")

# Sort
taxa = taxa.sort_values(["MajorGroup", "Subgroup", "Species"])

# Color mapping
colors = {
    "Streptophyte algae": "#d9ead3",
    "Land plants": "#b6d7a8",
    "Vascular plants": "#cfe2f3",
    "Seed plants": "#d9d2e9",
    "Angiosperm": "#fce5cd"
}

# Create figure
fig, ax = plt.subplots(figsize=(12, 20))
ax.axis('off')

# Create table
table = ax.table(
    cellText=taxa.values,
    colLabels=taxa.columns,
    loc='center',
    cellLoc='left'
)

# Adjust font size
table.auto_set_font_size(False)
table.set_fontsize(10)

# Scale table
table.scale(1, 1.5)

# Color rows
for i in range(len(taxa)):
    group = taxa.iloc[i]["MajorGroup"]
    color = colors.get(group, "white")
    for j in range(len(taxa.columns)):
        table[(i+1, j)].set_facecolor(color)

# Bold header
for j in range(len(taxa.columns)):
    table[(0, j)].set_text_props(weight='bold')

# Save
plt.savefig("taxa_table.png", dpi=300, bbox_inches='tight')
plt.savefig("taxa_table.pdf", bbox_inches='tight')

plt.show()