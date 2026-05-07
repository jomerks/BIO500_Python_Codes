import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Y = supported, N = not supported, ? = unclear/ambiguous
data = {
    "Tree": [
        "OG4257 ML", "OG4257 MP",
        "OG4414 ML", "OG4414 MP",
        "OG4771 ML", "OG4771 MP",
        "OG4947 ML", "OG4947 MP",
        "OG4958 ML", "OG4958 MP",
        "Concatenated"
    ],
    "Chlorophytes basal": ["Y","Y","Y","Y","Y","Y","Y","Y","Y","Y","Y"],
    "Streptophyte algae sister\nto land plants": ["?","?","?","?","Y","?","?","?","?","?","Y"],
    "Bryophytes\nclade recovered": ["?","?","?","N","?","?","?","?","?","?","Y"],
    "Vascular plants\nmonophyletic": ["Y","Y","Y","?","Y","?","N","Y","N","Y","Y"],
    "Seed plants\nmonophyletic": ["Y","Y","Y","N","Y","?","N","N","N","N","Y"],
    "Angiosperms\nmonophyletic": ["Y","?","Y","N","Y","Y","?","?","Y","?","Y"],
}

df = pd.DataFrame(data).set_index("Tree")

# Convert symbols to numeric scores for heatmap
score_map = {"N": 0, "?": 0.5, "Y": 1}
heatmap_data = df.replace(score_map).astype(float)

# Plot
plt.figure(figsize=(11, 6))

ax = sns.heatmap(
    heatmap_data,
    annot=df,
    fmt="",
    cmap=sns.color_palette(["#e89692", "#fcf1bd", "#93c2a8"], as_cmap=True),
    vmin=0,
    vmax=1,
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Topology support"}
)

# Customize colorbar labels
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0, 0.5, 1])
colorbar.set_ticklabels(["Not supported", "Unclear", "Supported"])

plt.title("Topology Concordance Across Single-Gene and Concatenated Trees", fontsize=14, pad=15)
plt.xlabel("Topology criterion")
plt.ylabel("Tree")
plt.xticks(rotation=35, ha="right")
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig("topology_concordance_heatmap.png", dpi=300)
plt.savefig("topology_concordance_heatmap.pdf")
plt.show()
