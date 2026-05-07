#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("missing_data_per_taxon.csv")
df = df.sort_values("Percent_missing", ascending=True)

plt.figure(figsize=(8, max(6, 0.22 * len(df))))
plt.barh(df["Taxon"], df["Percent_missing"])
plt.xlabel("Percent missing/gap characters")
plt.ylabel("Taxon")
plt.title("Missing data per taxon in concatenated alignment")
plt.tight_layout()
plt.savefig("missing_data_per_taxon.pdf", dpi=300)