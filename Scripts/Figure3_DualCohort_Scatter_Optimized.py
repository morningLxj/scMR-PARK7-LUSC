import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

df = pd.read_csv("My_MR_Project/Final_Classified_Results.csv")

plt.figure(figsize=(10, 8))
sns.set_style("ticks")

mask = (df["pp4_finngen"] > 0.1) | (df["pp4_ilcco"] > 0.1) | (df["gene"].isin(["PARK7","CTSW","TMEM50A"]))
df_h = df[mask].copy()
df_b = df[~mask].copy()
if df_h.shape[0] == 0:
    sc = (df[["pp4_finngen","pp4_ilcco"]].fillna(0).sum(axis=1)).values
    top = np.argsort(-sc)[:10]
    df_h = df.iloc[top].copy()
    df_b = df.drop(df.index[top]).copy()

plt.scatter(df_b["pp4_ilcco"], df_b["pp4_finngen"], color="lightgray", s=50, alpha=0.5)

palette = {"Strong": "#E64B35", "Moderate": "#4DBBD5", "Suggestive": "#00A087", "Weak": "#7F8C8D"}
sns.scatterplot(data=df_h, x="pp4_ilcco", y="pp4_finngen", hue="Evidence_Level", palette=palette, s=150, edgecolor="black", linewidth=1, zorder=10, legend="full")
ax = plt.gca()
ax.legend(loc="upper right", title="Evidence_Level")

mv = max(float(df["pp4_ilcco"].max()), float(df["pp4_finngen"].max())) * 1.1
plt.plot([0, mv], [0, mv], "--", color="gray", alpha=0.5)

texts = []
for _, r in df_h.iterrows():
    texts.append(plt.text(r["pp4_ilcco"], r["pp4_finngen"], f"{r['gene']}\n({r['cell']})", fontsize=9, fontweight="bold"))
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

plt.xlabel("Discovery Cohort (ILCCO) PP4", fontsize=12)
plt.ylabel("Validation Cohort (FinnGen) PP4", fontsize=12)
plt.xlim(-0.02, mv)
plt.ylim(-0.02, mv)

 
out_png = "My_MR_Project/Figure3_DualCohort_Scatter_Optimized.png"
out_pdf = "My_MR_Project/Figure3_DualCohort_Scatter_Optimized.pdf"
plt.tight_layout()
plt.savefig(out_png, dpi=600)
plt.savefig(out_pdf, dpi=600)
print(f"Optimized Figure 3 saved to: {out_png} and {out_pdf}")
