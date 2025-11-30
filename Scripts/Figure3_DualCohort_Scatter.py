import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text
from matplotlib.lines import Line2D

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

df = pd.read_csv("My_MR_Project/Final_Classified_Results.csv")
df = df.dropna(subset=["pp4_ilcco", "pp4_finngen"]) 

def p_to_size(p):
    if pd.isna(p):
        return 40
    if p == 0:
        return 200
    return float(np.clip(-np.log10(p), 0, 10)) * 20

size = df["p_ivw"].apply(p_to_size) if "p_ivw" in df.columns else df["mr_p"].apply(p_to_size)
sign = np.where(df["beta_ivw"] >= 0, "pos", "neg")
palette = {"pos": "#E64B35", "neg": "#4DBBD5"}

fig, ax = plt.subplots(figsize=(8, 8), dpi=150)
sns.set_style("whitegrid")
cells = df["cell"].astype(str).unique().tolist()
markers = ["o","s","^","D","P","X","v","<",">","h"]
marker_map = {c: markers[i % len(markers)] for i, c in enumerate(cells)}
for c in cells:
    msk = df["cell"].astype(str) == c
    col = np.where(df.loc[msk, "beta_ivw"] >= 0, palette["pos"], palette["neg"])
    ax.scatter(df.loc[msk, "pp4_ilcco"], df.loc[msk, "pp4_finngen"], s=size.loc[msk], c=col, alpha=0.75, edgecolors="black", linewidths=0.6, marker=marker_map[c])
ax.plot([0,1],[0,1], color="grey", linewidth=1, linestyle="--")
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.set_xlabel("PP4 (ILCCO)", fontsize=12)
ax.set_ylabel("PP4 (FinnGen)", fontsize=12)
sign_handles = [Line2D([],[],marker="o",color=palette["pos"],markerfacecolor=palette["pos"],markeredgecolor="black",linewidth=0,label="Positive"),
                Line2D([],[],marker="o",color=palette["neg"],markerfacecolor=palette["neg"],markeredgecolor="black",linewidth=0,label="Negative")]
leg1 = ax.legend(handles=sign_handles, loc="lower right", title="Effect")
ax.add_artist(leg1)
shape_handles = [Line2D([],[],marker=marker_map[c],color="black",markerfacecolor="none",markeredgecolor="black",linewidth=0,label=c) for c in cells[:10]]
ax.legend(handles=shape_handles, loc="upper left", bbox_to_anchor=(1.02,1), borderaxespad=0, title="Cell", frameon=True)

texts = []
ann = df[df["gene"].isin(["PARK7","CTSW","TMEM50A"])].copy()
for _, r in ann.iterrows():
    texts.append(ax.text(r["pp4_ilcco"], r["pp4_finngen"], f"{r['gene']} ({r['cell']})", fontsize=10))
adjust_text(texts, ax=ax, expand_points=(1.2,1.3), expand_text=(1.2,1.4), arrowprops=dict(arrowstyle="-", color="grey", lw=0.8))

plt.tight_layout()
pdf_out = "My_MR_Project/Figure3_DualCohort_Scatter.pdf"
png_out = "My_MR_Project/Figure3_DualCohort_Scatter.png"
fig.savefig(pdf_out)
fig.savefig(png_out, dpi=600)
print(f"Figure 3 saved to: {pdf_out} and {png_out}")
