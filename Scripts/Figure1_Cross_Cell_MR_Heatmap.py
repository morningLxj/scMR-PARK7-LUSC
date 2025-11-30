import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as pe

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

df = pd.read_csv("My_MR_Project/AllCells_MR_Results.csv")
mat = df.pivot(index="gene", columns="cell", values="beta_ivw")
pval = df.pivot(index="gene", columns="cell", values="p_ivw")
row_order = mat.abs().max(axis=1).sort_values(ascending=False).index
col_order = mat.abs().max(axis=0).sort_values(ascending=False).index
mat = mat.loc[row_order, col_order]
pval = pval.loc[row_order, col_order]

fig, ax = plt.subplots(figsize=(12, 9), dpi=150)
sns.heatmap(mat, ax=ax, cmap="RdBu_r", center=0, linewidths=0.5, linecolor="#e6e6e6", cbar_kws={"label": "IVW beta"})
ax.set_xlabel("Cell Type", fontsize=12)
ax.set_ylabel("Gene", fontsize=12)
ax.tick_params(axis="x", rotation=45, labelsize=10)
ax.tick_params(axis="y", rotation=0, labelsize=10)
xs = []
ys = []
ss = []
for i, g in enumerate(mat.index):
    for j, c in enumerate(mat.columns):
        pv = pval.loc[g, c]
        if pd.isna(pv):
            continue
        if pv == 0:
            size = 240
        else:
            size = float(np.clip(-np.log10(pv), 0, 10)) * 24
        if size >= 24:
            xs.append(j + 0.5)
            ys.append(i + 0.5)
            ss.append(size)
ax.scatter(xs, ys, s=ss, facecolors="none", edgecolors="black", linewidths=0.8)
plt.tight_layout()
out_pdf = "My_MR_Project/Figure1_Cross_Cell_MR_Heatmap.pdf"
fig.savefig(out_pdf, dpi=600)
out_png = "My_MR_Project/Figure1_Cross_Cell_MR_Heatmap.png"
fig.savefig(out_png, dpi=600)

fig2, ax2 = plt.subplots(figsize=(12, 9), dpi=150)
sns.heatmap(mat, ax=ax2, cmap="RdBu_r", center=0, linewidths=0.5, linecolor="#e6e6e6", cbar_kws={"label": "IVW beta"})
ax2.set_xlabel("Cell Type", fontsize=12)
ax2.set_ylabel("Gene", fontsize=12)
ax2.tick_params(axis="x", rotation=45, labelsize=10)
ax2.tick_params(axis="y", rotation=0, labelsize=10)
for i, g in enumerate(mat.index):
    for j, c in enumerate(mat.columns):
        pv = pval.loc[g, c]
        if pd.isna(pv):
            continue
        star = ""
        if pv == 0 or pv < 1e-6:
            star = "***"
        elif pv < 1e-3:
            star = "**"
        elif pv < 5e-2:
            star = "*"
        if star:
            fs = 11
            if star == "***":
                fs = 16
            elif star == "**":
                fs = 13
            txt = ax2.text(j + 0.5, i + 0.5, star, ha="center", va="center", color="white", fontsize=fs, fontweight="bold", zorder=10)
            txt.set_path_effects([pe.withStroke(linewidth=2.5, foreground="black")])
plt.tight_layout()
out_pdf2 = "My_MR_Project/Figure1_Cross_Cell_MR_Heatmap_stars.pdf"
fig2.savefig(out_pdf2, dpi=600)
out_png2 = "My_MR_Project/Figure1_Cross_Cell_MR_Heatmap_stars.png"
fig2.savefig(out_png2, dpi=600)
