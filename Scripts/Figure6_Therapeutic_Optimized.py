import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

drug_csv = "My_MR_Project/CMap_Results/Real_CMap_Top_Reversers_DrugLevel.csv"
mech_csv = "My_MR_Project/CMap_Results/DrugMechanism_Top15.csv"
df_drug = pd.read_csv(drug_csv)
df_drug["Connectivity Score"] = pd.to_numeric(df_drug["Connectivity Score"], errors="coerce")
df_drug = df_drug.dropna(subset=["Connectivity Score"]).copy()
df_mech = pd.read_csv(mech_csv)
df_mech["Connectivity Score"] = pd.to_numeric(df_mech["Connectivity Score"], errors="coerce")
df_mech = df_mech.dropna(subset=["Connectivity Score"]).copy()
df_mech = df_mech.sort_values("Connectivity Score")
df_mech = df_mech.groupby("Drug Name").apply(lambda g: g.loc[~g["Mechanism"].astype(str).str.contains("-666", na=False)].iloc[0] if (~g["Mechanism"].astype(str).str.contains("-666", na=False)).any() else g.iloc[0]).reset_index(drop=True)
df = df_drug.merge(df_mech[["Drug Name","Mechanism"]], on="Drug Name", how="left")
df = df.head(15).copy()
def to_group(m: str) -> str:
    s = str(m).lower()
    if "topoisomerase" in s:
        return "Topoisomerase"
    if "atm" in s or "atr" in s:
        return "ATM/ATR"
    return "Other"
df["Group"] = df["Mechanism"].apply(to_group)
order = {"Topoisomerase": 0, "ATM/ATR": 1, "Other": 2}
df = df.sort_values(["Group"], key=lambda col: col.map(order)).copy()
group_sizes = df.groupby("Group").size().reindex(["Topoisomerase","ATM/ATR","Other"], fill_value=0).tolist()
offsets = []
b = 0
for gsz in group_sizes:
    offsets.append(b)
    b += gsz + 1
grp_to_offset = {g:o for g,o in zip(["Topoisomerase","ATM/ATR","Other"], offsets)}
df["y"] = df.groupby("Group").cumcount() + df["Group"].map(grp_to_offset)
labels = df["Drug Name"].astype(str).tolist()
scores = df["Connectivity Score"].astype(float).values
ys = df["y"].values

fig, ax = plt.subplots(figsize=(10, 8), dpi=150)
sns.set_style("whitegrid")
ring = df["Group"].map({"Topoisomerase": "#b2182b", "ATM/ATR": "#4c1d95", "Other": "#bdbdbd"}).values
scatter = ax.scatter(x=scores, y=ys, c=scores, cmap='viridis', s=500, edgecolor=ring, linewidth=1.2, alpha=0.9)

for i in range(len(df)):
    score = float(df["Connectivity Score"].iloc[i])
    lab = labels[i]
    ax.annotate(lab,
                xy=(score, ys[i]), xycoords='data',
                xytext=(score + 0.08, ys[i]), textcoords='data',
                ha='left', va='center', fontsize=10, fontweight='bold',
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.6, shrinkA=5, shrinkB=5))

ax.set_yticks([])
ax.set_xlabel("Connectivity Score (Negative = Reversal)", fontsize=12)
ax.axvline(-1.00, linestyle='--', color='#3b4cc0', alpha=0.6)
for sep in offsets[1:]:
    ax.hlines(y=sep-0.5, xmin=min(scores)-0.05, xmax=max(scores)+0.05, color="#e0e0e0", linestyle="--", linewidth=0.8)
for t, off in grp_to_offset.items():
    n = df[df["Group"]==t].shape[0]
    if n > 0:
        ax.text(min(scores)-0.05, off + (n-1)/2, t, ha="right", va="center", fontsize=12)
sns.despine()

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Reversal Strength')
approved = {s.lower() for s in ["dexamethasone","dexamethasone-acetate","irinotecan","topotecan","teniposide","mitomycin-c","bosutinib","ruxolitinib","metformin"]}
ap_msk = df["Drug Name"].str.lower().isin(approved)
if ap_msk.any():
    ax.scatter(scores[ap_msk], ys[ap_msk], marker='*', s=320, color="#d97706", edgecolor="#000000", linewidth=0.9, zorder=5)

plt.tight_layout()
plt.tight_layout()
out_png = "My_MR_Project/Figure6_Therapeutic_Grouped_600dpi_v2.png"
out_pdf = "My_MR_Project/Figure6_Therapeutic_Grouped_600dpi_v2.pdf"
fig.savefig(out_png, dpi=600, bbox_inches='tight')
fig.savefig(out_pdf, dpi=600, bbox_inches='tight')
print(f"Figure 6 saved to: {out_png} and {out_pdf}")
