import pandas as pd
import tarfile
import io
import os
import matplotlib.pyplot as plt
import seaborn as sns

tar_path = "My_MR_Project/clue_result.tar.gz"
out_dir = "My_MR_Project/CMap_Results"
supp_dir = os.path.join(out_dir, "Supplement")
os.makedirs(out_dir, exist_ok=True)
os.makedirs(supp_dir, exist_ok=True)

with tarfile.open(tar_path, "r:gz") as tar:
    member = None
    for m in tar.getmembers():
        if m.name.endswith("query_result.gct"):
            member = m
            break
    if member is None:
        raise SystemExit(1)
    f = tar.extractfile(member)
    s = f.read().decode("utf-8", errors="replace")
    lines = s.splitlines()
    header_row = 0
    for i, line in enumerate(lines[:50]):
        if line.startswith("id") or line.startswith("cid") or line.startswith("Name\t"):
            header_row = i
            break
    df = pd.read_csv(io.StringIO(s), sep="\t", skiprows=header_row)

# columns: id, pert_id, pert_iname, cell_iname, pert_type, raw_cs, norm_cs, moa, target_name
score_col = "norm_cs" if "norm_cs" in df.columns else ("raw_cs" if "raw_cs" in df.columns else None)
if score_col is None:
    raise SystemExit(1)
df["Connectivity Score"] = pd.to_numeric(df[score_col], errors="coerce")

# filter compounds
type_col = "pert_type" if "pert_type" in df.columns else None
if type_col:
    df_cp = df[df[type_col].astype(str).str.contains("cp", case=False, na=False)].copy()
else:
    df_cp = df.copy()

# mechanism
mech = None
if "moa" in df_cp.columns:
    mech = df_cp["moa"].fillna("")
if (mech is None or (mech == "").all()) and "target_name" in df_cp.columns:
    mech = df_cp["target_name"].fillna("")
df_cp["Mechanism"] = mech if mech is not None else ""
df_mech = df_cp[df_cp["Mechanism"].astype(str) != ""].copy()

# mechanism top15
grp = df_mech.groupby("Mechanism", as_index=False)["Connectivity Score"].min()
top_moa = grp.sort_values("Connectivity Score", ascending=True).head(15)
top_moa["Example Drug"] = ""
if "pert_iname" in df_mech.columns:
    for i in range(len(top_moa)):
        mname = top_moa.iloc[i]["Mechanism"]
        sub = df_mech[df_mech["Mechanism"] == mname].sort_values("Connectivity Score", ascending=True)
        if len(sub) > 0:
            top_moa.at[top_moa.index[i], "Example Drug"] = str(sub.iloc[0]["pert_iname"]) if "pert_iname" in sub.columns else ""
moa_csv = os.path.join(out_dir, "Real_CMap_Top_Reversers_MOA.csv")
top_moa.to_csv(moa_csv, index=False)

# mechanism top3 examples and per-mechanism bar charts
moa_top3_rows = []
def sanitize(name: str):
    import re
    s = re.sub(r"[^a-zA-Z0-9]+", "_", str(name))
    return s[:60].strip("_") or "mechanism"
for i in range(len(top_moa)):
    mname = top_moa.iloc[i]["Mechanism"]
    best_score = float(top_moa.iloc[i]["Connectivity Score"]) if pd.notnull(top_moa.iloc[i]["Connectivity Score"]) else float('nan')
    sub = df_mech[df_mech["Mechanism"] == mname].copy()
    sub = sub.sort_values("Connectivity Score", ascending=True)
    drugs = list(sub["pert_iname"].astype(str)[:3]) if "pert_iname" in sub.columns else []
    row = {
        "Mechanism": mname,
        "BestScore": best_score,
        "Top1": drugs[0] if len(drugs) > 0 else "",
        "Top2": drugs[1] if len(drugs) > 1 else "",
        "Top3": drugs[2] if len(drugs) > 2 else "",
    }
    moa_top3_rows.append(row)
    # per-mechanism bar chart for up to 5 drugs
    take = sub[:5]
    if len(take) > 0:
        plt.figure(figsize=(8,5))
        plt.barh(take["pert_iname"], take["Connectivity Score"], color="steelblue")
        plt.xlabel("Connectivity Score")
        plt.title(f"Mechanism: {mname}")
        plt.tight_layout()
        fname = os.path.join(supp_dir, f"Mechanism_{sanitize(mname)}_Bar.png")
        plt.savefig(fname, dpi=300)

moa_top3 = pd.DataFrame(moa_top3_rows)
moa_top3_csv = os.path.join(out_dir, "Real_CMap_Top_Reversers_MOA_Top3.csv")
moa_top3.to_csv(moa_top3_csv, index=False)

# stratified by cell_iname
def cell_top(cell):
    if "cell_iname" not in df_cp.columns:
        return None
    sub = df_cp[df_cp["cell_iname"].astype(str).str.upper() == cell.upper()].copy()
    if "pert_iname" not in sub.columns:
        return None
    g = sub.groupby("pert_iname", as_index=False)["Connectivity Score"].min()
    return g.sort_values("Connectivity Score", ascending=True).head(15)

pc3 = cell_top("PC3")
hek = cell_top("HEK293")
if pc3 is not None:
    pc3.rename(columns={"pert_iname":"Drug Name"}, inplace=True)
    pc3.to_csv(os.path.join(out_dir, "DrugLevel_Top15_PC3.csv"), index=False)
if hek is not None:
    hek.rename(columns={"pert_iname":"Drug Name"}, inplace=True)
    hek.to_csv(os.path.join(out_dir, "DrugLevel_Top15_HEK293.csv"), index=False)

# supplemental figures
drug_csv = os.path.join(out_dir, "Real_CMap_Top_Reversers_DrugLevel.csv")
if os.path.exists(drug_csv):
    dtop = pd.read_csv(drug_csv)
    dtop["Connectivity Score"] = pd.to_numeric(dtop["Connectivity Score"], errors="coerce")
    plt.figure(figsize=(10,6))
    dtop_sorted = dtop.sort_values("Connectivity Score")
    plt.barh(dtop_sorted["Drug Name"], dtop_sorted["Connectivity Score"], color="purple")
    plt.xlabel("Connectivity Score")
    plt.title("Drug-level Top15 Reversal Strength")
    plt.tight_layout()
    plt.savefig(os.path.join(supp_dir, "Drug_Top15_Bar.png"), dpi=300)

if "cell_iname" in df_cp.columns:
    sub_v = df_cp[df_cp["cell_iname"].astype(str).str.upper().isin(["PC3","HEK293"])].copy()
    plt.figure(figsize=(8,6))
    sns.violinplot(data=sub_v, x="cell_iname", y="Connectivity Score")
    plt.title("Score Distribution by Cell Line")
    plt.tight_layout()
    plt.savefig(os.path.join(supp_dir, "Score_Violin_CellLine.png"), dpi=300)

# simple flowchart
plt.figure(figsize=(8,4))
plt.axis("off")
steps = [
    "Load CLUE GCT",
    "Filter compounds",
    "Aggregate (drug/mechanism)",
    "Sort by score",
    "Export tables & plots"
]
y = 0.8
for s in steps:
    plt.text(0.1, y, s, bbox=dict(boxstyle="round", fc="#e0f3f8"))
    y -= 0.18
plt.savefig(os.path.join(out_dir, "Method_Flowchart.png"), dpi=300)
