import pandas as pd
import tarfile
import io
import os
import matplotlib.pyplot as plt
import seaborn as sns

tar_path = "My_MR_Project/clue_result.tar.gz"
out_dir = "My_MR_Project/CMap_Results"
os.makedirs(out_dir, exist_ok=True)

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
    for i, line in enumerate(lines[:60]):
        if line.startswith("id") or line.startswith("cid") or line.startswith("Name\t"):
            header_row = i
            break
    df = pd.read_csv(io.StringIO(s), sep="\t", skiprows=header_row)

drug_csv = os.path.join(out_dir, "Real_CMap_Top_Reversers_DrugLevel.csv")
if not os.path.exists(drug_csv):
    raise SystemExit(1)
dtop = pd.read_csv(drug_csv)

df["pert_iname"] = df["pert_iname"].astype(str)
dtop["Drug Name"] = dtop["Drug Name"].astype(str)
map_cols = [c for c in ["moa","target_name","cell_iname"] if c in df.columns]
map_df = df[["pert_iname"] + map_cols].drop_duplicates()
merged = dtop.merge(map_df, left_on="Drug Name", right_on="pert_iname", how="left")
if "moa" in merged.columns:
    mech = merged["moa"].fillna("")
else:
    mech = pd.Series([""]*len(merged))
if (mech == "").all() and "target_name" in merged.columns:
    mech = merged["target_name"].fillna("")
merged["Mechanism"] = mech
out_csv = os.path.join(out_dir, "DrugMechanism_Top15.csv")
merged[["Drug Name","Mechanism","Connectivity Score"]].to_csv(out_csv, index=False)

plt.figure(figsize=(10,6))
sorted_m = merged.sort_values("Connectivity Score")
plt.barh(sorted_m["Drug Name"], sorted_m["Connectivity Score"], color="teal")
plt.xlabel("Connectivity Score")
plt.title("Drug-Mechanism Top15")
plt.tight_layout()
fig_path = os.path.join(out_dir, "DrugMechanism_Top15_Bar.png")
plt.savefig(fig_path, dpi=300)
print(out_csv)
print(fig_path)
