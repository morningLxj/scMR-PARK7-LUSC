import pandas as pd
import tarfile
import os
import io

tar_path = "My_MR_Project/clue_result.tar.gz"
output_dir = "My_MR_Project/CMap_Results"
os.makedirs(output_dir, exist_ok=True)

df_drug = None
with tarfile.open(tar_path, "r:gz") as tar:
    members = [m for m in tar.getmembers() if m.name.endswith(".gct")]
    for m in members:
        fobj = tar.extractfile(m)
        if not fobj:
            continue
        content = fobj.read().decode("utf-8", errors="replace")
        lines = content.splitlines()
        header_row = 0
        for i, line in enumerate(lines[:40]):
            if line.startswith("id") or line.startswith("cid") or line.startswith("Name\t"):
                header_row = i
                break
        df = pd.read_csv(io.StringIO(content), sep="\t", skiprows=header_row)
        cols = [str(c) for c in df.columns]
        has_drug = "pert_iname" in cols
        has_score = ("norm_cs" in cols) or ("raw_cs" in cols) or ("score" in cols)
        if has_drug and has_score:
            df_drug = df
            break

if df_drug is None:
    raise SystemExit(1)

if "pert_iname" in df_drug.columns:
    df_drug.rename(columns={"pert_iname": "Drug Name"}, inplace=True)
if "pert_type" in df_drug.columns:
    df_drug.rename(columns={"pert_type": "Type"}, inplace=True)

if "Connectivity Score" not in df_drug.columns:
    if "norm_cs" in df_drug.columns:
        df_drug["Connectivity Score"] = df_drug["norm_cs"]
    elif "raw_cs" in df_drug.columns:
        df_drug["Connectivity Score"] = df_drug["raw_cs"]
    elif "score" in df_drug.columns:
        df_drug["Connectivity Score"] = df_drug["score"]

if "Connectivity Score" not in df_drug.columns:
    raise SystemExit(1)
df_drug["Connectivity Score"] = pd.to_numeric(df_drug["Connectivity Score"], errors="coerce")
if "Type" in df_drug.columns:
    df_drug = df_drug[df_drug["Type"].astype(str).str.contains("cp", case=False, na=False)]

grouped = df_drug.groupby("Drug Name", as_index=False)["Connectivity Score"].min()
top15 = grouped.sort_values("Connectivity Score", ascending=True).head(15)

out_csv = os.path.join(output_dir, "Real_CMap_Top_Reversers_DrugLevel.csv")
top15.to_csv(out_csv, index=False)
print(out_csv)
