import pandas as pd
import tarfile
import os
import matplotlib.pyplot as plt
import seaborn as sns
import io

tar_path = "My_MR_Project/clue_result.tar.gz"
output_dir = "My_MR_Project/CMap_Results"
os.makedirs(output_dir, exist_ok=True)

print(f"Processing {tar_path}...")

target_df = None
try:
    with tarfile.open(tar_path, "r:gz") as tar:
        gct_files = [m for m in tar.getmembers() if m.name.endswith(".gct")]
        target_file = None
        for f in gct_files:
            name = f.name.lower()
            if ("knndown" in name) or ("result" in name) or ("ps_" in name):
                target_file = f
                break
        if target_file is None and gct_files:
            target_file = gct_files[0]
        if target_file:
            print(f"Found result file: {target_file.name}")
            fobj = tar.extractfile(target_file)
            content = fobj.read().decode('utf-8', errors='replace')
            lines = content.splitlines()
            header_row = 0
            for i, line in enumerate(lines[:30]):
                if line.startswith("id") or line.startswith("cid") or line.startswith("Name\t"):
                    header_row = i
                    break
            df = pd.read_csv(io.StringIO(content), sep='\t', skiprows=header_row)
            target_df = df
        else:
            print("Error: No .gct file found in the archive.")
            raise SystemExit(1)
except Exception as e:
    print(f"Error reading tar file: {e}")
    raise SystemExit(1)

if target_df is not None:
    col_map = {
        'pert_iname': 'Drug Name',
        'score': 'Connectivity Score',
        'pert_type': 'Type',
        'raw_cs': 'Raw Score'
    }
    renamed_cols = {}
    for k, v in col_map.items():
        for col in list(target_df.columns):
            if col == k:
                renamed_cols[col] = v
    target_df = target_df.rename(columns=renamed_cols)

    if 'Connectivity Score' not in target_df.columns:
        if 'raw_cs' in target_df.columns:
            target_df['Connectivity Score'] = target_df['raw_cs']
        elif 'score' in target_df.columns:
            target_df['Connectivity Score'] = target_df['score']
        else:
            for col in target_df.columns:
                if 'score' in str(col).lower():
                    target_df['Connectivity Score'] = target_df[col]
                    break

    if 'Type' in target_df.columns:
        df_cp = target_df[target_df['Type'].astype(str).str.contains("cp", case=False, na=False)].copy()
        if df_cp.empty:
            df_cp = target_df.copy()
    else:
        df_cp = target_df.copy()

    if 'Connectivity Score' not in df_cp.columns:
        print("Error: Connectivity Score column missing after parsing.")
        raise SystemExit(1)

    df_cp['Connectivity Score'] = pd.to_numeric(df_cp['Connectivity Score'], errors='coerce')
    top_reversers = df_cp.sort_values('Connectivity Score', ascending=True).head(15)

    csv_path = os.path.join(output_dir, "Real_CMap_Top_Reversers.csv")
    top_reversers.to_csv(csv_path, index=False)
    print(f"Top reversers saved to: {csv_path}")

    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    scores = top_reversers['Connectivity Score'].astype(float).values
    ys = list(range(len(top_reversers)))
    scatter = plt.scatter(
        x=scores,
        y=ys,
        c=scores,
        cmap='viridis',
        s=500,
        edgecolor='black',
        alpha=0.8
    )
    name_col = None
    for col in ['Drug Name','pert_iname']:
        if col in top_reversers.columns:
            name_col = col
            break
    for i in range(len(top_reversers)):
        score = float(top_reversers['Connectivity Score'].iloc[i])
        label = str(top_reversers[name_col].iloc[i]) if name_col else str(i)
        plt.text(score + 0.5, i, label, fontweight='bold', fontsize=11, va='center')
    plt.yticks([])
    plt.xlabel("Connectivity Score (Negative = Reversal)", fontsize=12)
    plt.title("Figure 6: Top Drug Candidates for Reversing PARK7 Signature\n(Based on Real CMap L1000 Analysis)", fontsize=14, fontweight='bold')
    plt.axvline(-90, linestyle='--', color='red', alpha=0.5)
    cbar = plt.colorbar(scatter)
    cbar.set_label('Reversal Strength')
    plot_path = os.path.join(output_dir, "Figure6_Real_Data.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    print(f"Real Data Plot saved to: {plot_path}")

    try:
        top_drug = top_reversers.iloc[0][name_col] if name_col else str(top_reversers.index[0])
        top_score = float(top_reversers.iloc[0]['Connectivity Score'])
        print("\n[Manuscript Text Suggestion]")
        print(f"To identify therapeutic agents, we queried the CMap L1000 database. The analysis identified **{top_drug}** (Score = {top_score:.2f}) as the top candidate for reversing the PARK7-associated transcriptional signature.")
    except Exception:
        pass
