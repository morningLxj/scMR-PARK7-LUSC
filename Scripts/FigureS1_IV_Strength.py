import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.font_manager as fm
from adjustText import adjust_text

# Settings
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'

# Paths
BASE_DIR = Path("My_MR_Project")
INPUT_CSV = BASE_DIR / "AllCells_MR_Results.csv"
EXP_DIR = BASE_DIR / "Exposure"
OUTPUT_PNG = BASE_DIR / "FigureS1_IV_Strength.png"
OUTPUT_PDF = BASE_DIR / "FigureS1_IV_Strength.pdf"

def calculate_f_stats(row):
    cell = row['cell']
    gene = row['gene']
    # Construct exposure filename
    fname = f"EXPOSURE_{cell}_{gene}_mr.tsv.gz"
    fpath = EXP_DIR / fname
    
    if not fpath.exists():
        print(f"Warning: Exposure file not found: {fname}")
        return np.nan

    try:
        # Read exposure data
        df_exp = pd.read_csv(fpath, sep='\t', compression='gzip')
        
        # Ensure numeric
        df_exp['beta'] = pd.to_numeric(df_exp['beta'], errors='coerce')
        df_exp['se'] = pd.to_numeric(df_exp['se'], errors='coerce')
        df_exp = df_exp.dropna(subset=['beta', 'se'])
        df_exp = df_exp[df_exp['se'] > 0]
        
        if df_exp.empty:
            return np.nan
            
        # Calculate F-statistic per SNP: F = beta^2 / se^2
        df_exp['F_stat'] = (df_exp['beta'] / df_exp['se']) ** 2
        
        # Return mean F-statistic
        return df_exp['F_stat'].mean()
        
    except Exception as e:
        print(f"Error processing {fname}: {e}")
        return np.nan

def main():
    # Load results
    if not INPUT_CSV.exists():
        print(f"Error: {INPUT_CSV} not found.")
        return
    
    df = pd.read_csv(INPUT_CSV)
    print(f"Loaded {len(df)} rows from {INPUT_CSV}")
    
    # Calculate F-statistics
    print("Calculating F-statistics from exposure files...")
    df['Mean_F_statistic'] = df.apply(calculate_f_stats, axis=1)
    
    # Drop rows where F calculation failed
    df_plot = df.dropna(subset=['Mean_F_statistic']).copy()
    print(f"Plotting {len(df_plot)} valid entries.")
    
    # SCI Optimization Settings
    plt.figure(figsize=(8, 6))  # Compact SCI size
    sns.set_style("whitegrid", {'axes.grid': True, 'grid.linestyle': '--', 'grid.alpha': 0.5})
    
    # Define color palette (using a high-contrast scientific palette)
    palette = sns.color_palette("Set2", n_colors=len(df_plot['gene'].unique()))
    
    # Sort by F-statistic or Cell Type for better visual layout
    df_plot = df_plot.sort_values('cell', ascending=False)
    
    # Scatter Plot with optimized aesthetics using stripplot for better categorical handling
    ax = sns.stripplot(
        data=df_plot,
        x='Mean_F_statistic',
        y='cell',
        hue='gene',
        jitter=0.15,  # Add slight vertical jitter to avoid perfect overlap
        size=9,       # Moderate size for clarity
        palette=palette,
        edgecolor='black',
        linewidth=0.8,
        alpha=0.9,
        zorder=3
    )
    
    # Add vertical line at F=10
    plt.axvline(x=10, color='#d62728', linestyle='--', linewidth=1.2, alpha=0.8, zorder=2)
    ymin, ymax = ax.get_ylim()
    plt.text(10.5, ymin, 'Threshold (F = 10)', color='#d62728', fontsize=10, va='bottom', fontweight='bold')
    
    # Add vertical line at F=30 if applicable
    max_f = df_plot['Mean_F_statistic'].max()
    if max_f > 32:
        plt.axvline(x=30, color='#2ca02c', linestyle=':', linewidth=1.2, alpha=0.8, zorder=2)
        plt.text(30.5, ymin, 'Strong IV (F > 30)', color='#2ca02c', fontsize=10, va='bottom', fontweight='bold')

    # Labels and Title
    plt.xlabel('Mean F-statistic', fontsize=12, fontweight='bold', labelpad=10)
    plt.ylabel('Cell Type', fontsize=12, fontweight='bold', labelpad=10)
    
    # Customize axis spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.2)
    ax.spines['bottom'].set_linewidth(1.2)
    ax.tick_params(width=1.2, labelsize=10)
    
    # Legend
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, title='Gene', frameon=False, fontsize=10, title_fontsize=11)
    
    plt.tight_layout()
    
    # Save
    plt.savefig(OUTPUT_PNG, dpi=600, bbox_inches='tight')
    plt.savefig(OUTPUT_PDF, dpi=600, bbox_inches='tight')
    print(f"Saved to {OUTPUT_PNG} and {OUTPUT_PDF}")

if __name__ == "__main__":
    main()
