import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.gridspec as gridspec
import textwrap

# Set style
sns.set_theme(style="whitegrid")
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10

# Create figure layout with MAXIMUM spacing
fig = plt.figure(figsize=(15, 7)) # Wider and taller
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.4], wspace=0.5) # More space for labels

# Panel A: Real Violin Plot
ax1 = plt.subplot(gs[0])

try:
    # Load real expression data
    expr_path = r"D:\2026YJ\My_MR_Project\PARK7_Expression_Data.csv"
    df_expr = pd.read_csv(expr_path)
    
    # Log normalization (log1p)
    df_expr['LogExpression'] = np.log1p(df_expr['Expression'])
    
    # Filter out rare cell types if too many
    top_cells = df_expr['Cell_Type'].value_counts().head(10).index
    df_expr = df_expr[df_expr['Cell_Type'].isin(top_cells)]
    
    # Plot
    sns.violinplot(data=df_expr, x='Cell_Type', y='LogExpression', ax=ax1, palette="Blues", linewidth=1)
    
    ax1.set_title('A. PARK7 Expression across Immune Cells', loc='left', fontweight='bold')
    ax1.set_ylabel('Log Normalized Expression')
    ax1.set_xlabel('')
    # Rotate x-labels to prevent overlap
    ax1.tick_params(axis='x', rotation=45)
    
except Exception as e:
    print(f"Error plotting Panel A: {e}")
    ax1.text(0.5, 0.5, f"Error: {e}", ha='center')

# Panel B: Functional Enrichment (GO)
ax2 = plt.subplot(gs[1])

try:
    # Path from previous LS
    csv_path = r"D:\2026YJ\My_MR_Project\FigureS3\B_cell_PARK7_High_GO_Results.csv" 
    df_go = pd.read_csv(csv_path)
    
    # Clean up column names
    df_go.columns = [c.replace('"', '') for c in df_go.columns]
    
    # Sort by p.adjust (ascending) and take top 10
    df_go = df_go.sort_values('p.adjust').head(10)
    
    # Calculate -log10(p.adjust)
    df_go['log_p'] = -np.log10(df_go['p.adjust'])
    
    # AGGRESSIVE WRAPPING: Wrap labels at 30 characters to force multi-line
    def wrap_labels(text, width=30):
        return "\n".join(textwrap.wrap(text, width))
    
    df_go['Description_Wrapped'] = df_go['Description'].apply(lambda x: wrap_labels(x, width=30))
    
    # Create scatter plot (Dot plot style)
    sc = ax2.scatter(x=df_go['FoldEnrichment'], y=df_go['Description_Wrapped'], 
                     s=df_go['Count']*10, c=df_go['log_p'], 
                     cmap='Reds', alpha=0.8, edgecolors='black', linewidth=0.5)
    
    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax2)
    cbar.set_label('-log10(Adj. P-value)')
    
    ax2.set_title('B. GO Enrichment in PARK7-High B Cells', loc='left', fontweight='bold')
    ax2.set_xlabel('Fold Enrichment')
    ax2.grid(True, linestyle='--', alpha=0.3)
    
    # Set smaller font for y-tick labels to help fit
    ax2.tick_params(axis='y', labelsize=9)
    
except Exception as e:
    print(f"Error plotting Panel B: {e}")
    ax2.text(0.5, 0.5, f"Error loading GO data: {e}", ha='center')

# Save figure as PNG 600 DPI
output_path = r"D:\2026YJ\论文写作\spring模板\FigureS6.png"
# Use tight_layout with EXTRA padding for the y-labels
plt.tight_layout(pad=4.0)
plt.savefig(output_path, dpi=600, bbox_inches='tight')
print(f"Figure S6 saved to {output_path}")
