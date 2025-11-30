import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

result_file = 'My_MR_Project/AllCells_MR_Results.csv'
if not os.path.exists(result_file):
    print(f'File not found: {result_file}')
    raise SystemExit(1)

df = pd.read_csv(result_file)
genes = ['PARK7', 'CTSW', 'TMEM50A']
df_sel = df[df['gene'].isin(genes)].copy()
beta_mat = df_sel.pivot(index='cell', columns='gene', values='beta_ivw')
pval_mat = df_sel.pivot(index='cell', columns='gene', values='p_ivw')
plt.figure(figsize=(10, 8))
ax = sns.heatmap(beta_mat, annot=True, cmap='RdBu_r', center=0, fmt='.2f', cbar_kws={'label': 'MR Beta (Effect Size)'})
for y in range(beta_mat.shape[0]):
    for x in range(beta_mat.shape[1]):
        pval = pval_mat.iloc[y, x]
        if not np.isnan(pval):
            if pval < 5e-8:
                ax.text(x + 0.5, y + 0.3, '★', ha='center', va='center', color='black', fontsize=18, fontweight='bold')
            elif pval < 1e-3:
                ax.text(x + 0.5, y + 0.3, '•', ha='center', va='center', color='black', fontsize=18)
plt.title('MR Effect of Genes on Lung Cancer (ILCCO) across Immune Cells')
plt.tight_layout()
plt.savefig('My_MR_Project/Cross_Cell_MR_Heatmap.png', dpi=300)
print('Heatmap saved to My_MR_Project/Cross_Cell_MR_Heatmap.png')
