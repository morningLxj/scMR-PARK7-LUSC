import pandas as pd
import numpy as np
import os

try:
    import seaborn as sns
    import matplotlib.pyplot as plt
except Exception:
    sns = None
    import matplotlib.pyplot as plt

df = pd.read_csv('My_MR_Project/AllCells_MR_Results.csv')
genes = ['PARK7','CTSW','TMEM50A']
df_sel = df[df['gene'].isin(genes)].copy()
heatmap_data = df_sel.pivot(index='cell', columns='gene', values='beta_ivw')
p_data = df_sel.pivot(index='cell', columns='gene', values='p_ivw')
plt.figure(figsize=(10,8))
if sns is not None:
    ax = sns.heatmap(heatmap_data, annot=True, cmap='RdBu_r', center=0, fmt='.2f')
else:
    im = plt.imshow(heatmap_data.values, cmap='RdBu_r', aspect='auto')
    for (i,j), val in np.ndenumerate(heatmap_data.values):
        plt.text(j, i, f'{val:.2f}', ha='center', va='center', color='black')
    plt.xticks(range(heatmap_data.shape[1]), heatmap_data.columns, rotation=0)
    plt.yticks(range(heatmap_data.shape[0]), heatmap_data.index)
for y in range(heatmap_data.shape[0]):
    for x in range(heatmap_data.shape[1]):
        pval = p_data.iloc[y, x]
        if pval < 5e-8:
            plt.text(x + 0.5, y + 0.7, '★', color='black', ha='center', va='center', fontsize=16)
plt.title('MR Effect Sizes (Beta) across Cell Types\n★ = p < 5e-8')
plt.tight_layout()
out_png = 'My_MR_Project/Cross_Cell_MR_Heatmap.png'
plt.savefig(out_png)
print('Heatmap saved to', out_png)
