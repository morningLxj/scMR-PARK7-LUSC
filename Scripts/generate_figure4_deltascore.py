import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 设置绘图参数
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

# 读取数据
data_path = 'd:\\2026YJ\\论文写作\\table\\Table 2. Histological subtype specificity of PARK7 colocalization across B-cell lineages.csv'
df = pd.read_csv(data_path)

# 如果数据中没有Delta Score，计算它
if 'Delta Score' not in df.columns:
    df['Delta Score'] = df['Coloc PP4 (LUSC)'] - df['Coloc PP4 (LUAD)']

# 按Delta Score降序排序
df = df.sort_values('Delta Score', ascending=False)

# 创建条形图
plt.figure(figsize=(8, 6))
colors = ['#2E9FDF' if gene == 'PARK7' else '#E7B800' for gene in df['Gene']]
sns.barplot(x='Gene', y='Delta Score', data=df, palette=colors)

# 添加标题和标签
plt.title('Delta Score (PP4_LUSC - PP4_LUAD) by Gene', fontsize=14, fontweight='bold')
plt.xlabel('Gene', fontsize=12)
plt.ylabel('Delta Score', fontsize=12)
plt.xticks(rotation=45, ha='right')

# 添加网格线
plt.grid(True, axis='y', linestyle='--', alpha=0.7)

# 添加数值标签
for i, v in enumerate(df['Delta Score']):
    plt.text(i, v + 0.01, f'{v:.4f}', ha='center', fontsize=10)

# 保存图表
plt.tight_layout()
plt.savefig('d:\\2026YJ\\My_MR_Project\\Figure4_DeltaScore_Barplot.pdf', dpi=300, bbox_inches='tight')
plt.savefig('d:\\2026YJ\\My_MR_Project\\Figure4_DeltaScore_Barplot.png', dpi=300, bbox_inches='tight')

print("Figure 4 (Delta Score Barplot) generated successfully!")
print(f"Data used: {data_path}")
print(f"Number of genes: {len(df)}")
print(f"Top gene: {df.iloc[0]['Gene']} with Delta Score: {df.iloc[0]['Delta Score']:.4f}")
