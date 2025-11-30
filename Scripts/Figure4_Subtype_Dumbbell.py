import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "Arial"

# 1. 读取数据
df = pd.read_csv('Final_Subtype_Specificity.csv')

# 筛选差异大的或者核心关注的基因
# 逻辑：Delta > 0.1 或者 是我们的主角 (PARK7, CTSW, TMEM50A)
plot_df = df[
    (df['Delta'] > 0.1) |
    (df['gene'].isin(['PARK7', 'CTSW', 'TMEM50A']))
].copy()

# 排序以便画图好看
plot_df = plot_df.sort_values('Delta', ascending=True)
plot_df['PP4_LUAD'] = pd.to_numeric(plot_df['PP4_LUAD'], errors='coerce')
plot_df['PP4_LUSC'] = pd.to_numeric(plot_df['PP4_LUSC'], errors='coerce')
plot_df['Delta'] = pd.to_numeric(plot_df['Delta'], errors='coerce')
plot_df = plot_df.dropna(subset=['PP4_LUAD','PP4_LUSC','Delta'])

# 2. 创建单个哑铃图
fig, ax = plt.subplots(figsize=(10, 8), dpi=150)

sns.set_style("whitegrid")

# 画连接线 (哑铃的杆)
ax.hlines(y=range(len(plot_df)),
         xmin=plot_df['PP4_LUAD'],
         xmax=plot_df['PP4_LUSC'],
         color='grey', alpha=0.5, linewidth=2, zorder=1)

# 画点 (哑铃的头)
# LUSC 点
ax.scatter(plot_df['PP4_LUSC'], range(len(plot_df)),
          color='#E64B35', s=150, label='Squamous (LUSC)', edgecolor='black', linewidth=0.6, alpha=0.9, zorder=3)
# LUAD 点
ax.scatter(plot_df['PP4_LUAD'], range(len(plot_df)),
          color='#4DBBD5', s=150, label='Adenocarcinoma (LUAD)', edgecolor='black', linewidth=0.6, alpha=0.9, zorder=3)

# 添加基因标签
texts = []
for i, row in enumerate(plot_df.itertuples()):
    mid_point = (row.PP4_LUAD + row.PP4_LUSC) / 2
    label = f"{row.gene} ({row.cell})"
    texts.append(ax.text(mid_point, i + 0.15, label,
                        ha='center', va='bottom', fontweight='bold', fontsize=10))
adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.6))

# 装饰哑铃图
ax.set_yticks([])
ax.set_xlabel("Colocalization Probability (PP4)", fontsize=12, fontweight='bold')
ax.legend(loc='lower right', frameon=True, fontsize=10)
ax.set_xlim(-0.05, 1.05)
sns.despine(ax=ax)

# 整体布局调整
plt.tight_layout()

# 3. 保存
output_path = 'Figure4_Subtype_Dumbbell.png'
output_path_png600 = 'Figure4_Subtype_Dumbbell_600dpi.png'
output_path_pdf600 = 'Figure4_Subtype_Dumbbell_600dpi.pdf'

fig.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Figure 4 saved to: {output_path}")

fig.savefig(output_path_png600, dpi=600, bbox_inches='tight')
print(f"Figure 4 saved to: {output_path_png600}")

fig.savefig(output_path_pdf600, dpi=600, bbox_inches='tight')
print(f"Figure 4 saved to: {output_path_pdf600}")
