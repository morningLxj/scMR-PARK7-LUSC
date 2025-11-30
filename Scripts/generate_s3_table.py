import pandas as pd

# 读取原始CSV文件
df = pd.read_csv('d:\\2026YJ\\My_MR_Project\\TCGA_Results\\TCGA_Batch_Survival_Summary.csv')

# 按队列分组并按P-value升序排序
df_sorted = df.groupby('Cohort').apply(lambda x: x.sort_values('PValue', ascending=True)).reset_index(drop=True)

# 格式化数据列
df_sorted['HR'] = df_sorted['HR'].round(4)
df_sorted['Cutoff'] = df_sorted['Cutoff'].round(4)
df_sorted['PValue'] = df_sorted['PValue'].apply(lambda x: '{:.2e}'.format(x))

# 保存优化后的结果
df_sorted.to_csv('d:\\2026YJ\\My_MR_Project\\TCGA_Results\\S3_TCGA_Survival_Detail.csv', index=False)

print("补充表S3已成功生成！")
print("文件路径：d:\\2026YJ\\My_MR_Project\\TCGA_Results\\S3_TCGA_Survival_Detail.csv")
