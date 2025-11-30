import pandas as pd

# 读取原始CSV文件
df = pd.read_csv('d:\\2026YJ\\My_MR_Project\\Final_Subtype_Specificity.csv')

# 按差异性(Delta)降序排列
df_sorted = df.sort_values('Delta', ascending=False)

# 格式化数值列，保留4位小数
df_sorted['PP4_LUAD'] = df_sorted['PP4_LUAD'].round(4)
df_sorted['PP4_LUSC'] = df_sorted['PP4_LUSC'].round(4)
df_sorted['Delta'] = df_sorted['Delta'].round(4)

# 保存优化后的结果
df_sorted.to_csv('d:\\2026YJ\\My_MR_Project\\S2_Subtype_Detail.csv', index=False)

print("补充表S2已成功生成！")
print("文件路径：d:\\2026YJ\\My_MR_Project\\S2_Subtype_Detail.csv")
