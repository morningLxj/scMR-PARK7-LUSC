import pandas as pd

# 读取生成的TableS1文件
df = pd.read_csv('TableS1_Full_MR_Statistics.csv')

# 显示文件的前5行
print("TableS1_Full_MR_Statistics.csv的前5行：")
print(df.head())

# 显示文件的列名
print("\n文件列名：")
print(df.columns.tolist())

# 显示文件的行数
print(f"\n文件行数：{len(df)}")
