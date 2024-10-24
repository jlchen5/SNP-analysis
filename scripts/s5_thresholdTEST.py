import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取数据
data = pd.read_csv('/home/jialechen/projects/snpDensity/4_result/diff_pvalues.tsv', sep='\t')

# 计算基本统计量
group1_stats = data['mean_group1_1KGp'].describe()
group2_stats = data['mean_group2_bk'].describe()

print("Group 1 (1KGp) Statistics:")
print(group1_stats)

print("\nGroup 2 (BK) Statistics:")
print(group2_stats)


# 绘制箱线图
# plt.figure(figsize=(10, 6))
# sns.boxplot(data=data[['mean_group1_1KGp', 'mean_group2_bk']])
# plt.xticks([0, 1], ['Group 1 (1KG)', 'Group 2 (BK)'])
# plt.ylabel('SNP Density')
# plt.title('SNP Density Distribution')
# plt.savefig('/home/jialechen/projects/snpDensity/4_result/boxplot.png')

# 阈值划定
threshold = data['mean_group2_bk'].quantile(0.25)
print("Threshold for SNP Density(Q1):", threshold)
