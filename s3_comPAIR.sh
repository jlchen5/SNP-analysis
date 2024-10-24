#!/bin/bash

# 定义文件路径
GROUP_FILE="/home/jialechen/projects/snpDensity/3_merged_data/group.tsv"
SNPDEN_FILE="/home/jialechen/projects/snpDensity/3_merged_data/snpden.txt"
OUTPUT_DIR="/home/jialechen/projects/snpDensity/4_result"
PYTHON_SCRIPT="$OUTPUT_DIR/plot_pairplot.py"

# 创建输出目录
mkdir -p $OUTPUT_DIR

# 提取分组信息
awk -v OFS='\t' '{print $1, $2}' $GROUP_FILE > $OUTPUT_DIR/group_info.tsv

# 提取样本ID和分组信息
awk '{print $1}' $GROUP_FILE | grep -v "ID" > $OUTPUT_DIR/sample_ids.txt
awk '{print $2}' $GROUP_FILE | grep -v "Group" > $OUTPUT_DIR/groups.txt

# 处理snpden.txt文件，提取需要的列
cut -f1-3,$(awk '{printf "%s,", NR+3}' $OUTPUT_DIR/sample_ids.txt | sed 's/,$//') $SNPDEN_FILE > $OUTPUT_DIR/snpden_filtered.txt

# 计算每个位点在两组之间的差异，计算p值
python3 - <<EOF
import pandas as pd
from scipy.stats import ttest_ind

# 读取数据
group_info = pd.read_csv('$OUTPUT_DIR/group_info.tsv', sep='\t', header=None, names=['SampleID', 'Group'])
snpden_data = pd.read_csv('$OUTPUT_DIR/snpden_filtered.txt', sep='\t')

# 分组
group1_samples = group_info[group_info['Group'] == '1KGp']['SampleID'].values
group2_samples = group_info[group_info['Group'] == 'bk']['SampleID'].values

# 计算均值和差异
group1_data = snpden_data[group1_samples]
group2_data = snpden_data[group2_samples]

mean_group1 = group1_data.mean(axis=1)
mean_group2 = group2_data.mean(axis=1)
mean_diff = mean_group1 - mean_group2
p_values = ttest_ind(group1_data, group2_data, axis=1).pvalue

# 计算整体样本的统计量
all_data = snpden_data.iloc[:, 3:]
overall_mean = all_data.mean(axis=1)
overall_median = all_data.median(axis=1)
overall_q1 = all_data.quantile(0.25, axis=1)
overall_q3 = all_data.quantile(0.75, axis=1)

# 保存结果
result = snpden_data.iloc[:, :3].copy()
result['mean_group1_1KGp'] = mean_group1
result['mean_group2_bk'] = mean_group2
result['mean_diff'] = mean_diff
result['p_value'] = p_values
result['overall_mean'] = overall_mean
result['overall_median'] = overall_median
result['overall_q1'] = overall_q1
result['overall_q3'] = overall_q3
result.to_csv('$OUTPUT_DIR/diff_pvalues.tsv', sep='\t', index=False)
EOF

# 生成绘图的Python脚本
cat <<EOF > $PYTHON_SCRIPT
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# 读取数据
data = pd.read_csv('$OUTPUT_DIR/diff_pvalues.tsv', sep='\t')

# 绘制pair图
sns.pairplot(data[['mean_group1_1KGp', 'mean_group2_bk', 'mean_diff', 'p_value']])
plt.savefig('$OUTPUT_DIR/pairplot.png')
EOF

# 运行Python脚本绘制图像
python3 $PYTHON_SCRIPT

echo "分析完成，结果保存在 $OUTPUT_DIR"
