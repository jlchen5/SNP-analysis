#!/bin/bash

# 检查输入参数
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <chromosome> <start> <end> <sample_file> <population_file>"
    exit 1
fi

CHROM=$1
START=$2
END=$3
SAMPLE_FILE=$4
POPULATION_FILE=$5

# 提取样本文件中的SNP密度
SAMPLE_SUM=$(awk -v chrom="$CHROM" -v start="$START" -v end="$END" '
BEGIN {sum=0; count=0}
NR>1 && $1==chrom && $2>=start && $2<end {
    sum+=$3; count++
}
END {if (count > 0) print sum/count; else print "NA"}
' "$SAMPLE_FILE")

# 提取总体文件中的SNP密度统计数据
POPULATION_STATS=$(awk -v chrom="$CHROM" -v start="$START" -v end="$END" '
NR>1 && $1==chrom && $2>=start && $3<=end {
    overall_mean=$8; overall_median=$9; overall_q1=$10; overall_q3=$11;
    print overall_mean, overall_median, overall_q1, overall_q3;
    exit
}
' "$POPULATION_FILE")

# 解析总体统计数据
OVERALL_MEAN=$(echo "$POPULATION_STATS" | awk '{print $1}')
OVERALL_MEDIAN=$(echo "$POPULATION_STATS" | awk '{print $2}')
OVERALL_Q1=$(echo "$POPULATION_STATS" | awk '{print $3}')
OVERALL_Q3=$(echo "$POPULATION_STATS" | awk '{print $4}')

# 计算差异
if [ "$SAMPLE_SUM" != "NA" ] && [ -n "$OVERALL_MEAN" ]; then
    DIFF_MEAN=$(echo "$SAMPLE_SUM - $OVERALL_MEAN" | bc -l)
    DIFF_MEDIAN=$(echo "$SAMPLE_SUM - $OVERALL_MEDIAN" | bc -l)
    DIFF_Q1=$(echo "$SAMPLE_SUM - $OVERALL_Q1" | bc -l)
    DIFF_Q3=$(echo "$SAMPLE_SUM - $OVERALL_Q3" | bc -l)
else
    DIFF_MEAN="NA"
    DIFF_MEDIAN="NA"
    DIFF_Q1="NA"
    DIFF_Q3="NA"
fi

# 计算GC含量
# 假设参考基因组文件为 hg19.fa
REFERENCE_GENOME="/home/jialechen/library/hg19/hg19.fa"

# 使用samtools和bedtools提取序列
SEQUENCE=$(samtools faidx "$REFERENCE_GENOME" "${CHROM}:${START}-${END}" | grep -v '^>')

# 计算GC含量
GC_CONTENT=$(echo "$SEQUENCE" | awk '
BEGIN {gc=0; total=0}
{
    for (i=1; i<=length($0); i++) {
        base=substr($0, i, 1)
        if (base == "G" || base == "C" || base == "g" || base == "c") {
            gc++
        }
        total++
    }
}
END {if (total > 0) print (gc/total)*100; else print "NA"}
')

# 输出结果
echo "Chromosome: $CHROM"
echo "Range: $START - $END"
echo "GC Content: $GC_CONTENT%"
echo "Sample SNP Density: $SAMPLE_SUM"
echo "Population Mean SNP Density: $OVERALL_MEAN"
echo "Population Median SNP Density: $OVERALL_MEDIAN"
echo "Population Q1 SNP Density: $OVERALL_Q1"
echo "Population Q3 SNP Density: $OVERALL_Q3"
echo "Difference with Population Mean: $DIFF_MEAN"
echo "Difference with Population Median: $DIFF_MEDIAN"
echo "Difference with Population Q1: $DIFF_Q1"
echo "Difference with Population Q3: $DIFF_Q3"
