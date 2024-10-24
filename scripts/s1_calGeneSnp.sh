#!/bin/bash

##################
##  SNP密度统计  ##
##################

# 输入文件
GENE_FILE="/home/jialechen/projects/snpDensity/3_merged_data/hg19.gene.sorted"
SNP_FILE="/home/jialechen/projects/snpDensity/3_merged_data/snpden.txt"

# 输出文件前缀
OUTPUT_PREFIX="snp_density_stats"

# 生成基因本身的区域文件
cp $GENE_FILE gene_0kb.bed

# 扩展基因区域并保存到临时文件，确保起始位置不小于0
awk 'BEGIN{OFS="\t"} {start = $2-50000; if (start < 0) start = 0; print $1, start, $3+50000, $4}' $GENE_FILE > gene_50kb.bed
awk 'BEGIN{OFS="\t"} {start = $2-100000; if (start < 0) start = 0; print $1, start, $3+100000, $4}' $GENE_FILE > gene_100kb.bed
awk 'BEGIN{OFS="\t"} {start = $2-200000; if (start < 0) start = 0; print $1, start, $3+200000, $4}' $GENE_FILE > gene_200kb.bed

# 函数：计算统计数据
calculate_stats() {
  local input_bed=$1
  local output_file=$2
  local total_samples=349  # 总样本数

  echo -e "Gene\tRegion\tMedian\tMean\tQ1\tQ3" > $output_file

  while read -r line; do
    gene=$(echo $line | awk '{print $4}')
    region=$(echo $line | awk '{print $1":"$2"-"$3}')
    snp_values=$(bedtools intersect -a <(echo "$line" | awk '{print $1"\t"$2"\t"$3"\t"$4}') -b $SNP_FILE -wa -wb | awk '{for(i=4;i<=NF;i++) print $i}' | grep -E '^[0-9]+(\.[0-9]+)?')

    if [ -z "$snp_values" ]; then
      echo -e "$gene\t$region\tNA\tNA\tNA\tNA" >> $output_file
    else
      median=$(echo "$snp_values" | awk '{a[NR]=$1} END{if (NR%2==1) {print a[int(NR/2)+1]} else {print (a[int(NR/2)]+a[int(NR/2)+1])/2}}')
      mean=$(echo "$snp_values" | awk -v total_samples=$total_samples '{sum+=$1} END{print sum/total_samples}')
      q1=$(echo "$snp_values" | awk '{a[NR]=$1} END{print a[int(NR/4)+1]}')
      q3=$(echo "$snp_values" | awk '{a[NR]=$1} END{print a[int(3*NR/4)+1]}')
      echo -e "$gene\t$region\t$median\t$mean\t$q1\t$q3" >> $output_file
    fi
  done < $input_bed
}

export -f calculate_stats
export SNP_FILE

# 计算基因本身区域的统计数据
echo "：：：：：正在统计基因区SNP密度"
calculate_stats gene_0kb.bed "${OUTPUT_PREFIX}_0kb.txt"
echo "：：：：：基因区SNP密度计算完成"

# # 计算50Kb上下游区域的统计数据
# echo "：：：：：正在统计基因区上下游50Kb区域SNP密度"
# calculate_stats gene_50kb.bed "${OUTPUT_PREFIX}_50kb.txt"
# echo "：：：：：基因区上下游50Kb区域SNP密度计算完成"

# # 计算100Kb上下游区域的统计数据
# echo "：：：：：正在统计基因区上下游100Kb区域SNP密度"
# calculate_stats gene_100kb.bed "${OUTPUT_PREFIX}_100kb.txt"
# echo "：：：：：基因区上下游100Kb区域SNP密度计算完成"

# # 计算200Kb上下游区域的统计数据
# echo "：：：：：正在统计基因区上下游200Kb区域SNP密度"
# calculate_stats gene_200kb.bed "${OUTPUT_PREFIX}_200kb.txt"
# echo "：：：：：基因区上下游200Kb区域SNP密度计算完成"


# # 清理临时文件
# echo "：：：：：正在清理临时文件"
# rm gene_0kb.bed gene_50kb.bed gene_100kb.bed gene_200kb.bed

echo "SNP密度统计已完成，结果保存在${OUTPUT_PREFIX}_0kb.txt"
