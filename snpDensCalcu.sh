#!/bin/bash

# 读取样本文件并赋值给数组
samples=($(cat sample.txt))

# 循环处理每个样本
for sample in "${samples[@]}"; do
    # 获取样本目录中的所有 VCF 文件
    vcf_files=(${sample}/*.vcf)

    
    # 循环处理每个 VCF 文件
    for vcf_file in "${vcf_files[@]}"; do
        # 提取染色体信息
        chr=$(basename "$vcf_file" | grep -oP 'chr([0-9]+|X|Y)')
        
        # 定义输出文件名
        output="${sample}/${sample}_${chr}_Het"
        
        # 运行 vcftools
        vcftools --vcf "$vcf_file" --SNPdensity 10000 --out "$output"
    done
    echo "Finished processing sample:" $sample
done

echo " All samples processed!"
