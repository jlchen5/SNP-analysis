# SNP-analysis

## 1 筛选vcf文件中纯合（homozygous）位点：

```bash
$ zcat U240707021F0.vcf.gz  |awk '$10 ~ /^0\/0/ || $10 ~ /^1\/1/ || $10 ~ /^2\/2/' > U240707021F0_Hom.vcf
```

## 2 筛选vcf文件中杂合（heterozygous）位点：

```bash
$ zcat U240707021F0.vcf.gz  |awk '$10 ~ /^0\/1/ || $10 ~ /^1\/2/ ' > U240707021F0_Het.vcf
```

## 3 计算vcf中SNP的密度

```
vcftools --vcf U240707021F0_Het.vcf  --SNPdensity 100000 --out U240707021F0_Het_Dens
```

## 4 输出SNP密度为0的结果

```bash
awk '{if ($3==0) print $0}' U240724001F0_chr10_Het.snpden
```

## 数据

- vcf 文件：`*hetsnp.vcf`

```
E150036213
E150028967
E150035507
E150035532
```

- steps

~~~
cp E150035532/result/*/02.align/*/alignsplit/*hetsnp.vcf ./

ls *vcf |sed 's/.chr/\t/'  |awk '{print $1}'| sort -u |while read line ;do mkdir -p $line ;done 

ls *vcf |sed 's/.chr/\t/'  |awk '{print $1}'| sort -u  > sample.txt

cat sample.txt  |while read line ;do mv  $line.chr*vcf $line/ ;done 

sh ~/scripts/snpDensCalcu.sh 
~~~


## 给每个样本的snp加名字 `sample_snp_NR`
```
cat sample.txt | while read sample; do
  awk -v sample="$sample" 'BEGIN{FS=OFS="\t"} 
  NR==1 {print $0, "SNP_NAME"} 
  NR>1 {print $0, sample"_SNP_" NR-1}' "${sample}_ALL_Het.snpden" > "${sample}_ALL_Het_2.snpden"
done
```
## 基因组拆分

```
bedtools makewindows -g ~/library/hg19/hg19.chrom.sizes_24  -w 10000  > hg19_10k.bed
```

## 添加注释信息

```
bedtools intersect -a hg19_gene.bed -b hg19_10k.bed  -wb  |awk '{print $9"\t"$10"\t"$11"\t"$8}'  |sed 's/"//g'   > hg19_10k_anno.bed

cat hg19_10k_anno.bed |sort -k1,1 -k2,2n > hg19_10k_anno_sorted.bed

bedtools map -a hg19_10k_sorted.bed -b hg19_10k_anno_sorted.bed -c 4 -o distinct > hg19_10k_with_genes.bed
```

## 计算每条染色体bins数目
```
for chr in {1..22} X Y; do
    echo -n "chr$chr: "
    awk -v chr="chr$chr" '$1 == chr' merged_snpden_anno.txt | wc -l
done
```

```
chr1: 128750
chr2: 40567
chr3: 19803
chr4: 19116
chr5: 18092
chr6: 17112
chr7: 15914
chr8: 14637
chr9: 14122
chr10: 13554
chr11: 13501
chr12: 13386
chr13: 11517
chr14: 10735
chr15: 10254
chr16: 9036
chr17: 8120
chr18: 7808
chr19: 5913
chr20: 6303
chr21: 4813
chr22: 5131
chrX: 15528
chrY: 5938

```





# 新增是否重复区，是否N区

### 计算与 hg19_gap.bed 的交集：

```
bedtools intersect -a hg19_10k_with_genes.bed -b hg19_gap.bed -c > hg19_10k_with_genes_with_gap.bed
```

### 计算与 hg19_seqdup.bed 的交集：

```
bedtools intersect -a hg19_10k_with_genes.bed -b hg19_seqdup.bed -c > hg19_10k_with_genes_with_dup.bed
```

### 合并交集结果并添加 gap 和 dup 列：

```
paste hg19_10k_with_genes_with_gap.bed hg19_10k_with_genes_with_dup.bed | \
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, ($5 > 0 ? "T" : "F"), ($10 > 0 ? "T" : "F")}' > hg19_10k_with_genes_with_gap_dup.bed
```

### 清理临时文件（可选）：

```
rm hg19_10k_with_genes_with_gap.bed hg19_10k_with_genes_with_dup.bed
```

### 添加标题

```bash
sed -i '1s/^/chrom\tstart\tend\tgene\tIF_GAP\tIF_DUP\n/' hg19_10k_with_genes_with_gap_dup.bed 

(base) jialechen ~/library/hg19 $ head hg19_10k_with_genes_with_gap_dup.bed 
chrom	start	end	gene	IF_GAP	IF_DUP
chr1	0	10000	.	T	F
chr1	10000	20000	DDX11L1,MIR6859-1,WASH7P	F	T
chr1	20000	30000	WASH7P	F	T
chr1	30000	40000	FAM138A,MIR1302-2	F	T
chr1	40000	50000	.	F	T
chr1	50000	60000	OR4G4P	F	T
chr1	60000	70000	OR4F5,OR4G11P	F	T
chr1	70000	80000	OR4F5	F	T
chr1	80000	90000	.	F	T
```



# 新增千人基因组样本信息

### 统计样本名
```bash
 bcftools query -l  ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz >      sample.txt
```

### 拆分多样本 VCF 文件 (/home/jialechen/projects/snpDensity/kiloGP/parallel_split_vcf.sh)

```bash
cat sample.txt | while read id; do
    mkdir -p $id
    for vcf in ALL.chr*.vcf.gz; do
        bcftools view -s $id -o $id/${id}_${vcf%.vcf.gz}.vcf.gz -O v $vcf
    done
done
```

# 绘制SNP密度结果
使用[plot_violin.R](https://github.com/jlchen5/SNP-analysis/blob/main/plot_violin.R)脚本，可以绘制位点上下游100kb的SNP密度。









