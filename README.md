# SNP Analysis Workflow

## File Processing

### 1. Filtering VCF Files

#### Homozygous Variants
```bash
zcat U240707021F0.vcf.gz | awk '$10 ~ /^0\/0/ || $10 ~ /^1\/1/ || $10 ~ /^2\/2/' > U240707021F0_Hom.vcf
```

#### Heterozygous Variants
```bash
zcat U240707021F0.vcf.gz | awk '$10 ~ /^0\/1/ || $10 ~ /^1\/2/' > U240707021F0_Het.vcf
```

### 2. Calculating SNP Density
```bash
vcftools --vcf U240707021F0_Het.vcf \
         --SNPdensity 100000 \
         --out U240707021F0_Het_Dens
```

### 3. Extracting Zero-Density Regions
```bash
awk '{if ($3==0) print $0}' U240724001F0_chr10_Het.snpden
```

## Data Organization

### Sample Processing Steps
```bash
# Copy files
cp E150035532/result/*/02.align/*/alignsplit/*hetsnp.vcf ./

# Create sample directories
ls *vcf | sed 's/.chr/\t/' | awk '{print $1}' | sort -u | while read line; do mkdir -p $line; done

# Create sample list
ls *vcf | sed 's/.chr/\t/' | awk '{print $1}' | sort -u > sample.txt

# Move files to sample directories
cat sample.txt | while read line; do mv $line.chr*vcf $line/; done

# Calculate SNP density
sh ~/scripts/snpDensCalcu.sh
```

### 4. Adding Sample Names to SNP IDs
```bash
cat sample.txt | while read sample; do
  awk -v sample="$sample" 'BEGIN{FS=OFS="\t"} 
  NR==1 {print $0, "SNP_NAME"} 
  NR>1 {print $0, sample"_SNP_" NR-1}' "${sample}_ALL_Het.snpden" > "${sample}_ALL_Het_2.snpden"
done
```

## Genomic Analysis

### 5. Genome Binning
```bash
bedtools makewindows -g ~/library/hg19/hg19.chrom.sizes_24 \
                     -w 10000 \
                     > hg19_10k.bed
```

### 6. Adding Gene Annotations
```bash
# Create annotated bed file
bedtools intersect -a hg19_gene.bed \
                   -b hg19_10k.bed \
                   -wb | awk '{print $9"\t"$10"\t"$11"\t"$8}' | sed 's/"//g' > hg19_10k_anno.bed

# Sort and map annotations
cat hg19_10k_anno.bed | sort -k1,1 -k2,2n > hg19_10k_anno_sorted.bed

bedtools map -a hg19_10k_sorted.bed \ 
             -b hg19_10k_anno_sorted.bed \ 
             -c 4 \ 
             -o distinct \ 
             > hg19_10k_with_genes.bed
```

### 7. Chromosome Bin Counts
```bash
for chr in {1..22} X Y; do
    echo -n "chr$chr: "
    awk -v chr="chr$chr" '$1 == chr' merged_snpden_anno.txt | wc -l
done
```

**Results:**
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

## Additional Annotations

### 8. Gap and Duplication Regions

#### Calculate overlaps:
```bash
# Gap regions
bedtools intersect -a hg19_10k_with_genes.bed -b hg19_gap.bed -c > hg19_10k_with_genes_with_gap.bed

# Duplication regions
bedtools intersect -a hg19_10k_with_genes.bed -b hg19_seqdup.bed -c > hg19_10k_with_genes_with_dup.bed
```

#### Merge results:
```bash
paste hg19_10k_with_genes_with_gap.bed hg19_10k_with_genes_with_dup.bed | \
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, ($5 > 0 ? "T" : "F"), ($10 > 0 ? "T" : "F")}' > hg19_10k_with_genes_with_gap_dup.bed
```

#### Add headers:
```bash
sed -i '1s/^/chrom\tstart\tend\tgene\tIF_GAP\tIF_DUP\n/' hg19_10k_with_genes_with_gap_dup.bed
```

**Example output:**
```
chrom	start	end	gene	IF_GAP	IF_DUP
chr1	0	10000	.	T	F
chr1	10000	20000	DDX11L1,MIR6859-1,WASH7P	F	T
chr1	20000	30000	WASH7P	F	T
...
```

## 1000 Genomes Project Data

### 9. Sample Processing

#### Extract sample names:
```bash
bcftools query -l ALL.chr13.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > sample.txt
```

#### Split multi-sample VCFs:
```bash
cat sample.txt | while read id; do
    mkdir -p $id
    for vcf in ALL.chr*.vcf.gz; do
        bcftools view -s $id -o $id/${id}_${vcf%.vcf.gz}.vcf.gz -O v $vcf
    done
done
```

## Visualization

### 10. SNP Density Plotting

Use the [plot_violin.R](https://github.com/jlchen5/SNP-analysis/blob/main/plot_violin.R) script to visualize SNP density in 100kb upstream/downstream regions.

## Supplementary Data

Additional required files available at:  
https://github.com/jlchen5/SNP-analysis/blob/main/data/1
