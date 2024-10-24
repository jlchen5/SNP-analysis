#!/bin/bash

#################
##  GC含量计算  ##
#################


# 定义输入文件和基因组文件路径
GENE_POSITIONS_FILE="/home/jialechen/projects/snpDensity/3_merged_data/hg19.gene.sorted"
GENOME_FILE="/home/jialechen/library/hg19/hg19.fa"
OUTPUT_FILE="hg19_gc_content.txt"

# Python script to calculate GC content
PYTHON_SCRIPT=$(cat <<'EOF'
import sys
from Bio import SeqIO

def calculate_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100 if len(seq) > 0 else 0

def main(genome_file, gene_positions_file, output_file):
    # Load the genome
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    
    with open(gene_positions_file, 'r') as pos_file, open(output_file, 'w') as out_file:
        for line in pos_file:
            parts = line.strip().split()
            if len(parts) != 4:
                continue
            chrom, start, end, gene_name = parts
            start, end = int(start), int(end)
            
            if chrom in genome:
                seq = genome[chrom].seq[start:end]
                gc_content = calculate_gc_content(seq)
                out_file.write("{}\t{:.2f}\n".format(gene_name, gc_content))
            else:
                sys.stderr.write("Chromosome {} not found in genome file.\n".format(chrom))

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: python script.py <genome_file> <gene_positions_file> <output_file>\n")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    gene_positions_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(genome_file, gene_positions_file, output_file)
EOF
)

# Write the Python script to a temporary file
TEMP_PYTHON_SCRIPT=$(mktemp)
echo "$PYTHON_SCRIPT" > "$TEMP_PYTHON_SCRIPT"

# Run the Python script
python3 "$TEMP_PYTHON_SCRIPT" "$GENOME_FILE" "$GENE_POSITIONS_FILE" "$OUTPUT_FILE"

# Clean up temporary Python script
rm "$TEMP_PYTHON_SCRIPT"

