#!/usr/bin/env Rscript

library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.table("/home/jialechen/projects/snpDensity/3_merged_data/snpden.txt",header = T)
group <- read.csv("/home/jialechen/projects/snpDensity/4_result/group_info.tsv",sep = '\t',header = T)

# 从命令行获取输入参数
args <- commandArgs(trailingOnly = TRUE)
chrom <- args[1]
pos <- as.numeric(args[2])
dir <- args[3]


# 读取平均值数据
chr_avg_data <- read.csv("/home/jialechen/projects/snpDensity/3_merged_data/chromosome_SNP_stats.csv", header = TRUE,sep="\t")
chr_avg <- chr_avg_data$mean[chr_avg_data$chrom == chrom]  # 获取对应染色体的mean值

# 如果找不到对应染色体的平均值，则终止
if (length(chr_avg) == 0) {
  stop("没有找到指定染色体的平均值，请检查文件路径和染色体名称。")
}


# 定义绘制函数
plot_snp_density <- function(chrom, pos, data, group, chr_avg, output_dir = "output") {
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # 找到对应的SNP区间（假设每行是一个10kbp的区间）
  snp_data <- data[data$chrom == chrom, ]
  
  # 选定目标区域：位置上下游100k的范围
  start_pos <- max(0, pos - 110000)
  end_pos <- pos + 110000
  
  # 获取该范围内的数据
  snp_data_range <- snp_data[snp_data$start >= start_pos & snp_data$end <= end_pos, ]
  
  # 将数据转化为长格式，便于绘图
  snp_data_long <- snp_data_range %>%
    pivot_longer(cols = -c(chrom, start, end), names_to = "sample", values_to = "snp_density") %>%
    left_join(group, by = "sample")  # 将分组信息合并到数据中

  # 调试输出：检查合并后的数据
  print("合并后的数据：")
  print(head(snp_data_long))
  
  # 如果数据为空，则停止
  if (nrow(snp_data_long) == 0) {
    stop("没有找到与输入样本相匹配的数据，请检查输入的染色体和位置。")
  }

  # 查看group列的唯一值，确保有两个组
  print("分组信息：")
  print(unique(snp_data_long$group))
  
  # 计算样本的平均值
  sample_avg <- mean(snp_data_long$snp_density, na.rm = TRUE)  # 样本内的均值
  
  # 绘制小提琴图
  p <- ggplot(snp_data_long, aes(x = factor(start), y = snp_density, fill = group)) +
    geom_violin(trim = FALSE,scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", color = "black", 
                 outlier.shape = NA, alpha = 0.5) +  # 添加箱线图，显示四分位数线
    geom_point(position = position_jitter(width = 0.1), alpha = 0.8,size=0.1) +             
    labs(title = paste("SNP Density around", chrom, ":", pos), 
         x = "Position (10kb intervals)", 
         y = "SNP Density") +
    theme_bw() +
    geom_hline(aes(yintercept = chr_avg, color = "CHR AVG"), linetype = "dashed") +  # 添加横线阈值，并指定颜色
    geom_hline(aes(yintercept = sample_avg, color = "AVG of ±100kb"), linetype = "solid") +  # 样本内均值线
    scale_fill_manual(values = c("1KGp" = "blue", "bk" = "red")) +
    scale_color_manual(values = c("CHR AVG" = "black", "AVG of ±100kb" = "green")) +  # 设置横线颜色和图例标签
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_wrap(~ group, scales = "free_y", ncol = 1)  # 分组小提琴图
  
  # 保存图像
  output_file <- file.path(dir, paste(chrom, pos, "snp_density_violin.png", sep = "_")) # 默认输出路径为当前路径:"./"
  ggsave(output_file, p, width = 10, height = 6)
  message("图像已保存:", output_file)
}

# 调用绘图函数
plot_snp_density(chrom, pos, data, group, chr_avg)

