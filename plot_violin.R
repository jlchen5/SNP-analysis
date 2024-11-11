library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.table("snpden.txt",header = T)
group <- read.csv("group_info.tsv",sep = '\t',header = T)

chrom <- "chr17" 
pos <- 29653271 
chr_avg <- 6.86


# 定义绘制函数
plot_snp_density <- function(chrom, pos, data, group, output_dir = "output") {
  # 确保输出目录存在
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  # 找到对应的SNP区间（假设每行是一个10kbp的区间）
  snp_data <- data[data$chrom == chrom, ]
  
  # 选定目标区域：位置上下游100k的范围
  start_pos <- max(0, pos - 100000)
  end_pos <- pos + 100000
  
  # 获取该范围内的数据
  snp_data_range <- snp_data[snp_data$start >= start_pos & snp_data$end <= end_pos, ]
  
  # 将数据转化为长格式，便于绘图
  snp_data_long <- snp_data_range %>%
    pivot_longer(cols =-c(chrom, start, end), names_to = "sample", values_to = "snp_density") %>%
    left_join(group, by = "sample")  # 将分组信息合并到数据中
  
  
  
  # 调试输出：检查合并后的数据
  print("合并后的数据：")
  print(head(snp_data_long))
  
  # 检查是否成功合并分组信息
  if (nrow(snp_data_long) == 0) {
    stop("没有找到与输入样本相匹配的数据，请检查输入的染色体和位置。")
  }
  
  # 查看group列的唯一值，确保有两个组
  print("分组信息：")
  print(unique(snp_data_long$group))
  
  sample_avg <- mean(snp_data_long$snp_density, na.rm = TRUE)  # 样本内的均值
  
  
  # 绘制小提琴图
  p <- ggplot(snp_data_long, aes(x = factor(start), y = snp_density, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", color = "black", 
                 outlier.shape = NA, alpha = 0.5) +  # 添加箱线图，显示四分位数线
    labs(title = paste("SNP Density around", chrom, ":", pos), 
         x = "Position (10kb intervals)", 
         y = "SNP Density") +
    theme_bw() +
    geom_hline(yintercept = chr_avg, color = "black", linetype = "dashed") +  # 添加横线阈值
    geom_hline(yintercept = sample_avg, color = "green", linetype = "solid") +  # 样本内均值线
    scale_fill_manual(values = c("1KGp" = "blue", "bk" = "red")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #scale_y_continuous(limits = c(0, NA)) +  # 设置y轴的最小值为0，最大值由数据自动决定
    facet_wrap(~ group, scales = "free_y", ncol = 1)  # 分组小提琴图
  
  # 保存图像
  output_file <- file.path("./", paste(chrom, pos, "snp_density_violin.png", sep = "_"))
  ggsave(output_file, p, width = 10, height = 6)
  message("图像已保存:", output_file)
}


# 调用绘图函数
plot_snp_density(chrom, pos, data, group)

