#!/bin/bash

Rscript - <<EOF
library(DNAcopy)

# 尝试读取并处理数据，如果失败则停止
tryCatch({
  # 读取输入数据
  data <- read.table('/home/jialechen/projects/snpDensity/3_merged_data/snpden.txt', header=TRUE, sep='\t')

  # 确保所有数值列都是数值类型
  data[, 4:ncol(data)] <- lapply(data[, 4:ncol(data)], function(x) as.numeric(as.character(x)))

  # 提取必要的列
  chrom <- as.character(data$chrom)
  start <- as.numeric(data$start)
  values <- as.matrix(data[, 4:ncol(data)])

  # 创建CNA对象
  cna_obj <- CNA(values, chrom, start, data.type="binary")

  # 进行CBS分析（跳过平滑步骤）
  segment_cna <- segment(cna_obj, verbose=1)

  # 保存结果
  write.table(segment_cna$output, file='CBS_merge_output.txt', sep='\t', row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  # 成功完成分析后输出信息
  cat("CBS分析完成，结果保存在output.txt中\n")
}, error = function(e) {
  # 捕获并打印错误信息
  cat("发生错误：", conditionMessage(e), "\n")
  quit(status = 1)
})
EOF