#!/usr/bin/env Rscript
library(SeuratObject)

# 从命令行参数获取输入和输出文件名
args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) < 2) {
  stop("用法: Rscript script_name.R <example.rds> <freq61_hbf2000_mono_example.rds> [output.csv]")
}

# 设置输入文件路径
example_file <- args[1]
freq_file <- args[2]

# 设置输出文件路径（可选第三个参数）
if (length(args) >= 3) {
  output_file <- args[3]
} else {
  # 如果没有提供输出文件名，使用第二个输入文件名并替换扩展名为.csv
  # 使用basename获取文件名（不包含路径）
  base_name <- tools::file_path_sans_ext(basename(freq_file))
  output_file <- paste0(base_name, ".csv")
  
  # 如果要保持完整路径，只替换扩展名：
  # output_file <- sub("\\.rds$", ".csv", freq_file, ignore.case = TRUE)
}

cat("正在读取文件...\n")
cat("读取文件:", example_file, "\n")
data <- readRDS(example_file)

cat("读取文件:", freq_file, "\n")
data_p <- readRDS(freq_file)

# 添加细胞类型信息
cat("正在添加细胞类型信息...\n")
data_p$type <- data$cell_type

# 保存为CSV文件
cat("正在保存到:", output_file, "\n")
write.csv(data_p, file = output_file, row.names = FALSE)

cat("处理完成！\n")
cat("输出文件:", output_file, "\n")

#Rscript freq61_type_transfer.R example.rds freq61_hbf2000_mono_example.rds