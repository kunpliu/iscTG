#!/usr/bin/env Rscript

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

# 检查是否有足够的参数
if (length(args) < 2) {
  stop("没有提供必要的标识符参数。用法：Rscript script.R <identifier>", call. = FALSE)
}

# 获取输入标识符
input_identifier <- args[1]
tissue <- args[2]
tissue_file <- paste0(tissue, '.rds')

SeuratSymbol <- readRDS(tissue_file)
data <- SeuratSymbol[['RNA']]$data
# 生成新的列名
new_col_names <- paste0("cell", 1:dim(data)[[2]])
# 替换列名
colnames(data) <- new_col_names

# 加载必要的库
library(Seurat)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(homologene)

# 使用input_identifier构造文件名
out_file_name <- paste0(input_identifier,".genes.out")
out <- fread(out_file_name)

# out$GENE from NCBI id to Symbol
res <- homologene(out$GENE, inTax=9606, outTax=9606)
common_genes <- intersect(res[["9606"]], rownames(data))
idx1 <- match(common_genes, res[["9606"]])
res_sub <- res[idx1, ]
idx2 <- match(res_sub[["9606_ID"]], out$GENE)
out_sub <- out[idx2, ]

# 使用match()函数获取sub_out$gene在sub_res$9606中的索引
idx <- match(out_sub$GENE, res_sub$'9606_ID')
# 将sub_res$9606_ID的对应值添加为新列
out_sub$SYMBOL <- res_sub$'9606'[idx]

# 保存结果，使用input_identifier作为文件名的一部分
magma_symbol_file_name <- paste0('MAGMA_SYMBOL_', input_identifier,'_',tissue ,'.rds')
saveRDS(out_sub, magma_symbol_file_name)

selected_rows <- data[out_sub$'SYMBOL', ]
data <- selected_rows

# 保存结果，使用input_identifier作为文件名的一部分
all_cell_exp_file_name <- paste0('AllCellexp_', input_identifier, '_',tissue ,'.rds')
saveRDS(data, all_cell_exp_file_name)

###############
data <- readRDS(all_cell_exp_file_name)
magma <- readRDS(magma_symbol_file_name)
magma_sorted <- magma[order(magma$P), ]
magma_sorted <- magma_sorted[1:2000,]
data <- as.matrix(data)
data <- data[magma_sorted$SYMBOL,]
data <- t(data)

gene <- colnames(data)
P <- magma_sorted$P[magma_sorted$SYMBOL %in% gene]

empty_cols <- matrix(nrow = nrow(data), ncol = 2000)
data <- cbind(data, empty_cols)

for (i in 1:length(P)) {
  col_name <- paste0("P_", i)  # 构造列名，如 P_1, P_2, ...
  data[, (2000 + i)] <- P[i]  # 将元素填充到对应的列
  colnames(data)[(2000 + i)] <- col_name  # 设置列名
}

# 保存结果，使用input_identifier作为文件名的一部分
gene_p_file_name <- paste0("2000gene2000p_", input_identifier, '_',tissue ,'.rds')
saveRDS(data, gene_p_file_name)
