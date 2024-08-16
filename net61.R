args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
data<-readRDS(paste0(tissue,".rds"))
mat<-data[['RNA']]$counts
mat<-as.matrix(mat)
# 直接计算相关系数矩
#print(class(mat))
cor_matrix <- cor(mat)
cor_matrix <- abs(cor_matrix)
column_sums <- colSums(cor_matrix)

# 对矩阵进行标准化
normalized_matrix <- cor_matrix / column_sums
saveRDS(normalized_matrix,paste0("Net_",tissue,".rds"))
