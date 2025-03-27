args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]

data<-readRDS(paste0(tissue,".rds"))
mat<-data[['RNA']]$counts
mat<-as.matrix(mat)

cor_matrix <- cor(mat)
cor_matrix <- abs(cor_matrix)
column_sums <- colSums(cor_matrix)

normalized_matrix <- cor_matrix / column_sums
saveRDS(normalized_matrix,paste0("Net_",tissue,".rds"))
