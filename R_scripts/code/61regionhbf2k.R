library(data.table)
library(MASS)
library(moments)
library(sampling)
args <- commandArgs(trailingOnly = TRUE)
disease <- args[1]
tissue <-  args[2]
matrix <- readRDS(paste0("2000gene2000p_",disease,"_",tissue,".rds"))
new_names <- paste0("cell", 1:dim(matrix)[[1]])
rownames(matrix) <- new_names
category1 <- matrix[,2001:4000]#p values
category2 <- matrix[,1:2000]#expression
Net <- readRDS(paste0('Net_',tissue,'.rds'))
rownames(Net) <- paste0("cell", 1:nrow(Net))
colnames(Net) <- paste0("cell", 1:ncol(Net))
cells <- rownames(matrix)

calculateCellWeights <- function(category1, category2) {
  #标准化-logp
  negative_log_category1 <- -log(category1)
  row_sums <- rowSums(negative_log_category1)
  normalized_negative_log_category1 <- t(t(negative_log_category1) / row_sums)
  score <- category2 * normalized_negative_log_category1
  # collect extra weight of candidate genes
  cat("Collecting extra weight of candidate cell...\n")
  # hotelling transform
  data <- score
  C <- cov(data)
  M <- colMeans(data)
  a.e <- eigen(C,symmetric = T)
  V <- a.e$vectors
  m <- round(0.6*ncol(category2))
  U <- V[1:m,]
  data_h <- U%*%(apply(data,1,function(x) x-M))
  data_h <- t(data_h)
  # box-cox transform
  for (i in 1:ncol(data_h)) {
    lam <- -floor(min(c(data_h[,i],0)))
    y <- data_h[,i] + lam
    if (all(y > 0)) {
      b <- boxcox(y ~ 1)
      lambda <- b$x[which.max(b$y)]
      if (lambda != 0) {
        data_h[, i] <- y^lambda
      } else {
        data_h[, i] <- log(y)
      }
    }
  }
  
  data_h <- scale(data_h)
  data_p <- apply(data_h,2,function(x) unlist(sapply(x,pnorm)))
  data_p <- data_p[, colSums(is.na(data_p)) == 0]
  extra_p <- data_p
  sub_cols <- 40
  
  num_groups <- 29
  
  sub_matrices <- lapply(1:num_groups, function(i) {
    start_col <- (i-1) * sub_cols + 1
    end_col <- i * sub_cols
    sub_matrix <- extra_p[, start_col:end_col]
    return(sub_matrix)
  })
  
  # 处理最后的子矩阵
  last_start_col <- num_groups * sub_cols + 1
  last_end_col <- ncol(extra_p)
  last_sub_matrix <- extra_p[, last_start_col:last_end_col]
  sub_matrices <- c(sub_matrices, list(last_sub_matrix))
  # 创建一个空的列表来存储p值
  p_values <- list()
  # 循环对split_matrices的每个元素求出p值
  for (i in 1:length(sub_matrices)) {
    p_values[[i]] <- apply(sub_matrices[[i]], 1, function(x) {
      pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
    })
  }
  # 为每个p值创建对应的变量名
  for (i in 1:length(p_values)) {
    assign(paste0("p", i), p_values[[i]])
  }
  p_matrix <- do.call(rbind, p_values)
  p_matrix <- t(p_matrix)
  # 找到非零最小值
  min_value <- min(p_matrix[p_matrix != 0])
  # 将零替换为非零最小值
  p_matrix[p_matrix == 0] <- min_value
  p_matrix[p_matrix < 1e-10] <- 1e-10
  extra_weight<-(-log(apply(p_matrix,1,function(x) { pchisq(-2*sum(log(x)),df=2*length(x),lower.tail=F)})))
  extra_weight <- data.frame(cell=cells,extra_weight=extra_weight)
  extra_weight <- extra_weight[match(unique(extra_weight$cell),extra_weight$cell),]
  #选择要标准化的列
  column_to_normalize <- extra_weight$extra_weight
  head(column_to_normalize)
  #最大值和第二大的值
  #max_value <- max(column_to_normalize)
  #second_max_value <- sort(column_to_normalize, decreasing = TRUE)[2]
  
  # 对列进行标准化，使得数据分布在0和1之间
  normalized_column <- (column_to_normalize - min(column_to_normalize)) / (max(column_to_normalize) - min(column_to_normalize))
  head(normalized_column)
  print(max(column_to_normalize))
  print(min(column_to_normalize))
  # 将标准化后的列赋值回DataFrame
  extra_weight$extra_weight <- normalized_column
  
  return(extra_weight)
}
extra_weight<-calculateCellWeights(category1,category2)
burnin_round <- 3000
after_burnin_round <- 3000
exclude_samegene <- T

####================= burn in step ======================####

#t0<-proc.time()
set.seed(1234)
thres <- 0.01; pickup <- 0;
num_riskcell <- 150
circle <- 1; chosen <- NULL
remaining <- sample(cells, num_riskcell)
num0 <- rep(0,length(cells))

dif <- thres + 1; dif_record <- NULL
while(dif > thres && circle<(burnin_round + 1))
{
  pickup <- pickup %% num_riskcell + 1
  if(pickup == 1)
    if(!is.null(chosen))
    {
      ###================================= calculate frequency =========================###
      num1 <- NULL
      num1 <- table(factor(chosen, levels=cells))
      num1 <- num1 + num0
      if(circle>1)
      {
        freq0 <- num0/(num_riskcell * (circle-1))
        freq1 <- num1/(num_riskcell * circle)
        dif <- (sum((freq0-freq1)^2))^0.5
        if( circle%%50==0 )
        {
          cat("Burnin sampling, sampling circle:",circle,"\n")
        }
        dif_record<-c(dif_record,dif)
      }
      num0<-num1; chosen<-NULL; circle<-circle+1
    }
  pickup_p<-Net[,is.element(colnames(Net),remaining[-pickup])]
  pickup_p<-pickup_p[!is.element(rownames(pickup_p),colnames(pickup_p)),]
  pickup_p <- apply(pickup_p,1,sum)
  pickup_p <- extra_weight[match(names(pickup_p),extra_weight$cell),]$extra_weight * pickup_p 

  if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i] <- 1/length(pickup_p)

  remaining[pickup] <- sample(names(pickup_p),1,replace=T,prob=pickup_p)
  chosen<-rbind(chosen,remaining)
}
###===================== end of burn in step ===============================###

###======================= post-burn in step ===================================###

pickup <- 0; num_riskcell <- 150; circle<-1; chosen <- NULL
num0 <- rep(0,length(cells))

thres <- 0.005; dif <- thres + 1

while(dif>thres && circle<(after_burnin_round + 1) )
{
  pickup <- pickup %% num_riskcell + 1
  if(pickup == 1)
    if(!is.null(chosen))
    {
      ###================================= calculate frequency =========================###
      num1 <- NULL
      num1 <- table(factor(chosen, levels=cells))
      num1 <- num1 + num0
      if(circle>1)
      {
        freq0 <- num0/(num_riskcell * (circle-1))
        freq1 <- num1/(num_riskcell * circle)
        dif <- (sum((freq0-freq1)^2))^0.5
        if( circle%%50==0 )
        {
          cat("PostBurnin sampling, sampling circle:",circle,"\n")
        }
        dif_record<-c(dif_record,dif)
      }
      num0<-num1; chosen<-NULL; circle<-circle+1
    }
      ###============================= end of calculating frequency =======================###

  pickup_p<-Net[,is.element(colnames(Net),remaining[-pickup])]
  pickup_p<-pickup_p[!is.element(rownames(pickup_p),colnames(pickup_p)),]
  pickup_p <- apply(pickup_p,1,sum)
  pickup_p <- extra_weight[match(names(pickup_p),extra_weight$cell),]$extra_weight * pickup_p 
  
  if(sum(pickup_p)==0) for(i in 1:length(pickup_p)) pickup_p[i] <- 1/length(pickup_p)
  
  remaining[pickup] <- sample(names(pickup_p),1,replace=T,prob=pickup_p)
  chosen<-rbind(chosen,remaining)
}

###=====================  end of post-burnin step  =======================================###

###===================== summarize and record the results ================================###
freq<-cbind(cells,freq1)

colnames(freq) <- c("cell","post_prob")
freq <- as.data.frame(freq,stringsAsFactors=F)
freq[,2] <- as.numeric(freq[,2])

set.seed(123)  # 设置随机数种子，以确保结果可复现

#计算Pvalue加到dataframe中
add_p_norm <- function(df) {
  if (!all(c("post_prob") %in% names(df))) {
    stop("Input resframe must contain a column named 'post_prob'")
  }
  # 生成随机样本
  random_sample <- rnorm(n = length(df$post_prob), 
                         mean = mean(df$post_prob), 
                         sd = sd(df$post_prob))
  # 计算经验P值
  empirical_p_values <- sapply(df$post_prob, function(p) {
    (sum(random_sample >= p) + 1) / length(random_sample)
  })
  df$p_norm <- empirical_p_values
  return(df)
}
freq <- add_p_norm(freq)
saveRDS(freq,paste0('freq61_hbf2000_',disease,'_',tissue,'.rds'))
