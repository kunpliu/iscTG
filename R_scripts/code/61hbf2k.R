library(data.table)
library(MASS)
library(moments)
library(sampling)

run_freq61_hbf2000 <- function(disease,
                               tissue,
                               matrix_prefix = "2000gene2000p_",
                               net_prefix = "Net_",
                               burnin_round = 3000,
                               after_burnin_round = 3000,
                               exclude_samegene = TRUE,
                               num_riskcell = 150,
                               burnin_thres = 0.01,
                               post_thres = 0.005,
                               sub_cols = 40,
                               num_groups = 29,
                               eigen_keep_ratio = 0.6,
                               p_floor = 1e-10,
                               out_prefix = "freq61_hbf2000_",
                               save_dir = ".") {

  matrix <- readRDS(paste0(matrix_prefix, disease, "_", tissue, ".rds"))
  new_names <- paste0("cell", 1:dim(matrix)[[1]])
  rownames(matrix) <- new_names
  category1 <- matrix[, 2001:4000] #p values
  category2 <- matrix[, 1:2000]    #expression

  Net <- readRDS(paste0(net_prefix, tissue, ".rds"))
  rownames(Net) <- paste0("cell", 1:nrow(Net))
  colnames(Net) <- paste0("cell", 1:ncol(Net))

  cells <- rownames(matrix)

  calculateCellWeights <- function(category1, category2) {
    negative_log_category1 <- -log(category1)
    row_sums <- rowSums(negative_log_category1)
    normalized_negative_log_category1 <- t(t(negative_log_category1) / row_sums)
    score <- category2 * normalized_negative_log_category1

    cat("Collecting extra weight of candidate cell...\n")

    data <- score
    C <- cov(data)
    M <- colMeans(data)
    a.e <- eigen(C, symmetric = TRUE)
    V <- a.e$vectors

    m <- round(eigen_keep_ratio * ncol(category2))
    U <- V[1:m, ]

    data_h <- U %*% (apply(data, 1, function(x) x - M))
    data_h <- t(data_h)

    for (i in 1:ncol(data_h)) {
      lam <- -floor(min(c(data_h[, i], 0)))
      y <- data_h[, i] + lam
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
    data_p <- apply(data_h, 2, function(x) unlist(sapply(x, pnorm)))
    data_p <- data_p[, colSums(is.na(data_p)) == 0]
    extra_p <- data_p

    sub_matrices <- lapply(1:num_groups, function(i) {
      start_col <- (i - 1) * sub_cols + 1
      end_col <- i * sub_cols
      sub_matrix <- extra_p[, start_col:end_col]
      return(sub_matrix)
    })

    last_start_col <- num_groups * sub_cols + 1
    last_end_col <- ncol(extra_p)
    last_sub_matrix <- extra_p[, last_start_col:last_end_col]
    sub_matrices <- c(sub_matrices, list(last_sub_matrix))

    p_values <- list()
    for (i in 1:length(sub_matrices)) {
      p_values[[i]] <- apply(sub_matrices[[i]], 1, function(x) {
        pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
      })
    }

    for (i in 1:length(p_values)) {
      assign(paste0("p", i), p_values[[i]])
    }

    p_matrix <- do.call(rbind, p_values)
    p_matrix <- t(p_matrix)

    min_value <- min(p_matrix[p_matrix != 0])
    p_matrix[p_matrix == 0] <- min_value
    p_matrix[p_matrix < p_floor] <- p_floor

    extra_weight <- (-log(apply(p_matrix, 1, function(x) {
      pchisq(-2 * sum(log(x)), df = 2 * length(x), lower.tail = FALSE)
    })))

    extra_weight <- data.frame(cell = cells, extra_weight = extra_weight)
    extra_weight <- extra_weight[match(unique(extra_weight$cell), extra_weight$cell), ]

    column_to_normalize <- extra_weight$extra_weight
    head(column_to_normalize)

    normalized_column <- (column_to_normalize - min(column_to_normalize)) /
      (max(column_to_normalize) - min(column_to_normalize))

    head(normalized_column)
    print(max(column_to_normalize))
    print(min(column_to_normalize))

    extra_weight$extra_weight <- normalized_column
    return(extra_weight)
  }

  extra_weight <- calculateCellWeights(category1, category2)

  ####================= burn in step ======================####
  set.seed(1234)
  thres <- burnin_thres
  pickup <- 0

  circle <- 1
  chosen <- NULL
  remaining <- sample(cells, num_riskcell)
  num0 <- rep(0, length(cells))

  dif <- thres + 1
  dif_record <- NULL

  while (dif > thres && circle < (burnin_round + 1)) {
    pickup <- pickup %% num_riskcell + 1
    if (pickup == 1)
      if (!is.null(chosen)) {
        num1 <- NULL
        num1 <- table(factor(chosen, levels = cells))
        num1 <- num1 + num0
        if (circle > 1) {
          freq0 <- num0 / (num_riskcell * (circle - 1))
          freq1 <- num1 / (num_riskcell * circle)
          dif <- (sum((freq0 - freq1)^2))^0.5
          if (circle %% 50 == 0) {
            cat("Burnin sampling, sampling circle:", circle, "\n")
          }
          dif_record <- c(dif_record, dif)
        }
        num0 <- num1
        chosen <- NULL
        circle <- circle + 1
      }

    pickup_p <- Net[, is.element(colnames(Net), remaining[-pickup])]
    pickup_p <- pickup_p[!is.element(rownames(pickup_p), colnames(pickup_p)), ]
    pickup_p <- apply(pickup_p, 1, sum)
    pickup_p <- extra_weight[match(names(pickup_p), extra_weight$cell), ]$extra_weight * pickup_p

    if (sum(pickup_p) == 0) for (i in 1:length(pickup_p)) pickup_p[i] <- 1 / length(pickup_p)

    remaining[pickup] <- sample(names(pickup_p), 1, replace = TRUE, prob = pickup_p)
    chosen <- rbind(chosen, remaining)
  }
  ###===================== end of burn in step ===============================###

  ###======================= post-burn in step ===================================###
  pickup <- 0
  circle <- 1
  chosen <- NULL
  num0 <- rep(0, length(cells))

  thres <- post_thres
  dif <- thres + 1

  while (dif > thres && circle < (after_burnin_round + 1)) {
    pickup <- pickup %% num_riskcell + 1
    if (pickup == 1)
      if (!is.null(chosen)) {
        num1 <- NULL
        num1 <- table(factor(chosen, levels = cells))
        num1 <- num1 + num0
        if (circle > 1) {
          freq0 <- num0 / (num_riskcell * (circle - 1))
          freq1 <- num1 / (num_riskcell * circle)
          dif <- (sum((freq0 - freq1)^2))^0.5
          if (circle %% 50 == 0) {
            cat("PostBurnin sampling, sampling circle:", circle, "\n")
          }
          dif_record <- c(dif_record, dif)
        }
        num0 <- num1
        chosen <- NULL
        circle <- circle + 1
      }

    pickup_p <- Net[, is.element(colnames(Net), remaining[-pickup])]
    pickup_p <- pickup_p[!is.element(rownames(pickup_p), colnames(pickup_p)), ]
    pickup_p <- apply(pickup_p, 1, sum)
    pickup_p <- extra_weight[match(names(pickup_p), extra_weight$cell), ]$extra_weight * pickup_p

    if (sum(pickup_p) == 0) for (i in 1:length(pickup_p)) pickup_p[i] <- 1 / length(pickup_p)

    remaining[pickup] <- sample(names(pickup_p), 1, replace = TRUE, prob = pickup_p)
    chosen <- rbind(chosen, remaining)
  }
  ###=====================  end of post-burnin step  =======================================###

  ###===================== summarize and record the results ================================###
  freq <- cbind(cells, freq1)
  colnames(freq) <- c("cell", "post_prob")
  freq <- as.data.frame(freq, stringsAsFactors = FALSE)
  freq[, 2] <- as.numeric(freq[, 2])

  set.seed(123)


  add_p_gamma <- function(df) {
    if (!all(c("post_prob") %in% names(df))) {
      stop("Input resframe must contain a column named 'post_prob'")
    }  
    prob_values <- df$post_prob
    nonzero_values <- prob_values[prob_values > 0]
    stopifnot(length(nonzero_values) >= 10)  
    fit_gamma <- fitdist(nonzero_values, "gamma")
    shape_hat <- fit_gamma$estimate["shape"]
    rate_hat  <- fit_gamma$estimate["rate"]  
    random_sample <- rgamma(length(prob_values), shape = shape_hat, rate = rate_hat)  
    df$p_gamma <- sapply(prob_values, function(p) {
      (sum(random_sample >= p) + 1) / (length(random_sample) + 1)
    })  
    return(df)
  }


  freq <- add_p_gamma(freq)

  out_path <- file.path(save_dir, paste0(out_prefix, disease, "_", tissue, ".rds"))
  saveRDS(freq, out_path)

  return(freq)
}


run_freq61_hbf2000('mono', 'example')
