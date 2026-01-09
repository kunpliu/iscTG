#!/usr/bin/env Rscript
library(SeuratObject)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("用法: Rscript script_name.R <example.rds> <freq61_hbf2000_mono_example.rds> [output.csv]")
}

example_file <- args[1]
freq_file <- args[2]

if (length(args) >= 3) {
  output_file <- args[3]
} else {
  base_name <- tools::file_path_sans_ext(basename(freq_file))
  output_file <- paste0(base_name, ".csv")
  # output_file <- sub("\\.rds$", ".csv", freq_file, ignore.case = TRUE)
}
cat("read...\n")
cat("read ", example_file, "\n")
data <- readRDS(example_file)

cat("read ", freq_file, "\n")
data_p <- readRDS(freq_file)

cat("annotating...\n")
data_p$type <- data$cell_type

cat("Save to ", output_file, "\n")
write.csv(data_p, file = output_file, row.names = FALSE)

cat("Over\n")
cat("Output:", output_file, "\n")

#Rscript freq61_type_transfer.R example.rds freq61_hbf2000_mono_example.rds
