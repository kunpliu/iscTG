# iscTG: integrate single cell Transcriptome and Genome-wide association studies

iscTG is a comprehensive workflow designed to integrate single-cell transcriptome data with genome-wide association studies (GWAS) to uncover the genetic basis of complex traits at the cellular level. This repository provides a set of tools and scripts to facilitate this integration, along with a workflow to guide users through the process.

![Flowchart](workflow.png)

## Workflow

The iscTG workflow is designed to streamline the integration of single-cell RNA-seq data with GWAS summary statistics. The key steps include:

1. **Data Preparation**: Ensure single-cell RNA-seq data is in Seurat format and GWAS summary statistics are properly formatted.
2. **Feature Extraction**: Use `extra_features.R` to extract relevant features from gene expression matrices and MAGMA output files.
3. **Testing**: Use `test.R` to run a quick test of the workflow using the provided `mono500.rds` dataset.
4. **Visualization**: Use `out_fig.R` to visualize the results in the form of heatmaps or other relevant plots.
## main function

net61.R
### To get PCC net

61regionhbf2k.R
### 描述
computeP.R
### 描述

## Quick Start

To get started with iscTG, follow these steps:

1. **Install dependencies**:
   ```bash
   Rscript -e "install.packages(c('Seurat', 'ggplot2'))"
   ```

2. **Run the test script**:
   ```bash
   Rscript test.R
   ```

3. **Visualize the results**:
   ```bash
   Rscript out_fig.R
   ```

## Small Tools

### extra_features.R

This script is used to extract features required for iscTG from gene expression matrices and MAGMA output files. It assumes that the single-cell data is in Seurat format.

**Usage**:
```R
source("extra_features.R")
extract_features(seurat_object, magma_output_file)
```

### test.R

This script is used to test the iscTG workflow using the provided `mono500.rds` dataset. It demonstrates the basic functionality of the workflow.

**Usage**:
```R
source("test.R")
run_test("mono500.rds")
```

### mono500.rds

This is a test dataset containing 500 monocytes and 500 other cells (e.g., red cells). It is designed to provide a quick test of the workflow due to its small size.

### out_fig.R

This script is used to visualize the results of the iscTG workflow. It can generate heatmaps or other relevant plots from the Seurat-formatted data.

**Usage**:
```R
source("out_fig.R")
generate_heatmap(seurat_object)
```
