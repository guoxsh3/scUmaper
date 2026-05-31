# **scUmaper: Single Cell Utility Matrices Processing Engine in R**

## **Introduction**

scUmaper is an R package designed to streamline the processing of single-cell transcriptomics sequencing data. It automatically removes heterotypic doublets and annotates major cell lineages with high reliability by supervising biological markers. It allows researchers to process raw 10x Genomics gene expression matrices directly, bypassing complex intermediate steps and avoiding potential analytical errors.

## **Requirements and Installation**

scUmaper is implemented by R. Make sure your R version \>= 4.3.0. The packages below are required for scUmaper:

-   dplyr

-   grDevices

-   ggplot2

-   Seurat \>= 5.0.0

-   SeuratObject \>= 5.0.0

-   presto

Download and install scUmaper:

``` r
devtools::install_github("guoxsh3/scUmaper")
```

## **Before you start**

Please review the following requirements before running the pipeline:

-   **Supported Species:** Currently, only tissue samples from *Homo sapiens* (human) or *Mus musculus* (mouse) can be analyzed directly from raw matrices using the `run_scumaper` function.

-   **Input Format:** The input format is standard 10x Genomics Cell Ranger output. Three files, `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz` (or their unzipped versions), must be included in a subfolder named after each sample.

-   **Directory Structure:** The subfolders for each sample must be placed within a single parent directory. The path to this parent directory will be the `input_dir` argument. The file structure should look like this:

    ``` bash
    /path/to/input/
    ├── Sample1/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    ├── Sample2/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    └── ...
    ```

-   **Gene Naming**: Gene symbols (e.g., `"PTPRC"` for human or `"Ptprc"` for mouse) must be used in your `features.tsv.gz` files.

-   **Custom Mode**: Samples from other species or PBMC samples are supported via the `run_custom_scumaper` function. In this mode, users must provide a pre-processed Seurat object and a customized list of marker genes.

## **Usage**

### 1. For Raw Human or Mouse Tissue Samples

``` r
library(scUmaper)

# Quick Start
run_scumaper(input_dir = "path/to/input", species = "human")
```

**Full Arguments:**

``` r
run_scumaper(input_dir = NULL,
             output_dir = NULL,
             species = "human",
             output_plot = TRUE,
             output_plot_format = "png",
             n_features = 3000,
             qc_all = 25,
             qc_immune = 5,
             doublet_reserve = FALSE)
```

-   `input_dir` *(Required)*: Path to the parent directory containing sample subfolders.

-   `output_dir` *(Optional)*: Path to store output files. Defaults to creating a `scUmaper_output` folder in the current working directory.

-   `species` *(Optional)*: Species of the data, either `"human"` or `"mouse"`. Defaults to `"human"`.

-   `output_plot` *(Optional)*: Logical. Whether to output UMAP plots. Defaults to `TRUE`.

-   `output_plot_format` *(Optional)*: Format for UMAP plots, either `"png"` or `"pdf"`. Defaults to `"png"`.

-   `n_features` *(Optional)*: Number of highly variable genes (HVGs) used for clustering. Defaults to `3000`. Recommended range: 1000 to 4000.

-   `qc_all` *(Optional)*: Mitochondrial percentage cutoff for all cells. Cells exceeding this threshold are removed. Defaults to `25`. Recommended range: 10 to 25.

-   `qc_immune` *(Optional)*: Mitochondrial percentage cutoff for immune cells. Defaults to `5`. Recommended range: 5 to 10.

-   `doublet_reserve` *(Optional)*: Logical. Whether to save the removed doublets as a separate `.rds` Seurat object. Defaults to `FALSE`. Recommended to set as `TRUE` for experienced analysts.

### 2. For Other Samples with Customized Marker Genes

``` r
library(scUmaper)

# Quick Start Example

genelist <- list(
  "B/Plasma_cell" = c('CD79A', 'CD79B', 'CD19', 'MS4A1', 'MZB1', 'IGHG1', 'IGHA1'),
  "Endothelial_cell" = c('PECAM1', 'VWF'),
  "Myeloid" = c('CD68', 'CSF3R', 'TPSAB1'),
  "T/NK_cell" = c('CD3D', 'CD3E'),
  "Fibroblast" = c('COL1A1', 'ACTA2', 'ACTG2', 'COL6A3'),
  "Epithelial_cell" = c('EPCAM', 'KRT19', 'KRT8')
)

run_custom_scumaper(input_seurat = scRNA_object, input_genelist = genelist)
```

**Full Arguments:**

``` r
run_custom_scumaper(input_seurat = NULL,
                    input_genelist = NULL,
                    output_dir = NULL,
                    output_plot = TRUE,
                    output_plot_format = "png",
                    n_features = 3000,
                    doublet_reserve = FALSE)
```

-   `input_seurat` *(Required)*: A pre-built Seurat object to be processed.

-   `input_genelist` *(Required)*: A named list of character vectors, where names are the target cell clusters and vectors contain the corresponding marker genes.

## **Hints & Best Practices**

-   **Computational Resources**: Running scUmaper on datasets exceeding 100,000 cells can be memory-intensive and time-consuming. We highly recommend running scUmaper on a Linux-based server with at least 64 GB of RAM. Executing the script as a background job (e.g., using `nohup Rscript script.R &`) is also advised.

-   **Strict Doublet Removal**: scUmaper employs a highly strict standard to preserve pure singlets. Consequently, some edge-case singlets might be mislabeled as doublets and removed. For experienced analysts, we recommend setting `doublet_reserve = TRUE`. You can manually inspect the removed doublets and merge any mislabeled singlets back into your final dataset.

-   **Custom Mode Input**: When using `run_custom_scumaper`, inputting a Seurat object that has already been quality-controlled (e.g., filtering out high-mitochondrial cells) will significantly improve the accuracy of doublet removal and cell type annotation.

## **Citation**

If you find scUmaper helpful in your research, please cite our publication:

Guo X, Zhang M, Xie Z, Wang Y, Zhang S, Zhou G. **scUmaper: An automated framework for doublet removal and cell-type annotation in single-cell transcriptomics**. *iScience*. 2026;29(6):115850. <doi:10.1016/j.isci.2026.115850>. PMID: 42164849.

## **Contact & Support**

If you encounter any code-related issues or have questions regarding the software, please contact:

-   **Dr. Xushun Guo** (guoxsh3\@mail3.sysu.edu.cn)

    MD, Resident Doctor,

    Guangxi Hospital Division of The First Affiliated Hospital, Sun Yat-sen University, Nanning, China.

For other inquiries or collaboration requests, please contact the Principal Investigator of our bioinformatics lab:

-   **Dr. Gaoshi Zhou** (zhougshi\@mail3.sysu.edu.cn)

    MD, Research Associate, Resident Doctor,

    The First Affiliated Hospital, Sun Yat-sen University, Guangzhou, China.
