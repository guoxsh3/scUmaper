# **scUmaper: Single Cell Utility Matrices Processing Engine in R**

## **Requirements and Installation**

scUmaper is implemented by R. Make sure your R version \>= 4.3.0. The packages below are required for scUmaper:

-   tibble

-   dplyr

-   grDevices

-   ggplot2

-   Seurat \>= 5.0.0

-   SeuratObject \>= 5.0.0

Download and install scUmaper:

``` r
devtools::install_github("guoxsh3/scUmaper")
```

## **Before you start**

Please make sure:

-   Currently, only tissue samples from *Homo sapiens* (human) can be analyzed directly from raw matrices through *scUmaper::run_scumaper* function.

-   Input format of *scUmaper::run_scumaper* is similar to the output of CellRanger. Three files, named "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz", should be included in a folder named by each sample.

-   Gene symbols, i.e. "PTPRC", are avaliable in "features.tsv.gz" files.

-   The folders of each sample should be included in same parent folder, whose directory will be the input of scUmaper.

-   Samples from other species, i.e. *Mus musculus* (mouse), or PBMC samples, are only supported by processed Seurat objects through *scUmaper::run_custom_scumaper* function. Meanwhile, users should provide a customized list of marker genes to use.

-   Your Seurat package has a version \>= 5.0.0.

## **Usage**

### For samples of human tissue

``` r
library(scUmaper)
run_scumaper(input_dir = NULL,
             output_dir = NULL,
             output_plot = TRUE,
             output_plot_format = "png",
             n_features = 3000,
             qc_all = 25,
             qc_immune = 5,
             doublet_reserve = FALSE)
```

Required arguments:

-   input_dir: Path of the folder storing raw single cell gene expression matrices.

Optional arguments:

-   output_dir: Path to store output files, defaults to create a new folder named "scUmaper_output" under working directory.

-   output_plot: Whether outputting UMAP plots or not, defaults to TRUE.

-   output_plot_format: Format of outputting UMAP plots, should be "png" or "pdf", defaults to "png".

-   n_features: Number of feature genes when clustering, defaults to 3000. Recommend 1000 to 3000.

-   qc_all: Quality control standard of all cells, the cells above how many percentage of mitochondrial gene expression will be removed, defaults to 25. Recommend up to 25.

-   qc_immune: Quality control standard of immune cells, the cells above how many percentage of mitochondrial gene expression will be removed, defaults to 5. Recommend up to 10.

-   doublet_reserve: Whether reserving doublet cells as another Seurat object, defaults to FALSE. Recommend TRUE if you are experienced in single cell analysis.

Quick start:

``` r
library(scUmaper)
run_scumaper('path/to/input')
```

### For other samples with customized marker genes

``` r
library(scUmaper)
run_custom_scumaper(input_seurat = NULL,
                    input_genelist = NULL,
                    output_dir = NULL,
                    output_plot = TRUE,
                    output_plot_format = "png",
                    n_features = 3000,
                    doublet_reserve = FALSE)
```

Required arguments:

-   input_seurat: Seurat object to remove doublets and annotate.
-   input_genelist: A list of several vectors, which named by target clusters and composed of corresponding marker genes.

Optional arguments:

-   output_dir: Path to store output files, defaults to create a new folder named "scUmaper_output" under working directory.
-   output_plot: Whether outputting UMAP plots or not, defaults to TRUE.
-   output_plot_format: Format of outputting UMAP plots, should be "png" or "pdf", defaults to "png".
-   n_features: Number of feature genes when clustering, defaults to 3000.
-   doublet_reserve: Whether reserving doublet cells as another Seurat object, defaults to FALSE.

Quick start:

``` r
library(scUmaper)
genelist = list("B/Plasma_cell" = c('CD79A', 'CD79B', 'CD19', 'MS4A1', 'MZB1', 'IGHG1', 'IGHA1'),
                "Endothelial_cell" = c('PECAM1', 'VWF'),
                "Myeloid" = c('CD68', 'CSF3R', 'TPSAB1'),
                "T/NK_cell" = c('CD3D', 'CD3E'),
                "Fibroblast" = c('COL1A1', 'ACTA2', 'ACTG2', 'COL6A3'),
                "Epithelial_cell" = c('EPCAM', 'KRT19', 'KRT8'))
run_custom_scumaper(scRNA, 
                    genelist)
```

## **Hints**

-   Running scUmaper might require a long duration, especially when you input cells over 100,000. Executing scUmaper as a background R job through *nohup* and *Rscript* commands could save your time.

-   scUmaper might need a large RAM space over 32 GBs, especially when you input cells over 100,000. Make sure your server or PC have enough RAM. Recommend run scUmaper in a linux-based server with RAM larger than 64 GBs.

-   scUmaper follows a extremely strict standard to keep singlets. Therefore, some singlets might be mis-labelled as doublets and removed. For experienced analyst, *doublet_reserve* argument is recommended to TRUE. After automatic analysis, the operator could exterminate the doublets manually, and merge mis-labelled singlets into auto-labelled singlets for further analysis.

-   When using *scUmaper::run_custom_scumaper* function, inputting a Seurat object quality-controlled by expression of mitochondrial genes might improve the effectiveness of doublet-removal and annotation.

## **Citation**

If you find our work helpful in your research or work, please cite us.

X. Guo, G. Zhou. scUmaper: An Automated Tool for Doublet Removal and Cell Type Annotation in scRNA-Seq Data. Under review, 2025.

## **Questions & Problems**

If you have any questions or problems for code questions, please contact Xushun Guo.

-   Xushun Guo (guoxsh3\@mail2.sysu.edu.cn)\
    MD, Resident Doctor,\
    Zhongshan School of Medicine, Sun Yat-sen University, Canton, China.

For any other further questions or requests, please contact Gaoshi Zhou, the Principle Investigator of our bioinfomatic lab.

-   Gaoshi Zhou (zhougshi\@mail2.sysu.edu.cn)\
    MD, Research Associate, Resident Doctor,\
    The First Affiliated Hospital, Sun Yat-sen University, Canton, China.
