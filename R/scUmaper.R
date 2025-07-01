.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\u001b[31mWelcome to use Single Cell Utility Matrices Processing Engine in R (scUmaper) package.\u001b[0m")
  packageStartupMessage("\u001b[31mscUmaper package can remove doublets and annotate single cell clusters automatically.\u001b[0m")
  packageStartupMessage("\u001b[31mPlease notice that only tissue samples from Homo sapiens can be used by scUmaper, meanwhile the gene names used should be Gene Symbol.\u001b[0m\n")
  packageStartupMessage("\u001b[31mThis package was written by Guo X. et al. from Zhoulab.\u001b[0m")
  packageStartupMessage("\u001b[31mAll copyrights reserved. \u00A9 2025 Zhoulab.\u001b[0m\n")
}
### gobal variables ####
utils::globalVariables(c("percent.mt", "nFeature_RNA", "cluster", "seurat_clusters"))

#' run_scumaper
#' @title run_scumaper
#'
#' @import dplyr
#' @import tibble
#' @import grDevices
#' @import ggplot2
#' @import Seurat
#' @import SeuratObject
#'
#' @description run_scumaper() is a function for executing scUmaper from raw single cell gene expression matrices.
#'
#' @param input_dir Path of the folder storing raw single cell gene expression matrices. This path should include several folders named by the samples, which should include 3 files named "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz".
#' @param output_dir Path to store output files, defaults to create a new folder named "scUmaper_output" under working directory.
#' @param output_plot Whether outputting UMAP plots or not, defaults to TRUE.
#' @param output_plot_format Format of outputting UMAP plots, should be "png" or "pdf", defaults to "png".
#' @param n_features Number of feature genes when clustering, defaults to 3000.
#' @param qc_all Quality control standard of all cells, the cells above how many percentage of mitochondrial gene expression will be removed, defaults to 25.
#' @param qc_immune Quality control standard of immune cells, the cells above how many percentage of mitochondrial gene expression will be removed, defaults to 5.
#' @param doublet_reserve Whether reserving doublet cells as another Seurat object, defaults to FALSE.
#'
#' @return NULL
#'
#' @examplesIf interactive()
#' run_scumaper(input_dir = 'path/to/input', output_dir = 'path/to/output')
#'
#' @export
run_scumaper = function(input_dir = NULL,
                        output_dir = NULL,
                        output_plot = TRUE,
                        output_plot_format = "png",
                        n_features = 3000,
                        qc_all = 25,
                        qc_immune = 5,
                        doublet_reserve = FALSE) {
  #####reading parameters##########################################################
  if(is.null(input_dir)) {cat("Error: `input_dir` parameter cannot be empty. Please type in a valid path.\n"); return()}
  if(!dir.exists(input_dir)) {cat("Error: Please type in a valid `input_dir` path.\n"); return()}
  if(is.null(output_dir)) {
    dir.create(path = paste0(getwd(), "/scUmaper_output"))
    output_dir = paste0(getwd(), "/scUmaper_output")
  }
  if(all(!is.null(output_dir), !dir.exists(output_dir))) {cat("Error: Please type in a valid `output_dir` path.\n"); return()}
  cat("================================================================================\n\n")
  cat("Single cell matrices are reading from:", input_dir, "\n")
  cat("Output files of this scUmaper run will be located in:", output_dir, "\n")
  cat("Feature number of clustering major groups is", n_features, "\n")
  cat("Quality control standard of all cells is less than", qc_all, "percent mitochondrial gene expression.", "\n")
  cat("Quality control standard of immune cells is less than", qc_immune, "percent mitochondrial gene expression.", "\n")
  if (!(output_plot_format %in% c("png", "pdf"))) {output_plot = F; cat("Warning: format", output_plot_format, "is not supported by scUmaper, only `pdf` or `png` are supported.\n")}
  if (!output_plot) cat("UMAP plots will not be created in this run.\n")
  if (doublet_reserve) cat("The doublets deleted will be kept in:", output_dir, "\n")
  cat("\n", "\u001b[31mStart running scUmaper.\u001b[0m\n\n", sep = "")

  cat("\u001b[31m================================================================================\u001b[0m\n\n")
  ref.df = data.frame("cell_types" = c("B/Plasma_cell", "Endothelial_cell", "Myeloid", "T/NK_cell", "Fibroblast", "Epithelial_cell"),
                      "markers" = c("CD79A;CD79B;CD19;MS4A1;MZB1;IGHG1;IGHA1",
                                    "PECAM1;VWF",
                                    "CD68;CSF3R;TPSAB1",
                                    "CD3D;CD3E",
                                    "COL1A1;ACTA2;ACTG2;COL6A3",
                                    "EPCAM;KRT19;KRT8"),
                      "mt_threshold" = c(qc_immune, qc_all, qc_immune, qc_immune, qc_all, qc_all))

  #####read seurats##########################################################
  input.list = list.dirs(path = input_dir, recursive = F, full.names = F) %>% unlist()

  sc = NULL
  for (each in input.list) {
    cat("Reading", each, "\n")
    temp = Seurat::Read10X(paste(input_dir, each, sep = "/")) %>%
      Seurat::CreateSeuratObject(project = each, min.cells = 3, min.features = 200) %>%
      SeuratObject::RenameCells(new.names = paste(each, rownames(.@meta.data), sep = "_"))
    temp[["percent.mt"]] = Seurat::PercentageFeatureSet(temp, pattern = "^MT-")
    temp = temp %>% subset(percent.mt < qc_all) %>% subset(nFeature_RNA > 200 & nFeature_RNA < 6000) ###Qualification
    if (is.null(sc)) {
      sc = temp
    } else {
      sc = merge(sc, temp)
    }
    rm(temp)
    gc()
  }

  rm(input.list, input_dir, each)
  gc()

  cat("Single cell matrices are successfully loaded.\n")

  #####add preliminary annotation##############################################
  cat("================================================================================\n\n")
  cat("\u001b[31mStart preliminarily annotating.\u001b[0m\n\n")

  sc = sc %>% SeuratObject::JoinLayers() %>%
    Seurat::NormalizeData(verbose = F) %>%
    Seurat::FindVariableFeatures(nfeatures = n_features, verbose = F) %>%
    Seurat::ScaleData(verbose = F) %>%
    Seurat::RunPCA(features = setdiff(x = VariableFeatures(.), y = grep("^MT-|^RP[SL]", VariableFeatures(.), value = T)), verbose = F) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = F) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

  run_res = 0.05
  repeat {
    sc = sc %>% Seurat::FindClusters(resolution = run_res, verbose = F)
    if (sc@meta.data$seurat_clusters %>% levels() %>% length() < 8 & run_res < 5) {
      run_res = run_res + 0.05
      gc()
      next
    }

    run_markers = Seurat::FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
    sc$cluster = NA

    for (temp.ident in levels(sc$seurat_clusters)) {
      for (i in 1:nrow(ref.df)) {
        temp.cluster = ref.df[i, "cell_types"]
        temp.markers = strsplit(ref.df[i, "markers"], split = ";") %>% unlist()
        if(any(run_markers[run_markers$cluster == temp.ident, "gene"] %in% temp.markers)) {
          sc$cluster[sc$seurat_clusters == temp.ident] = temp.cluster
        }
      }
    }

    sc$cluster[is.na(sc$cluster)] = "undefined"

    if (unique(sc$cluster) %>% length() < 5 & run_res < 3) {
      run_res = run_res + 0.05
      gc()
      next
    }

    sc@meta.data = sc@meta.data %>% .[,!grepl(pattern = "^RNA_snn_res|seurat_clusters", x = colnames(.))]
    cat("Final clustering resolution for preliminary annotation is", run_res, "\n\n")
    break
  }

  rm(list = c(ls(pattern = "^qc|^run_|^temp"), "i"))
  gc()

  cat("================================================================================\n\n")

  #####remove doublets########################################################
  cat("\u001b[31mStart seeking for doublets\u001b[0m.\n\n")

  for (each in unique(sc$cluster)) {
    temp = subset(sc, cluster == each)
    gc()
    each = each %>% gsub(pattern = "/", replacement = "&", x = .) #replace '/'
    saveRDS(temp, paste0(output_dir, "/scUmaper_temp_", each, ".rds"))
  }

  rm(temp, sc, each)
  gc()

  ref.df = ref.df %>% tibble::column_to_rownames(var = "cell_types")

  sc = NULL
  dbl = NULL
  for (each in list.files(path = output_dir, pattern = "scUmaper_temp_.*\\.rds")) {
    temp = readRDS(paste0(output_dir, "/", each))
    each = each %>% gsub(pattern = "scUmaper_temp_|\\.rds", replacement = "", x = .) %>% gsub(pattern = "&", replacement = "/", x = .) #replace '&'
    cat("\nStart finding doublet of", each, "\n\n")

    if (ncol(temp) < 200) { ###skip too small subclusters
      cat(each, "has too few cells. No cells removed.\n\n")
      if (is.null(sc)) {
        sc = temp
      } else {
        sc = merge(sc, temp)
      }
      gc()
      next
    }

    if (each != "undefined") temp = temp %>% subset(percent.mt < ref.df[each, "mt_threshold"])
    temp = temp %>% Seurat::NormalizeData() %>%
      Seurat::FindVariableFeatures(nfeatures = 2000, verbose = F) %>%
      Seurat::ScaleData(verbose = F) %>%
      Seurat::RunPCA(features = setdiff(x = VariableFeatures(.), y = grep("^MT-|^RP[SL]", VariableFeatures(.), value = T)), verbose = F) %>%
      Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = F) %>%
      Seurat::FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)

    ###find doublets then remove
    run_res = 1
    repeat {
      temp = temp %>% Seurat::FindClusters(resolution = run_res, verbose = F)
      gc()
      if (temp@meta.data$seurat_clusters %>% levels() %>% length() < 15 & run_res < 5) {
        run_res = run_res + 0.1
        next
      }

      if (temp@meta.data$seurat_clusters %>% levels() %>% length() < 3) {
        cat(each, "has too few cells. No cells removed.\n\n")
        break
      }

      run_markers = Seurat::FindAllMarkers(temp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F) %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_head(n = 50)

      if (each %in% c("B/Plasma_cell", "Myeloid", "T/NK_cell")) {
        impossb.genes = ref.df[rownames(ref.df) %>% setdiff(each), "markers"] %>%
          paste(collapse = ";") %>%
          strsplit(split = ";") %>%
          unlist()
        for (i in levels(temp$seurat_clusters)) {
          if (any(run_markers[run_markers$cluster == i,]$gene %in% impossb.genes)) {
            temp.dbl = subset(temp, seurat_clusters == i)
            if (is.null(dbl)) {
              dbl = temp.dbl
            } else {
              dbl = merge(dbl, temp.dbl)
            }
            temp = temp %>% subset(seurat_clusters != i) #remove doublets
            cat("cluster", i, "of", each, "is removed.\n")
          }
        }
      } else if (each %in% c("Fibroblast", "Endothelial_cell", "Epithelial_cell")){
        impossb.genes = ref.df[rownames(ref.df) %>% setdiff(each), "markers"] %>%
          paste(collapse = ";") %>%
          strsplit(split = ";") %>%
          unlist() %>%
          c(., "PTPRC")
        for (i in levels(temp$seurat_clusters)) {
          if (any(run_markers[run_markers$cluster == i,]$gene %in% impossb.genes)) {
            temp.dbl = subset(temp, seurat_clusters == i)
            if (is.null(dbl)) {
              dbl = temp.dbl
            } else {
              dbl = merge(dbl, temp.dbl)
            }
            temp = temp %>% subset(seurat_clusters != i) #remove doublets
            cat("cluster", i, "of", each, "is removed.\n")
          }
        }
      } else {
        impossb.genes = ref.df[c("B/Plasma_cell", "Myeloid", "T/NK_cell"), "markers"] %>%
          paste(collapse = ";") %>%
          strsplit(split = ";") %>%
          unlist() %>%
          c(., "PTPRC")
        impossb.genes2 = ref.df[c("Fibroblast", "Endothelial_cell", "Epithelial_cell"), "markers"] %>%
          paste(collapse = ";") %>%
          strsplit(split = ";") %>%
          unlist()
        for (i in levels(temp$seurat_clusters)) {
          if (any(run_markers[run_markers$cluster == i,]$gene %in% impossb.genes) & any(run_markers[run_markers$cluster == i,]$gene %in% impossb.genes2)) {
            temp.dbl = subset(temp, seurat_clusters == i)
            if (is.null(dbl)) {
              dbl = temp.dbl
            } else {
              dbl = merge(dbl, temp.dbl)
            }
            temp = temp %>% subset(seurat_clusters != i) #remove doublets
            cat("cluster", i, "of", each, "is removed.\n")
          }
        }
      }

      cat("\nFinal clustering resolution of", each, "is", run_res, "\n\n")
      break
    }

    if (is.null(sc)) {
      sc = temp
    } else {
      sc = merge(sc, temp)
    }
    gc()
  }

  rm(list = c(ls(pattern = "^temp|^run_"), "each", "i", "impossb.genes"))
  gc()

  if (!is.null(dbl)) {
    dbl = dbl %>% SeuratObject::JoinLayers()
    cat("\nTotally", ncol(dbl), "possible doublet cells have been removed.\n")
    if (doublet_reserve) {
      saveRDS(dbl, paste0(output_dir, "/scUmaper_removed_doublet_seurat.rds"))
      cat(paste("Removed doublets are saved at:", paste0(output_dir, "/final_removed_doublet_seurat.rds")), "\n\n")
    }
  } else {
    cat("No doublet was found or removed.\n\n")
  }

  rm(dbl)
  gc()

  #####get final results########################################################
  cat("================================================================================\n\n")
  cat("\u001b[31mStart final annotating.\u001b[0m\n\n")

  sc = sc %>% SeuratObject::JoinLayers() %>%
    Seurat::NormalizeData(verbose = F) %>%
    Seurat::FindVariableFeatures(nfeatures = n_features, verbose = F) %>%
    Seurat::ScaleData(verbose = F) %>%
    Seurat::RunPCA(features = setdiff(x = VariableFeatures(.), y = grep("^MT-|^RP[SL]", VariableFeatures(.), value = T)), verbose = F) %>%
    Seurat::RunUMAP(reduction = "pca", dims = 1:20, verbose = F) %>%
    Seurat::FindNeighbors(reduction = "pca", dims = 1:20, verbose = F)
  gc()

  ref.df = ref.df %>% tibble::rownames_to_column(var = "cell_types")
  run_res = 0.05
  repeat {
    sc = sc %>% Seurat::FindClusters(resolution = run_res, verbose = F)
    if (sc@meta.data$seurat_clusters %>% levels() %>% length() < 8 & run_res < 5) {
      run_res = run_res + 0.05
      next
    }

    run_markers = Seurat::FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = F)
    sc$cluster = NA

    for (temp.ident in levels(sc$seurat_clusters)) {
      for (i in 1:nrow(ref.df)) {
        temp.cluster = ref.df[i, "cell_types"]
        temp.markers = strsplit(ref.df[i, "markers"], split = ";") %>% unlist()
        if(any(run_markers[run_markers$cluster == temp.ident, "gene"] %in% temp.markers)) {
          sc$cluster[sc$seurat_clusters == temp.ident] = temp.cluster
        }
      }
    }

    sc$cluster[is.na(sc$cluster)] = "undefined"

    if (unique(sc$cluster) %>% length() < 5 & run_res < 3) {
      run_res = run_res + 0.05
      next
    }

    sc@meta.data = sc@meta.data %>% .[,!grepl(pattern = "^RNA_snn_res|seurat_clusters", x = colnames(.))]
    cat("\nLast clustering resolution for final annotation is", run_res, "\n\n")
    cat("================================================================================\n\n")
    break
  }

  sc$cluster = sc$cluster %>% factor(levels = c("T/NK_cell", "B/Plasma_cell", "Myeloid", "Fibroblast", "Endothelial_cell", "Epithelial_cell", "undefined"))
  Seurat::Idents(sc) = sc$cluster

  file.remove(list.files(path = output_dir, pattern = "scUmaper_temp_.*\\.rds", full.names = T))
  rm(list = c(ls(pattern = "^run_|^temp"), "i"))
  gc()

  saveRDS(sc, paste0(output_dir, "/scUmaper_annotated_singlet_seurat.rds"))
  cat("\nEventually,", ncol(sc), "singlet cells have been kept and annotated.\n")
  cat("Final processed seurat object has been successfully saved as:\n")
  cat(output_dir, "/final_annotated_singlet_seurat.rds\n\n", sep = "")

  if (output_plot) {
    dir.create(path = paste0(output_dir, "/markers_plot"))
    if (output_plot_format == "pdf") {
      temp.p = Seurat::DimPlot(sc,
                               reduction = "umap",
                               group.by = "cluster",
                               label = F,
                               pt.size = .1, raster = F,
                               cols = c("T/NK_cell" = "#E56F5E",
                                        "B/Plasma_cell" = "#F6C957",
                                        "Myeloid" = "#FBE8D5",
                                        "Fibroblast" = "#ABD0F1",
                                        "Endothelial_cell" = "#79CEED",
                                        "Epithelial_cell" = "#43978F",
                                        "undefined" = "#B5B5B5"))
      grDevices::pdf(paste0(output_dir, "/final_umap.pdf"), height = 6, width = 7.5)
      print(temp.p)
      grDevices::dev.off()

      for (each in c("PTPRC", "CD3E", "CD79A", "MS4A1", "MZB1", "IGHG1", "IGHA1", "CD68", "CSF3R", "KIT", "TPSAB1",
                     "EPCAM", "KRT8", "KRT19", "PECAM1", "COL1A1", "ACTA2")) {
        temp.p = Seurat::FeaturePlot(sc,
                                     features = each,
                                     cols = c("lightgrey","red"),
                                     pt.size = 0.1,
                                     raster = F,
                                     ncol = 1)
        grDevices::pdf(paste0(output_dir, "/markers_plot/", each, "_umap.pdf"), height = 6, width = 6.5)
        print(temp.p)
        grDevices::dev.off()
      }
    } else if(output_plot_format == "png") {
      temp.p = Seurat::DimPlot(sc,
                               reduction = "umap",
                               group.by = "cluster",
                               label = F,
                               pt.size = .1, raster = F,
                               cols = c("T/NK_cell" = "#E56F5E",
                                        "B/Plasma_cell" = "#F6C957",
                                        "Myeloid" = "#FBE8D5",
                                        "Fibroblast" = "#ABD0F1",
                                        "Endothelial_cell" = "#79CEED",
                                        "Epithelial_cell" = "#43978F",
                                        "undefined" = "#B5B5B5"))
      ggplot2::ggsave(filename = paste0(output_dir, "/final_umap.png"), plot = temp.p, height = 6, width = 7.5)

      for (each in c("PTPRC", "CD3E", "CD79A", "MS4A1", "MZB1", "IGHG1", "IGHA1", "CD68", "CSF3R", "KIT", "TPSAB1",
                     "EPCAM", "KRT8", "KRT19", "PECAM1", "COL1A1", "ACTA2")) {
        temp.p = Seurat::FeaturePlot(sc,
                                     features = each,
                                     cols = c("lightgrey","red"),
                                     pt.size = 0.1,
                                     raster = F,
                                     ncol = 1)
        ggplot2::ggsave(paste0(output_dir, "/markers_plot/", each, "_umap.png"), plot = temp.p, height = 6, width = 6.5)
      }
    }
    cat("Final umaps have been successfully saved at:", paste0(output_dir, "/markers_plot"), "\n\n")
  }

  rm(list = ls())
  gc()

  cat("\u001b[31m================================================================================\u001b[0m\n\n")
  cat("\u001b[31mAll progressions have been completed. Thank you for using scUmaper.\u001b[0m\n\n")
  cat("\u001b[31mPlease cite scUmaper while publishing your papers.\u001b[0m\n")
  cat("\u001b[31m\u00A9 2025 Zhoulab. All rights reserved.\u001b[0m\n")
}
