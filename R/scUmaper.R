utils::globalVariables(c(".", "percent.mt", "nFeature_RNA", "cluster", "seurat_clusters", "gene"))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to scUmaper: Single Cell Utility Matrices Processing Engine in R.")
  packageStartupMessage("scUmaper automatically removes heterotypic doublets and annotates major cell lineages from scRNA-seq data.")
  packageStartupMessage("Use run_scumaper() for human/mouse tissues with Gene Symbol input.")
  packageStartupMessage("Use run_custom_scumaper() with a user-supplied marker list otherwise.")
  packageStartupMessage("Citation: Guo X, et al. iScience. 2026;29(6):115850. doi:10.1016/j.isci.2026.115850. PMID: 42164849.")
  packageStartupMessage("Type cite_scumaper() to view the citation. (c) 2026 Zhoulab. All rights reserved.")
}

#' Print citation for scUmaper
#'
#' @return Invisibly returns the citation string.
#' @export
cite_scumaper <- function() {
  msg <- paste0(
    "Please cite scUmaper:\n",
    "Guo X, Zhang M, Xie Z, Wang Y, Zhang S, Zhou G. ",
    "scUmaper: An automated framework for doublet removal and cell-type annotation ",
    "in single-cell transcriptomics. iScience. 2026;29(6):115850. ",
    "doi:10.1016/j.isci.2026.115850. PMID: 42164849."
  )
  message(msg)
  invisible(msg)
}

#' Built-in marker gene libraries used by scUmaper
#'
#' @param species Either "human" or "mouse".
#' @return A data.frame with two columns: cell types and markers.
#' @export
scumaper_marker_genes <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    data.frame(
      cell_types = c("B/Plasma_cell", "Endothelial_cell", "Myeloid", "T/NK_cell", "Fibroblast", "Epithelial_cell"),
      markers = c("CD79A;CD79B;CD19;MS4A1;MZB1;IGHG1;IGHA1",
                  "PECAM1;VWF",
                  "CD68;CSF3R;TPSAB1",
                  "CD3D;CD3E",
                  "COL1A1;ACTA2;ACTG2;COL6A3",
                  "EPCAM;KRT19;KRT8"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      cell_types = c("B/Plasma_cell", "Endothelial_cell", "Myeloid", "T/NK_cell", "Fibroblast", "Epithelial_cell"),
      markers = c("Cd79a;Cd19;Ms4a1;Mzb1;Sdc1;Ighg1;Igha;Jchain",
                  "Pecam1;Eng;Cdh5",
                  "Cd68;Csf3r;S100a8;Mcpt4",
                  "CD3d;Cd3e",
                  "Col1a1;Col1a2;Acta2;Actg2",
                  "Epcam;Krt8;Krt18"),
      stringsAsFactors = FALSE
    )
  }
}

.scumaper_section <- function() {
  message(strrep("=", 80))
}

.scumaper_calc_dims <- function(n_features) {
  d <- round(n_features / 100)
  if (d < 10) d <- 10L
  if (d > 40) d <- 40L
  as.integer(d)
}

.scumaper_check_presto <- function() {
  if (!requireNamespace("presto", quietly = TRUE)) {
    stop(
      "Package 'presto' is required by scUmaper but not installed.\n",
      "scUmaper relies on presto to accelerate marker gene detection within Seurat v5.\n",
      "Please install presto first via one of the following commands:\n",
      "  remotes::install_github('immunogenomics/presto')\n",
      "or:\n",
      "  install.packages('presto')\n",
      "Then re-run scUmaper.",
      call. = FALSE
    )
  }
  message("Package 'presto' is available; FindAllMarkers will be accelerated.")
  invisible(TRUE)
}

.scumaper_find_markers <- function(obj, only.pos = TRUE, min.pct = 0.25,
                                   logfc.threshold = 0.25, top_n = NULL) {
  m <- Seurat::FindAllMarkers(obj, only.pos = only.pos, min.pct = min.pct,
                              logfc.threshold = logfc.threshold, verbose = FALSE)
  if (!is.null(top_n) && nrow(m) > 0) {
    m <- dplyr::ungroup(dplyr::slice_head(dplyr::group_by(m, cluster), n = top_n))
  }
  as.data.frame(m)
}

.scumaper_preprocess <- function(obj, n_features, dims,
                                 exclude_pattern = "^MT-|^RP[SL]|^mt-|^Rp[sl]") {
  obj <- SeuratObject::JoinLayers(obj)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_features, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  vf <- Seurat::VariableFeatures(obj)
  if (!is.null(exclude_pattern) && length(vf) > 0) {
    vf <- setdiff(vf, grep(exclude_pattern, vf, value = TRUE))
  }
  obj <- Seurat::RunPCA(obj, features = vf, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = seq_len(dims), verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = seq_len(dims), verbose = FALSE)
  obj
}

.scumaper_marker_list_from_df <- function(ref_df) {
  out <- lapply(ref_df$markers, function(x) strsplit(x, ";", fixed = TRUE)[[1]])
  names(out) <- ref_df$cell_types
  out
}

.scumaper_iterative_annotate <- function(obj, marker_list, res_step = 0.05,
                                         res_init = 0.05, min_clusters = 8L,
                                         max_res_phase1 = 2, min_types = 5L,
                                         max_res_phase2 = 2) {
  run_res <- res_init
  repeat {
    obj <- Seurat::FindClusters(obj, resolution = run_res, verbose = FALSE)
    n_clust <- length(levels(obj@meta.data$seurat_clusters))
    if (n_clust < min_clusters && run_res < max_res_phase1) {
      run_res <- run_res + res_step
      gc()
      next
    }
    run_markers <- .scumaper_find_markers(obj)
    obj$cluster <- NA_character_
    for (cl_id in levels(obj$seurat_clusters)) {
      cl_genes <- run_markers[run_markers$cluster == cl_id, "gene"]
      for (j in seq_along(marker_list)) {
        if (any(cl_genes %in% marker_list[[j]])) {
          obj$cluster[obj$seurat_clusters == cl_id] <- names(marker_list)[j]
        }
      }
    }
    obj$cluster[is.na(obj$cluster)] <- "undefined"
    if (length(unique(obj$cluster)) < min_types && run_res < max_res_phase2) {
      run_res <- run_res + res_step
      gc()
      next
    }
    keep <- !grepl("^RNA_snn_res|seurat_clusters", colnames(obj@meta.data))
    obj@meta.data <- obj@meta.data[, keep, drop = FALSE]
    message("Final clustering resolution = ", round(run_res, 3))
    break
  }
  obj
}

.scumaper_split_save_groups <- function(sc, output_dir) {
  for (each in unique(sc$cluster)) {
    sub <- subset(sc, cluster == each)
    gc()
    safe_name <- gsub("/", "&", each, fixed = TRUE)
    saveRDS(sub, file.path(output_dir, paste0("_tmp_scumaper_", safe_name, ".rds")))
  }
  invisible(NULL)
}

.scumaper_cleanup_tmp <- function(output_dir) {
  fs <- list.files(output_dir, pattern = "^_tmp_scumaper_.*\\.rds$", full.names = TRUE)
  if (length(fs) > 0) file.remove(fs)
  invisible(NULL)
}

.scumaper_detect_doublets_default <- function(temp, group_name, marker_list,
                                              n_features = 2000,
                                              mt_threshold = NULL,
                                              immune_groups = c("B/Plasma_cell", "Myeloid", "T/NK_cell"),
                                              stromal_groups = c("Fibroblast", "Endothelial_cell", "Epithelial_cell")) {
  dims <- .scumaper_calc_dims(n_features)
  if (!is.null(mt_threshold) && group_name != "undefined") {
    temp <- subset(temp, percent.mt < mt_threshold)
  }
  if (ncol(temp) < 200) {
    message(group_name, " has too few cells; no cells removed.")
    return(list(singlet = temp, doublet = NULL))
  }
  temp <- .scumaper_preprocess(temp, n_features = n_features, dims = dims)

  dbl <- NULL
  run_res <- 1
  repeat {
    temp <- Seurat::FindClusters(temp, resolution = run_res, verbose = FALSE)
    gc()
    n_clust <- length(levels(temp@meta.data$seurat_clusters))
    if (n_clust < 15 && run_res < 5) {
      run_res <- run_res + 0.1
      next
    }
    if (n_clust < 3) {
      message(group_name, " has too few subclusters; no cells removed.")
      break
    }
    run_markers <- .scumaper_find_markers(temp, top_n = 50)

    doublet_clusters <- character()
    if (group_name %in% immune_groups) {
      others <- setdiff(c(immune_groups, stromal_groups), group_name)
      impossb <- unique(unlist(marker_list[others]))
      for (i in levels(temp$seurat_clusters)) {
        cl_genes <- run_markers[run_markers$cluster == i, "gene"]
        if (any(cl_genes %in% impossb)) doublet_clusters <- c(doublet_clusters, i)
      }
    } else if (group_name %in% stromal_groups) {
      others <- setdiff(c(immune_groups, stromal_groups), group_name)
      impossb <- unique(c(unlist(marker_list[others]), "PTPRC", "Ptprc"))
      for (i in levels(temp$seurat_clusters)) {
        cl_genes <- run_markers[run_markers$cluster == i, "gene"]
        if (any(cl_genes %in% impossb)) doublet_clusters <- c(doublet_clusters, i)
      }
    } else {
      ptprc_markers <- c("PTPRC", "Ptprc")
      stromal_markers <- unique(unlist(marker_list[stromal_groups]))
      for (i in levels(temp$seurat_clusters)) {
        cl_genes <- run_markers[run_markers$cluster == i, "gene"]
        group_hits <- vapply(marker_list, function(mk) any(cl_genes %in% mk), logical(1))
        cond_a <- sum(group_hits) >= 2
        has_ptprc <- any(cl_genes %in% ptprc_markers)
        has_stromal <- any(cl_genes %in% stromal_markers)
        cond_b <- has_ptprc && has_stromal
        if (cond_a || cond_b) doublet_clusters <- c(doublet_clusters, i)
      }
    }

    for (i in doublet_clusters) {
      tmp_dbl <- subset(temp, seurat_clusters == i)
      dbl <- if (is.null(dbl)) tmp_dbl else merge(dbl, tmp_dbl)
      temp <- subset(temp, seurat_clusters != i)
      message("cluster ", i, " of ", group_name, " is removed.")
    }
    message("Final subclustering resolution of ", group_name, ": ", round(run_res, 3))
    break
  }
  list(singlet = temp, doublet = dbl)
}

.scumaper_detect_doublets_custom <- function(temp, group_name, marker_list, n_features = 2000) {
  dims <- .scumaper_calc_dims(n_features)
  if (ncol(temp) < 200) {
    message(group_name, " has too few cells; no cells removed.")
    return(list(singlet = temp, doublet = NULL))
  }
  temp <- .scumaper_preprocess(temp, n_features = n_features, dims = dims, exclude_pattern = NULL)

  dbl <- NULL
  run_res <- 1
  repeat {
    temp <- Seurat::FindClusters(temp, resolution = run_res, verbose = FALSE)
    gc()
    n_clust <- length(levels(temp@meta.data$seurat_clusters))
    if (n_clust < 15 && run_res < 5) {
      run_res <- run_res + 0.1
      next
    }
    if (n_clust < 3) {
      message(group_name, " has too few subclusters; no cells removed.")
      break
    }
    run_markers <- .scumaper_find_markers(temp, top_n = 50)

    if (group_name != "undefined") {
      impossb <- setdiff(unique(unlist(marker_list)), marker_list[[group_name]])
      for (i in levels(temp$seurat_clusters)) {
        cl_genes <- run_markers[run_markers$cluster == i, "gene"]
        if (any(cl_genes %in% impossb)) {
          tmp_dbl <- subset(temp, seurat_clusters == i)
          dbl <- if (is.null(dbl)) tmp_dbl else merge(dbl, tmp_dbl)
          temp <- subset(temp, seurat_clusters != i)
          message("cluster ", i, " of ", group_name, " is removed.")
        }
      }
    } else {
      for (i in levels(temp$seurat_clusters)) {
        cl_genes <- run_markers[run_markers$cluster == i, "gene"]
        group_hits <- vapply(marker_list, function(mk) any(cl_genes %in% mk), logical(1))
        if (sum(group_hits) >= 2) {
          tmp_dbl <- subset(temp, seurat_clusters == i)
          dbl <- if (is.null(dbl)) tmp_dbl else merge(dbl, tmp_dbl)
          temp <- subset(temp, seurat_clusters != i)
          message("cluster ", i, " of ", group_name, " is removed.")
        }
      }
    }
    message("Final subclustering resolution of ", group_name, ": ", round(run_res, 3))
    break
  }
  list(singlet = temp, doublet = dbl)
}

.scumaper_loop_groups <- function(output_dir, marker_list, mode = c("default", "custom"),
                                  ref_df_indexed = NULL, n_features_sub = 2000) {
  mode <- match.arg(mode)
  files <- list.files(output_dir, pattern = "^_tmp_scumaper_.*\\.rds$", full.names = FALSE)
  sc <- NULL
  dbl <- NULL
  for (f in files) {
    temp <- readRDS(file.path(output_dir, f))
    each <- gsub("^_tmp_scumaper_|\\.rds$", "", f)
    each <- gsub("&", "/", each, fixed = TRUE)
    message("\nProcessing major group: ", each)

    if (mode == "default") {
      mt_thr <- if (!is.null(ref_df_indexed) && each %in% rownames(ref_df_indexed))
        ref_df_indexed[each, "mt_threshold"] else NULL
      res <- .scumaper_detect_doublets_default(temp, each, marker_list,
                                               n_features = n_features_sub,
                                               mt_threshold = mt_thr)
    } else {
      res <- .scumaper_detect_doublets_custom(temp, each, marker_list,
                                              n_features = n_features_sub)
    }

    if (!is.null(res$doublet)) {
      dbl <- if (is.null(dbl)) res$doublet else merge(dbl, res$doublet)
    }
    sc <- if (is.null(sc)) res$singlet else merge(sc, res$singlet)
    gc()
  }
  list(singlet = sc, doublet = dbl)
}

.scumaper_default_palette <- function() {
  c("T/NK_cell" = "#E56F5E",
    "B/Plasma_cell" = "#F6C957",
    "Myeloid" = "#FBE8D5",
    "Fibroblast" = "#ABD0F1",
    "Endothelial_cell" = "#79CEED",
    "Epithelial_cell" = "#43978F",
    "undefined" = "#B5B5B5")
}

.scumaper_custom_palette <- function() {
  c("#EE3743","#33a02c","#0B7FAB","#6a3d9a","#ff5084","#ff7f00","#01ceff","#F6A9BD","#fdbf6f","#c51b7d",
    "#6F6C9E","#01665e","#a6cee3","#b2df8a","#ffff99","#bf812d","#999999","#BB7DB2","#824615","#ffff33","#86DBD4",
    "#BFE2E3","#A1CFFA","#78BDAD","#D45651","#397A7F","#F0918E","#EEE8DA","#1F5392","#A0BFAF",
    "#AE98D6","#ECCBDC","#54BAD3","#8b4a4b","#DB896C","#AABAC2","#ffae3b",
    "#CCCCCC","#B5B5B5","#A092C3","#DDA0DD","#03A4C6","#7AC5CD","#A3D9ED","#00E5EE","#F9EBDA","#F5DEB3",
    "#98689E","#E84115","#FFB5C5","#00FF7F","#F9D01C","#B03060","#00ABDC","#D2691E","#03A464","#FF7F00",
    "#8968CD","#1C5B75")
}

.scumaper_save_dimplot <- function(sc, output_dir, format, group_by = "cluster",
                                   cols = NULL, filename = "final_umap",
                                   width = 7.5, height = 6) {
  p <- Seurat::DimPlot(sc, reduction = "umap", group.by = group_by, label = FALSE,
                       pt.size = 0.1, raster = FALSE, cols = cols) +
    ggplot2::theme(aspect.ratio = 1)
  out_path <- file.path(output_dir, paste0(filename, ".", format))
  if (format == "pdf") {
    grDevices::pdf(out_path, height = height, width = width); print(p); grDevices::dev.off()
  } else {
    ggplot2::ggsave(filename = out_path, plot = p, height = height, width = width)
  }
  invisible(out_path)
}

.scumaper_save_featureplot <- function(sc, gene, output_dir, format,
                                       width = 6.5, height = 6) {
  p <- Seurat::FeaturePlot(sc, features = gene, cols = c("lightgrey", "red"),
                           pt.size = 0.1, raster = FALSE, ncol = 1) +
    ggplot2::theme(aspect.ratio = 1)
  out_path <- file.path(output_dir, paste0(gene, "_umap.", format))
  if (format == "pdf") {
    grDevices::pdf(out_path, height = height, width = width); print(p); grDevices::dev.off()
  } else {
    ggplot2::ggsave(filename = out_path, plot = p, height = height, width = width)
  }
  invisible(out_path)
}

.scumaper_default_marker_plotgenes <- function(species) {
  if (species == "human") {
    c("PTPRC","CD3E","CD79A","MS4A1","MZB1","IGHG1","IGHA1","CD68","CSF3R","KIT","TPSAB1",
      "EPCAM","KRT8","KRT19","PECAM1","COL1A1","ACTA2")
  } else {
    c("Ptprc","Cd3e","Cd79a","Ms4a1","Mzb1","Ighg1","Igha","Cd68","Csf3r","Kit","Mcpt4",
      "Epcam","Krt8","Krt18","Pecam1","Col1a1","Acta2")
  }
}

.scumaper_read_inputs <- function(input_dir, qc_all) {
  input_list <- list.dirs(path = input_dir, recursive = FALSE, full.names = FALSE)
  if (length(input_list) == 0) {
    stop("No sample subfolders found inside input_dir: ", input_dir)
  }
  sc <- NULL
  for (each in input_list) {
    message("Reading sample: ", each)
    temp <- Seurat::Read10X(file.path(input_dir, each))
    temp <- Seurat::CreateSeuratObject(temp, project = each, min.cells = 3, min.features = 200)
    temp <- SeuratObject::RenameCells(temp, new.names = paste(each, rownames(temp@meta.data), sep = "_"))
    temp[["percent.mt"]] <- Seurat::PercentageFeatureSet(temp, pattern = "^MT-|^mt-")
    temp <- subset(temp, percent.mt < qc_all)
    temp <- subset(temp, nFeature_RNA > 200 & nFeature_RNA < 6000)
    sc <- if (is.null(sc)) temp else merge(sc, temp)
    gc()
  }
  sc
}

#' Run scUmaper on raw 10x scRNA-seq data of human or mouse tissues
#'
#' @param input_dir Folder containing per-sample subfolders with 10x output format.
#' @param output_dir Output folder. Default creates "scUmaper_output" under getwd().
#' @param species Either "human" or "mouse".
#' @param output_plot Whether to output UMAP plots.
#' @param output_plot_format Either "png" or "pdf".
#' @param n_features Number of HVGs (default 3000). PCA dims auto-set to round(n_features/100), clamped to [10,40].
#' @param qc_all Mitochondrial percent cutoff for all cells (default 25).
#' @param qc_immune Mitochondrial percent cutoff for immune cells (default 5).
#' @param doublet_reserve Whether to also save the removed doublets as an RDS.
#'
#' @return Invisibly returns the singlet Seurat object; results also saved to output_dir.
#'
#' @examplesIf interactive()
#' run_scumaper(input_dir = "path/to/input", output_dir = "path/to/output")
#'
#' @import dplyr
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import presto
#' @export
run_scumaper <- function(input_dir = NULL, output_dir = NULL, species = "human",
                         output_plot = TRUE, output_plot_format = "png",
                         n_features = 3000, qc_all = 25, qc_immune = 5,
                         doublet_reserve = FALSE) {
  .scumaper_check_presto()

  if (is.null(input_dir)) { message("Error: 'input_dir' is required."); return(invisible(NULL)) }
  if (!dir.exists(input_dir)) { message("Error: invalid 'input_dir'."); return(invisible(NULL)) }
  if (!(species %in% c("human", "mouse"))) {
    message("Error: 'species' must be 'human' or 'mouse'. Use run_custom_scumaper() otherwise.")
    return(invisible(NULL))
  }
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "scUmaper_output")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }
  if (!dir.exists(output_dir)) { message("Error: invalid 'output_dir'."); return(invisible(NULL)) }
  if (!(output_plot_format %in% c("png", "pdf"))) {
    output_plot <- FALSE
    message("Warning: format '", output_plot_format, "' not supported; UMAP plots disabled.")
  }

  .scumaper_section()
  message("Input directory : ", input_dir)
  message("Output directory: ", output_dir)
  message("Species         : ", species)
  message("HVGs            : ", n_features)
  message("QC: non-immune mito% < ", qc_all, "; immune mito% < ", qc_immune)
  if (doublet_reserve) message("Removed doublets will be saved as a Seurat object.")
  if (!output_plot) message("UMAP plots will not be created.")
  .scumaper_section()

  ref_df <- scumaper_marker_genes(species)
  ref_df$mt_threshold <- ifelse(ref_df$cell_types %in% c("B/Plasma_cell", "Myeloid", "T/NK_cell"),
                                qc_immune, qc_all)
  ref_df_indexed <- ref_df
  rownames(ref_df_indexed) <- ref_df_indexed$cell_types

  message("Loading input matrices...")
  sc <- .scumaper_read_inputs(input_dir, qc_all = qc_all)
  message("Input matrices loaded. Cells = ", ncol(sc))

  .scumaper_section()
  dims_main <- .scumaper_calc_dims(n_features)
  message("Step 1/4: Preliminary clustering (dims = ", dims_main, ").")
  sc <- .scumaper_preprocess(sc, n_features = n_features, dims = dims_main)
  marker_list <- .scumaper_marker_list_from_df(ref_df)
  sc <- .scumaper_iterative_annotate(sc, marker_list = marker_list)

  .scumaper_section()
  message("Step 2/4: Splitting major groups for doublet detection.")
  .scumaper_split_save_groups(sc, output_dir)
  rm(sc); gc()

  message("Step 3/4: Per-group doublet detection.")
  result <- .scumaper_loop_groups(output_dir, marker_list = marker_list,
                                  mode = "default", ref_df_indexed = ref_df_indexed,
                                  n_features_sub = 2000)
  sc <- result$singlet; dbl <- result$doublet; rm(result); gc()

  if (!is.null(dbl)) {
    dbl <- SeuratObject::JoinLayers(dbl)
    message("Total doublets removed: ", ncol(dbl))
    if (doublet_reserve) {
      saveRDS(dbl, file.path(output_dir, "scUmaper_removed_doublet_seurat.rds"))
      message("Removed doublets saved at: ", file.path(output_dir, "scUmaper_removed_doublet_seurat.rds"))
    }
  } else {
    message("No doublets were identified or removed.")
  }
  rm(dbl); gc()

  .scumaper_section()
  message("Step 4/4: Final reclustering and annotation of singlets (dims = ", dims_main, ").")
  sc <- .scumaper_preprocess(sc, n_features = n_features, dims = dims_main)
  sc <- .scumaper_iterative_annotate(sc, marker_list = marker_list)

  type_levels <- c("T/NK_cell", "B/Plasma_cell", "Myeloid", "Fibroblast",
                   "Endothelial_cell", "Epithelial_cell", "undefined")
  sc$cluster <- factor(sc$cluster, levels = type_levels)
  Seurat::Idents(sc) <- sc$cluster

  .scumaper_cleanup_tmp(output_dir)
  saveRDS(sc, file.path(output_dir, "scUmaper_annotated_singlet_seurat.rds"))
  message("Singlet Seurat object retained: ", ncol(sc), " cells.")
  message("Saved to: ", file.path(output_dir, "scUmaper_annotated_singlet_seurat.rds"))

  if (output_plot) {
    plot_dir <- file.path(output_dir, "markers_plot")
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
    .scumaper_save_dimplot(sc, output_dir, output_plot_format, cols = .scumaper_default_palette())
    plotgenes <- .scumaper_default_marker_plotgenes(species)
    for (g in plotgenes) {
      if (g %in% rownames(sc)) .scumaper_save_featureplot(sc, g, plot_dir, output_plot_format)
    }
    message("UMAP plots saved at: ", plot_dir)
  }

  .scumaper_section()
  message("scUmaper run completed. Thank you for using scUmaper.")
  cite_scumaper()
  invisible(sc)
}

#' Run scUmaper with user-defined marker libraries (custom mode)
#'
#' @param input_seurat A pre-built Seurat object as input.
#' @param input_genelist Named list of marker gene character vectors per cell type.
#' @param output_dir Output folder. Default creates "scUmaper_output" under getwd().
#' @param output_plot Whether to output UMAP plots.
#' @param output_plot_format Either "png" or "pdf".
#' @param n_features Number of HVGs (default 3000). PCA dims auto-set to round(n_features/100), clamped to [10,40].
#' @param doublet_reserve Whether to save the removed doublets as an RDS.
#'
#' @return Invisibly returns the singlet Seurat object.
#'
#' @examplesIf interactive()
#' run_custom_scumaper(input_seurat = scRNA,
#'                     input_genelist = list(T_cell = c("CD3D","CD3E"),
#'                                           Plasma  = c("MZB1","JCHAIN")))
#'
#' @import dplyr
#' @import Seurat
#' @import SeuratObject
#' @import ggplot2
#' @import presto
#' @export
run_custom_scumaper <- function(input_seurat = NULL, input_genelist = NULL,
                                output_dir = NULL, output_plot = TRUE,
                                output_plot_format = "png", n_features = 3000,
                                doublet_reserve = FALSE) {
  .scumaper_check_presto()

  if (is.null(input_seurat) || !inherits(input_seurat, "Seurat")) {
    message("Error: 'input_seurat' must be a Seurat object."); return(invisible(NULL))
  }
  if (is.null(input_genelist) || !is.list(input_genelist) || is.null(names(input_genelist))) {
    message("Error: 'input_genelist' must be a named list of marker gene vectors.")
    return(invisible(NULL))
  }
  if (is.null(output_dir)) {
    output_dir <- file.path(getwd(), "scUmaper_output")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }
  if (!dir.exists(output_dir)) { message("Error: invalid 'output_dir'."); return(invisible(NULL)) }
  if (!(output_plot_format %in% c("png", "pdf"))) {
    output_plot <- FALSE
    message("Warning: format '", output_plot_format, "' not supported; UMAP plots disabled.")
  }

  .scumaper_section()
  message("Output directory  : ", output_dir)
  message("HVGs              : ", n_features)
  message("Cell types defined: ", paste(names(input_genelist), collapse = ", "))
  if (doublet_reserve) message("Removed doublets will be saved as a Seurat object.")
  if (!output_plot) message("UMAP plots will not be created.")
  .scumaper_section()

  dims_main <- .scumaper_calc_dims(n_features)
  message("Step 1/4: Preliminary clustering (dims = ", dims_main, ").")
  sc <- input_seurat
  sc <- SeuratObject::JoinLayers(sc)
  sc <- subset(sc, nFeature_RNA > 200 & nFeature_RNA < 6000)
  sc <- .scumaper_preprocess(sc, n_features = n_features, dims = dims_main, exclude_pattern = NULL)

  min_types <- max(1L, ceiling(length(input_genelist) * 0.7))
  sc <- .scumaper_iterative_annotate(sc, marker_list = input_genelist, min_types = min_types)

  .scumaper_section()
  message("Step 2/4: Splitting major groups for doublet detection.")
  .scumaper_split_save_groups(sc, output_dir)
  rm(sc); gc()

  message("Step 3/4: Per-group doublet detection.")
  result <- .scumaper_loop_groups(output_dir, marker_list = input_genelist,
                                  mode = "custom", n_features_sub = 2000)
  sc <- result$singlet; dbl <- result$doublet; rm(result); gc()

  if (!is.null(dbl)) {
    dbl <- SeuratObject::JoinLayers(dbl)
    message("Total doublets removed: ", ncol(dbl))
    if (doublet_reserve) {
      saveRDS(dbl, file.path(output_dir, "scUmaper_removed_doublet_seurat.rds"))
      message("Removed doublets saved at: ", file.path(output_dir, "scUmaper_removed_doublet_seurat.rds"))
    }
  } else {
    message("No doublets were identified or removed.")
  }
  rm(dbl); gc()

  .scumaper_section()
  message("Step 4/4: Final reclustering and annotation of singlets (dims = ", dims_main, ").")
  sc <- .scumaper_preprocess(sc, n_features = n_features, dims = dims_main, exclude_pattern = NULL)
  sc <- .scumaper_iterative_annotate(sc, marker_list = input_genelist, min_types = min_types)

  sc$cluster <- factor(sc$cluster, levels = c(names(input_genelist), "undefined"))
  Seurat::Idents(sc) <- sc$cluster

  .scumaper_cleanup_tmp(output_dir)
  saveRDS(sc, file.path(output_dir, "scUmaper_annotated_singlet_seurat.rds"))
  message("Singlet Seurat object retained: ", ncol(sc), " cells.")
  message("Saved to: ", file.path(output_dir, "scUmaper_annotated_singlet_seurat.rds"))

  if (output_plot) {
    n_types_total <- length(input_genelist) + 1L
    palette <- if (length(input_genelist) <= 58) .scumaper_custom_palette()[seq_len(min(n_types_total, 59L))] else NULL
    .scumaper_save_dimplot(sc, output_dir, output_plot_format, cols = palette)
    message("Final UMAP saved at: ", output_dir)
  }

  .scumaper_section()
  message("scUmaper (custom) run completed. Thank you for using scUmaper.")
  cite_scumaper()
  invisible(sc)
}
