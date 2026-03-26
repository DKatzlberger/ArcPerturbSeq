## max. code width ============================================================



.script_startup <- function(){

  t0 <- Sys.time()

  # renv
  renv::load()
  t1 <- Sys.time()

  # Libraries
  suppressPackageStartupMessages({
    library(data.table, quietly = TRUE)
    library(anndataR, quietly = TRUE)
    library(SingleCellExperiment, quietly = TRUE)
    library(dream, quietly = TRUE)
    library(edgeR, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    library(variancePartition, quietly = TRUE)
    library(scater, quietly = TRUE)
    library(scran, quietly = TRUE)
    library(patchwork, quietly = TRUE)
  })

  t2 <- Sys.time()
  
  # Report
  cat(sprintf(
    "Startup completed in %.2f seconds\n",
    as.numeric(difftime(t2, t0, units = "secs"))
  ))

  cat(sprintf(
    "- renv activation:    %.2f seconds\n",
    as.numeric(difftime(t1, t0, units = "secs"))
  ))

  cat(sprintf(
    "- library attachment: %.2f seconds\n",
    as.numeric(difftime(t2, t1, units = "secs"))
  ))

  # Return
  return(invisible(TRUE))
}



.load_pseudobulk <- function(
  base_dir, 
  files,
  tissue,
  x_mapping = "counts",
  verbose = TRUE
) {

  stopifnot(
    is.character(files),
    !is.null(names(files)),
    is.character(tissue),
    length(tissue) == 1
  )

  # Load anndata obj. and convert to SCE obj.
  sce_list <- lapply(names(files), function(ct) {
    ad <- read_h5ad(file.path(base_dir, files[[ct]]))
    ad$obs$cell_type <- ct
    ad$as_SingleCellExperiment(x_mapping = x_mapping)
  })

  # Add celltype and tissue
  names(sce_list) <- names(files)
  sce <- do.call(cbind, sce_list)
  colData(sce)$tissue <- tissue

  # Verbose
  if (verbose) {
    cat(sprintf("%s:\n", tissue))

    ct <- colData(sce)$cell_type
    split_cols <- split(seq_len(ncol(sce)), ct)

    for (type in names(split_cols)) {
      idx <- split_cols[[type]]
      dim <- sprintf("%d x %d (sample x feature)",length(idx), nrow(sce))
      cat(sprintf("- %s %s\n", type, dim))
    }
  }

  # Return
  return(sce)
}



.combine_pseudobulk <- function(
  specs,
  exp_name,
  x_mapping = "counts",
  verbose = TRUE,
  output = c("combined", "list")
) {

  output <- match.arg(output)

  stopifnot(
    is.list(specs), length(specs) > 0,
    is.character(exp_name), length(exp_name) == 1
  )

  # Load per-tissue SCEs
  sce_list <- lapply(names(specs), function(tissue) {
    spec <- specs[[tissue]]

    .load_pseudobulk(
      base_dir  = spec$base_dir,
      files     = spec$files,
      tissue    = tissue,
      x_mapping = x_mapping,
      verbose   = FALSE
    )
  })
  names(sce_list) <- names(specs)

  # Helper: shared verbose printer
  .print_summary <- function(header, sce_list, combined_sce = NULL) {
    cat(sprintf("\n%s | %s:\n", header, exp_name))

    for (t in names(sce_list)) {
      cat(sprintf("• %s:\n", t))
      sce_t <- sce_list[[t]]

      ct_tab <- table(colData(sce_t)$cell_type)

      for (ct in names(ct_tab)) {
        idx_n <- ct_tab[[ct]]
        n_feat <- nrow(sce_t)

        dim_str <- sprintf(
          "%d x %d (sample x feature)",
          idx_n,
          n_feat
        )

        cat(sprintf("  - %s %s\n", ct, dim_str))
      }
    }
  }

  # Return list
  if (output == "list") {
    if (verbose) {
      .print_summary(
        header = "Pseudobulk SCE list",
        sce_list = sce_list
      )
    }
    return(sce_list)
  }

  # Combine all tissues
  sce <- do.call(cbind, sce_list)
  colData(sce)$experiment <- exp_name

  # Combined verbose
  if (verbose) {
    .print_summary(
      header = "Pseudobulk SCE object",
      sce_list = sce_list,
      combined_sce = sce
    )
  }

  # Return
  return(sce)
}



.normalize_pseudobulk <- function(sce) {
  
  stopifnot(inherits(sce, "SingleCellExperiment"))

  # scuttle logNormCounts
  sce <- scuttle::logNormCounts(sce)

  # Return
  return(sce)
}



.plot_mean_expr <- function(
  sce,
  assay_name
) {

  stopifnot(
    assay_name %in% assayNames(sce)
  )

  logc  <- assay(sce, assay_name)
  yvals <- Matrix::colMeans(logc)

  df <- data.frame(
    sample = colnames(sce),
    cell_type = colData(sce)$cell_type,
    tissue = colData(sce)$tissue,
    value = yvals
  )

  # ggplot
  gg <- ggplot(
    data = df, 
    mapping = aes(
      x = cell_type, 
      y = value
      )
    ) +
    geom_violin(
      mapping = aes(
        fill = cell_type
      ),
      trim = FALSE, 
      scale = "width", 
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      linewidth = 0.1,
      outlier.shape = NA
    ) +
    scale_fill_manual(
      values = .celltype_pal
    ) +
    labs(
      x = "Cell type",
      y = sprintf("Mean %s", assay_name),
      fill = NULL
    ) +
    .theme_DK(
      rotate = 45,
      small_lgd = TRUE
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical"
    )
  
  # Return
  return(gg)
}


.plot_mean_expr <- function(
  sce,
  assay_name
) {

  stopifnot(
    assay_name %in% assayNames(sce)
  )

  logc  <- assay(sce, assay_name)
  yvals <- Matrix::colMeans(logc)

  df <- data.frame(
    sample = colnames(sce),
    cell_type = colData(sce)$cell_type,
    tissue = colData(sce)$tissue,
    value = yvals
  )

  # ggplot
  gg <- ggplot(
    data = df, 
    mapping = aes(
      x = cell_type, 
      y = value
      )
    ) +
    geom_violin(
      mapping = aes(
        fill = cell_type
      ),
      trim = FALSE, 
      scale = "width", 
      alpha = 0.5,
      show.legend = FALSE
    ) +
    geom_boxplot(
      width = 0.15,
      linewidth = 0.1,
      outlier.shape = NA
    ) +
    scale_fill_manual(
      values = .celltype_pal
    ) +
    labs(
      x = "Cell type",
      y = sprintf("Mean %s", assay_name),
      fill = NULL
    ) +
    .theme_DK(
      rotate = 45,
      small_lgd = TRUE
    ) +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical"
    )
  
  # Return
  return(gg)
}



.plot_smpl_density <- function(
  sce,
  assay_name,
  n_samples = 10
) {

  stopifnot(assay_name %in% assayNames(sce))

  mat <- assay(sce, assay_name)

  n_samples <- min(n_samples, ncol(mat))
  mat <- mat[, seq_len(n_samples), drop = FALSE]

  df <- data.frame(
    value  = as.vector(mat),
    sample = rep(colnames(mat), each = nrow(mat))
  )

  # ggplot
  gg <- ggplot(
    data = df,
    mapping = aes(
      x = value, 
      color = sample
    )
  ) +
  geom_density(
    linewidth = 0.1,
    show.legend = FALSE
  ) +
  labs(
    x = assay_name,
    y = "Density"
  ) +
  .theme_DK()

  # Return
  return(gg)
}