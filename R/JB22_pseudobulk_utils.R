## max. code width ============================================================

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
  output = c("sce", "list")
) {

  output <- match.arg(output)

  stopifnot(
    is.list(specs), length(specs) > 0,
    is.character(exp_name), length(exp_name) == 1
  )

  # Load per-tissue SCEs
  sce_list <- lapply(names(specs), function(tissue) {
    spec <- specs[[tissue]]

    sce_obj = .load_pseudobulk(
      base_dir = spec$base_dir,
      files = spec$files,
      tissue = tissue,
      x_mapping = x_mapping,
      verbose = FALSE
    )
    
    # Add exp. description
    colData(sce_obj)$experiment = exp_name
    sce_obj
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

    # List sce verbose
    if (verbose) {
      .print_summary(
        header = "Pseudobulk SCE list",
        sce_list = sce_list
      )
    }

    # Return
    return(sce_list)

  } else {

    # Combine all tissues
    sce <- do.call(cbind, sce_list)

    # Combined sce verbose
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
    value = as.vector(mat),
    sample = rep(
      colnames(mat), 
      each = nrow(mat)
    )
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

.plot_PCA <- function(
  sce
){

  # Coordinates
  df <- data.frame(
    PC1 = reducedDim(sce, "PCA")[,1],
    PC2 = reducedDim(sce, "PCA")[,2],
    tissue = colData(sce)$tissue,
    cell_type = colData(sce)$cell_type
  )

  # Plot
  gg <- ggplot(
    data = df, 
    mapping = aes(
      x = PC1, 
      y = PC2,
      color = cell_type,
      shape = tissue
      )
    ) +
    geom_point(
      size = 0.5
    ) +
    scale_color_manual(
      values = .celltype_pal
    ) +
    labs(
      color = "Cell type",
      shape = "Tissue"
    ) +
    .theme_DK(
      small_lgd = TRUE
    )

  # Return
  return(gg)
}

.prepare_metadata <- function(
  sce_obj,
  formula,
  ref_level = NULL,
  shuffle = FALSE,
  verbose = TRUE
) {

  meta <- as.data.frame(colData(sce_obj))

  # First variable = covariate of interest
  vars <- all.vars(formula)
  covariate <- vars[1]

  missing <- setdiff(vars, colnames(meta))
  if (length(missing) > 0) {
    stop("Missing variables in colData: ", paste(missing, collapse = ", "))
  }

  # Helper: clean factor levels
  .clean_factor_levels <- function(x) {
    levels(x) <- make.names(levels(x))
    x
  }

  # Ensure vars are factors + clean coef names
  for (v in vars) {
    if (!is.factor(meta[[v]])) {
      meta[[v]] <- factor(meta[[v]])
    }
    meta[[v]] <- .clean_factor_levels(meta[[v]])
  }

  # Handle reference level
  if (!is.null(ref_level)) {

    ref_level_clean <- make.names(ref_level)

    # Only inform if it actually changed
    if (!identical(ref_level, ref_level_clean)) {
      cat(sprintf(
        " - Reference level cleaned: %s → %s\n",
        ref_level, ref_level_clean
      ))
    }

    if (!ref_level_clean %in% levels(meta[[covariate]])) {
      stop(
        "Reference level not found after cleaning in ",
        covariate, ": ",
        ref_level, " → ", ref_level_clean
      )
    }

    meta[[covariate]] <- relevel(
      droplevels(meta[[covariate]]),
      ref = ref_level_clean
    )

    ref_level <- ref_level_clean
  }

  # Conditional permutation of the covariate of interest
  #
  # If `shuffle = TRUE`, the covariate of interest (the first variable in the
  # model formula) is randomly permuted across samples in order to destroy any
  # real association between gene expression and that covariate. This is used
  # to generate null data for benchmarking or permutation-based calibration.
  #
  # When additional variables are present in the model formula (e.g. batch, sex,
  # subject), permutations are performed *within strata defined by those
  # variables*. This is known as a conditional or restricted permutation.
  #
  # Example: formula = ~ condition + sex + batch
  #
  # In this case, `condition` will be shuffled only within groups defined by the
  # combination of `(sex, batch)`. This preserves the design structure while
  # removing the association between expression and the covariate of interest.
  #
  # This procedure ensures the permutation corresponds to the correct null
  # hypothesis tested by the model:
  #
  #     expression ⟂ covariate | nuisance variables
  #
  # In other words, gene expression is independent of the covariate of interest
  # conditional on the remaining variables in the design.
  #
  # If the model contains only the covariate of interest, a global permutation
  # across all samples is performed.
  #
  # Note:
  # Blocks that contain only a single level of the covariate cannot be permuted
  # and will remain unchanged. This is expected and preserves statistical
  # validity.
  if (shuffle) {
    shuffle_blocks <- setdiff(vars, covariate)

    if (length(shuffle_blocks) == 0) {
      perm <- sample.int(nrow(meta))
      meta[[covariate]] <- meta[[covariate]][perm]

    } else {
      blocks <- interaction(meta[shuffle_blocks], drop = TRUE)
      meta[[covariate]] <- ave(
        meta[[covariate]],
        blocks,
        FUN = sample
      )
    }
  }

  # meta[[covariate]] <- relevel(
  #   droplevels(meta[[covariate]]),
  #   ref = ref_level
  # )
  colData(sce_obj) <- S4Vectors::DataFrame(meta)

  # Verbose
  if (verbose){
    cat("Model design\n")
    cat(sprintf(" - Formula: %s\n", deparse(formula)))
    cat(sprintf(" - Intercept: %s\n", if (!is.null(ref_level)) ref_level else levels(meta[[covariate]])[1]))
    cat(sprintf(" - Nr. coefs: %d\n", ncol(model.matrix(formula, meta))))
    cat(sprintf(" - Shuffled: %s\n", if (shuffle) "TRUE" else "FALSE"))
  }

  # Return
  return(sce_obj)
}

.edgeR_filterByExpr <- function(
  sce_obj,
  formula
){

  meta   <- as.data.frame(colData(sce_obj))
  design <- model.matrix(formula, meta)
  counts <- assay(sce_obj, "counts")

  # edgeR: filter approach
  dge  <- edgeR::DGEList(counts)

  keep <- edgeR::filterByExpr(dge, design = design)
  sce_obj <- sce_obj[keep, , drop = FALSE]

  # Numbers for reporting
  n_total <- nrow(counts)
  n_keep  <- sum(keep)

  # Verbose
  cat(sprintf(" - Nr. features: %d (%d)\n", n_keep, n_total))

  # Return
  return(
    list(
      sce_obj = sce_obj,
      method = "edgeR_filterByExpr"
    )
  )
}

.dge_fit_limma_voom <- function(
  sce_obj,
  formula
){

  meta <- as.data.frame(colData(sce_obj))

  # Model design
  design <- model.matrix(formula, meta)

  # edgeR obj.
  counts <- assay(sce_obj, "counts")

  dge <- edgeR::DGEList(counts)
  dge <- edgeR::calcNormFactors(dge)

  # Verbose
  cat(" - DGE approach: limma_voom\n")

  # limma: fit & shrinkage
  vobj <- limma::voom(dge, design, plot = FALSE)
  fit  <- limma::lmFit(vobj, design)
  fit  <- limma::eBayes(fit)

  # Return
  return(
    list(
      method = "limma_voom",
      formula = formula,
      fit_obj = fit,
      design = design,
      meta = meta
    )
  )
}

.dge_fit_edgeR_QLF <- function(
  sce_obj, 
  formula
) {

  meta <- as.data.frame(colData(sce_obj))

  # Model design
  design <- model.matrix(formula, meta)

  # edgeR obj.
  counts <- assay(sce_obj, "counts")

  dge <- edgeR::DGEList(counts)
  dge <- edgeR::calcNormFactors(dge)

  # Verbose
  cat(" - DGE approach: edgeR_QLF\n")

  # edgeR: fit & shrinkage
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design)

  # Return
  return(
    list(
      method = "edgeR_QLF",
      formula = formula,
      fit_obj = fit,
      design = design,
      meta = meta
    )
  )
}

.dge_fit_DESeq2_Wald <- function(
  sce_obj,
  formula
) {

  meta <- as.data.frame(colData(sce_obj))

  # Counts
  counts <- assay(sce_obj, "counts")

  # DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = formula
  )

  # Verbose
  cat(" - DGE approach: DESeq_Wald\n")

  # DESeq2: fit & shrink
  fit_obj <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
  design  <- model.matrix(DESeq2::design(fit_obj), colData(fit_obj))

  # Return
  return(
    list(
      method = "DESeq2_Wald",
      formula = formula,
      fit_obj = fit_obj,
      design = design,
      meta = meta
    )
  )
}

.harmonize_tt <- function(tt) {

  rename_map <- c(
    log2FoldChange = "logFC",
    P.Value = "pval",
    PValue = "pval",
    pvalue = "pval",
    adj.P.Val = "padj",
    FDR = "padj",
    t = "stat",
    F = "stat"
  )

  for (old in names(rename_map)) {
    if (old %in% names(tt)) {
      names(tt)[names(tt) == old] <- rename_map[[old]]
    }
  }

  return(tt)
}

.dge_extract <- function(x) {

  fit    <- x$fit_obj
  method <- x$method

  # Coefs from fit obj.
  coefs <- switch(
    method,
    limma_voom = colnames(fit$coefficients),
    edgeR_QLF = colnames(fit$coefficients),
    DESeq2_Wald = DESeq2::resultsNames(fit)
  )

  # Remove intercept if present
  coefs <- coefs[!grepl("Intercept", coefs)]

  if (length(coefs) == 0) {
    stop("No coefficients found in fit object")
  }

  n_tests <- nrow(fit) * length(coefs)

  # Verbose
  cat("\nExtracted results\n")
  cat(sprintf(" - Nr. coefs: %-8d (without intercept)\n", length(coefs)))
  cat(sprintf(" - Nr. tests: %-8d (features × coefs)\n", n_tests))

  # Function per coef
  get_tt <- function(cf) {

    tt <- switch(
      method,

      # limma
      limma_voom = limma::topTable(
        fit,
        coef = cf,
        number = Inf,
        sort.by = "none"
      ),

      # edgeR
      edgeR_QLF = {
        qlf <- edgeR::glmQLFTest(
          fit,
          coef = match(cf, colnames(fit$coefficients))
        )
        edgeR::topTags(qlf, n = Inf)$table
      },

      # DESeq2
      DESeq2_Wald = as.data.frame(
        DESeq2::results(fit, name = cf)
      )
    )

    tt$feature <- rownames(tt)
    tt$coef <- cf
    rownames(tt) <- NULL
    tt
  }

  # Harmonize
  res <- do.call(rbind, lapply(coefs, get_tt))
  res <- .harmonize_tt(res)

  # Return
  return(res)
}

.run_dge <- function(
  sce_obj,
  formula,
  ref_level = NULL,
  shuffle = FALSE,
  filter_fun = NULL,
  dge_fun
) {

  if(!is.function(dge_fun)){
    stop("dge_fun must a function")
  }

  if (!is.null(filter_fun) && !is.function(filter_fun)) {
    stop("filter_fun must be NULL or a function")
  }

  # Prepare meta + counts
  sce_obj <- .prepare_metadata(
    sce_obj = sce_obj,
    formula = formula,
    ref_level = ref_level,
    shuffle = shuffle
  )

  # Filter features
  if (is.null(filter_fun)) {
    cat("Feature selection\n")
    cat(sprintf(
      " - Nr. features: %d (no filter)\n", nrow(sce_obj)
    ))
  } else {
    filter_obj <- filter_fun(
      sce_obj = sce_obj,
      formula = formula
    )
    sce_obj <- filter_obj$sce_obj
  }

  # Model fit
  fit <- dge_fun(
    sce_obj = sce_obj,
    formula = formula
  )

  # Extract coef as tt
  tt_obj <- .dge_extract(fit)
  tt_obj$method    <- fit$method
  tt_obj$formula   <- deparse(fit$formula)
  tt_obj$ref_level <- ref_level
  tt_obj$shuffle   <- shuffle

  # Return
  return(
    list(
      method = fit$method,
      formula = fit$formula,
      ref_level = ref_level,
      meta = fit$meta,
      filter_fun = if (is.null(filter_fun)) "none" else filter_obj$method,
      design = fit$design,
      fit_obj = fit$fit_obj,
      tt_obj = tt_obj
    )
  )
}

.safe_run <- function(...) {
  tryCatch(
    {
      res <- .run_dge(...)
      res$success <- TRUE
      res
    },
    error = function(e) {
      list(
        success = FALSE,
        error   = conditionMessage(e)
      )
    }
  )
}