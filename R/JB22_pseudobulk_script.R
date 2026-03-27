## max. code width ============================================================


## Experimental description ===================================================

## Title: Arc Institute in vivo pertubation screens
## Author: dkatzlberger
## Date: 21.01.2026

## Sample            (biological replicate)
## └── Genotype      (knockout)
##     └── Guide     (technical replication of knockout)
##         └── Cells (cells per guide)

## Single cell pseudobulks per tissue
## Tissues: Peritoneum, Liver, Spleen

## Available cell types per tissue
## Peritoneum: Macrophages
## Liver:      Bcells, Monocytes, Neutrophiles
## Spleen:     Bcells


## Libraries =================================================================

# Scource R scripts
source("R/utils.R")
source("R/JB22_pseudobulk_utils.R")

# Activate renv + load libs.
.script_startup()

## Directories ================================================================

# Define paths
res_path <- file.path("results", "JB22_pseudobulk")
fig_path <- file.path(res_path, "fig")
dge_path <- file.path(res_path, "dge")
tab_path <- file.path(res_path, "tab")

# Create output dirs
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)
dir.create(dge_path, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_path, recursive = TRUE, showWarnings = FALSE)

# Shared directory (location: vscratch)
JB22_pb_dir  <- "/vscratch/wes/arc/share"

## Color schemes ==============================================================

.celltype_pal <- c(
  Bcells = "#fbb4ae",
  Macrophages = "#b3cde3",
  Neutrophils = "#ccebc5"
)

## Dataset ====================================================================

# Pseudobulk specs
JBB_pb_specs <- list(
  Peritoneum = list(
    base_dir = file.path(JB22_pb_dir, "JB22_PM"),
    files = c(
      Macrophages = "pseudobulk_Macrophages.h5ad"
    )
  ),
  Liver = list(
    base_dir = file.path(JB22_pb_dir, "JB22_Liver"),
    files = c(
      Bcells = "pseudobulk_B_cells.h5ad",
      Macrophages = "pseudobulk_Monocytes.h5ad",
      Neutrophils = "pseudobulk_Neutrophils.h5ad"
    )
  ),
  Spleen = list(
    base_dir = file.path(JB22_pb_dir, "JB22_Spleen"),
    files = c(
      Bcells = "pseudobulk_B_cells.h5ad"
    )
  )
)

# Experiment JB22 | SingleCellExperiment obj.
JB22_pb_sce_list <- .combine_pseudobulk(
  JBB_pb_specs, 
  exp_name = "JB22_pseudobulk",
  output = "list"
)

## Normalization ==============================================================

# Normalize
JB22_pb_norm_sce_list <- lapply(JB22_pb_sce_list, .normalize_pseudobulk)

# Plot: mean expr
JB22_mean_expr_p <- lapply(
  names(JB22_pb_norm_sce_list), 
  function(tissue) {
    .plot_mean_expr(
      sce = JB22_pb_norm_sce_list[[tissue]],
      assay_name = "logcounts"
    ) + 
    ggtitle(tissue)
})
JB22_mean_expr_p <- wrap_plots(JB22_mean_expr_p, nrow = 1)

# Save the plot
ggsave(
  filename = file.path(fig_path, "mean_expr.svg"),
  plot = JB22_mean_expr_p,
  width = 10,
  height = 5,
  units = "cm"
)

# Plot: sample density
JB22_smpl_density_p <- lapply(
  names(JB22_pb_norm_sce_list),
  function(tissue) {
    .plot_smpl_density(
      sce = JB22_pb_norm_sce_list[[tissue]],
      assay_name = "logcounts",
      n_samples = 50
    ) +
    ggtitle(tissue)
})
JB22_smpl_density_p <- wrap_plots(JB22_smpl_density_p, nrow = 1)

# Save the plot
ggsave(
  filename = file.path(fig_path, "smpl_density.svg"),
  plot = JB22_smpl_density_p,
  width = 10,
  height = 5,
  units = "cm"
)

## Cluster ====================================================================
### fun. ----------------------------------------------------------------------

# PCA: PC1 vs PC2
.plot_PCA <- function(
  sce
){

  # Coordinates
  df <- data.frame(
  PC1 = reducedDim(JB22_pb_sce, "PCA")[,1],
  PC2 = reducedDim(JB22_pb_sce, "PCA")[,2],
  tissue = colData(JB22_pb_sce)$tissue,
  cell_type = colData(JB22_pb_sce)$cell_type
)

# Plot
p <- ggplot(
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
}

### JB22 clusters -------------------------------------------------------------

JB22_pb_sce <- do.call(cbind, JB22_pb_norm_sce_list)

# PCA
JB22_pb_sce <- runPCA(
  JB22_pb_sce,
  assay.type = "logcounts",
  ncomponents = 10,
  subset_row = getTopHVGs(
    modelGeneVar(
      JB22_pb_sce, 
      assay.type = "logcounts"
    ),
    n = 2000
  )
)

# Plot: PC1 vs PC2
JB22_PCA_p <- .plot_PCA(
  sce = sce
)

# Save the plot
ggsave(
  filename = file.path(fig_path, "PCA.svg"),
  plot = JB22_PCA_p,
  width = 10,
  height = 5,
  units = "cm"
)

## DGE analysis ===============================================================
### fun. ----------------------------------------------------------------------
#### preparation --------------------------------------------------------------

.prepare_metadata <- function(
  sce_obj,
  formula,
  ref_level = NULL,
  shuffle = FALSE
) {

  meta <- as.data.frame(colData(sce_obj))

  # First variable = covariate of interest
  vars <- all.vars(formula)
  covariate <- vars[1]

  missing <- setdiff(vars, colnames(meta))
  if (length(missing) > 0) {
    stop("Missing variables in colData: ", paste(missing, collapse = ", "))
  }

  # Ensure vars are factors
  for (v in vars) {
    if (!is.factor(meta[[v]])) {
      meta[[v]] <- factor(meta[[v]])
    }
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

  # Select reference
  if (!is.null(ref_level)) {
    if (!ref_level %in% levels(meta[[covariate]])) {
      stop(
        "Reference level not found in ",
        covariate, ": ",
        ref_level
      )
    }
    meta[[covariate]] <- relevel(
      droplevels(meta[[covariate]]),
      ref = ref_level
    )
  }
  colData(sce_obj) <- S4Vectors::DataFrame(meta)

  # Verbose
  cat("Model design\n")
  cat(sprintf(" - Formula:   %s\n", deparse(formula)))
  cat(sprintf(" - Intercept: %s\n", if (!is.null(ref_level)) ref_level else levels(meta[[covariate]])[1]))
  cat(sprintf(" - Nr. coefs: %d\n", ncol(model.matrix(formula, meta))))
  cat(sprintf(" - Shuffled:  %s\n", if (shuffle) "TRUE" else "FALSE"))

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
  dge  <- edgeR::calcNormFactors(dge)

  keep <- edgeR::filterByExpr(dge, design = design)
  sce_obj <- sce_obj[keep, , drop = FALSE]

  # Numbers for reporting
  n_total <- nrow(counts)
  n_keep  <- sum(keep)

  # Verbose
  cat("Feature selection\n")
  cat(sprintf(" - Nr. features: %d (%d)\n", n_keep, n_total))

  # Return
  return(
    list(
      sce_obj = sce_obj,
      method = "edgeR_filterByExpr"
    )
  )
}

#### dge method ---------------------------------------------------------------

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
  cat("Model fit\n")
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
  cat("Model fit\n")
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
  cat("Model fit\n")
  cat(" - DGE approach: DESeq_Wald\n")

  # DESeq2: fit & shrink
  fit_obj <- DESeq2::DESeq(dds, test = "Wald", quiet = TRUE)
  design  <- stats::model.matrix(design(fit_obj), colData(fit_obj))

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

#### tt extraction method -----------------------------------------------------

.get_coefs <- function(
  formula, 
  design, 
  n_features
) {

  covariate <- all.vars(formula)[1]
  coefs <- colnames(design)[grepl(paste0("^", covariate), colnames(design))]

  if (length(coefs) == 0) {
    stop("No coefficients found for covariate: ", covariate)
  }

  n_tests <- n_features * length(coefs)

  # Verbose
  cat("Extracted results\n")
  cat(sprintf(" - Nr. coefs: %-8d (without intercept)\n", length(coefs)))
  cat(sprintf(" - Nr. tests: %-8d (features × coefs)\n", n_tests))

  # Return
  return(coefs)
}

.harmonize_tt <- function(tt) {

  rename_map <- c(
    log2FoldChange = "logFC",
    P.Value        = "pval",
    PValue         = "pval",
    pvalue         = "pval",
    adj.P.Val      = "padj",
    FDR            = "padj",
    t              = "stat",
    F              = "stat"
  )

  common <- intersect(names(rename_map), names(tt))
  names(tt)[match(common, names(tt))] <- rename_map[common]

  tt
}

.dge_extract <- function(x) {

  fit     <- x$fit_obj
  design  <- x$design
  formula <- x$formula
  method  <- x$method

  coefs <- .get_coefs(formula, design, nrow(fit))

  get_tt <- function(cf) {

    tt <- switch(
      method,

      limma_voom =
        limma::topTable(fit, coef = cf, number = Inf, sort.by = "none"),

      edgeR_QLF = {
        qlf <- edgeR::glmQLFTest(fit, coef = match(cf, colnames(design)))
        edgeR::topTags(qlf, n = Inf)$table
      },

      DESeq2_Wald =
        as.data.frame(DESeq2::results(fit, name = cf))
    )

    tt$feature <- rownames(tt)
    tt$coef <- cf
    rownames(tt) <- NULL

    .harmonize_tt(tt)
  }

  do.call(rbind, lapply(coefs, get_tt))
}

#### core runner --------------------------------------------------------------

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

  # Prepare meta
  sce_obj <- .prepare_metadata(
    sce_obj = sce_obj,
    formula = formula,
    ref_level = ref_level,
    shuffle = shuffle
  )

  # Filter
  if (is.null(filter_fun)) {
    cat("Feature selection\n")
    cat(sprintf(" - Nr. features: %d (no filter)\n", nrow(sce_obj)))
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
      fit_obj = fit$fit_obj,
      design = fit$design,
      meta = fit$meta,
      ref_level = ref_level,
      filter_fun = if (is.null(filter_fun)) "none" else filter_obj$method,
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

#### pipelines ----------------------------------------------------------------

tissue_cell_type_pipeline <- function(
  sce_obj,
  tissue,
  cell_type,
  formula,
  ref_level,
  filter_fun,
  dge_fun,
  out_path
) {

  outfile <- paste0(out_path, ".rds")

  if (file.exists(outfile)) {
    return(invisible(NULL))
  }

  cat("------------------------------------------------------------------------\n")
  cat(sprintf("Run:     %s  %s\n", tissue, cell_type))
  cat(sprintf("Out dir: %s\n", outfile))

  # Filter for cell type
  meta   <- as.data.frame(colData(sce_obj))
  ct_sce <- sce_obj[, meta$cell_type == cell_type]

  if (ncol(ct_sce) == 0) {
    stop("No cells found for cell_type: ", cell_type)
  }

  # Output is a list obj.
  res <- .safe_run(
    sce_obj = ct_sce,
    formula = formula,
    ref_level = ref_level,
    filter_fun = filter_fun,
    dge_fun = dge_fun
  )

  res$tt_obj$tissue    <- tissue
  res$tt_obj$cell_type <- cell_type
  res$timestamp <- Sys.time()

  if (isFALSE(res$success)) {

    error_file <- paste0(out_path, ".error.txt")
    
    cat(
      "TIME: ", Sys.time(), "\n",
      "TISSUE: ", tissue, "\n",
      "CELL_TYPE: ", cell_type, "\n",
      "ERROR:\n", res$error, "\n",
      sep = "",
      file = error_file
    )

    cat(
      sprintf(
        "\nERROR! Check error file:\n%s\n",
        error_file
      )
    )

  } else {
    # Save
    saveRDS(res, outfile)
  }

  # Return
  cat("------------------------------------------------------------------------\n")
  return(NULL)
}

### JB22 limma bulk RNA-seq ---------------------------------------------------
#### Within tissue within cell type -------------------------------------------
dge_methods <- list(
  edgeR_QLF = .dge_fit_edgeR_QLF,
  limma_voom = .dge_fit_limma_voom,
  DESeq2_Wald = .dge_fit_DESeq2_Wald
)

dge_formulas <- list(
  genotype = ~ genotype,
  genotype_sample = ~ genotype + sample
)

dge_bios <- do.call(
  rbind,
  lapply(names(JB22_pb_sce_list), function(tissue) {
    sce <- JB22_pb_sce_list[[tissue]]
    cts <- unique(as.data.frame(colData(sce))$cell_type)
    data.frame(tissue = tissue, cell_type = cts)
  })
)

tasks <- merge(
  dge_bios,
  expand.grid(
    method  = names(dge_methods),
    formula = names(dge_formulas),
    stringsAsFactors = FALSE
  ),
  by = NULL
)
#tasks <- tasks[7,]

# Set the number of used cores
n_cores <- as.integer(Sys.getenv("NSLOTS", unset = 1))
n_workers <- min(n_cores, nrow(tasks))
register(MulticoreParam(n_workers))

m <- bplapply(seq_len(nrow(tasks)), function(i) {

  tissue <- tasks$tissue[i]
  ct     <- tasks$cell_type[i]
  method <- tasks$method[i]
  form   <- tasks$formula[i]

  # Outfile
  out_path <- file.path(
    "results/JB22_pseudobulk/dge",
    paste0(method, "_", form, "_", tissue, "_", ct)
  )

  tissue_cell_type_pipeline(
    sce_obj = JB22_pb_sce_list[[tissue]],
    tissue = tissue,
    cell_type = ct,
    formula = dge_formulas[[form]],
    ref_level = "NTC",
    filter_fun = .edgeR_filterByExpr,
    dge_fun = dge_methods[[method]],
    out_path = out_path
  )
})

# limma_res <- readRDS("results/JB22_pseudobulk/dge/limma_voom_genotype_Liver_Bcells.rds")
# as.data.table(limma_res$tt_obj)
