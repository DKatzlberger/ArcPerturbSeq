## max. code width ============================================================


## Experimental description ===================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Libraries =================================================================

# Scource R scripts
source("R/utils.R")
source("R/JB22_pseudobulk_utils.R")

# Activate renv + load libs.
.script_startup(
  c(
    "anndataR",
    "SingleCellExperiment"
  )
)

## Directories ================================================================

# Define paths
res_path <- file.path("results", "JB22_pseudobulk_orig")
fig_path <- file.path(res_path, "fig")
dge_path <- file.path(res_path, "dge")
tab_path <- file.path(res_path, "tab")

# Create output dirs
dir.create(res_path, recursive = TRUE, showWarnings = FALSE)
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
JB22_pb_norm_sce_list <- lapply(
  JB22_pb_sce_list, 
  .normalize_pseudobulk
)

# ggplot: mean expr
JB22_mean_expr_gg_list <- lapply(
  names(JB22_pb_norm_sce_list), 
  function(tissue) {
    .plot_mean_expr(
      sce = JB22_pb_norm_sce_list[[tissue]],
      assay_name = "logcounts"
    ) + 
    ggtitle(tissue)
})

# Save the plot
ggsave(
  plot = wrap_plots(
    JB22_mean_expr_gg_list, 
    nrow = 1
  ),
  filename = file.path(
    fig_path, 
    "mean_expr.svg"
  ),
  width = 10,
  height = 5,
  units = "cm"
)

# ggpllot: sample density
JB22_smpl_density_gg_list <- lapply(
  names(JB22_pb_norm_sce_list),
  function(tissue) {
    .plot_smpl_density(
      sce = JB22_pb_norm_sce_list[[tissue]],
      assay_name = "logcounts",
      n_samples = 50
    ) +
    ggtitle(tissue)
})

# Save the plot
ggsave(
  plot = wrap_plots(
    JB22_smpl_density_gg_list, 
    nrow = 1
  ),
  filename = file.path(
    fig_path, 
    "smpl_density.svg"
  ),
  width = 10,
  height = 5,
  units = "cm"
)

## Cluster ====================================================================

# Join list
JB22_pb_sce <- do.call(
  cbind, 
  JB22_pb_norm_sce_list
)

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

# Save the plot
ggsave(
  plot = .plot_PCA(
    sce = JB22_pb_sce
  ),
  filename = file.path(
    fig_path, 
    "PCA.svg"
  ),
  width = 10,
  height = 5,
  units = "cm"
)

## JB22 bulk RNA-seq ==========================================================

pipeline <- function(
  sce_obj,
  tissue,
  cell_type,
  formula,
  ref_level,
  shuffle,
  filter_fun,
  dge_fun,
  out_path
) {

  outfile <- paste0(out_path, ".rds")

  if (file.exists(outfile)) {
    return(invisible(NULL))
  }

  cat("\n------------------------------------------------------------------\n")
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
    shuffle = shuffle,
    filter_fun = filter_fun,
    dge_fun = dge_fun
  )

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

    # Add meta
    res$tt_obj$tissue    <- tissue
    res$tt_obj$cell_type <- cell_type
    res$timestamp <- Sys.time()
    
    # Save
    saveRDS(res, outfile)
  }

  # Return
  cat("------------------------------------------------------------------\n")
  return(NULL)
}

dge_methods <- list(
  edgeR_QLF = .dge_fit_edgeR_QLF,
  limma_voom = .dge_fit_limma_voom,
  DESeq2_Wald = .dge_fit_DESeq2_Wald
)

dge_formulas <- list(
  genotype = ~ genotype,
  genotype_sample = ~ genotype + sample
)

dge_datasets <- do.call(
  rbind,
  lapply(names(JB22_pb_sce_list), function(tissue) {
    sce <- JB22_pb_sce_list[[tissue]]
    cts <- unique(as.data.frame(colData(sce))$cell_type)
    data.frame(tissue = tissue, cell_type = cts)
  })
)

dge_tasks <- merge(
  dge_datasets,
  expand.grid(
    method  = names(dge_methods),
    formula = names(dge_formulas),
    stringsAsFactors = FALSE
  ),
  by = NULL
)
# dge_tasks <- dge_tasks[7,]

# Set the number of used cores
n_cores <- as.integer(Sys.getenv("NSLOTS", unset = 1))
n_workers <- min(n_cores, nrow(dge_tasks))
register(MulticoreParam(n_workers))

m <- bplapply(seq_len(nrow(dge_tasks)), function(i) {

  tissue <- dge_tasks$tissue[i]
  ct     <- dge_tasks$cell_type[i]
  method <- dge_tasks$method[i]
  form   <- dge_tasks$formula[i]

  # Outfile
  out_path <- file.path(
    res_path, "dge",
    paste0(method, "_", form, "_", tissue, "_", ct)
  )

  pipeline(
    sce_obj = JB22_pb_sce_list[[tissue]],
    tissue = tissue,
    cell_type = ct,
    formula = dge_formulas[[form]],
    ref_level = "NTC",
    shuffle = FALSE,
    filter_fun = .edgeR_filterByExpr,
    dge_fun = dge_methods[[method]],
    out_path = out_path
  )
})

# limma_res <- readRDS("results/JB22_pseudobulk_orig/dge/limma_voom_genotype_Liver_Bcells.rds")
# as.data.table(limma_res$tt_obj)

