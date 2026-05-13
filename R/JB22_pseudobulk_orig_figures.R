## max. code width ============================================================


## Experimental description ===================================================

## Title: Arc Institute in vivo pertubation screens
## Author: DKatzlberger
## Date: 21.01.2026

## Libraries =================================================================

# Scource R scripts
source("R/utils.R")

# Activate renv + load libs.
.script_startup(
  c(
    "parallel",
    "data.table", 
    "ggplot2",
    "ComplexHeatmap",
    "circlize",
    "RColorBrewer",
    "wesanderson"
  )
)

## Palettes ===================================================================

.library_pal <- c(
  Br    = "#1b9e77",
  Pool1 = "#d95f02",
  Pool2 = "#7570b3"
)

.cluster_pal <- c(
  "1"  = "#1f78b4",
  "2"  = "#33a02c",
  "3"  = "#e31a1c",
  "4"  = "#ff7f00",
  "5"  = "#6a3d9a",
  "6"  = "#b15928",
  "7"  = "#a6cee3",
  "8"  = "#b2df8a",
  "9"  = "#fb9a99",
  "10" = "#fdbf6f",
  "11" = "#cab2d6",
  "12" = "#ffff99"
)

## Directories ================================================================

# Result directories
dge_dir <- file.path("results", "JB22_pseudobulk_orig", "dge")
lib_dir <- file.path("results", "JB22_gRNA_library", "tab")

# Figures
fig_dir <- file.path("results", "JB22_pseudobulk_orig", "fig")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

## Results ====================================================================

# gRNA library
gRNA_lib <- fread(file.path(lib_dir, "JB22_gRNA_library.csv"))

# DGE res 
dge_files <- list.files(dge_dir)
# DESeq2_files <- dge_files[grep("^DESeq2_Wald", dge_files)]
# edgeR_files <- dge_files[grep("^edgeR_QLF", dge_files)]
limma_files <- dge_files[grep("^limma_voom", dge_files)]

# Load each method
load_tt_obj <- function(
  file_vec, 
  cores = parallel::detectCores() - 1
) {
  t0 <- Sys.time()

  res <- mclapply(file_vec, function(f) {
    tt_obj <- readRDS(file.path(dge_dir, f))$tt_obj

    # Convert to data.table (by reference if already DT)
    setDT(tt_obj)

    # Required stats columns
    req_stats <- c("feature", "logFC", "pval", "padj")
    miss_stats <- setdiff(req_stats, names(tt_obj))
    if (length(miss_stats)) {
      stop(sprintf(
        "File %s is missing required columns: %s",
        f, paste(miss_stats, collapse = ", ")
      ))
    }

    return(tt_obj)
  }, mc.cores = cores)

  # Size + timing
  res_size_mb <- as.numeric(object.size(res)) / (1024^2)
  t1 <- Sys.time()

  cat(sprintf(
    paste0(
      "Files loaded:  %d dge res\n",
      "- load time:   %.2f sec\n",
      "- object size: %.2f MB\n"
    ),
    length(res),
    as.numeric(difftime(t1, t0, units = "secs")),
    res_size_mb
  ))

  return(res)
}

# List of (tt_obj)
# DESeq2_list <- load_tt_obj(DESeq2_files)
# edgeR_list  <- load_tt_obj(edgeR_files)
limma_list  <- load_tt_obj(limma_files)

# Combine -> data.table
cols <- c("method", "formula", "tissue", "cell_type", "coef", "feature", "logFC", "pval", "padj")
# DESeq2_dt <- rbindlist(DESeq2_list)[, ..cols][grepl("^genotype", coef)][, coef := sub("^genotype_(.*?)_vs_.*$", "\\1", coef)]
# edgeR_dt  <- rbindlist(edgeR_list)[, ..cols][grepl("^genotype", coef)][, coef := sub("^genotype(.*?)", "\\1", coef)]
limma_dt  <- rbindlist(limma_list)[, ..cols][grepl("^genotype", coef)][, coef := sub("^genotype(.*?)", "\\1", coef)]

# DGE_res <- rbindlist(list(DESeq2_dt, edgeR_dt, limma_dt), use.names = TRUE, fill = TRUE)

## Heatmap ====================================================================

limma_logFC <- limma_dt[formula == "~genotype"]

cor_list <- limma_logFC[, {
  
  wide <- dcast(.SD, feature ~ coef, value.var = "logFC")
  mat <- cor(
    wide[, -1],  
    use = "pairwise.complete.obs",
    method = "pearson"
  )
  
  .(cor_matrix = list(mat))
  
}, by = .(tissue, cell_type)]

# Correlation heatmap
mat <- cor_list$cor_matrix[[1]]
diag(mat) <- NA

# Top column annonation
top_col_anno <- unique(
  gRNA_lib[
    .(target = colnames(mat)),
    on = "target",
    .(
      target,
      library,
      cluster = factor(cluster)
    )
  ],
  by = "target"
)

top_col_anno <- .cmplx_anno(
  Library = top_col_anno$library,
  Cluster = top_col_anno$cluster,
  annotation_name_side = "left",
  which = "column",
  col = list(
    Library = .library_pal,
    Cluster = .cluster_pal
  )
)

# Heatmap
ht <- .cmplx_ht(
  mat,
  name = "Correlation logFC\n(Pearson)",
  top_annotation = top_col_anno,
  show_row_names = FALSE,
  show_column_names = TRUE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "Target",
  column_title_side = "bottom",
  na_col = "grey",
  cell_size = 1.55
)

# Save
.cmplx_save(
  plot = ht,
  filename = file.path(
    fig_dir,
    "effect_size_correlation.svg"
  ),
  height = 0.155 * nrow(mat) + 3,
  width = 0.155 * nrow(mat) + 3,
  units = "cm",
)

for (i in seq_len(nrow(cor_list))) {
  
  mat <- cor_list$cor_matrix[[i]]
  diag(mat) <- NA
  
  # Name for file
  nm <- paste(cor_list$tissue[i], cor_list$cell_type[i], sep = "_")



  mf_levels <- c(
    "Transcription", "Chromatin", "Immune", "Signaling",
    "Protein Binding", "Enzymatic", "RNA Binding", "Other", "None"
  )
  col_anno_df[, MF := factor(MF, levels = mf_levels)]
  ord <- order(col_anno_df$MF)

  mat <- mat[, ord]
  col_anno_df <- col_anno_df[ord]

  # Colors
  col_anno_cols <- c(
    "Transcription" = "#1b9e77",
    "Chromatin" = "#d95f02",
    "Immune" = "#7570b3",
    "Signaling" = "#e7298a",
    "Protein Binding" = "#66a61e",
    "Enzymatic" = "#e6ab02",
    "RNA Binding" = "#a6761d",
    "Other" = "grey80"
  )

  # Top annonations
  top_ha <- .cmplx_anno(
    `Molecular fun.` = col_anno_df$MF,              
    col = list(`Molecular fun.` = col_anno_cols),
    which = "column",
    annotation_name_side = "left",   
    base_size = base_size
  )
    
  # color function (fixed scale for comparability)
  col_fun <- colorRamp2(
      c(-1, 0, 1),
      c(wes_palette("Zissou1")[1], "white", wes_palette("BottleRocket2")[1])
    )
  
  ht <- .cmplx_ht(
    mat,
    name = "Pearson (logFC)",
    col = col_fun,
    top_annotation = top_ha,
    show_row_names = FALSE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = nm,
    row_title = "Genotype",
    column_title_side = "top",
    row_title_side = "right",
    na_col = "grey",
    cell_size = 0.2,
    base_size = base_size
  )

  # Save
  for (ext in formats) {
    .cmplx_save(
      plot = ht,
      filename = file.path(
        fig_dir,
        paste0("cor_logFC_", nm, ".", ext)
      ),
      height = 0.02 * nrow(mat) + 2,
      width = 0.02 * nrow(mat) + 2,
      units = "cm",
    )
  }
}


### Nr. DEGs ------------------------------------------------------------------

# Nr. significant genes
nr_degs <- DGE_res[
  , .(
    nr_degs   = sum(padj < alpha, na.rm = TRUE),
    nr_total = .N
  ),
  by = .(method, formula, tissue, cell_type, coef)
]

# Heatmap
nr_degs[, row_id := paste(method, formula, tissue, cell_type, sep = " | ")]

# Matrix
mat_dt <- dcast(
  nr_degs,
  row_id ~ coef,
  value.var = "nr_degs",
  fill = NA
)

mat_raw <- as.matrix(mat_dt[, -1])
rownames(mat_raw) <- mat_dt$row_id 

# Log transform
mat_log <- log10(mat_raw + 1)

row_totals <- rowSums(mat_raw, na.rm = TRUE)

# Column Annotations
mf_anno <- emf_res[, .(SYMBOL, MF = top_level)]
col_anno_df <- mf_anno[match(colnames(mat_raw), SYMBOL)]
col_anno_df[is.na(MF), MF := "Other"]

mf_levels <- c(
  "Transcription", "Chromatin", "Immune", "Signaling",
  "Protein Binding", "Enzymatic", "RNA Binding", "Other", "None"
)
col_anno_df[, MF := factor(MF, levels = mf_levels)]
ord <- order(col_anno_df$MF)

mat_raw <- mat_raw[, ord]
mat_log <- mat_log[, ord]
col_anno_df <- col_anno_df[ord]


# Colors
col_anno_cols <- c(
  "Transcription" = "#1b9e77",
  "Chromatin" = "#d95f02",
  "Immune" = "#7570b3",
  "Signaling" = "#e7298a",
  "Protein Binding" = "#66a61e",
  "Enzymatic" = "#e6ab02",
  "RNA Binding" = "#a6761d",
  "Other" = "grey80"
)

# Top annonations
top_ha <- .cmplx_anno(
  `Molecular fun.` = col_anno_df$MF,              
  col = list(`Molecular fun.` = col_anno_cols),
  which = "column",
  annotation_name_side = "left",   
  base_size = base_size
)

# Row Annotations
row_anno_split <- tstrsplit(rownames(mat_raw), " \\| ")
row_anno_df <- data.table(
  Method = trimws(row_anno_split[[1]]),
  Formula = trimws(row_anno_split[[2]]),
  Tissue = trimws(row_anno_split[[3]]),
  `Cell type` = trimws(row_anno_split[[4]])
)

ord <- order(
  row_anno_df$`Cell type`,
  row_anno_df$Tissue,
  row_anno_df$Formula,
  row_anno_df$Method
)

mat_raw <- mat_raw[ord, ]
mat_log <- mat_log[ord, ]
row_anno_df <- row_anno_df[ord]
row_totals  <- row_totals[ord]

# Colors
row_anno_cols <- list(
  Method = structure(
    .cmplx_colors(length(unique(anno_df$Method)), function() {
      RColorBrewer::brewer.pal(8, "Set1")
    }),
    names = unique(anno_df$Method)
  ),
  
  Formula = structure(
    .cmplx_colors(length(unique(anno_df$Formula)), function() {
      RColorBrewer::brewer.pal(8, "Dark2")
    }),
    names = unique(anno_df$Formula)
  ),
  
  Tissue = structure(
    .cmplx_colors(length(unique(anno_df$Tissue)), function() {
      RColorBrewer::brewer.pal(8, "Set2")
    }),
    names = unique(anno_df$Tissue)
  ),
  
  `Cell type` = structure(
    .cmplx_colors(length(unique(anno_df$`Cell type`)), function() {
      RColorBrewer::brewer.pal(12, "Paired")
    }),
    names = unique(anno_df$`Cell type`)
  )
)

# Left annotatios
left_ha <- .cmplx_anno(
  df = row_anno_df,
  col = row_anno_cols,
  annotation_name_side = "bottom",
  which = "row",
  base_size = base_size
)

# Right annotation
right_ha <- .cmplx_anno(
  `Total DEGs` = row_anno_barplot(
    row_totals,
    border = FALSE,
    gp = gpar(col = NA, fill = "#4c72b0"),
    axis_param = list(
      gp = gpar(fontsize = base_size)
    )  
  ),
  annotation_name_side = "bottom",
  which = "row",
  base_size = base_size
)

ht <- .cmplx_ht(
  mat_log,
  name = "# DEGs\n(log10p)",
  col_fun <- colorRamp2(
    c(0, max(mat_log, na.rm = TRUE)),
    c("white", wes_palette("BottleRocket2")[1])
  ),
  row_split = anno_df$Method,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_title = "Genotype",
  row_title = NULL,
  column_title_side = "bottom",
  left_annotation = left_ha,
  right_annotation = right_ha,
  top_annotation = top_ha,
  na_col = "grey",
  cell_size = NULL,
  base_size = base_size
)

# Save
.cmplx_save(
  plot = ht,
  filename = file.path(
    fig_dir,
    "nr_degs.svg"
  ),
  height = 4.2,
  width = 15,
  units = "cm"
)


### DESeq2 (full) -------------------------------------------------------------
DESeq2_dt <- rbindlist(DESeq2_list, use.names = TRUE, fill = TRUE)
DESeq2_dt <- DESeq2_dt[grepl("^genotype_", coef)]
DESeq2_dt[
  , coef_clean := sub("^genotype_(.*?)_vs_.*$", "\\1", coef)
][
  , condition := paste(method, formula, tissue, cell_type, coef_clean, sep = ".")
][]

# Total number of conditions
n_cond <- uniqueN(DESeq2_dt$condition)
DESeq2_dt_complete <- DESeq2_dt[
  , .SD[uniqueN(condition) == n_cond], 
  by = feature
]

DESeq2_dcast <- dcast(
  DESeq2_dt_complete,
  feature ~ condition,
  value.var = "logFC"
)

DESeq2_mat <- as.matrix(DESeq2_dcast[, -1])
rownames(DESeq2_mat) <- DESeq2_dcast$feature

# Correlation across features
DESeq2_cor_mat <- cor(DESeq2_mat, method = "pearson")
diag(DESeq2_cor_mat) <- NA

# Annotation
anno_dt <- data.table(condition = colnames(DESeq2_mat))

anno_dt[, c("Method", "Formula", "Tissue", "Cell type") :=
  tstrsplit(condition, ".", fixed = TRUE, keep = 1:4)
]

anno_df <- as.data.frame(anno_dt[, .(Method, Formula, Tissue, `Cell type`)])
rownames(anno_df) <- anno_dt$condition


# Colors
anno_colors <- list(
  Method = structure(
    get_colors2(length(unique(anno_df$Method)), function() {
      RColorBrewer::brewer.pal(8, "Set1")
    }),
    names = unique(anno_df$Method)
  ),
  
  Formula = structure(
    get_colors2(length(unique(anno_df$Formula)), function() {
      RColorBrewer::brewer.pal(8, "Dark2")
    }),
    names = unique(anno_df$Formula)
  ),
  
  Tissue = structure(
    get_colors2(length(unique(anno_df$Tissue)), function() {
      RColorBrewer::brewer.pal(8, "Set2")
    }),
    names = unique(anno_df$Tissue)
  ),
  
  `Cell type` = structure(
    get_colors2(length(unique(anno_df$`Cell type`)), function() {
      RColorBrewer::brewer.pal(12, "Paired")
    }),
    names = unique(anno_df$`Cell type`)
  )
)

target_size_cm <- 8
nr <- nrow(DESeq2_cor_mat)
nc <- ncol(DESeq2_cor_mat)
cell_cm <- target_size_cm / max(nr, nc)

top_ha <- .cmplx_anno(
  df = anno_df,
  col = anno_colors,
  annotation_name_side = "left",
  which = "column"
)

ht <- .cmplx_ht(
  DESeq2_cor_mat,
  name = "Pearson",
  col = colorRamp2(
    c(-1, 0, 1),
    c(wes_palette("Zissou1")[1], "white", wes_palette("BottleRocket2")[1])
  ),
  top_annotation = top_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "gRNAs",
  column_title_side = "bottom",
  na_col = "grey",
  width  = unit(nc * cell_cm, "cm"),
  height = unit(nr * cell_cm, "cm")
)

.cmplx_save(
  plot = ht,
  filename = file.path(
    "results/JB22_pseudobulk_orig/fig", 
    "DESeq2_cor.svg"
  ),
  height = 11.5,
  width = 11.5,
  units = "cm"
)


### Limma ---------------------------------------------------------------------

limma_dt <- rbindlist(limma_list, use.names = TRUE, fill = TRUE)
limma_dt <- limma_dt[grepl("^genotype", coef)]

# Filter (condition)
limma_dt <- limma_dt[formula == "~genotype" & tissue == "Liver"]

limma_dt[
  , coef_clean := sub("^genotype(.*?)", "\\1", coef)
]
limma_dt[
  , condition := paste(method, formula, tissue, cell_type, coef_clean, sep = ".")
]

# Total number of conditions
n_cond <- uniqueN(limma_dt$condition)
limma_dt_complete <- limma_dt[
  , .SD[uniqueN(condition) == n_cond], 
  by = feature
]

limma_dcast <- dcast(
  limma_dt_complete,
  feature ~ condition,
  value.var = "logFC"
)

limma_mat <- as.matrix(limma_dcast[, -1])
rownames(limma_mat) <- limma_dcast$feature

# Correlation across features
limma_cor_mat <- cor(limma_mat, method = "pearson")
diag(limma_cor_mat) <- NA


# Annotation
anno_dt <- data.table(condition = colnames(limma_mat))

anno_dt[, c("Method", "Formula", "Tissue", "Cell type") :=
  tstrsplit(condition, ".", fixed = TRUE, keep = 1:4)
]

anno_df <- as.data.frame(anno_dt[, .(Method, Formula, Tissue, `Cell type`)])
rownames(anno_df) <- anno_dt$condition

# Colors
anno_colors <- list(
  Method = structure(
    .cmplx_colors(length(unique(anno_df$Method)), function() {
      RColorBrewer::brewer.pal(8, "Set1")
    }),
    names = unique(anno_df$Method)
  ),
  
  Formula = structure(
    .cmplx_colors(length(unique(anno_df$Formula)), function() {
      RColorBrewer::brewer.pal(8, "Dark2")
    }),
    names = unique(anno_df$Formula)
  ),
  
  Tissue = structure(
    .cmplx_colors(length(unique(anno_df$Tissue)), function() {
      RColorBrewer::brewer.pal(8, "Set2")
    }),
    names = unique(anno_df$Tissue)
  ),
  
  `Cell type` = structure(
    .cmplx_colors(length(unique(anno_df$`Cell type`)), function() {
      RColorBrewer::brewer.pal(12, "Paired")
    }),
    names = unique(anno_df$`Cell type`)
  )
)


target_size_cm <- 12
nr <- nrow(limma_cor_mat)
nc <- ncol(limma_cor_mat)
cell_cm <- target_size_cm / max(nr, nc)

top_ha <- .cmplx_anno(
  df = anno_df,
  col = anno_colors,
  annotation_name_side = "left",
  which = "column"
)

ht <- .cmplx_ht(
  limma_cor_mat,
  name = "Pearson",
  col = colorRamp2(
    c(-1, 0, 1),
    c(wes_palette("Zissou1")[1], "white", wes_palette("BottleRocket2")[1])
  ),
  top_annotation = top_ha,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_title = "gRNAs",
  column_title_side = "bottom",
  na_col = "grey",
  width  = unit(nc * cell_cm, "cm"),
  height = unit(nr * cell_cm, "cm")
)

.cmplx_save(
  plot = ht,
  filename = file.path(
    "results/JB22_pseudobulk_orig/fig", 
    "limma_cor.png"
  ),
  height = 20,
  width = 20,
  units = "cm"
)
