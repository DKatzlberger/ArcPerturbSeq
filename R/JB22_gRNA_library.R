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
    "SingleCellExperiment",
    "data.table", 
    "clusterProfiler",
    "enrichplot",
    "igraph",
    "org.Mm.eg.db",
    "ggplot2",
    "ggraph",
    "patchwork"
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

# Shared directory (location: vscratch)
JB22_pb_dir <- "/vscratch/wes/arc/share"

# Define paths
res_dir <- file.path("results", "JB22_gRNA_library")
fig_dir <- file.path(res_dir, "fig")
tab_dir <- file.path(res_dir, "tab")

# Create output dirs
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

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

## gRNA library ===============================================================

# Pool of target genes
gRNA_lib <- rbindlist(lapply(JB22_pb_sce_list, function(x) {
  as.data.table(colData(x))
}), use.names = TRUE, fill = TRUE)

# Extract lib. from sample name
gRNA_lib[
  , library := fifelse(
    grepl("-Br-", sample), "Br",
    fifelse(
      grepl("-Pool1-", sample), "Pool1",
      fifelse(
        grepl("-Pool2-", sample), "Pool2", 
        NA_character_
      )
    )
  )
]

# Exlude NTC guides
gRNA_lib <- gRNA_lib[
  genotype != "NTC"
]

gRNA_lib <- gRNA_lib[, .(library, genotype, SYMBOL = genotype)]
gRNA_lib <- unique(gRNA_lib, by = "genotype")

## GO enrichment ==============================================================

gRNA_lib_list <- split(gRNA_lib$SYMBOL, gRNA_lib$library)
gRNA_lib_list <- lapply(gRNA_lib_list, unique)

# Gene ontology (Biological processes)
ego_list <- lapply(gRNA_lib_list, function(x) {
  clusterProfiler::enrichGO(
    gene = x,
    OrgDb = org.Mm.eg.db,
    keyType  = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    readable = TRUE
  )
})

# Simplify terms
ego_list <- lapply(ego_list, function(x) {
  
  # Keep sig. terms (< 0.01)
  df <- x@result
  df <- df[df$p.adjust < 0.01, ]
  
  if (nrow(df) == 0) {
    x@result <- df
    return(x)
  }
  
  # Keep top 300 terms
  df <- df[order(df$p.adjust), ][1:min(300, nrow(df)), ]
  x@result <- df
  
  clusterProfiler::simplify(
    x,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )
})

## Clusters ===================================================================

ego_cluster_list <- Map(function(x, nm) {

  # Generate graph
  ego <- pairwise_termsim(x)
  sim <- ego@termsim

  sim[is.na(sim)] <- 0
  rownames(sim) <- colnames(sim) <- ego@result$ID

  # Cluster GO graph
  set.seed(42)
  g <- graph_from_adjacency_matrix(
    sim,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # Louvain algo.
  cl <- cluster_louvain(g, weights = E(g)$weight)
  mm <- membership(cl)

  V(g)$cluster <- mm[V(g)$name]
  V(g)$label   <- ego@result$Description[match(V(g)$name, ego@result$ID)]
  V(g)$library <- nm

  # GO res. clusters
  dt <- as.data.table(ego@result)
  dt[, cluster := mm[match(ID, names(mm))]]
  dt[, library := nm]

  # Cluster stats.
  dt[, n_terms_cluster := .N, by = cluster]                                      # Nr. terms per cluster
  dt[, n_genes_cluster := uniqueN(unlist(strsplit(geneID, "/"))), by = cluster]  # Nr. target genes per cluster
  dt[, n_term_genes := lengths(strsplit(geneID, "/"))]                           # Nr. target genes per GO-term
  dt[, term_coverage := n_term_genes / n_genes_cluster]                          # Frac. target genes covered by GO-term per cluster

  # Cluster summaries
  dt_top <- dt[
    order(p.adjust),
    {
      top_dt        <- head(.SD, 5)
      top_genes     <- unique(unlist(strsplit(top_dt$geneID, "/")))
      cluster_genes <- unique(unlist(strsplit(geneID, "/")))

      list(
        n_top_genes  = length(top_genes),
        top_coverage = length(top_genes) /length(cluster_genes),
        top_terms    = paste(top_dt$Description, collapse = " | ")
      )
    },
    by = cluster
  ]

  dt <- merge(
    dt,
    dt_top,
    by = "cluster",
    all.x = TRUE
  )

  # Return
  list(ego_res = dt, graph = g)

}, ego_list, names(ego_list))

## Gene cluster assignment ====================================================
ego_score_list <- Map(function(x, nm) {

  # Gene frequency within cluster
  gene_scores <- x$ego_res[
    , .(
      gene = unlist(strsplit(geneID, "/")),
      p.adjust
    ), 
    by = .(cluster, ID)
  ][
    , .(
      gene_count  = .N,
      gene_weight = sum(-log10(p.adjust))
    ),
    by = .(gene, cluster)
  ]

  # Cluster metadata
  cluster_meta <- unique(
    x$ego_res[, .(cluster, n_terms_cluster, top_terms)]
  )

  gene_scores <- merge(
    gene_scores,
    cluster_meta,
    by = "cluster",
    all.x = TRUE
  )

  # TF score (specificity for cluster)  
  gene_scores[, library := nm]    
  gene_scores[, tf := gene_weight / n_terms_cluster]                  

  # Priority:
  # 1. Highest score
  # 2. Highest cumulative enrichment
  # 3. Highest GO-term frequency
  # 4. Lowest cluster ID
  setorder(
      gene_scores,
      gene,
      -tf,
      -gene_weight,
      -gene_count,
      cluster
    )

  gene_scores[, is_best := seq_len(.N) == 1, by = gene]
  setorder(gene_scores, gene, -tf)
  gene_scores[]

}, ego_cluster_list, names(ego_cluster_list))

## Save gRNA library ==========================================================

# Enriched targets (contain all genes in GO-terms)
gRNA_lib_enriched <- rbindlist(ego_score_list)[
  is_best == TRUE,
  .(
    library,
    target = gene,
    cluster,
    terms = top_terms
  )
]

# Reduce to original library
gRNA_lib_final <- merge(

  # Original library
  unique(
    gRNA_lib[
      ,
      .(
        library,
        target = as.character(genotype)
      )
    ]
  ),

  # Enrichment annotations
  gRNA_lib_enriched,

  by = c("library", "target"),
  all.x = TRUE,
  sort = FALSE
)

# Save
fwrite(
  gRNA_lib_final, 
  file = file.path(
    tab_dir, 
    "JB22_gRNA_library.csv"
  )
)

## Visualize ==================================================================
d_graph_layout <- Filter(Negate(is.null), lapply(ego_cluster_list, `[[`, "graph"))
d_graph_layout <- do.call(disjoint_union, d_graph_layout)

# Plot
p_graph_layout <- ggraph(
  graph = d_graph_layout, 
  layout = "fr"
) +
geom_node_point(
  mapping = aes(color = factor(cluster)), 
  size = 0.5
) +
scale_color_manual(
  values = .cluster_pal
) +
facet_nodes(
  ~library, 
  scales = "free"
) +
labs(
  x = NULL,
  y = NULL,
  color = "Cluster",
  caption = "Each node = GO-term"
) +
.theme_DK(
  show_grid = FALSE,
  show_ticks_x = FALSE,
  show_ticks_y = FALSE
)

# Grap stats
p_graph_stats <- wrap_plots(
  # Nr. terms per cluster
  ggplot(
    data = unique(rbindlist(lapply(ego_cluster_list, `[[`, "ego_res"))[
      , .(library, cluster, n_terms_cluster)
    ]), 
    mapping = aes(
      x = factor(cluster), 
      y = n_terms_cluster
    )
  ) +
  geom_col(
    mapping = aes(fill = library),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = .library_pal
  ) +
  facet_grid(
    cols = vars(library),
    space = "free",
    scales = "free"
  ) +
  labs(
    x = "Cluster",
    y = "Nr. GO-terms per cluster"
  ) +
  .theme_DK(
  ),
  # Nr. genes per cluster
  ggplot(
    data = unique(rbindlist(lapply(ego_cluster_list, `[[`, "ego_res"))[
      , .(library, cluster, n_genes_cluster)
    ]), 
    mapping = aes(
      x = factor(cluster), 
      y = n_genes_cluster
    )
  ) +
  geom_col(
    mapping = aes(fill = library),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = .library_pal
  ) +
  facet_grid(
    cols = vars(library),
    space = "free",
    scales = "free"
  ) +
  labs(
    x = "Cluster",
    y = "Nr. target genes per cluster"
  ) +
  .theme_DK(),
  ncol = 2,
  nrow = 1
)

# Cluster terms
d_graph_terms <- rbindlist(lapply(ego_cluster_list, `[[`, "ego_res"))[
  order(library, cluster, p.adjust),
  head(.SD, 5),
  by = .(library, cluster)
]

# Background
bg_graph_terms <- unique(d_graph_terms[, .(library, cluster)])
bg_graph_terms[, xmin := seq_len(.N) - 0.5, by = library]
bg_graph_terms[, xmax := seq_len(.N) + 0.5, by = library]

p_graph_terms <- ggplot(
  data = d_graph_terms, 
  mapping = aes(
    x = factor(cluster),
    y = reorder(Description, -p.adjust),
    size = -log10(p.adjust),
    color = term_coverage
  )
) +
geom_rect(
  data = bg_graph_terms,
  mapping = aes(
    xmin = xmin,
    xmax = xmax,
    ymin = -Inf,
    ymax = Inf,
    fill = factor(cluster)
  ),
  inherit.aes = FALSE,
  show.legend = FALSE,
  alpha = 0.15
) +
scale_fill_manual(
  values = .cluster_pal
) +
geom_point() +
scale_color_gradient(
  limits = c(0, 1),
  low ="gold",
  high = "brown"
) +
scale_size_continuous(
  range = c(0.1, 1.5),
) +
facet_grid(
  cols = vars(library),
  space = "free",
  scales = "free"
) +
labs(
  x = "Cluster",
  y = "Term",
  color = "Target genes covered\nby term within cluster"
) +
.theme_DK(
  show_grid = FALSE
)

# Cluster stats top terms
p_cluster_stats <- ggplot(
  data = unique(rbindlist(lapply(ego_cluster_list, `[[`, "ego_res"))[
    , .(library, cluster, top_coverage)
  ]), 
  mapping = aes(
    x = factor(cluster), 
    y = top_coverage
  )
) +
geom_col(
  mapping = aes(fill = library),
  show.legend = FALSE
) +
scale_fill_manual(
  values = .library_pal
) +
facet_grid(
  cols = vars(library),
  space = "free",
  scales = "free"
) +
labs(
  x = "Cluster",
  y = "Target genes covered\nby top terms"
) +
.theme_DK()

p_assignments <- wrap_plots(
  ggplot(
    data = gRNA_lib_final[
      ,
      .(n_targets = uniqueN(target)),
      by = .(library, cluster)
    ],
    mapping = aes(
      x = cluster,
      y = n_targets
    ) 
  ) +
  geom_col(
    mapping = aes(fill = library),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = .library_pal
  ) +
  facet_grid(
    cols = vars(library), 
    scales = "free",
    space = "free"
  ) +
  labs(
    x = "Cluster",
    y = "Nr. target genes\nassigned to cluster"
  ) +
  .theme_DK(),
  ggplot(
    data = gRNA_lib_final[
      ,
      .(n_targets = uniqueN(target)),
      by = .(library, cluster)
    ][
      ,
      pct_targets := 100 * n_targets / sum(n_targets),
      by = library
    ],
    mapping = aes(
      x = library,
      y = pct_targets,
      fill = factor(cluster)
    )
  ) +
  geom_col() +
  scale_fill_manual(
    values = .cluster_pal
  ) +
  labs(
    x = "Library",
    y = "% target genes",
    fill = "Cluster",
  ) +
  .theme_DK(),
  ncol = 2,
  nrow = 1
)

# Patchwork
p <- (
  (
  free(p_graph_layout) +
  p_graph_stats +
  plot_layout(widths = c(0.1, 1))
  ) /
  p_graph_terms /
  (
  free(p_cluster_stats) +
  p_assignments +
  plot_layout(widths = c(0.1, 1))
  )
) +
plot_layout(
  heights = c(0.2, 1, 0.2)
)

# Save
ggsave(
  plot = p,
  file = file.path(
    fig_dir,
    "JB22_gRNA_lib_clusters.svg"
  ),
  width = 15,
  height = 20,
  unit = "cm"
)
