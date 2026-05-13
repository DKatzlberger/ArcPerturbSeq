## max. code width ============================================================

## script ---------------------------------------------------------------------

.script_startup <- function(
  packages = NULL,
  status = FALSE,
  detach = FALSE
){

  # Helper fun.
  format_path <- function(path) {
    home <- normalizePath(path.expand("~"), winslash = "/")
    path <- normalizePath(path, winslash = "/")
    sub(paste0("^", home), "~", path)
  }

  cat("\n")
  t0 <- Sys.time()

  # renv (only when needed)
  current_project <- tryCatch(renv::project(), error = function(e) NULL)

  if (is.null(current_project) || detach) {
    renv::load()
    current_project <- renv::project()
  } else {

    renv_version <- as.character(utils::packageVersion("renv"))
    cat(sprintf(
      "- Project '%s' active. [renv %s]\n",
      format_path(current_project),
      renv_version
    ))
  }

  if (status){
    renv::status()
  }

  t1 <- Sys.time()

  # Libraries (only if provided)
  if (!is.null(packages)) {

    if (!is.character(packages)) {
      stop("`packages` must be a character vector.")
    }

    suppressPackageStartupMessages({
      for (pkg in unique(packages)) {
        library(pkg, character.only = TRUE, quietly = TRUE)
      }
    })
  }

  t2 <- Sys.time()
  
  # Report
  cat(sprintf(
    "Startup completed in %.2f seconds\n",
    as.numeric(difftime(t2, t0, units = "secs"))
  ))

  cat(sprintf(
    "- renv activation: %.2f seconds\n",
    as.numeric(difftime(t1, t0, units = "secs"))
  ))

  cat(sprintf(
    "- lib. attachment: %.2f seconds\n",
    as.numeric(difftime(t2, t1, units = "secs"))
  ))

  # Return
  return(invisible(TRUE))
}

## datasets -------------------------------------------------------------------

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

## ggplot2 theme  -------------------------------------------------------------

.theme_DK <- function(
  base_size = 5,
  lgd_key_prop = 1,
  show_axis = TRUE,
  show_grid = TRUE,
  show_facet_x = TRUE,
  show_facet_y = TRUE,
  show_ticks_x = TRUE,
  show_ticks_y = TRUE,
  show_borders = TRUE,
  zero_margins = TRUE,
  rotate = NULL,
  base_family = "Sans",
  ...
) {

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Base theme
  p <- ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    strip.background = ggplot2::element_rect(fill = "white", color = NA),

    # Facets
    strip.placement = "inside",
    strip.text.x = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    strip.text.y = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, t = 0, unit = "mm")),
    panel.spacing.x = grid::unit(0.5, "mm"),
    panel.spacing.y = grid::unit(0.5, "mm"),

    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(1, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey80", linewidth = 0.3) else ggplot2::element_blank(),

    # Text
    axis.text  = ggplot2::element_text(size = base_size, family = base_family),
    axis.title = ggplot2::element_text(size = base_size, family = base_family),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.caption = element_text(size = base_size, hjust = 0, family = base_family),
    ...
  )

  # Legend sizing
  lgd_text_size <- base_size * lgd_key_prop
  lgd_key_size  <- base_size * lgd_key_prop

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 1, unit = "mm")),
    legend.text = ggplot2::element_text(size = base_size, family = base_family),

    legend.key.height = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.width = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.spacing = ggplot2::unit(0.5, "mm"),
    legend.box.spacing = ggplot2::unit(1, "mm"),
    legend.background = ggplot2::element_rect(fill = NA, color = NA)
  )

  # Adaptive margins (zero margins)
  if (zero_margins) {
    p <- p + ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
      plot.margin = ggplot2::margin(0.25, 0.25, 0.25, 0.25, unit = "mm")
    )
  }

  # Adaptive facets (x and y facet)
  if (!show_facet_x) p <- p + ggplot2::theme(strip.text.x = ggplot2::element_blank())
  if (!show_facet_y) p <- p + ggplot2::theme(strip.text.y = ggplot2::element_blank())

  # Adaptive boarders
  if (show_borders) { p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line = ggplot2::element_blank()
    )
  }

  # Adaptive axis
  if (!show_axis) {
    p <- p + ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.5, b = 0, unit = "mm")), # pull x labels closer
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, l = 0, unit = "mm"))  # pull y labels closer
    )
  }

  if (!show_ticks_x) {
    p <- p + ggplot2::theme(
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank()
    )
  }
  if (!show_ticks_y) {
    p <- p + ggplot2::theme(
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank()
    )
  }

  # Adaptive rotation
  if (!is.null(rotate)) {
    angle <- rotate %% 360
    # Determine justification heuristically
    if (angle == 0) {
      hjust <- 0.5; vjust <- 0.5
    } else if (angle == 90) {
      hjust <- 1; vjust <- 0.5
    } else if (angle == 270) {
      hjust <- 0; vjust <- 0.5
    } else if (angle > 0 && angle < 90) {
      hjust <- 1; vjust <- 1
    } else if (angle > 90 && angle < 180) {
      hjust <- 1; vjust <- 0
    } else if (angle > 180 && angle < 270) {
      hjust <- 0; vjust <- 0
    } else { # e.g., 315°
      hjust <- 0; vjust <- 1
    }

    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = angle,
        hjust = hjust,
        vjust = vjust,
        family = base_family
      )
    )
  }

  # Return
  return(p)
}

## ComplexHeatmap -------------------------------------------------------------
.cmplx_gp <- function(
  base_size = 5, 
  base_family = "sans", 
  ...
) {
  grid::gpar(
    fontsize = base_size,
    fontfamily = base_family,
    fontface = "plain",
    ...
  )
}

.cmplx_colors <- function(n, palette_fun) {
  base_cols <- palette_fun()
  
  if (length(base_cols) >= n) {
    return(base_cols[seq_len(n)])
  } else {
    return(colorRampPalette(base_cols)(n))
  }
}

.cmplx_anno <- function(
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  which = c("column", "row")
) {

  which <- match.arg(which)

  dots <- list(...)
  name_side <- dots$annotation_name_side

  # Only validate if provided
  if (!is.null(name_side)) {
    if (which == "column" && 
        !name_side %in% c("left", "right")) {
      stop("For column annotations, annotation_name_side must be 'left' or 'right'")
    }

    if (which == "row" && 
        !name_side %in% c("top", "bottom")) {
      stop("For row annotations, annotation_name_side must be 'top' or 'bottom'")
    }
  }

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Annotation scaling
  gp_base <- .cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")

  # Legend scaling
  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  anno <- ComplexHeatmap::HeatmapAnnotation(
    ...,

    which = which,

    # Annotation scaling
    simple_anno_size = anno_height,

    # Annotation names
    annotation_name_gp = gp_base,

    # Legend styling
    annotation_legend_param = list(
      title_position = "topleft",
      title_gp = gp_base,
      labels_gp = gp_base,
      gap = grid::unit(0.5, "mm"),
      grid_width = lgd_unit,
      grid_height = lgd_unit,
      legend_height = lgd_unit
    )
  )

  # Return
  return(anno)
}

.cmplx_ht <- function(
  mat,
  ...,
  base_size = 5,
  lgd_key_prop = 1,
  base_family = "sans",
  cell_size = NULL
) {

  stopifnot(is.matrix(mat) || is.data.frame(mat))
  mat <- as.matrix(mat)

  # Dimensions
  n_row <- nrow(mat)
  n_col <- ncol(mat)

  # Unit conversion helpers 
  mm_to_pt <- function(mm) mm / 0.352778
  pt_to_mm <- function(pt) pt * 0.352778

  # Annotation scaling
  gp_base <- .cmplx_gp(base_size, base_family)
  anno_height <- grid::unit(pt_to_mm(base_size), "mm")

  # Legend scaling
  lgd_key_size <- base_size * lgd_key_prop
  lgd_unit <- grid::unit(pt_to_mm(lgd_key_size), "mm")

  # Conditionally enforce square tiles
  size_args <- list()
  if (!is.null(cell_size)) {
    if (is.numeric(cell_size)) {
      cell_size <- grid::unit(cell_size, "mm")
    }
    size_args$width  <- cell_size * n_col
    size_args$height <- cell_size * n_row
  }

  ht <- do.call(ComplexHeatmap::Heatmap, c(list(
    mat,
    ...,

    # Heatmap text
    row_names_gp = gp_base,
    column_names_gp = gp_base,
    column_title_gp = gp_base,
    row_title_gp = gp_base,

    # Heatmap legend
    heatmap_legend_param = list(
      title_gp = gp_base,
      labels_gp = gp_base,
      grid_width = lgd_unit,
      legend_height = lgd_unit
    ),

    # Dendrogram styling
    row_dend_gp = grid::gpar(lwd = 0.3),
    column_dend_gp = grid::gpar(lwd = 0.3)
  ), size_args))

  # Return
  return(ht)
}

.cmplx_save <- function(
  plot,
  filename,
  width,
  height,
  units = "in",
  device = NULL,
  res = 300,
  ...
) {

  to_inches <- function(x, units) {
    switch(
      units,
      "in" = x,
      "cm" = x / 2.54,
      "mm" = x / 25.4,
      "px" = x / res,
      stop("Unsupported units: ", units)
    )
  }

  width_in  <- to_inches(width, units)
  height_in <- to_inches(height, units)

  if (is.null(device)) {
    ext <- tools::file_ext(filename)
    device <- switch(
      tolower(ext),
      pdf = grDevices::pdf,
      svg = grDevices::svg,
      png = grDevices::png,
      jpeg = grDevices::jpeg,
      jpg = grDevices::jpeg,
      tiff = grDevices::tiff,
      stop("Unsupported file type: ", ext)
    )
  }

  # raster vs vector handling
  if (identical(device, grDevices::png) ||
      identical(device, grDevices::jpeg) ||
      identical(device, grDevices::tiff)) {

    device(
      filename,
      width = width_in,
      height = height_in,
      units = "in",
      res  = res,
      ...
    )

  } else {

    device(
      filename,
      width = width_in,
      height = height_in,
      ...
    )
  }

  ComplexHeatmap::ht_opt(
    legend_gap = grid::unit(1, "mm")
  )

  grid::grid.newpage()
  ComplexHeatmap::draw(
    plot,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    legend_grouping = "original",
    align_heatmap_legend = "heatmap_top",
    align_annotation_legend = "heatmap_top",
    merge_legend = TRUE
  )

  invisible(grDevices::dev.off())
}
