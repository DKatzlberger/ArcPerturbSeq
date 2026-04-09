## max. code width ============================================================
## ggplot2 themes  ------------------------------------------------------------

.theme_DK <- function(
  base_size = 6,
  lgd_text_prop = 1,
  lgd_key_prop = 1,
  show_axis = TRUE,
  show_facets = TRUE,
  show_grid = FALSE,
  show_borders = FALSE,
  small_lgd = FALSE,
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
    strip.placement  = "inside",

    axis.line  = ggplot2::element_line(color = "black", linewidth = 0.5),
    axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = ggplot2::unit(1, "mm"),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "grey80", linewidth = 0.3) else ggplot2::element_blank(),

    axis.text  = ggplot2::element_text(size = base_size, family = base_family),
    axis.title = ggplot2::element_text(size = base_size, family = base_family),
    plot.title = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0.5, family = base_family),
    strip.text = ggplot2::element_text(size = base_size, family = base_family),
    plot.caption = element_text(size = base_size, hjust = 0, family = base_family),
    ...
  )

  # Legend sizing
  lgd_text_size <- base_size * lgd_key_prop
  lgd_key_size  <- base_size * lgd_key_prop

  p <- p + ggplot2::theme(
    legend.title = ggplot2::element_text(size = base_size, family = base_family, margin = ggplot2::margin(b = 2, unit = "mm")),
    legend.text = ggplot2::element_text(size = base_size, family = base_family),

    legend.key.height = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.width = ggplot2::unit(pt_to_mm(lgd_key_size), "mm"),
    legend.key.spacing = ggplot2::unit(0.5, "mm"),
    legend.box.spacing = ggplot2::unit(1, "mm"),
    legend.background = ggplot2::element_rect(fill = NA, color = NA),
    legend.margin = if (small_lgd) ggplot2::margin(t = 0, b = 0, l = 0, r = 0, unit = "mm")
  )

  # Adaptive facets
  if (!show_facets) p <- p + ggplot2::theme(strip.text = ggplot2::element_blank())

  # Adaptive boarders
  if (show_borders) { p <- p + ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line    = ggplot2::element_blank()
    )
  }

  # Adaptive axis
  if (!show_axis) {
    p <- p + ggplot2::theme(
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(), 
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0.5, b = 0, unit = "mm")), # pull x labels closer
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0.5, l = 0, unit = "mm"))  # pull y labels closer
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

  return(p)
}

## SCE objects ================================================================

