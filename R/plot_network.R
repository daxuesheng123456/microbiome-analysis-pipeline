plot_one_network <- function(net, grp, cfg_plot = list()) {

  check_packages(c("ggplot2", "ggnewscale"))
  suppressMessages({
    library(ggplot2)
    library(ggnewscale)
  })

  width      <- cfg_plot$width       %||% 12
  height     <- cfg_plot$height      %||% 10
  pt_size    <- cfg_plot$point_size  %||% 3
  ln_width   <- cfg_plot$line_width  %||% 0.5
  pos_col    <- cfg_plot$pos_col     %||% "#D95F02"
  neg_col    <- cfg_plot$neg_col     %||% "#1B9E77"

  edge <- net$edges
  dat  <- net$layout

  mod_tbl <- unique(data.frame(mod = dat$model, col = dat$color,
                               stringsAsFactors = FALSE))
  col_map <- stats::setNames(mod_tbl$col, mod_tbl$mod)

  edge$cor <- factor(edge$cor, levels = c("+", "-"))

  p <- ggplot() +
    geom_segment(
      data = edge,
      aes(x = X1, y = Y1, xend = X2, yend = Y2, colour = cor),
      linewidth = ln_width,
      show.legend = TRUE
    ) +
    scale_colour_manual(
      name   = "Association",
      values = c("+" = pos_col, "-" = neg_col),
      drop   = FALSE
    ) +
    ggnewscale::new_scale_colour() +
    geom_point(
      data = dat,
      aes(x = X1, y = X2, colour = model),
      size = pt_size
    ) +
    scale_colour_manual(name = "Module", values = col_map) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_void() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", hjust = 0.5, size = 14)) +
    labs(title = sprintf("%s  (nodes: %d, edges: %d, +:%d / -:%d)",
                         grp, igraph::vcount(net$igraph), igraph::ecount(net$igraph),
                         net$n_pos, net$n_neg))

  p
}

save_network_plot <- function(p, path, cfg_plot = list()) {
  check_packages(c("ggplot2"))
  width  <- cfg_plot$width  %||% 12
  height <- cfg_plot$height %||% 10
  ggplot2::ggsave(path, p, width = width, height = height)
  log_msg("  figure saved: ", basename(path))
  invisible(path)
}
