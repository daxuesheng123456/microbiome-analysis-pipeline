relative_abundance <- function(ps) {
  phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
}

run_pcoa <- function(ps, dist = "bray") {
  suppressMessages(library(phyloseq))

  d <- phyloseq::distance(ps, method = dist, type = "samples")
  pcoa <- stats::cmdscale(d, k = 2, eig = TRUE)

  points <- as.data.frame(pcoa$points)
  colnames(points) <- c("x", "y")
  points$ID <- rownames(points)

  map <- plain_df(phyloseq::sample_data(ps))
  points$Group <- map[match(points$ID, rownames(map)), "Group"]

  list(points = points, eig = pcoa$eig, distance = d)
}

pcoa_variance <- function(eig, n_axes = 2) {
  pos <- pmax(eig, 0)
  total <- sum(pos)
  if (total <= 0) return(rep(NA_real_, n_axes))
  100 * pos[seq_len(n_axes)] / total
}

group_centroids <- function(points, group_column = "Group") {
  stats::aggregate(points[, c("x", "y")],
                   by = list(Group = points[[group_column]]),
                   FUN = mean)
}

plot_pcoa_spider <- function(points, var_pct = c(NA, NA), title = NULL,
                             pal = NULL, ellipse_level = 0.68,
                             show_labels = FALSE) {
  check_packages(c("ggplot2", "ggrepel"))
  suppressMessages({
    library(ggplot2)
  })

  points$Group <- as.factor(points$Group)
  if (is.null(pal)) {
    pal <- default_group_palette(levels(points$Group))
    names(pal) <- levels(points$Group)
  }

  cent <- group_centroids(points)
  segs <- merge(points, stats::setNames(cent, c("Group", "cx", "cy")),
                by = "Group", sort = FALSE)

  xlab <- if (is.na(var_pct[1])) "PCoA 1" else
    sprintf("PCoA 1 (%.2f%%)", var_pct[1])
  ylab <- if (is.na(var_pct[2])) "PCoA 2" else
    sprintf("PCoA 2 (%.2f%%)", var_pct[2])

  p <- ggplot(points, aes(x = x, y = y)) +

    geom_segment(data = segs,
                 aes(xend = cx, yend = cy, colour = Group),
                 linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +

    stat_ellipse(aes(group = Group, colour = Group),
                 type = "norm", level = ellipse_level,
                 linetype = 2, linewidth = 0.5, show.legend = FALSE) +

    geom_point(aes(fill = Group), shape = 21, size = 4, stroke = 0.6,
               colour = "grey20", alpha = 0.85) +

    geom_point(data = cent, aes(x = x, y = y),
               shape = 24, size = 5, colour = "black", fill = "yellow") +
    scale_fill_manual(values = pal, name = "Group") +
    scale_colour_manual(values = pal) +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = c(0.18, 0.18),
          legend.background = element_rect(fill = "white", colour = "grey70"),
          legend.key = element_blank(),
          plot.title = element_text(size = 10, hjust = 0))

  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(aes(label = ID), size = 3,
                                      max.overlaps = 20, seed = 1)
  }
  p
}

default_group_palette <- function(n_or_levels) {
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
           "#66A61E", "#E6AB02", "#A6761D", "#666666",
           "#1F78B4", "#B2DF8A")
  n <- if (is.numeric(n_or_levels)) n_or_levels else length(n_or_levels)
  if (n > length(pal)) {
    pal <- grDevices::colorRampPalette(pal)(n)
  }
  pal[seq_len(n)]
}
