ps_to_microtable <- function(ps) {
  check_packages(c("microeco", "phyloseq"))
  suppressMessages(library(phyloseq))

  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  otu <- as.data.frame(otu, stringsAsFactors = FALSE)

  tax <- as.data.frame(as(phyloseq::tax_table(ps), "matrix"),
                       stringsAsFactors = FALSE)
  tax <- tax[rownames(otu), , drop = FALSE]

  samp <- plain_df(phyloseq::sample_data(ps))
  samp <- samp[colnames(otu), , drop = FALSE]

  if (!identical(rownames(tax), rownames(otu))) {
    die("tax_table rownames do not match otu_table rownames after alignment")
  }
  if (!identical(rownames(samp), colnames(otu))) {
    die("sample_table rownames do not match otu_table colnames after alignment")
  }

  microeco::microtable$new(otu_table = otu, tax_table = tax,
                           sample_table = samp)
}

taxa_composition <- function(mt, taxrank = "Species", ntaxa = 8,
                             groupmean = "Group", ...) {
  check_packages(c("microeco"))
  suppressMessages(
    microeco::trans_abund$new(dataset = mt, taxrank = taxrank,
                              ntaxa = ntaxa, groupmean = groupmean, ...)
  )
}

composition_palette <- function(n) {
  pal <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
           "#66A61E", "#E6AB02", "#A6761D", "#666666",
           "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
           "#A65628", "#F781BF")
  if (n > length(pal)) pal <- grDevices::colorRampPalette(pal)(n)
  pal[seq_len(n)]
}

composition_shown_taxa <- function(ta) {
  nm <- ta$data_taxanames
  if (is.null(nm)) nm <- unique(ta$data_abund$Taxonomy)
  nm
}

plot_composition_donut <- function(ta, label = FALSE, colors = NULL,
                                   facet_nrow = 1) {
  check_packages(c("microeco", "ggplot2"))
  if (is.null(colors)) {
    n <- length(composition_shown_taxa(ta)) + 1L
    colors <- composition_palette(n)
  }
  ta$plot_donut(color_values = colors, label = label, facet_nrow = facet_nrow)
}

plot_composition_bar <- function(ta, others_color = "grey90",
                                 coord_flip = FALSE, ...) {
  check_packages(c("microeco", "ggplot2"))
  ta$plot_bar(others_color = others_color, coord_flip = coord_flip, ...)
}

composition_table <- function(ta) {
  as.data.frame(ta$data_abund, stringsAsFactors = FALSE)
}

composition_plot_table <- function(ta) {
  dat   <- composition_table(ta)
  shown <- composition_shown_taxa(ta)

  dat$Taxonomy[!dat$Taxonomy %in% shown] <- "Others"
  out <- stats::aggregate(dat$Abundance,
                          by = list(Taxonomy = dat$Taxonomy, Sample = dat$Sample),
                          FUN = sum)
  colnames(out)[3] <- "Abundance"
  out$Taxonomy <- factor(out$Taxonomy, levels = c(shown, "Others"))
  out <- out[order(out$Sample, out$Taxonomy), ]
  out$Taxonomy <- as.character(out$Taxonomy)
  rownames(out) <- NULL
  out
}
