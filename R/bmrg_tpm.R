# ---------------------------------------------------------------------------
# BMRG (BacMet metal-resistance gene) TPM utilities
#
# Gene-level TPM values of the 42-sample metagenome run were summed over all
# genes annotated to each BacMet BMRG category (data/BMRG_gene2BMRG.csv), and
# the resulting BMRG-level TPM matrix for the 21 rhizosphere samples is stored
# in data/BMRG_TPM_rhizosphere.csv (22 genes x 21 samples).
# ---------------------------------------------------------------------------

BMRG_GROUPS <- c("CK", "RSB", "GRSB", "FRSB", "RHB", "GRHB", "FRHB")

read_bmrg_tpm <- function(path = "data/BMRG_TPM_rhizosphere.csv") {
  if (!file.exists(path)) die("Input file not found: ", path)
  utils::read.csv(path, fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  stringsAsFactors = FALSE)
}

bmrg_matrix <- function(df) {
  sample_cols <- setdiff(colnames(df), c("Gene", "BMRG"))
  m <- as.matrix(df[, sample_cols, drop = FALSE])
  rownames(m) <- df$Gene
  m
}

bmrg_group_means <- function(df, groups = BMRG_GROUPS) {
  m   <- bmrg_matrix(df)
  grp <- sub("_[0-9]+$", "", colnames(m))
  out <- vapply(groups, function(g) {
    cols <- which(grp == g)
    if (!length(cols)) die("No samples found for group: ", g)
    rowMeans(m[, cols, drop = FALSE])
  }, numeric(nrow(m)))
  rownames(out) <- rownames(m)
  out
}

bmrg_log2fc <- function(means, ref = "CK") {
  if (!ref %in% colnames(means)) die("Reference group not found: ", ref)
  log2(means / means[, ref])
}

plot_bmrg_abundance <- function(means) {
  long <- as.data.frame(as.table(means), stringsAsFactors = FALSE)
  colnames(long) <- c("Gene", "Treatment", "TPM")
  long$Gene      <- factor(long$Gene, levels = rev(rownames(means)))
  long$Treatment <- factor(long$Treatment, levels = colnames(means))
  long$Label     <- ifelse(long$TPM >= 10,
                           sprintf("%.0f", long$TPM),
                           sprintf("%.2f", long$TPM))
  ggplot2::ggplot(long, ggplot2::aes(.data$Treatment, .data$Gene,
                                     fill = log10(.data$TPM))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), size = 2.4) +
    ggplot2::scale_fill_gradientn(
      colours = c("#08306B", "#2171B5", "#6BAED6", "#C7E9B4", "#FFFFCC"),
      name    = expression(log[10] * "(TPM)")) +
    ggplot2::labs(x = NULL, y = NULL, title = "BMRG abundance (TPM)") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid   = ggplot2::element_blank())
}

plot_bmrg_log2fc <- function(lfc, threshold = 0.2) {
  long <- as.data.frame(as.table(lfc), stringsAsFactors = FALSE)
  colnames(long) <- c("Gene", "Treatment", "log2FC")
  long$Gene      <- factor(long$Gene, levels = rev(rownames(lfc)))
  long$Treatment <- factor(long$Treatment, levels = colnames(lfc))
  long$Label     <- ifelse(abs(long$log2FC) >= threshold,
                           sprintf("%.2f", long$log2FC), "")
  lim <- max(0.5, ceiling(max(abs(lfc)) * 10) / 10)
  ggplot2::ggplot(long, ggplot2::aes(.data$Treatment, .data$Gene,
                                     fill = .data$log2FC)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), size = 2.2) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-lim, lim),
                                  name = "log2FC") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "log2 fold change (treatment/CK)") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid   = ggplot2::element_blank())
}
