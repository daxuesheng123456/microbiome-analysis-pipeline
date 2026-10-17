top_features <- function(otu, n = NULL) {
  m <- as.matrix(otu)
  storage.mode(m) <- "double"
  rel <- sweep(m, 2, colSums(m), "/")
  if (!is.null(n) && n < nrow(m)) {
    ord <- order(rowMeans(rel, na.rm = TRUE), decreasing = TRUE)
    m <- m[ord[seq_len(n)], , drop = FALSE]
  }
  m
}

matcor_plot <- function(env, otu, group = "bac", top_n = NULL,
                        method_cor = "spearman", cor_p = 0.05, p_thur = 0.05,
                        onlysig = TRUE,
                        edge_cols = c("blue", "red"),
                        fill_palette = "PRGn", ...) {

  check_packages(c("ggClusterNet", "ggplot2", "dplyr"))

  env <- as.data.frame(env)
  otu <- top_features(otu, top_n)

  if (is.null(rownames(env))) die("`env` must have sample names as rownames")
  miss <- setdiff(rownames(env), colnames(otu))
  if (length(miss)) {
    die("Samples in `env` missing from the OTU table: ",
        paste(utils::head(miss, 10), collapse = ", "))
  }
  otu <- otu[, rownames(env), drop = FALSE]

  suppressMessages({ library(dplyr); library(ggplot2) })

  tabOTU <- list(otu)
  names(tabOTU) <- group

  p <- ggClusterNet::MatCorPlot(
    env.dat    = env,
    tabOTU     = tabOTU,
    method.cor = method_cor,
    cor.p      = cor_p,
    p.thur     = p_thur,
    onlysig    = onlysig,
    ...
  )

  p <- p +
    ggplot2::scale_colour_manual(values = edge_cols) +
    ggplot2::scale_fill_distiller(palette = fill_palette)

  attr(p, "n_features") <- nrow(otu)
  attr(p, "n_env")      <- ncol(env)
  attr(p, "n_samples")  <- nrow(env)
  p
}

read_env_table <- function(path) {
  df <- read_csv_safe(path, row.names = NULL)
  ids <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  if (any(vapply(df, function(x) all(is.na(x)), logical(1)))) {
    die("Non-numeric columns in the environmental table")
  }
  as.data.frame(df)
}

read_otu_table <- function(path) {
  df <- read_csv_safe(path, row.names = NULL)
  ids <- as.character(df[[1]])
  df <- df[, -1, drop = FALSE]
  rownames(df) <- ids
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  df[is.na(df)] <- 0
  as.matrix(df)
}
