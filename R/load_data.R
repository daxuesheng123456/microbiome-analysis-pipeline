TAX_RANKS <- c("Domain", "Kingdom", "Phylum", "Class",
               "Order", "Family", "Genus", "Species")

read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) die("Input file not found: ", path)
  utils::read.csv(path, fileEncoding = "UTF-8-BOM",
                  check.names = FALSE, stringsAsFactors = FALSE, ...)
}

split_combined <- function(df) {
  present <- intersect(TAX_RANKS, colnames(df))
  if (length(present) < 4) {
    die("Combined table must contain taxonomy columns (",
        paste(TAX_RANKS, collapse = ", "), "). Found: ",
        paste(present, collapse = ", "))
  }
  tax_cols <- TAX_RANKS[TAX_RANKS %in% colnames(df)]
  tax_df   <- df[, tax_cols, drop = FALSE]

  rest     <- setdiff(colnames(df), tax_cols)

  is_num <- vapply(df[, rest, drop = FALSE],
                   function(x) is.numeric(x) ||
                     all(is.na(suppressWarnings(as.numeric(x))) == FALSE),
                   logical(1))
  id_candidates <- rest[!is_num]

  if (length(id_candidates) == 1) {
    ids <- as.character(df[[id_candidates]])
    rest <- setdiff(rest, id_candidates)
    log_msg("  feature ID column detected: ", id_candidates)
  } else if (length(id_candidates) > 1) {
    die("More than one non-numeric column besides taxonomy: ",
        paste(id_candidates, collapse = ", "))
  } else {
    ids <- sprintf("ASV%05d", seq_len(nrow(df)))
    log_msg("  no feature ID column; generated sequential IDs (ASV00001 ...)")
  }

  if (anyDuplicated(ids)) {
    die("Duplicated feature IDs detected (", sum(duplicated(ids)),
        " duplicates). Ensure the ID column is unique.")
  }

  abund <- df[, rest, drop = FALSE]
  rownames(abund) <- ids
  rownames(tax_df) <- ids
  list(abund = abund, tax = tax_df, ids = ids)
}

load_phyloseq <- function(otu_path,
                          tax_path   = NULL,
                          meta_path  = NULL,
                          id_col     = NULL,
                          sample_col = NULL) {

  check_packages(c("phyloseq"))
  suppressMessages(library(phyloseq))

  log_msg("Reading abundance table: ", basename(otu_path))
  otu_raw <- read_csv_safe(otu_path, row.names = NULL)

  if (is.null(tax_path)) {

    parts  <- split_combined(otu_raw)
    abund  <- parts$abund
    tax_df <- parts$tax
  } else {
    log_msg("Reading taxonomy table: ", basename(tax_path))
    tax_raw <- read_csv_safe(tax_path, row.names = NULL)
    id_o <- id_col %||% colnames(otu_raw)[1]
    id_t <- colnames(tax_raw)[1]
    abund  <- otu_raw[, setdiff(colnames(otu_raw), id_o), drop = FALSE]
    rownames(abund) <- otu_raw[[id_o]]
    tax_df <- tax_raw[match(rownames(abund), tax_raw[[id_t]]),
                      TAX_RANKS[TAX_RANKS %in% colnames(tax_raw)], drop = FALSE]
    rownames(tax_df) <- rownames(abund)
    tax_df[is.na(tax_df)] <- ""
  }

  if (is.null(rownames(abund))) {
    rownames(abund) <- sprintf("ASV%05d", seq_len(nrow(abund)))
    rownames(tax_df) <- rownames(abund)
  }

  abund[] <- lapply(abund, function(x) suppressWarnings(as.numeric(x)))
  bad <- vapply(abund, function(x) all(is.na(x)), logical(1))
  if (any(bad)) {
    die("Non-numeric abundance columns: ", paste(colnames(abund)[bad], collapse = ", "))
  }
  abund[is.na(abund)] <- 0
  abund <- as.matrix(abund)

  log_msg(sprintf("  %d features x %d samples", nrow(abund), ncol(abund)))

  tax_df <- as.matrix(tax_df)
  rownames(tax_df) <- rownames(abund)
  colnames(tax_df) <- TAX_RANKS[seq_len(ncol(tax_df))]

  otu <- phyloseq::otu_table(abund, taxa_are_rows = TRUE)
  tax <- phyloseq::tax_table(tax_df)

  if (is.null(meta_path)) {
    log_msg("No metadata supplied; using sample names as a single group", level = "WARN")
    sam <- phyloseq::sample_data(data.frame(
      Group = rep("all", ncol(abund)), row.names = colnames(abund)))
  } else {
    log_msg("Reading metadata: ", basename(meta_path))
    meta <- read_csv_safe(meta_path, row.names = NULL)
    sc   <- sample_col %||% colnames(meta)[1]
    rownames(meta) <- meta[[sc]]
    meta[[sc]] <- NULL
    missing_samples <- setdiff(colnames(abund), rownames(meta))
    if (length(missing_samples)) {
      die("Samples present in abundance table but missing from metadata: ",
          paste(utils::head(missing_samples, 10), collapse = ", "),
          if (length(missing_samples) > 10) " ..." else "")
    }
    sam <- phyloseq::sample_data(meta[colnames(abund), , drop = FALSE])
  }

  ps <- phyloseq::phyloseq(otu, tax, sam)
  ps
}

`%||%` <- function(a, b) if (is.null(a)) b else a
