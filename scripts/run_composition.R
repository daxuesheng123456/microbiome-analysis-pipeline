script_dir <- local({
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else getwd()
})
repo_root <- if (basename(script_dir) == "scripts") dirname(script_dir) else script_dir
if (!dir.exists(file.path(repo_root, "R"))) {
  stop("Cannot locate the R/ directory. Run this script from inside the repository.")
}

for (f in list.files(file.path(repo_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- parse_args()

if (isTRUE(args$help) || is.null(args$config)) {
  cat("
Taxonomic composition module (microeco::trans_abund)

Options:
  --config <path>   path to the config file (required)
  --out <dir>       override the composition output directory
  --taxrank <rank>  Phylum | Class | Order | Family | Genus | Species (default Species)
  --ntaxa <n>       number of most abundant taxa shown separately (default 8)
  --per-sample      plot individual samples instead of group means
  --help            show this message
")
  quit(status = if (isTRUE(args$help)) 0 else 1)
}

setwd(repo_root)
log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              level, paste0(...)))
}

cfg <- load_config(args$config)

taxrank   <- args$taxrank %||% cfg$composition$taxrank %||% "Species"
ntaxa     <- as.integer(args$ntaxa %||% cfg$composition$ntaxa %||% 8)
out_dir   <- args$out %||% cfg$composition$output_dir %||%
             file.path(cfg$output_dir, "composition")
group_col <- cfg$group_column %||% "Group"
groupmean <- if (isTRUE(args[["per-sample"]])) NULL else group_col

ensure_dir(out_dir)
logfile <- file.path(out_dir, "composition.log")
if (file.exists(logfile)) invisible(file.remove(logfile))
log_msg <- log_to(logfile)

log_msg("=== Taxonomic composition ===", level = "STEP")
log_msg("taxrank  : ", taxrank)
log_msg("ntaxa    : ", ntaxa)
log_msg("groupmean: ", if (is.null(groupmean)) "(per sample)" else groupmean)
log_msg("output   : ", normalizePath(out_dir, mustWork = FALSE))

check_packages(c("phyloseq", "microeco", "ggplot2"))

inp <- cfg$input
ps <- tryCatch(
  load_phyloseq(otu_path   = inp$otu_table,
                tax_path   = inp$taxonomy,
                meta_path  = inp$metadata,
                id_col     = inp$id_col,
                sample_col = inp$sample_col),
  error = function(e) die("Failed to build phyloseq object: ", conditionMessage(e))
)

mt <- tryCatch(
  ps_to_microtable(ps),
  error = function(e) die("Could not build the microtable: ", conditionMessage(e))
)
log_msg(sprintf("microtable: %d samples x %d features",
                nrow(mt$sample_table), nrow(mt$otu_table)))

ta <- tryCatch(
  taxa_composition(mt, taxrank = taxrank, ntaxa = ntaxa, groupmean = groupmean),
  error = function(e) die("trans_abund failed for rank '", taxrank, "': ",
                          conditionMessage(e))
)
shown  <- composition_shown_taxa(ta)
log_msg("taxa drawn individually: ", length(shown),
        " (+ Others)")
log_msg("  ", paste(shown, collapse = ", "))

colors <- composition_palette(length(shown) + 1L)

figs <- list(
  donut_plain = list(fn = plot_composition_donut, label = FALSE),
  donut_label = list(fn = plot_composition_donut, label = TRUE)
)

for (nm in names(figs)) {
  spec <- figs[[nm]]
  p <- tryCatch(
    plot_composition_donut(ta, label = spec$label, colors = colors),
    error = function(e) {
      log_msg(nm, " failed: ", conditionMessage(e), level = "WARN")
      NULL
    }
  )
  if (!is.null(p)) {
    f <- file.path(out_dir, sprintf("composition_%s_%s.pdf", taxrank, nm))
    ggplot2::ggsave(f, p, width = 9, height = 7)
    ggplot2::ggsave(sub("\\.pdf$", ".png", f), p, width = 9, height = 7, dpi = 300)
    log_msg("  figure saved: ", basename(f))
  }
}

p_bar <- tryCatch(
  plot_composition_bar(ta, others_color = "grey90"),
  error = function(e) {
    log_msg("bar plot failed: ", conditionMessage(e), level = "WARN")
    NULL
  }
)
if (!is.null(p_bar)) {
  f <- file.path(out_dir, sprintf("composition_%s_bar.pdf", taxrank))
  ggplot2::ggsave(f, p_bar, width = 9, height = 6)
  ggplot2::ggsave(sub("\\.pdf$", ".png", f), p_bar, width = 9, height = 6, dpi = 300)
  log_msg("  figure saved: ", basename(f))
}

tab <- composition_table(ta)
utils::write.csv(tab, file.path(out_dir, sprintf("composition_%s_abundance_all.csv", taxrank)),
                 row.names = FALSE)
plot_tab <- composition_plot_table(ta)
utils::write.csv(plot_tab,
                 file.path(out_dir, sprintf("composition_%s_abundance_top.csv", taxrank)),
                 row.names = FALSE)
log_msg("  tables saved: composition_", taxrank, "_abundance_all.csv (", nrow(tab),
        " rows) and _top.csv (", nrow(plot_tab), " rows)")

write_session_info(file.path(out_dir, "sessionInfo.txt"))
log_msg("=== finished ===", level = "STEP")
log_msg("results: ", normalizePath(out_dir, mustWork = FALSE))
