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
Beta diversity module

Options:
  --config <path>   path to the config file (required)
  --out <dir>       override the beta output directory
  --dist <name>     distance metric, default bray
  --pairwise        also run pairwise PERMANOVA
  --labels          label individual samples on the figure
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

dist_metric  <- args$dist %||% cfg$beta$dist %||% "bray"
out_dir      <- args$out  %||% cfg$beta$output_dir %||% file.path(cfg$output_dir, "beta")
group_column <- cfg$group_column %||% "Group"
pal          <- cfg$beta$palette
show_labels  <- isTRUE(args$labels)

ensure_dir(out_dir)
logfile <- file.path(out_dir, "beta.log")
if (file.exists(logfile)) invisible(file.remove(logfile))
log_msg <- log_to(logfile)

log_msg("=== Beta diversity (PCoA) ===", level = "STEP")
log_msg("config : ", args$config)
log_msg("dist   : ", dist_metric)
log_msg("output : ", normalizePath(out_dir, mustWork = FALSE))

check_packages(c("phyloseq", "vegan", "ggplot2", "ggrepel"))

inp <- cfg$input
ps <- tryCatch(
  load_phyloseq(otu_path   = inp$otu_table,
                tax_path   = inp$taxonomy,
                meta_path  = inp$metadata,
                id_col     = inp$id_col,
                sample_col = inp$sample_col),
  error = function(e) die("Failed to build phyloseq object: ", conditionMessage(e))
)

sd <- plain_df(phyloseq::sample_data(ps))
if (!group_column %in% colnames(sd)) {
  die("Group column '", group_column, "' not found. Columns: ",
      paste(colnames(sd), collapse = ", "))
}
n_per_group <- table(sd[[group_column]])
log_msg("samples per group: ",
        paste(names(n_per_group), n_per_group, sep = "=", collapse = ", "))

ps_rel <- relative_abundance(ps)
ord    <- run_pcoa(ps_rel, dist = dist_metric)
varpct <- pcoa_variance(ord$eig, 2)

log_msg(sprintf("PCoA 1 = %.2f%%, PCoA 2 = %.2f%% (positive eigenvalues only)",
                varpct[1], varpct[2]))

perm     <- cfg$beta$permutations %||% 999
micromet <- cfg$beta$micromet    %||% "adonis"
test <- tryCatch(
  beta_test(ps_rel, group_column = group_column, Micromet = micromet,
            dist = dist_metric, permutations = perm),
  error = function(e) {
    log_msg(micromet, " test failed: ", conditionMessage(e), level = "WARN")
    list(label = NULL, stat = NA_real_, p = NA_real_, method = micromet)
  }
)
if (!is.null(test$label)) log_msg("  ", test$label)

pairwise <- NULL
if (isTRUE(args$pairwise)) {
  log_msg("Running pairwise ", micromet, " ...")
  pairwise <- tryCatch(
    pair_beta_test(ps_rel, group_column = group_column, Micromet = micromet,
                   dist = dist_metric, permutations = perm),
    error = function(e) {
      log_msg("pairwise test failed: ", conditionMessage(e), level = "WARN")
      NULL
    }
  )
  if (!is.null(pairwise)) {
    log_msg(sprintf("  %d pairs; %d with p < 0.05 before correction",
                    nrow(pairwise), sum(pairwise$p < 0.05, na.rm = TRUE)))
  }
}

p <- plot_pcoa_spider(ord$points, var_pct = varpct, title = test$label,
                      pal = pal, show_labels = show_labels)

pdf_path <- file.path(out_dir, sprintf("beta_PCoA_%s.pdf", dist_metric))
png_path <- file.path(out_dir, sprintf("beta_PCoA_%s.png", dist_metric))
ggplot2::ggsave(pdf_path, p, width = 7, height = 6)

if (isTRUE(args$labels)) {
  pdf_path <- file.path(out_dir, sprintf("beta_PCoA_%s_labelled.pdf", dist_metric))
  ggplot2::ggsave(pdf_path, p, width = 7, height = 6)
}

ggplot2::ggsave(sub("\\.pdf$", ".png", pdf_path), p, width = 7, height = 6, dpi = 300)
log_msg("  figure saved: ", basename(pdf_path))

utils::write.csv(ord$points, file.path(out_dir, "beta_PCoA_scores.csv"), row.names = FALSE)
utils::write.csv(group_centroids(ord$points),
                 file.path(out_dir, "beta_PCoA_centroids.csv"), row.names = FALSE)
utils::write.csv(
  data.frame(axis = c("PCoA1", "PCoA2"),
             eigenvalue = ord$eig[1:2],
             variance_pct = varpct),
  file.path(out_dir, "beta_PCoA_eigenvalues.csv"), row.names = FALSE
)
if (!is.null(test$label)) {
  utils::write.csv(
    data.frame(test = paste0(toupper(micromet), " (overall)"),
               statistic = test$stat, p = test$p,
               distance = dist_metric, permutations = perm,
               label = test$label, stringsAsFactors = FALSE),
    file.path(out_dir, sprintf("beta_%s_overall.csv", micromet)), row.names = FALSE
  )
}
if (!is.null(pairwise)) {
  utils::write.csv(pairwise,
                   file.path(out_dir, sprintf("beta_%s_pairwise.csv", micromet)),
                   row.names = FALSE)
}

write_session_info(file.path(out_dir, "sessionInfo.txt"))
log_msg("=== finished ===", level = "STEP")
log_msg("results: ", normalizePath(out_dir, mustWork = FALSE))
