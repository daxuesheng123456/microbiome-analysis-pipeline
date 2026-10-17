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
Microbiome co-occurrence network pipeline

Options:
  --config <path>   path to the config file (required)
  --out <dir>       override output directory
  --groups <a,b,c>  override group list (comma separated)
  --no-plot         skip figure generation
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

out_dir  <- args$out %||% cfg$output_dir %||% "results"
data_dir <- file.path(out_dir, "network_data")
plot_dir <- file.path(out_dir, "network_plots")
groups   <- if (!is.null(args$groups)) strsplit(args$groups, ",")[[1]] else cfg$groups
make_plots <- !isTRUE(args$no.plot)

ensure_dir(out_dir); ensure_dir(data_dir); ensure_dir(plot_dir)
logfile <- file.path(out_dir, "run.log")
if (file.exists(logfile)) invisible(file.remove(logfile))
log_msg <- log_to(logfile)

log_msg("=== Microbiome network pipeline ===", level = "STEP")
log_msg("config : ", args$config)
log_msg("output : ", normalizePath(out_dir, mustWork = FALSE))
log_msg("groups : ", paste(groups, collapse = ", "))

check_packages(c("phyloseq", "ggClusterNet", "igraph",
                 "ggplot2", "ggnewscale", "dplyr"))

inp <- cfg$input
ps <- tryCatch(
  load_phyloseq(otu_path   = inp$otu_table,
                tax_path   = inp$taxonomy,
                meta_path  = inp$metadata,
                id_col     = inp$id_col,
                sample_col = inp$sample_col),
  error = function(e) die("Failed to build phyloseq object: ", conditionMessage(e))
)

group_column <- cfg$group_column %||% "Group"
sd <- plain_df(phyloseq::sample_data(ps))
absent <- setdiff(groups, unique(sd[[group_column]]))
if (length(absent)) {
  die("Groups not present in metadata: ", paste(absent, collapse = ", "),
      "\n  Available: ", paste(unique(sd[[group_column]]), collapse = ", "))
}

net_cfg  <- cfg$network
plot_cfg <- cfg$plot %||% list()
results  <- list()
summary_rows <- list()

for (grp in groups) {
  log_msg("Processing group: ", grp, level = "STEP")

  net <- try_group(paste0("group ", grp),
                   build_one_network(ps, grp, net_cfg, group_column))
  if (is.null(net)) {
    log_msg("  group ", grp, " skipped", level = "WARN")
    next
  }

  results[[grp]] <- net

  utils::write.csv(net$nodes, file.path(data_dir, paste0(grp, "_nodes.csv")), row.names = FALSE)
  utils::write.csv(net$edges, file.path(data_dir, paste0(grp, "_edges.csv")), row.names = FALSE)
  utils::write.csv(net$properties, file.path(data_dir, paste0(grp, "_net_properties.csv")), row.names = FALSE)

  summary_rows[[grp]] <- data.frame(
    Group      = grp,
    n_samples  = phyloseq::nsamples(phyloseq::prune_samples(
                    rownames(sd)[sd[[group_column]] == grp], ps)),
    n_nodes    = igraph::vcount(net$igraph),
    n_edges    = igraph::ecount(net$igraph),
    n_positive = net$n_pos,
    n_negative = net$n_neg,
    modularity = tryCatch(igraph::modularity(
                    igraph::cluster_fast_greedy(
                      as_undir(net$igraph))),
                  error = function(e) NA_real_),
    stringsAsFactors = FALSE
  )

  if (make_plots) {
    p <- try_group(paste0("plot ", grp),
                   plot_one_network(net, grp, plot_cfg))
    if (!is.null(p)) {
      save_network_plot(p, file.path(plot_dir, paste0(grp, "_network.pdf")), plot_cfg)
    }
  }

  log_msg(sprintf("  done: %d nodes, %d edges (+%d / -%d)",
                  igraph::vcount(net$igraph), igraph::ecount(net$igraph),
                  net$n_pos, net$n_neg))
}

if (length(results) == 0) {
  die("No group produced a network. Check thresholds or group labels.")
}

all_props <- dplyr::bind_rows(lapply(results, `[[`, "properties"))
utils::write.csv(all_props, file.path(data_dir, "All_groups_net_properties.csv"), row.names = FALSE)

summary_df <- dplyr::bind_rows(summary_rows)
utils::write.csv(summary_df, file.path(out_dir, "network_summary.csv"), row.names = FALSE)

write_session_info(file.path(out_dir, "sessionInfo.txt"))

log_msg(sprintf("=== finished: %d/%d groups produced a network ===",
                length(results), length(groups)), level = "STEP")
log_msg("summary: ", file.path(out_dir, "network_summary.csv"))
