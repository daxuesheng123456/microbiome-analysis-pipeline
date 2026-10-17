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

if (isTRUE(args$help) || (is.null(args$config) && is.null(args$env))) {
  cat("
Environment-correlation / OTU-environment network plot

Options:
  --config <path>   config file (or use --env/--otu directly)
  --env <path>      environmental table: first column = sample ID
  --otu <path>      abundance table: first column = feature ID
  --top <n>         keep the n most abundant features (0 = all)
  --out <dir>       output directory
  --help            show this message
")
  quit(status = 0)
}

setwd(repo_root)
log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              level, paste0(...)))
}

cfg <- if (!is.null(args$config)) load_config(args$config) else list()
ec  <- cfg$envcor %||% list()

env_path <- args$env %||% ec$env_table
otu_path <- args$otu %||% ec$otu_table
top_n    <- as.integer(args$top %||% ec$top_otus %||% 10)
if (!is.na(top_n) && top_n <= 0) top_n <- NULL
out_dir  <- args$out %||% ec$output_dir %||%
            file.path(cfg$output_dir %||% "results", "envcor")

ensure_dir(out_dir)
logfile <- file.path(out_dir, "envcor.log")
if (file.exists(logfile)) invisible(file.remove(logfile))
log_msg <- log_to(logfile)

log_msg("=== Environment correlation / network ===", level = "STEP")
log_msg("env : ", env_path)
log_msg("otu : ", otu_path)
log_msg("top : ", if (is.null(top_n)) "all features" else top_n)
log_msg("out : ", normalizePath(out_dir, mustWork = FALSE))

check_packages(c("ggClusterNet", "ggplot2", "dplyr"))

env <- read_env_table(env_path)
otu <- read_otu_table(otu_path)

log_msg(sprintf("env: %d samples x %d variables", nrow(env), ncol(env)))
log_msg(sprintf("otu: %d features x %d samples", nrow(otu), ncol(otu)))
log_msg("variables: ", paste(colnames(env), collapse = ", "))

p <- tryCatch(
  matcor_plot(
    env        = env,
    otu        = otu,
    group      = ec$group %||% "bac",
    top_n      = top_n,
    method_cor = ec$method_cor %||% "spearman",
    cor_p      = ec$cor_p      %||% 0.05,
    p_thur     = ec$p_thur     %||% 0.05,
    onlysig    = isTRUE(ec$onlysig %||% TRUE),
    edge_cols  = ec$edge_cols  %||% c("blue", "red"),
    fill_palette = ec$fill_palette %||% "PRGn",
    diag       = isFALSE(ec$diag %||% FALSE),
    range      = ec$range      %||% 0.1,
    numpoint   = ec$numpoint   %||% 21,
    numpoint2  = ec$numpoint2  %||% 21,
    curvature  = ec$curvature  %||% 0.1,
    shownum    = isTRUE(ec$shownum %||% FALSE),
    sig        = TRUE,
    siglabel   = FALSE,
    lacx       = ec$lacx %||% "left",
    lacy       = ec$lacy %||% "bottom"
  ),
  error = function(e) die("matcor_plot failed: ", conditionMessage(e))
)

w <- ec$width  %||% 12
h <- ec$height %||% 8
pdf_f <- file.path(out_dir, "envcor_matrix_network.pdf")
png_f <- file.path(out_dir, "envcor_matrix_network.png")
ggplot2::ggsave(pdf_f, p, width = w, height = h)
ggplot2::ggsave(png_f, p, width = w, height = h, dpi = 300)
log_msg("figure: ", basename(pdf_f), " / ", basename(png_f),
        sprintf("  (%g x %g in)", w, h))

utils::write.csv(env, file.path(out_dir, "env_variables_used.csv"))
write_session_info(file.path(out_dir, "sessionInfo.txt"))

log_msg("=== finished ===", level = "STEP")
