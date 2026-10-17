log_msg <- function(..., level = "INFO") {
  cat(sprintf("[%s] %-5s %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              level,
              paste0(...)))
  invisible(NULL)
}

log_to <- function(logfile) {
  force(logfile)
  function(..., level = "INFO") {
    msg <- sprintf("[%s] %-5s %s\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   level, paste0(...))
    cat(msg)
    cat(msg, file = logfile, append = TRUE)
    invisible(NULL)
  }
}

die <- function(...) {
  log_msg(..., level = "ERROR")
  quit(status = 1)
}

check_packages <- function(pkgs, strict = TRUE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    msg <- paste0("Missing package(s): ", paste(missing, collapse = ", "))
    if (strict) {
      die(msg, "\n  Install with: install.packages(c(",
          paste0('"', missing, '"', collapse = ", "), "))")
    } else {
      log_msg(msg, level = "WARN")
    }
  }
  invisible(missing)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

load_config <- function(path) {
  if (!file.exists(path)) die("Config file not found: ", path)
  env <- new.env()
  sys.source(path, envir = env)
  if (!exists("config", envir = env)) {
    die("Config file must define an object named `config`: ", path)
  }
  cfg <- get("config", envir = env)
  cfg
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (grepl("^--[A-Za-z0-9_.-]+$", a)) {
      key <- sub("^--", "", a)
      if (i < length(args) && !grepl("^--", args[[i + 1]])) {
        out[[key]] <- args[[i + 1]]
        i <- i + 2
      } else {
        out[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  out
}

write_session_info <- function(path) {
  ensure_dir(dirname(path))
  writeLines(capture.output(sessionInfo()), path)
  log_msg("Session info written to ", path)
}

plain_df <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x) && !inherits(x, "sample_data")) return(x)
  out <- as.data.frame(lapply(seq_len(ncol(x)), function(i) x[[i]]),
                       stringsAsFactors = FALSE)
  names(out) <- colnames(x)
  rownames(out) <- rownames(x)
  out
}

as_undir <- function(g, mode = "collapse") {
  if (exists("as_undirected", where = asNamespace("igraph"), inherits = FALSE)) {
    igraph::as_undirected(g, mode = mode)
  } else {
    igraph::as.undirected(g, mode = mode)
  }
}

try_group <- function(label, expr) {
  res <- tryCatch(expr, error = function(e) e)
  if (inherits(res, "error")) {
    log_msg(label, " failed: ", conditionMessage(res), level = "ERROR")
    return(NULL)
  }
  res
}
