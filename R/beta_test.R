BETA_METHODS <- c("adonis", "MRPP", "anosim")

beta_stat_prefix <- function(method) {
  switch(method,
         adonis = "Adonis:R",
         MRPP   = "MRPP.delta",
         anosim = "ANOSIM.r",
         method)
}

beta_test <- function(ps, group_column = "Group", Micromet = "adonis",
                      dist = "bray", permutations = 999, seed = 1) {

  check_packages(c("phyloseq", "vegan"))
  suppressMessages(library(phyloseq))
  Micromet <- match.arg(Micromet, BETA_METHODS)

  ps_rel <- relative_abundance(ps)
  d      <- phyloseq::distance(ps_rel, method = dist, type = "samples")

  map  <- plain_df(phyloseq::sample_data(ps_rel))
  map  <- map[labels(d), , drop = FALSE]
  grp  <- as.factor(map[[group_column]])
  if (nlevels(grp) < 2) {
    die("Need at least two groups for a beta-diversity test")
  }

  set.seed(seed)

  if (Micromet == "adonis") {
    fit  <- vegan::adonis2(d ~ grp, data = data.frame(grp = grp),
                           permutations = permutations, by = NULL)
    stat <- fit$R2[1]
    p    <- fit$`Pr(>F)`[1]
    label <- sprintf("Adonis:R %s p: %s", round(stat, 3), round(p, 3))

  } else if (Micromet == "MRPP") {
    fit  <- vegan::mrpp(d, grp)
    stat <- fit$delta
    p    <- fit$Pvalue
    label <- sprintf("MRPP.delta %s p: %s", round(stat, 3), round(p, 3))

  } else {
    fit  <- vegan::anosim(d, grp)
    stat <- fit$statistic
    p    <- fit$signif
    label <- sprintf("ANOSIM.r %s p: %s", round(stat, 3), round(p, 3))
  }

  list(label = label, stat = as.numeric(stat), p = as.numeric(p),
       method = Micromet, dist = dist, n_groups = nlevels(grp),
       n_samples = length(grp), fit = fit)
}

pair_beta_test <- function(ps, group_column = "Group", Micromet = "adonis",
                           dist = "bray", permutations = 999, seed = 1,
                           p.adjust.method = "BH") {

  check_packages(c("phyloseq", "vegan"))
  suppressMessages(library(phyloseq))
  Micromet <- match.arg(Micromet, BETA_METHODS)

  ps_rel <- relative_abundance(ps)
  map    <- plain_df(phyloseq::sample_data(ps_rel))
  grps   <- levels(as.factor(map[[group_column]]))
  if (length(grps) < 2) die("Need at least two groups for pairwise tests")

  pairs <- utils::combn(grps, 2)
  rows  <- vector("list", ncol(pairs))

  for (i in seq_len(ncol(pairs))) {
    g1 <- pairs[1, i]; g2 <- pairs[2, i]
    keep <- rownames(map)[map[[group_column]] %in% c(g1, g2)]
    sub  <- phyloseq::prune_samples(keep, ps_rel)
    sub  <- phyloseq::prune_taxa(phyloseq::taxa_sums(sub) > 0, sub)

    res <- tryCatch(
      beta_test(sub, group_column = group_column, Micromet = Micromet,
                dist = dist, permutations = permutations, seed = seed),
      error = function(e) {
        log_msg("pair ", g1, " vs ", g2, " failed: ",
                conditionMessage(e), level = "WARN")
        list(label = NA_character_, stat = NA_real_, p = NA_real_)
      }
    )

    rows[[i]] <- data.frame(
      group1 = g1, group2 = g2, pair = paste0(g1, "_VS_", g2),
      method = Micromet, n1 = sum(map[[group_column]] == g1),
      n2 = sum(map[[group_column]] == g2),
      stat = res$stat, p = res$p,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p, method = p.adjust.method)
  out[order(out$p), , drop = FALSE]
}

MicroTest <- function(otu = NULL, tax = NULL, map = NULL, ps = NULL,
                      group = "Group", Micromet = "adonis",
                      dist = "bray", permutations = 999) {
  if (is.null(ps)) {
    if (is.null(otu) || is.null(map)) {
      die("MicroTest() needs either `ps`, or both `otu` and `map`")
    }
    check_packages(c("phyloseq", "ggClusterNet"))
    ps <- ggClusterNet::inputMicro(otu, tax, map, NULL, NULL, group = group)
  }
  beta_test(ps, group_column = group, Micromet = Micromet,
            dist = dist, permutations = permutations)$label
}

pairMicroTest <- function(ps = NULL, Micromet = "adonis", dist = "bray",
                          group = "Group", permutations = 999) {
  if (is.null(ps)) die("pairMicroTest() needs a phyloseq object `ps`")
  res <- pair_beta_test(ps, group_column = group, Micromet = Micromet,
                        dist = dist, permutations = permutations)

  data.frame(
    ID   = res$pair,
    stat = sprintf("%s %s", beta_stat_prefix(Micromet),
                   round(res$stat, 3)),
    p    = sprintf("p: %s", round(res$p, 3)),
    stat_numeric = res$stat,
    p_numeric    = res$p,
    p_adj        = res$p_adj,
    stringsAsFactors = FALSE
  )
}
