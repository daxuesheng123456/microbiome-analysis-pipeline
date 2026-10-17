build_one_network <- function(ps, grp, cfg, group_column = "Group") {

  suppressMessages({
    library(phyloseq)
    library(igraph)
    library(ggClusterNet)
  })

  sd <- plain_df(phyloseq::sample_data(ps))
  if (!group_column %in% colnames(sd)) {
    die("Group column '", group_column, "' not found in metadata")
  }
  keep <- rownames(sd)[sd[[group_column]] == grp]
  if (length(keep) < 3) {
    log_msg("Group ", grp, ": only ", length(keep),
            " sample(s) - need >= 3, skipping", level = "WARN")
    return(NULL)
  }

  ps_grp <- phyloseq::prune_samples(keep, ps)
  ps_grp <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_grp) > 0, ps_grp)

  top_n <- cfg$top_otus %||% 2000
  if (length(top_n) != 1 || is.na(top_n) || !is.numeric(top_n) || top_n <= 0) {
    top_n <- 2000
    log_msg("  group ", grp, ": top_otus not set, defaulting to ", top_n,
            level = "WARN")
  }
  if (phyloseq::ntaxa(ps_grp) > top_n) {
    ps_grp <- ggClusterNet::filter_OTU_ps(ps_grp, Top = top_n)
  }

  log_msg(sprintf("  group %s: %d samples, %d features (top_otus = %d)",
                  grp, phyloseq::nsamples(ps_grp), phyloseq::ntaxa(ps_grp), top_n))

  result <- ggClusterNet::cor_Big_micro2(
    ps           = ps_grp,
    N            = 0,
    r.threshold  = cfg$r_threshold,
    p.threshold  = cfg$p_threshold,
    method       = cfg$method,
    scale        = cfg$scale
  )
  cor_mat <- result[[1]]

  result2 <- ggClusterNet::model_igraph2(
    cor    = cor_mat,
    method = cfg$cluster_method,
    seed   = cfg$seed
  )
  node <- result2[[1]]
  dat  <- result2[[2]]

  edge <- ggClusterNet::edgeBuild(cor = cor_mat, node = node)
  if (nrow(edge) == 0) {
    log_msg("Group ", grp, ": no edges passed thresholds, skipping", level = "WARN")
    return(NULL)
  }
  if (ncol(edge) >= 8) colnames(edge)[8] <- "cor"

  otu_table_df <- as.data.frame(t(vegan_otu(ps_grp)))
  tax_table_df <- as.data.frame(vegan_tax(ps_grp))
  nodes <- ggClusterNet::nodeadd(plotcord  = node,
                                 otu_table = otu_table_df,
                                 tax_table = tax_table_df)

  edge_for_graph <- data.frame(
    from   = edge$OTU_1,
    to     = edge$OTU_2,
    sign   = ifelse(edge$cor == "+", 1, -1),
    stringsAsFactors = FALSE
  )
  g <- igraph::graph_from_data_frame(edge_for_graph, directed = FALSE)
  g <- igraph::simplify(g, edge.attr.comb = "first")

  node_info <- as.data.frame(node)
  node_info$name <- rownames(node_info)
  idx <- match(igraph::V(g)$name, node_info$name)
  g  <- igraph::set_vertex_attr(g, "X1", value = node_info$X1[idx])
  g  <- igraph::set_vertex_attr(g, "X2", value = node_info$X2[idx])
  g  <- igraph::set_vertex_attr(g, "model", value = node_info$model[idx])

  n_pos <- sum(igraph::E(g)$sign > 0)
  n_neg <- sum(igraph::E(g)$sign < 0)

  g_prop <- g
  igraph::E(g_prop)$weight <- 1
  props <- tryCatch(
    as.data.frame(ggClusterNet::net_properties(g_prop)),
    error = function(e) {
      log_msg("net_properties failed for ", grp, ": ",
              conditionMessage(e), level = "WARN")
      data.frame(value = NA_real_, row.names = "unavailable")
    }
  )
  props$Property <- rownames(props)
  props$Group    <- grp
  rownames(props) <- NULL

  props <- rbind(
    props,
    data.frame(value = n_pos, Property = "num.pos.edges.raw", Group = grp),
    data.frame(value = n_neg, Property = "num.neg.edges.raw", Group = grp)
  )

  list(igraph = g, nodes = nodes, edges = edge,
       layout = dat, properties = props,
       n_pos = n_pos, n_neg = n_neg)
}
