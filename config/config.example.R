config <- list(
  input = list(
    otu_table = "data/taxonomy.csv",
    taxonomy  = NULL,
    metadata  = "data/metadata.csv",
    id_col     = NULL,
    sample_col = NULL
  ),
  group_column = "Group",
  groups = c("CK", "RSB", "GRSB", "FRSB", "RHB", "GRHB", "FRHB"),
  network = list(
    top_otus       = 2000,
    r_threshold    = 0.9,
    p_threshold    = 0.05,
    method         = "pearson",
    scale          = FALSE,
    cluster_method = "cluster_fast_greedy",
    seed           = 12
  ),
  plot = list(
    width      = 12,
    height     = 10,
    point_size = 3,
    line_width = 0.5,
    pos_col    = "#D95F02",
    neg_col    = "#1B9E77"
  ),
  beta = list(
    dist         = "bray",
    micromet     = "adonis",
    permutations = 999,
    output_dir   = NULL,
    palette      = NULL
  ),
  composition = list(
    taxrank    = "Species",
    ntaxa      = 8,
    output_dir = NULL
  ),
  envcor = list(
    env_table    = NULL,
    otu_table    = NULL,
    group        = "bac",
    top_otus     = 10,
    method_cor   = "spearman",
    cor_p        = 0.05,
    p_thur       = 0.05,
    onlysig      = TRUE,
    edge_cols    = c("blue", "red"),
    fill_palette = "PRGn",
    diag         = FALSE,
    range        = 0.1,
    numpoint     = 21,
    numpoint2    = 21,
    curvature    = 0.1,
    shownum      = FALSE,
    lacx         = "left",
    lacy         = "bottom",
    width        = 12,
    height       = 8,
    output_dir   = NULL
  ),
  output_dir = "results"
)
