#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# BMRG gene abundance (TPM) heatmaps - Supplementary Fig. S14 / S15
#
#   Rscript scripts/run_bmrg_tpm.R
#
# Input : data/BMRG_TPM_rhizosphere.csv
# Output: results/bmrg_tpm/
#           BMRG_group_mean_TPM.csv
#           BMRG_log2FC_vs_CK.csv
#           FigS14_BMRG_abundance.pdf
#           FigS15_BMRG_log2FC.pdf
# ---------------------------------------------------------------------------

source("R/utils.R")
source("R/bmrg_tpm.R")

check_packages(c("ggplot2"))

out_dir <- file.path("results", "bmrg_tpm")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_msg("reading BMRG TPM matrix")
tpm <- read_bmrg_tpm()

means <- bmrg_group_means(tpm)
lfc   <- bmrg_log2fc(means, ref = "CK")

utils::write.csv(round(means, 2), file.path(out_dir, "BMRG_group_mean_TPM.csv"))
utils::write.csv(round(lfc, 2), file.path(out_dir, "BMRG_log2FC_vs_CK.csv"))

p1 <- plot_bmrg_abundance(means)
ggplot2::ggsave(file.path(out_dir, "FigS14_BMRG_abundance.pdf"),
                p1, width = 6.7, height = 5.0)

p2 <- plot_bmrg_log2fc(lfc, threshold = 0.2)
ggplot2::ggsave(file.path(out_dir, "FigS15_BMRG_log2FC.pdf"),
                p2, width = 6.7, height = 5.5)

log_msg("done: ", out_dir)
