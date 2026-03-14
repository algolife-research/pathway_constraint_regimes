library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(scales)

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, args)
  if (length(match) > 0) return(dirname(normalizePath(sub(needle, "", args[match]))))
  if (exists("ofile", envir = sys.frame(1))) return(dirname(sys.frame(1)$ofile))
  return(getwd())
}
BASE    <- normalizePath(file.path(get_script_dir(), ".."))
DATA    <- file.path(BASE, "data")
OUTPUT  <- file.path(BASE, "output")
FIG_DIR <- file.path(BASE, "figures")

REGIME_COLS <- c(
  "Buffered"       = "#f4903d",
  "Equilibrium"       = "#b0b0b0",
  "Decoupled"         = "#43b1bb",
  "Topology-Amplified" = "#4dba3f"
)
REGIME_ORDER <- c("Buffered", "Equilibrium", "Decoupled", "Topology-Amplified")

theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position  = "bottom"
  )
theme_set(theme_pub)

cat("Loading data...\n")

# Pathway-level TCM model (1,275 pathways with kappa fit)
pw_kappa <- read_csv(file.path(OUTPUT, "pathway_kappa_model.csv"),
                     show_col_types = FALSE)

# Pathway-level full metrics with regime labels (1,708 pathways)
pw_regimes <- read_csv(file.path(OUTPUT, "pathway_regimes_classified.csv"),
                       show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

# Gene-level topology metrics
gene_topo <- read_csv(file.path(OUTPUT, "gene_topology_metrics.csv"),
                      show_col_types = FALSE)

# Cassa et al. 2017 s_het
shet <- read_excel(file.path(DATA,
  "Cassa_2017_41588_2017_BFng3831_MOESM71_ESM.xlsx")) %>%
  select(gene = gene_symbol, s_het) %>%
  filter(!is.na(s_het))

# Reactome GMT (gene sets)
gmt_path <- file.path(DATA, "ReactomePathways.gmt.zip")
gmt_con  <- unz(gmt_path, unzip(gmt_path, list = TRUE)$Name[1])
gmt_lines <- readLines(gmt_con)
close(gmt_con)
pw_genes <- setNames(
  lapply(gmt_lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    parts[-(1:2)]
  }),
  sapply(gmt_lines, function(x) strsplit(x, "\t")[[1]][1])
)

# GSEA results
gsea_alpha <- read_csv(file.path(OUTPUT, "gsea_alpha_rank_results.csv"),
                       show_col_types = FALSE)
gsea_kappa <- read_csv(file.path(OUTPUT, "gsea_kappa_rank_results.csv"),
                       show_col_types = FALSE)

# HPO enrichment
hpo <- read_csv(file.path(OUTPUT, "hpo_enrichment_results.csv"),
                show_col_types = FALSE)

# ClinGen overlap
clingen <- read_csv(file.path(OUTPUT, "clingen_overlap_results.csv"),
                    show_col_types = FALSE)

# Each gene is assigned to the regime of its most constrained pathway
# (lowest mean_loeuf), matching the Python pipeline logic.
gene_regime <- pw_regimes %>%
  select(pathway, regime, mean_loeuf) %>%
  inner_join(
    enframe(pw_genes, name = "pathway", value = "genes") %>%
      unnest(genes) %>%
      rename(gene = genes),
    by = "pathway"
  ) %>%
  arrange(mean_loeuf) %>%
  distinct(gene, .keep_all = TRUE) %>%
  select(gene, regime) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

# Merge gene data
gene_df <- gene_topo %>%
  mutate(degree = round(degree_centrality * (nrow(gene_topo) - 1))) %>%
  inner_join(shet, by = "gene") %>%
  inner_join(gene_regime, by = "gene")

cat(sprintf("  %d genes with degree, s_het, and regime\n", nrow(gene_df)))

# ===========================================================================
# FIGURE 1: Background cost of complexity + TCM deviation space
# ===========================================================================
cat("Figure 1: Cost of complexity + TCM deviation space\n")

pk <- pw_kappa %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

# Panel A: density vs mean s_het per pathway (background cost of complexity)
fig1a <- ggplot(pk %>% filter(density > 0), aes(density, mean_shet)) +
  geom_point(size = 1, alpha = 0.5, colour = "#b0b0b0") +
  geom_smooth(method = "lm", formula = y ~ log(x),
              se = TRUE, colour = "grey30", linewidth = 0.5) +
  labs(x = "Network density",
       y = expression("Mean " * s[het]))

# Panel B: TCM deviation space with rectangular threshold boundaries
# Compute marginal percentile thresholds matching prepare_tables.R
OUTLIER_PROB <- 0.90
resid_upper <- quantile(pk$mean_residual, OUTLIER_PROB)
resid_lower <- quantile(pk$mean_residual, 1 - OUTLIER_PROB)
kappa_upper <- quantile(pk$deviation_kappa, OUTLIER_PROB)
kappa_lower <- quantile(pk$deviation_kappa, 1 - OUTLIER_PROB)

REGIME_SHAPES <- c(
  "Buffered"            = 15,   # square
  "Equilibrium"         = 16,   # filled circle
  "Decoupled"           = 17,   # triangle
  "Topology-Amplified"  = 18    # diamond
)

fig1b <- ggplot(pk, aes(mean_residual, deviation_kappa)) +
  # Rectangular threshold lines (dashed)
  geom_hline(yintercept = c(kappa_lower, kappa_upper),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(resid_lower, resid_upper),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  # 95% ellipse for Equilibrium pathways
  stat_ellipse(data = pk %>% filter(regime == "Equilibrium"),
               level = 0.95, type = "norm",
               colour = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_point(aes(colour = regime, shape = regime), size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  scale_shape_manual(values = REGIME_SHAPES, name = "Regime") +
  labs(x = "Mean residual (pathway constraint vs. degree expectation)",
       y = expression("Deviation " * italic(b) * " (excess topology-dependence)")) +
  guides(colour = guide_legend(override.aes = list(size = 3)),
         shape  = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = "right")

fig1 <- fig1a + fig1b +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "figure1_kappa_model.png"),
       fig1, width = 10, height = 4.5, dpi = 300)

# ===========================================================================
# FIGURE 3: Topology–fitness decoupling
# ===========================================================================
cat("Figure 3: Topology-fitness decoupling\n")

# Panel A: degree vs s_het with LOESS per regime (95% CI ribbon)
fig2a <- ggplot(gene_df, aes(degree, s_het, colour = regime, fill = regime)) +
  geom_point(size = 0.3, alpha = 0.15) +
  geom_smooth(method = "loess", se = TRUE, level = 0.95,
              linewidth = 1, span = 0.5, alpha = 0.2) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = REGIME_COLS, guide = "none") +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  labs(x = "Degree (PPI partners)", y = expression(s[het]))

# Panel B: degree-binned mean s_het with 95% CI (degree tertiles)
deg_breaks <- quantile(gene_df$degree, probs = c(0, 1/3, 2/3, 1))
deg_labels <- c("Low", "Mid", "High")

gene_df_binned <- gene_df %>%
  mutate(degree_bin = cut(degree, breaks = deg_breaks,
                          labels = deg_labels, include.lowest = TRUE)) %>%
  filter(!is.na(degree_bin)) %>%
  group_by(regime, degree_bin) %>%
  summarise(
    mean_shet = mean(s_het, na.rm = TRUE),
    se_shet   = sd(s_het, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(ci_lo = pmax(0, mean_shet - 1.96 * se_shet),
         ci_hi = mean_shet + 1.96 * se_shet) %>%
  complete(regime, degree_bin)

dodge <- position_dodge(width = 0.8)

fig2b <- ggplot(gene_df_binned, aes(degree_bin, mean_shet, fill = regime)) +
  geom_col(position = dodge, width = 0.7) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = dodge, width = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = REGIME_COLS, name = "Regime") +
  labs(x = "Degree bin (log₁₀)", y = expression("Mean " * s[het])) +
  theme(legend.position = "right")

fig2 <- fig2a + fig2b +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "figure3_hub_fallacy.png"),
       fig2b + theme(legend.position = "bottom"),
       width = 7, height = 3, dpi = 300)

# ===========================================================================
# FIGURE 2: Dosage sensitivity & structural context
# ===========================================================================
cat("Figure 2: Structural context\n")
library(ggsignif)

pw_all <- pw_regimes

# Pairwise comparisons: each outlier regime vs Equilibrium
regime_pairs <- list(
  c("Buffered", "Equilibrium"),
  c("Equilibrium", "Decoupled"),
  c("Equilibrium", "Topology-Amplified"),
  c("Decoupled", "Topology-Amplified")
)

# Panel A: genomic clustering by regime (violin + boxplot, log10 y)
# Add small offset to handle zeros
fig3a <- ggplot(pw_all, aes(regime, clustering_score + 1e-5, fill = regime)) +
  geom_violin(alpha = 0.6, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, linewidth = 0.3) +
  geom_signif(comparisons = regime_pairs, map_signif_level = TRUE,
              test = "wilcox.test", textsize = 3, step_increase = 0.08,
              tip_length = 0.01) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2, 1e-1, 1),
                labels = expression(10^-4, 10^-3, 10^-2, 10^-1, 10^0)) +
  labs(x = NULL, y = "Genomic clustering score") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel B: pLI by regime with significance
fig3b <- ggplot(pw_all, aes(regime, mean_pli, fill = regime)) +
  geom_violin(alpha = 0.6, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, linewidth = 0.3) +
  geom_signif(comparisons = regime_pairs, map_signif_level = TRUE,
              test = "wilcox.test", textsize = 3, step_increase = 0.08,
              tip_length = 0.01) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  labs(x = NULL, y = "Mean pLI") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Panel C: CORUM complex fraction by regime with significance
fig3c <- ggplot(pw_all, aes(regime, complex_fraction, fill = regime)) +
  geom_violin(alpha = 0.6, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.5, linewidth = 0.3) +
  geom_signif(comparisons = regime_pairs, map_signif_level = TRUE,
              test = "wilcox.test", textsize = 3, step_increase = 0.08,
              tip_length = 0.01) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  labs(x = NULL, y = "CORUM complex fraction") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig3 <- fig3b / fig3c / fig3a +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "figure2_structural_context.png"),
       fig3, width = 6, height = 10, dpi = 300)

# ===========================================================================
# FIGURE 4: Functional & phenotypic validation
# ===========================================================================
cat("Figure 4: Functional validation\n")

# Panel A: GSEA on residual axis — top terms per direction
gsea_a_top <- gsea_alpha %>%
  filter(`FDR q-val` < 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Constrained", "Buffered")) %>%
  group_by(direction) %>%
  slice_max(abs(NES), n = 5) %>%
  ungroup() %>%
  mutate(Term = str_remove(Term, " \\(GO:\\d+\\)"),
         Term = str_trunc(Term, 45),
         Term = fct_reorder(Term, NES))

fig4a <- ggplot(gsea_a_top, aes(NES, Term,
                                 fill = ifelse(NES > 0, "Constrained", "Buffered"))) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Constrained" = "#43b1bb", "Buffered" = "#f4903d"),
                    name = "Direction") +
  labs(x = "Normalised Enrichment Score", y = NULL,
       subtitle = "Residual axis") +
  theme(legend.position = "none")

# Panel B: GSEA on kappa axis — top terms per direction
gsea_k_top <- gsea_kappa %>%
  filter(`FDR q-val` < 0.05) %>%
  mutate(direction = ifelse(NES > 0, "Topology-Amplified", "Topology-Independent")) %>%
  group_by(direction) %>%
  slice_max(abs(NES), n = 5) %>%
  ungroup() %>%
  mutate(Term = str_remove(Term, " \\(GO:\\d+\\)"),
         Term = str_trunc(Term, 45),
         Term = fct_reorder(Term, NES))

fig4b <- ggplot(gsea_k_top, aes(NES, Term,
                                 fill = ifelse(NES > 0, "Topology-Amplified",
                                               "Topology-Independent"))) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = c("Topology-Amplified" = "#4dba3f",
                               "Topology-Independent" = "#f4903d"),
                    name = "Direction") +
  labs(x = "Normalised Enrichment Score", y = NULL,
       subtitle = expression(italic(b) * " axis")) +
  theme(legend.position = "none")

# Panel C: HPO enrichment dot plot — top terms per regime with 95% CI
hpo_top <- hpo %>%
  filter(fdr < 0.1, odds_ratio > 0) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER)) %>%
  group_by(regime) %>%
  slice_min(fdr, n = 8, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    # Log OR CI from Fisher's test contingency table (Woolf method)
    log_or = log(odds_ratio),
    se_log_or = sqrt(1/pmax(regime_hits,1) + 1/pmax(regime_total - regime_hits,1) +
                     1/pmax(equil_hits,1) + 1/pmax(equil_total - equil_hits,1)),
    log2_or = log2(pmax(odds_ratio, 0.01)),
    log2_or_lo = log2(exp(log_or - 1.96 * se_log_or)),
    log2_or_hi = log2(exp(log_or + 1.96 * se_log_or)),
    hpo_short = str_trunc(hpo_name, 50),
    hpo_short = fct_reorder(hpo_short, log2_or)
  )

dodge_v <- position_dodge(width = 0.6)

fig4c <- ggplot(hpo_top, aes(log2_or, hpo_short,
                              colour = regime, size = -log10(fdr))) +
  geom_errorbarh(aes(xmin = log2_or_lo, xmax = log2_or_hi),
                 height = 0.3, linewidth = 0.3, size = NULL,
                 position = dodge_v) +
  geom_point(alpha = 0.8, position = dodge_v) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  scale_size_continuous(range = c(1.5, 4), name = expression(-log[10](FDR))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = expression(log[2](OR)), y = NULL,
       subtitle = "HPO enrichment") +
  theme(legend.position = "right")

# Panel D: ClinGen haploinsufficiency bar plot with 95% binomial CI
clingen_plot <- clingen %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER)) %>%
  filter(!is.na(regime)) %>%
  mutate(
    pct_lo = 100 * qbeta(0.025, hi_count + 0.5, total - hi_count + 0.5),
    pct_hi = 100 * qbeta(0.975, hi_count + 0.5, total - hi_count + 0.5)
  )

fig4d <- ggplot(clingen_plot, aes(regime, pct, fill = regime)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = pct_lo, ymax = pct_hi),
                width = 0.2, linewidth = 0.3) +
  geom_text(aes(y = pct_hi, label = sprintf("%.1f%%", pct)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% ClinGen HI Score 3",
       subtitle = "ClinGen haploinsufficiency") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Main figure: 2-panel (HPO + ClinGen only), single-column width
fig4 <- fig4c / fig4d +
  plot_layout(heights = c(3, 1.2)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(FIG_DIR, "figure4_functional_validation.png"),
       fig4, width = 6, height = 7, dpi = 300)
