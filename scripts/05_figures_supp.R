library(tidyverse)
library(readxl)
library(patchwork)
library(scales)
library(ggcorrplot)
library(data.table)

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
dir.create(FIG_DIR, showWarnings = FALSE)

REGIME_COLS <- c(
  "Buffered"            = "#f4903d",
  "Equilibrium"         = "#b0b0b0",
  "Decoupled"           = "#43b1bb",
  "Topology-Amplified"  = "#4dba3f"
)
REGIME_ORDER <- c("Buffered", "Equilibrium", "Decoupled", "Topology-Amplified")
REGIME_SHAPES <- c(
  "Buffered"            = 15,   # square
  "Equilibrium"         = 16,   # filled circle
  "Decoupled"           = 17,   # triangle
  "Topology-Amplified"  = 18    # diamond
)

theme_pub <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position  = "bottom"
  )
theme_set(theme_pub)

cat("Loading shared data...\n")

gene_topo <- read_csv(file.path(OUTPUT, "gene_topology_metrics.csv"),
                      show_col_types = FALSE)

# Load gene-level s_het from Cassa et al. 2017
shet_dt <- as.data.table(read_excel(
  file.path(DATA, "Cassa_2017_41588_2017_BFng3831_MOESM71_ESM.xlsx")))
shet_dt <- shet_dt[!is.na(s_het), .(gene = gene_symbol, s_het)]
gene_topo <- gene_topo %>% left_join(as_tibble(shet_dt), by = "gene")

pw_kappa <- read_csv(file.path(OUTPUT, "pathway_kappa_model.csv"),
                     show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

pw_regimes <- read_csv(file.path(OUTPUT, "pathway_regimes_classified.csv"),
                       show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

pw_regime_metrics <- read_csv(file.path(OUTPUT, "pathway_regime_metrics.csv"),
                              show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))


# ===========================================================================
# SUPPLEMENTARY FIGURE S1: The Topological Landscape
# ===========================================================================
cat("Figure S1: Topological landscape\n")

gt <- gene_topo %>% filter(!is.na(s_het))
gt <- gt %>% mutate(degree = round(degree_centrality * (nrow(gene_topo) - 1)))

# Panel A: Degree centrality vs s_het (log scale, loess)
rho_deg <- cor.test(gt$degree_centrality, gt$s_het, method = "spearman")
s1a <- ggplot(gt, aes(degree_centrality, s_het)) +
  geom_point(size = 0.3, alpha = 0.2, colour = "#4a7cb5") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "red", linewidth = 0.6) +
  scale_x_log10() +
  annotate("text", x = max(gt$degree_centrality[gt$degree_centrality > 0]) * 0.3,
           y = max(gt$s_het) * 0.95,
           label = sprintf("rho == %.2f~~p < 10^{-99}", rho_deg$estimate),
           parse = TRUE, size = 3, hjust = 1) +
  labs(x = "Degree centrality", y = expression(s[het]),
       subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel B: Clustering coefficient vs s_het (loess)
rho_cc <- cor.test(gt$clustering_coeff, gt$s_het, method = "spearman")
s1b <- ggplot(gt, aes(clustering_coeff, s_het)) +
  geom_point(size = 0.3, alpha = 0.2, colour = "#4a7cb5") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE,
              colour = "red", linewidth = 0.6) +
  annotate("text", x = max(gt$clustering_coeff, na.rm = TRUE) * 0.9,
           y = max(gt$s_het) * 0.95,
           label = sprintf("rho == %.2f", rho_cc$estimate),
           parse = TRUE, size = 3, hjust = 1) +
  labs(x = "Clustering coefficient", y = expression(s[het]),
       subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel C: K-core shell index vs s_het (boxplot, binned)
gt_kcore <- gt %>%
  mutate(kcore_bin = cut(kcore, breaks = quantile(kcore, probs = seq(0, 1, length.out = 7)),
                         include.lowest = TRUE, dig.lab = 4))
rho_kc <- cor.test(gt$kcore, gt$s_het, method = "spearman")

s1c <- ggplot(gt_kcore %>% filter(!is.na(kcore_bin)), aes(kcore_bin, s_het)) +
  geom_boxplot(aes(fill = kcore_bin), outlier.size = 0.3, linewidth = 0.3) +
  scale_fill_viridis_d(direction = -1, guide = "none") +
  annotate("text", x = 1, y = max(gt$s_het) * 0.95,
           label = sprintf("rho == %.2f", rho_kc$estimate),
           parse = TRUE, size = 3, hjust = 0) +
  labs(x = "K-core shell index", y = expression(s[het]),
       subtitle = "C") +
  theme(plot.subtitle = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 30, hjust = 1))

# Panel D: Spearman correlation matrix
corr_vars <- gt %>%
  select(Degree = degree_centrality, Clustering = clustering_coeff,
         Betweenness = betweenness, `K-Core` = kcore, s_het = s_het) %>%
  drop_na()
corr_mat <- cor(corr_vars, method = "spearman")

s1d <- ggcorrplot(corr_mat, type = "full", lab = TRUE, lab_size = 3,
                  colors = c("#2166ac", "white", "#b2182b"),
                  outline.color = "grey80",
                  title = "D") +
  theme(plot.title = element_text(face = "bold", size = 12))

fig_s1 <- (s1a | s1b) / (s1c | s1d)

ggsave(file.path(FIG_DIR, "supp_s1_topological_landscape.png"),
       fig_s1, width = 10, height = 8, dpi = 300)
cat("  Saved supp_s1_topological_landscape.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S2: Pathway Landscape & Bias Exclusion
# ===========================================================================
cat("Figure S2: Pathway landscape & bias exclusion\n")

pw_land <- pw_regimes %>%
  left_join(pw_kappa %>% select(pathway, mean_shet), by = "pathway") %>%
  filter(density > 0)

# Panel A: Density vs mean LOEUF with regime highlights
rho_pw <- cor.test(pw_land$density, pw_land$mean_loeuf, method = "spearman")

s3a <- ggplot(pw_land, aes(density, mean_loeuf)) +
  geom_point(size = 1, alpha = 0.3, colour = "grey60") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linetype = "dashed", linewidth = 0.5) +
  scale_x_log10() +
  annotate("text", x = max(pw_land$density) * 0.5,
           y = max(pw_land$mean_loeuf) * 0.95,
           label = sprintf("rho == %.2f~~p == %.1e", rho_pw$estimate, rho_pw$p.value),
           parse = TRUE, size = 3, hjust = 1) +
  labs(x = "Pathway network density (intra-module PPI)",
       y = "Mean LOEUF", subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel B: Partial residual plot (density after removing confounders)
ols_data <- pw_land %>%
  filter(!is.na(mean_paralog), !is.na(mean_cds_length),
         !is.na(clustering_score), !is.na(mean_shet)) %>%
  mutate(across(c(density, mean_paralog, mean_cds_length, clustering_score),
                ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE),
                .names = "{.col}_z"))

fit_partial <- lm(mean_shet ~ density_z + mean_paralog_z + mean_cds_length_z +
                    clustering_score_z, data = ols_data)
ols_data$partial_resid <- residuals(fit_partial) +
  coef(fit_partial)["density_z"] * ols_data$density_z

r_partial_dens <- cor.test(ols_data$density_z, ols_data$partial_resid)

s3b <- ggplot(ols_data, aes(density_z, partial_resid)) +
  geom_point(size = 0.5, alpha = 0.4, colour = "#4a7cb5") +
  geom_smooth(method = "lm", se = FALSE, colour = "black",
              linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = min(ols_data$density_z) * 0.5,
           y = max(ols_data$partial_resid) * 0.95,
           label = sprintf("Partial~r == %.2f", r_partial_dens$estimate),
           parse = TRUE, size = 3, hjust = 0) +
  labs(x = "Density (standardised, confounders removed)",
       y = expression("Partial residual (mean " * s[het] * ")"),
       subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel C: Landscape coloured by mean paralog count
s3c <- ggplot(pw_land %>% filter(!is.na(mean_paralog)),
              aes(density, mean_shet, colour = mean_paralog)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = "black", linetype = "dashed", linewidth = 0.5) +
  scale_colour_distiller(palette = "YlOrRd", direction = 1,
                         name = "Mean paralog\ncount") +
  scale_x_log10() +
  labs(x = "Pathway network density",
       y = expression("Mean " * s[het]), subtitle = "C") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

fig_s3 <- s3a + s3b + s3c + plot_layout(ncol = 3)

ggsave(file.path(FIG_DIR, "supp_s2_pathway_landscape.png"),
       fig_s3, width = 14, height = 4.5, dpi = 300)
cat("  Saved supp_s2_pathway_landscape.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S3: DFE-Parameter Validation
# ===========================================================================
cat("Figure S3: DFE-parameter validation\n")

pw_dfe <- read_csv(file.path(OUTPUT, "pathway_dfe_validation.csv"),
                   show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

# Panel A: Density vs mean s_het
pw_dfe_shet <- pw_dfe %>% filter(!is.na(mean_shet))
rho_shet <- cor.test(pw_dfe_shet$density, pw_dfe_shet$mean_shet, method = "spearman")

s5a <- ggplot(pw_dfe_shet, aes(density, mean_shet, colour = regime)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              colour = "black", linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  annotate("text", x = max(pw_dfe_shet$density) * 0.6,
           y = max(pw_dfe_shet$mean_shet, na.rm = TRUE) * 0.95,
           label = sprintf("rho == %.2f~~p == %.1e",
                           rho_shet$estimate, rho_shet$p.value),
           parse = TRUE, size = 3, hjust = 1) +
  labs(x = "Network density",
       y = expression("Mean " * s[het]), subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12),
        legend.position = "none")

# Panel B: Density vs mean missense Z
pw_dfe_misz <- pw_dfe %>% filter(!is.na(mean_missense_z))
rho_misz <- cor.test(pw_dfe_misz$density, pw_dfe_misz$mean_missense_z, method = "spearman")

s5b <- ggplot(pw_dfe_misz, aes(density, mean_missense_z, colour = regime)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              colour = "black", linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  annotate("text", x = max(pw_dfe_misz$density) * 0.6,
           y = max(pw_dfe_misz$mean_missense_z, na.rm = TRUE) * 0.95,
           label = sprintf("rho == %.2f~~p == %.1e",
                           rho_misz$estimate, rho_misz$p.value),
           parse = TRUE, size = 3, hjust = 1) +
  labs(x = "Network density",
       y = "Mean missense Z-score", subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12),
        legend.position = "none")

# Panel C: Mean missense Z by regime (boxplot)
kw_misz <- kruskal.test(mean_missense_z ~ regime, data = pw_dfe_misz)

s5c <- ggplot(pw_dfe_misz, aes(regime, mean_missense_z, fill = regime)) +
  geom_boxplot(alpha = 0.6, outlier.size = 0.5, linewidth = 0.3) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  annotate("text", x = 2, y = max(pw_dfe_misz$mean_missense_z, na.rm = TRUE) * 1.05,
           label = sprintf("Kruskal-Wallis~p == %.1e", kw_misz$p.value),
           parse = TRUE, size = 3) +
  labs(x = NULL, y = "Mean missense Z-score", subtitle = "C") +
  theme(plot.subtitle = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 30, hjust = 1))

fig_s5 <- s5a + s5b + s5c + plot_layout(ncol = 3)

ggsave(file.path(FIG_DIR, "supp_s3_dfe_validation.png"),
       fig_s5, width = 14, height = 4.5, dpi = 300)
cat("  Saved supp_s3_dfe_validation.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S4: Permutation Null Model
# ===========================================================================
cat("Figure S4: Permutation null\n")

perm <- read_csv(file.path(OUTPUT, "supp_permutation_null.csv"),
                 show_col_types = FALSE)

# Observed values (from pw_kappa)
obs_n_outlier <- sum(pw_kappa$regime != "Equilibrium")
obs_spread <- sd(pw_kappa$mean_residual, na.rm = TRUE)

p_outlier <- (sum(perm$n_outlier >= obs_n_outlier) + 1) / (nrow(perm) + 1)
p_spread <- (sum(perm$spread_mr >= obs_spread) + 1) / (nrow(perm) + 1)

# Panel A: Distribution of outlier counts
s8a <- ggplot(perm, aes(n_outlier)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "grey40", alpha = 0.7) +
  geom_vline(xintercept = obs_n_outlier, colour = "red",
             linewidth = 0.8, linetype = "solid") +
  annotate("text", x = obs_n_outlier * 0.95,
           y = Inf, vjust = 1.5,
           label = sprintf("Observed = %d\np < %.3f", obs_n_outlier, p_outlier),
           colour = "red", size = 3, hjust = 1) +
  labs(x = "Number of outlier pathways",
       y = "Permutation count", subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel B: Distribution of mean_residual spread
s8b <- ggplot(perm, aes(spread_mr)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "grey40", alpha = 0.7) +
  geom_vline(xintercept = obs_spread, colour = "red",
             linewidth = 0.8, linetype = "solid") +
  annotate("text", x = obs_spread * 0.98,
           y = Inf, vjust = 1.5,
           label = sprintf("Observed = %.4f\np < %.3f", obs_spread, p_spread),
           colour = "red", size = 3, hjust = 1) +
  labs(x = "SD of pathway mean residual",
       y = "Permutation count", subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

fig_s8 <- s8a | s8b

ggsave(file.path(FIG_DIR, "supp_s4_permutation_null.png"),
       fig_s8, width = 10, height = 4, dpi = 300)
cat("  Saved supp_s4_permutation_null.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S5: Jackknife Stability
# ===========================================================================
cat("Figure S5: Jackknife stability\n")

jk <- read_csv(file.path(OUTPUT, "supp_jackknife_stability.csv"),
               show_col_types = FALSE) %>%
  mutate(original_regime = factor(original_regime, levels = REGIME_ORDER))

# Panel A: Stability distribution by regime (overlapping histograms)
s9a <- ggplot(jk, aes(stability, fill = original_regime)) +
  geom_histogram(bins = 20, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0.9, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  scale_fill_manual(values = REGIME_COLS, name = "Regime") +
  labs(x = "Jackknife stability (fraction retaining label)",
       y = "Number of pathways", subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel B: Stability vs pathway size
rho_jk <- cor.test(jk$n_genes, jk$stability, method = "spearman")

s9b <- ggplot(jk, aes(n_genes, stability, colour = original_regime)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  scale_x_log10() +
  labs(x = "Number of genes in pathway",
       y = "Jackknife stability", subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

fig_s9 <- s9a | s9b

ggsave(file.path(FIG_DIR, "supp_s5_jackknife_stability.png"),
       fig_s9, width = 10, height = 4, dpi = 300)
cat("  Saved supp_s5_jackknife_stability.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S6: STRING Network Replication
# ===========================================================================
cat("Figure S6: STRING network replication\n")

string_rep <- read_csv(file.path(OUTPUT, "supp_string_replication.csv"),
                       show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

string_conc <- read_csv(file.path(OUTPUT, "supp_string_concordance.csv"),
                        show_col_types = FALSE) %>%
  mutate(regime_biogrid = factor(regime_biogrid, levels = REGIME_ORDER),
         regime_string = factor(regime_string, levels = REGIME_ORDER))

# Panel A: STRING deviation space with rectangular thresholds
string_mr <- string_rep$mean_residual[!is.na(string_rep$mean_residual)]
string_dk <- string_rep$deviation_kappa[!is.na(string_rep$deviation_kappa)]
s_resid_upper <- quantile(string_mr, 0.90)
s_resid_lower <- quantile(string_mr, 0.10)
s_kappa_upper <- quantile(string_dk, 0.90)
s_kappa_lower <- quantile(string_dk, 0.10)

s11a <- ggplot(string_rep, aes(mean_residual, deviation_kappa)) +
  geom_hline(yintercept = c(s_kappa_lower, s_kappa_upper),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(s_resid_lower, s_resid_upper),
             linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  stat_ellipse(data = string_rep %>% filter(regime == "Equilibrium"),
               level = 0.95, type = "norm",
               colour = "grey50", linetype = "solid", linewidth = 0.4) +
  geom_point(aes(colour = regime, shape = regime), size = 1.2, alpha = 0.7) +
  scale_colour_manual(values = REGIME_COLS, name = "Regime") +
  scale_shape_manual(values = REGIME_SHAPES, name = "Regime") +
  labs(x = "Mean residual", y = expression("Deviation " * italic(b)),
       subtitle = "A") +
  guides(colour = guide_legend(override.aes = list(size = 3)),
         shape  = guide_legend(override.aes = list(size = 3))) +
  theme(plot.subtitle = element_text(face = "bold", size = 12),
        legend.position = "right")

# Panel B: Mean residual concordance (BioGRID vs STRING)
r_mr <- cor.test(string_conc$mean_residual_biogrid,
                 string_conc$mean_residual_string)

p_mr_label <- if (r_mr$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("%.1e", r_mr$p.value)

s11b <- ggplot(string_conc, aes(mean_residual_biogrid, mean_residual_string)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.4) +
  geom_point(size = 1, alpha = 0.4, colour = "grey50") +
  annotate("text", x = min(string_conc$mean_residual_biogrid, na.rm = TRUE) * 0.5,
           y = max(string_conc$mean_residual_string, na.rm = TRUE) * 0.9,
           label = sprintf("r == %.2f~~%s", r_mr$estimate, p_mr_label),
           parse = TRUE, size = 3, hjust = 0) +
  labs(x = "Mean residual (BioGRID)", y = "Mean residual (STRING)",
       subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel C: Deviation kappa concordance
r_dk <- cor.test(string_conc$deviation_kappa_biogrid,
                 string_conc$deviation_kappa_string)

p_dk_label <- if (r_dk$p.value < 2.2e-16) "p < 2.2e-16" else sprintf("%.1e", r_dk$p.value)

s11c <- ggplot(string_conc, aes(deviation_kappa_biogrid, deviation_kappa_string)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.4) +
  geom_point(size = 1, alpha = 0.4, colour = "grey50") +
  annotate("text", x = min(string_conc$deviation_kappa_biogrid, na.rm = TRUE) * 0.5,
           y = max(string_conc$deviation_kappa_string, na.rm = TRUE) * 0.9,
           label = sprintf("r == %.2f~~%s", r_dk$estimate, p_dk_label),
           parse = TRUE, size = 3, hjust = 0) +
  labs(x = expression("Deviation " * italic(b) * " (BioGRID)"),
       y = expression("Deviation " * italic(b) * " (STRING)"),
       subtitle = "C") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

fig_s11 <- s11a + s11b + s11c + plot_layout(ncol = 3, widths = c(1.3, 1, 1))

ggsave(file.path(FIG_DIR, "supp_s6_string_replication.png"),
       fig_s11, width = 14, height = 4.5, dpi = 300)
cat("  Saved supp_s6_string_replication.png\n")


# ===========================================================================
# SUPPLEMENTARY FIGURE S7: Pathway Independence (Jaccard Clustering)
# ===========================================================================
cat("Figure S7: Pathway independence\n")

pw_indep <- read_csv(file.path(OUTPUT, "supp_pathway_independence.csv"),
                     show_col_types = FALSE) %>%
  mutate(regime = factor(regime, levels = REGIME_ORDER))

# Panel A: Raw vs effective pathway counts per regime
regime_counts <- pw_indep %>%
  group_by(regime) %>%
  summarise(
    raw = n(),
    effective = n_distinct(cluster_id),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(raw, effective), names_to = "count_type", values_to = "n") %>%
  mutate(count_type = factor(count_type, levels = c("raw", "effective")),
         regime_short = recode(regime,
                               "Buffered" = "Buff.",
                               "Equilibrium" = "Equil.",
                               "Decoupled" = "Decoup.",
                               "Topology-Amplified" = "Topo-Amp."))

s12a <- ggplot(regime_counts, aes(regime_short, n, fill = regime,
                                   alpha = count_type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = REGIME_COLS, guide = "none") +
  scale_alpha_manual(values = c("raw" = 0.4, "effective" = 0.9),
                     name = NULL,
                     labels = c("Raw", "Effective (clustered)")) +
  labs(x = NULL, y = "Number of pathways", subtitle = "A") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

# Panel B: Cluster size distribution (outlier regimes only)
outlier_clusters <- pw_indep %>%
  filter(regime != "Equilibrium") %>%
  group_by(regime, cluster_id) %>%
  summarise(cluster_size = n(), .groups = "drop")

s12b <- ggplot(outlier_clusters, aes(factor(cluster_size), fill = regime)) +
  geom_bar(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  scale_fill_manual(values = REGIME_COLS, name = "Regime") +
  labs(x = "Cluster size (pathways per cluster)",
       y = "Number of clusters", subtitle = "B") +
  theme(plot.subtitle = element_text(face = "bold", size = 12))

fig_s12 <- s12a | s12b

ggsave(file.path(FIG_DIR, "supp_s7_pathway_independence.png"),
       fig_s12, width = 10, height = 4.5, dpi = 300)
cat("  Saved supp_s7_pathway_independence.png\n")

