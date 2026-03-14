suppressPackageStartupMessages({
  library(data.table)
})

# ---- paths ---------------------------------------------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, args)
  if (length(match) > 0) return(dirname(normalizePath(sub(needle, "", args[match]))))
  if (exists("ofile", envir = sys.frame(1))) return(dirname(sys.frame(1)$ofile))
  return(getwd())
}
BASE      <- normalizePath(file.path(get_script_dir(), ".."))
OUTPUT    <- file.path(BASE, "output")
MANU_DIR  <- file.path(BASE, "manuscript")

REGIME_ORDER <- c("Buffered", "Equilibrium", "Decoupled", "Topology-Amplified")

# ---- load data -----------------------------------------------------------
cat("Loading data...\n")
pw_kappa      <- fread(file.path(OUTPUT, "pathway_kappa_model.csv"))
pw_classified <- fread(file.path(OUTPUT, "pathway_regimes_classified.csv"))
pw_indep      <- fread(file.path(OUTPUT, "supp_pathway_independence.csv"))

# Merge independence cluster info onto classified pathways
pw_all <- merge(pw_classified, pw_kappa[, .(pathway, mean_residual, deviation_kappa,
                                             mean_shet, mean_degree, regime)],
                by = "pathway", suffixes = c("", ".kappa"))
# Use regime from kappa (the TCM classification)
pw_all[, regime := regime.kappa]
pw_all[is.na(regime), regime := "Equilibrium"]

# Escape LaTeX special characters
escape_latex <- function(x) {
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("#", "\\\\#", x)
  x
}

# Format numbers compactly
fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "---", formatC(x, format = "f", digits = digits))
}


# ==========================================================================
# TABLE S1: Regime summary statistics
# ==========================================================================
cat("=== Table S1: Regime summary statistics ===\n")

summary_dt <- pw_all[regime %in% REGIME_ORDER, {
  list(
    n_pathways = .N,
    n_clusters = uniqueN(pw_indep[regime == .BY$regime]$cluster_id),
    med_mr     = median(mean_residual, na.rm = TRUE),
    iqr_mr     = IQR(mean_residual, na.rm = TRUE),
    med_dk     = median(deviation_kappa, na.rm = TRUE),
    iqr_dk     = IQR(deviation_kappa, na.rm = TRUE),
    med_shet   = median(mean_shet, na.rm = TRUE),
    iqr_shet   = IQR(mean_shet, na.rm = TRUE),
    med_pli    = median(mean_pli, na.rm = TRUE),
    iqr_pli    = IQR(mean_pli, na.rm = TRUE),
    med_cf     = median(complex_fraction, na.rm = TRUE),
    iqr_cf     = IQR(complex_fraction, na.rm = TRUE),
    med_ngenes = median(as.numeric(n_genes), na.rm = TRUE)
  )
}, by = regime]

# Order
summary_dt[, regime := factor(regime, levels = REGIME_ORDER)]
setorder(summary_dt, regime)

# Build LaTeX
tex_s1 <- c(
  "\\begin{table}[!ht]",
  "\\centering",
  paste0("\\caption[Summary statistics by regime]{\\textbf{Summary statistics by regime.} ",
         "Median (IQR) for key metrics across the four regimes. ",
         "$n$: number of pathways; clusters: independent pathway clusters ",
         "after collapsing by Jaccard similarity $\\geq 0.50$.}\\label{tab:regime_summary}"),
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{lcccccccc}",
  "\\hline",
  paste0("Regime & $n$ & Clusters & Mean residual & Deviation $b$ & ",
         "Mean $s_\\mathrm{het}$ & Mean pLI & Complex frac. & Med.\\ genes \\\\"),
  "\\hline"
)

for (i in seq_len(nrow(summary_dt))) {
  r <- summary_dt[i]
  row <- sprintf(
    "%s & %d & %d & %s (%s) & %s (%s) & %s (%s) & %s (%s) & %s (%s) & %.0f \\\\",
    escape_latex(as.character(r$regime)),
    r$n_pathways, r$n_clusters,
    fmt(r$med_mr, 3), fmt(r$iqr_mr, 3),
    fmt(r$med_dk, 3), fmt(r$iqr_dk, 3),
    fmt(r$med_shet, 3), fmt(r$iqr_shet, 3),
    fmt(r$med_pli, 2), fmt(r$iqr_pli, 2),
    fmt(r$med_cf, 2), fmt(r$iqr_cf, 2),
    r$med_ngenes
  )
  tex_s1 <- c(tex_s1, row)
}

tex_s1 <- c(tex_s1, "\\hline", "\\end{tabular}}", "\\end{table}")

writeLines(tex_s1, file.path(MANU_DIR, "supp_table_regime_summary.tex"))
cat("  Saved supp_table_regime_summary.tex\n")


# ==========================================================================
# TABLE S2: Outlier pathway catalogue
# ==========================================================================
cat("=== Table S2: Outlier pathway catalogue ===\n")

outliers <- pw_all[regime != "Equilibrium"]
outliers[, regime := factor(regime, levels = REGIME_ORDER)]
setorder(outliers, regime, -mean_residual)

# Merge cluster info (sort=FALSE to preserve order, then re-sort)
outliers <- merge(outliers,
                  pw_indep[, .(pathway, cluster_id, is_cluster_representative)],
                  by = "pathway", all.x = TRUE, sort = FALSE)
outliers[, regime := factor(regime, levels = REGIME_ORDER)]
setorder(outliers, regime, -mean_residual)

# Build longtable (spans multiple pages)
tex_s2 <- c(
  "\\begin{footnotesize}",
  "\\setlength{\\LTcapwidth}{\\textwidth}",
  paste0("\\begin{longtable}{p{7.5cm}lrrrrrrr}"),
  paste0("\\caption[Outlier pathways]{\\textbf{Outlier pathways.} ",
         "All non-Equilibrium pathways with regime assignment, key metrics, ",
         "and independence cluster ID. ",
         "MR: mean residual; D$b$: deviation $b$; ",
         "$\\bar{s}$: mean $s_\\mathrm{het}$; ",
         "pLI: mean pLI; CF: complex fraction; ",
         "$n$: number of genes; Cl.: cluster ID (pathways sharing the same cluster ",
         "have Jaccard similarity $\\geq 0.50$).}\\label{tab:outlier_catalogue} \\\\"),
  "\\hline",
  "Pathway & Regime & $n$ & MR & D$b$ & $\\bar{s}$ & pLI & CF & Cl. \\\\",
  "\\hline",
  "\\endfirsthead",
  "\\hline",
  "Pathway & Regime & $n$ & MR & D$b$ & $\\bar{s}$ & pLI & CF & Cl. \\\\",
  "\\hline",
  "\\endhead",
  "\\hline",
  "\\endfoot"
)

# Short regime labels for compactness
regime_short <- c(
  "Buffered" = "Buf.",
  "Decoupled" = "Dec.",
  "Topology-Amplified" = "T-Amp."
)

current_regime <- ""
for (i in seq_len(nrow(outliers))) {
  r <- outliers[i]
  reg <- as.character(r$regime)

  # Add a midrule between regimes
  if (reg != current_regime) {
    if (current_regime != "") tex_s2 <- c(tex_s2, "\\hline")
    current_regime <- reg
  }

  # Full pathway name (no truncation; column width handles wrapping)
  pw_name <- as.character(r$pathway)

  # Extract cluster number
  cl_num <- sub(".*_", "", as.character(r$cluster_id))

  row <- sprintf(
    "%s & %s & %d & %s & %s & %s & %s & %s & %s \\\\",
    escape_latex(pw_name),
    regime_short[reg],
    r$n_genes,
    fmt(r$mean_residual, 3),
    fmt(r$deviation_kappa, 3),
    fmt(r$mean_shet, 3),
    fmt(r$mean_pli, 2),
    fmt(r$complex_fraction, 2),
    cl_num
  )
  tex_s2 <- c(tex_s2, row)
}

tex_s2 <- c(tex_s2, "\\hline", "\\end{longtable}", "\\end{footnotesize}")

writeLines(tex_s2, file.path(MANU_DIR, "supp_table_outlier_pathways.tex"))
cat("  Saved supp_table_outlier_pathways.tex\n")

# ==========================================================================
# TABLE 1 (main text): OLS regression of mean s_het
# ==========================================================================
cat("=== Table 1: OLS regression (main text) ===\n")

ols <- fread(file.path(OUTPUT, "table1_ols_shet.csv"))
mstats <- fread(file.path(OUTPUT, "table1_ols_shet_summary.csv"))

# Pretty predictor names
pred_labels <- c(
  density           = "Network density",
  mean_paralog      = "Mean paralog count",
  mean_cds_length   = "Mean CDS length",
  clustering_score  = "Genomic clustering",
  mean_pubmed       = "Mean PubMed count"
)

# Format p-value as LaTeX scientific notation or plain decimal
fmt_p <- function(p) {
  if (is.na(p)) return("---")
  if (p >= 0.01) return(sprintf("%.2f", p))
  exp <- floor(log10(p))
  mant <- p / 10^exp
  sprintf("$%.1f \\times 10^{%d}$", mant, exp)
}

# Format beta/SE
fmt_coef <- function(x) {
  if (is.na(x)) return("---")
  if (x < 0) return(sprintf("$-$%s", formatC(abs(x), format = "f", digits = 4)))
  sprintf("%s", formatC(x, format = "f", digits = 4))
}

n_primary <- mstats[model == "primary"]$n

tex_t1 <- c(
  "\\begin{table*}[!t]",
  "\\centering",
  sprintf(paste0(
    "\\caption{\\textbf{Multivariate linear regression of ",
    "pathway-level mean $s_\\mathrm{het}$ on network density ",
    "and confounders} ($n = %s$ Reactome pathways). All ",
    "predictors are standardized to zero mean and unit variance. ",
    "The primary model controls for paralogy, CDS length, and ",
    "genomic clustering. The secondary model additionally ",
    "includes mean PubMed publication count per gene. $\\beta$: ",
    "standardized regression coefficient; SE: standard error.}%%"),
    format(n_primary, big.mark = "{,}")),
  "\\label{tab:linreg}",
  "\\begin{tabular}{lrrrrrrr}",
  "\\hline",
  " & \\multicolumn{3}{c}{Primary model}",
  " & \\multicolumn{3}{c}{Sensitivity model} \\\\",
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
  "Predictor & $\\beta$ & SE & $p$",
  "  & $\\beta$ & SE & $p$ \\\\"
)

tex_t1 <- c(tex_t1, "\\hline")

for (i in seq_len(nrow(ols))) {
  r <- ols[i]
  pname <- pred_labels[r$predictor]
  row <- sprintf(
    "%s\n  & %s & %s & %s\n  & %s & %s & %s \\\\",
    pname,
    fmt_coef(r$beta_primary), fmt_coef(r$se_primary), fmt_p(r$p_primary),
    fmt_coef(r$beta_sensitivity), fmt_coef(r$se_sensitivity), fmt_p(r$p_sensitivity)
  )
  tex_t1 <- c(tex_t1, row)
}

# R² row
mp <- mstats[model == "primary"]
ms <- mstats[model == "sensitivity"]
tex_t1 <- c(tex_t1,
  "\\hline",
  sprintf("$R^2$ / Adj.\\ $R^2$\n  & \\multicolumn{3}{c}{%.3f / %.3f}\n  & \\multicolumn{3}{c}{%.3f / %.3f} \\\\",
          mp$r_squared, mp$adj_r_squared, ms$r_squared, ms$adj_r_squared),
  "\\hline",
  "\\end{tabular}",
  "\\end{table*}"
)

writeLines(tex_t1, file.path(MANU_DIR, "table1_ols_regression.tex"))
cat("  Saved table1_ols_regression.tex\n")


cat("\nDone.\n")
