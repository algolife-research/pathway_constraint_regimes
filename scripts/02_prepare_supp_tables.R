suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(readxl)
  library(future.apply)
})

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
dir.create(OUTPUT, showWarnings = FALSE)

REGIME_ORDER     <- c("Buffered", "Equilibrium", "Decoupled", "Topology-Amplified")
MIN_GENES_KAPPA  <- 15
OUTLIER_PROB     <- 0.90
N_PERMUTATIONS   <- 500
N_BOOTSTRAP      <- 2000
STRING_EXP_THR   <- 300
JACCARD_THR      <- 0.50
RNG_SEED         <- 42

set.seed(RNG_SEED)

# ==========================================================================
# SECTION 0: Shared data loading
# ==========================================================================
cat("=== Loading shared data ===\n")

gene_topo     <- fread(file.path(OUTPUT, "gene_topology_metrics.csv"))
pw_kappa      <- fread(file.path(OUTPUT, "pathway_kappa_model.csv"))
pw_classified <- fread(file.path(OUTPUT, "pathway_regimes_classified.csv"))

cat(sprintf("  %d genes, %d pathways with TCM, %d pathways classified\n",
            nrow(gene_topo), nrow(pw_kappa), nrow(pw_classified)))

gmt_path <- file.path(DATA, "ReactomePathways.gmt.zip")
gmt_file <- unzip(gmt_path, list = TRUE)$Name[1]
gmt_con <- unz(gmt_path, gmt_file)
gmt_lines <- readLines(gmt_con)
close(gmt_con)

pw_genes_list <- setNames(
  lapply(gmt_lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    list(name = parts[1], id = parts[2], genes = parts[-(1:2)])
  }),
  sapply(gmt_lines, function(x) strsplit(x, "\t")[[1]][1])
)
pw_genes_list <- Filter(function(x) length(x$genes) >= 10 & length(x$genes) <= 300,
                         pw_genes_list)
cat(sprintf("  %d Reactome pathways (10-300 genes)\n", length(pw_genes_list)))

hierarchy <- fread(file.path(DATA, "ReactomePathwaysRelation.txt"),
                   header = FALSE, col.names = c("parent", "child"))
pw_names <- fread(file.path(DATA, "ReactomePathways.txt"),
                  header = FALSE, col.names = c("id", "name", "species"))
pw_names <- pw_names[species == "Homo sapiens"]

get_top_ancestor <- function(rid) {
  visited <- character(0)
  current <- rid
  while (TRUE) {
    parents <- hierarchy[child == current]$parent
    parents <- parents[grepl("^R-HSA-", parents)]
    if (length(parents) == 0 || current %in% visited) break
    visited <- c(visited, current)
    current <- parents[1]
  }
  current
}

SUPER_MAP <- c(
  "Immune System" = "Immune System",
  "Metabolism" = "Metabolism",
  "Signal Transduction" = "Signal Transduction",
  "Gene Expression" = "Gene Expression",
  "Cell Cycle" = "Cell Cycle & DNA",
  "DNA Repair" = "Cell Cycle & DNA",
  "DNA Replication" = "Cell Cycle & DNA",
  "Cellular responses to stimuli" = "Signal Transduction",
  "Transport of small molecules" = "Transport & Vesicles",
  "Vesicle-mediated transport" = "Transport & Vesicles",
  "Membrane Trafficking" = "Transport & Vesicles"
)

# 0d. Auxiliary data
shet <- as.data.table(read_excel(
  file.path(DATA, "Cassa_2017_41588_2017_BFng3831_MOESM71_ESM.xlsx")))
shet <- shet[!is.na(s_het), .(gene = gene_symbol, s_het)]

gnomad <- fread(file.path(DATA, "gnomad_constraint.txt.bgz"),
                select = c("gene", "transcript", "oe_lof_upper", "pLI", "mis_z"))
gnomad <- gnomad[!is.na(oe_lof_upper)][order(oe_lof_upper)][, .SD[1], by = gene]
setnames(gnomad, "oe_lof_upper", "loeuf")

paralogs <- fread(file.path(DATA, "paralog_counts.tsv"))
cds <- fread(file.path(DATA, "ensembl_cds_length.tsv"))
if ("cds_length" %in% names(cds)) {
  cds_max <- cds[, .(cds_length = max(cds_length, na.rm = TRUE)), by = gene]
} else {
  cds_alt <- fread(file.path(DATA, "gencode_cds_lengths.tsv"))
  cds_max <- cds_alt[, .(cds_length = max(cds_length, na.rm = TRUE)), by = gene]
}
coords <- fread(file.path(DATA, "gene_coordinates.tsv"))

# 0e. Build BioGRID graph
cat("  Building BioGRID PPI graph...\n")
biogrid_zip <- file.path(DATA, "biogrid_latest.tab3.zip")
biogrid_file <- unzip(biogrid_zip, list = TRUE)$Name[1]
bg <- fread(cmd = sprintf('unzip -p "%s" "%s"', biogrid_zip, biogrid_file),
            sep = "\t", header = TRUE, select = c(8, 9, 13, 16, 17),
            col.names = c("geneA", "geneB", "system_type", "taxA", "taxB"))
bg <- bg[system_type == "physical" & taxA == 9606 & taxB == 9606 & geneA != geneB]
bg <- unique(bg[, .(geneA, geneB)])

g <- graph_from_data_frame(bg, directed = FALSE)
g <- simplify(g)
comps <- components(g)
giant_nodes <- names(comps$membership[comps$membership == which.max(comps$csize)])
g <- induced_subgraph(g, giant_nodes)
net_genes <- V(g)$name
cat(sprintf("  Giant component: %d nodes, %d edges\n", vcount(g), ecount(g)))

# 0f. Reconstruct gene_for_tcm
gene_degree <- data.table(gene = V(g)$name, degree = degree(g))
gene_for_tcm <- merge(gene_degree, shet, by = "gene")
gene_for_tcm <- gene_for_tcm[degree > 0]
gene_for_tcm[, log_degree := log10(degree)]
baseline <- lm(s_het ~ log_degree, data = gene_for_tcm)
gene_for_tcm[, residual := residuals(baseline)]
cat(sprintf("  %d genes with degree + s_het for TCM\n", nrow(gene_for_tcm)))

# 0g. Rectangular thresholds (from original data)
kappa_complete <- pw_kappa[!is.na(mean_residual) & !is.na(deviation_kappa)]
resid_upper <- quantile(kappa_complete$mean_residual, OUTLIER_PROB)
resid_lower <- quantile(kappa_complete$mean_residual, 1 - OUTLIER_PROB)
kappa_upper <- quantile(kappa_complete$deviation_kappa, OUTLIER_PROB)
kappa_lower <- quantile(kappa_complete$deviation_kappa, 1 - OUTLIER_PROB)
cat(sprintf("  Thresholds: resid [%.4f, %.4f], kappa [%.4f, %.4f]\n",
            resid_lower, resid_upper, kappa_lower, kappa_upper))

# 0h. Helper functions
classify_regimes_rect <- function(mr, dk, ru = resid_upper, rl = resid_lower,
                                   ku = kappa_upper, kl = kappa_lower) {
  fifelse(
    mr > ru & dk > ku, "Topology-Amplified",
    fifelse(mr > ru & dk < kl, "Decoupled",
    fifelse(mr < rl & dk < kl, "Buffered",
    "Equilibrium")))
}

compute_pathway_tcm <- function(genes, gene_data, min_genes = MIN_GENES_KAPPA) {
  gd <- gene_data[gene %in% genes]
  gd <- gd[!is.na(residual) & !is.na(log_degree)]
  if (nrow(gd) < min_genes) return(list(mr = NA_real_, dk = NA_real_, n = nrow(gd)))
  mr <- mean(gd$residual)
  if (nrow(gd) >= 3 && var(gd$log_degree) > 0) {
    X <- cbind(1, gd$log_degree)
    fit <- .lm.fit(X, gd$residual)
    dk <- fit$coefficients[2]
  } else {
    dk <- NA_real_
  }
  list(mr = unname(mr), dk = unname(dk), n = nrow(gd))
}


# ==========================================================================
# SECTION 1: pathway_regime_metrics.csv
# ==========================================================================
cat("\n=== [1/9] Building pathway_regime_metrics.csv ===\n")

pw_regime_metrics <- rbindlist(lapply(pw_genes_list, function(pw) {
  genes <- pw$genes
  genes_in_net <- intersect(genes, net_genes)
  n <- length(genes)
  n_net <- length(genes_in_net)

  if (n_net >= 2) {
    sg <- induced_subgraph(g, genes_in_net)
    e_int <- ecount(sg)
    density <- 2 * e_int / (n_net * (n_net - 1))
    total_deg <- sum(degree(g, genes_in_net))
    insularity <- if (total_deg > 0) 2 * e_int / total_deg else 0
    n_edges <- e_int
  } else {
    density <- 0
    insularity <- 0
    n_edges <- 0L
  }

  loeuf_vals <- gnomad[gene %in% genes]$loeuf
  mean_loeuf <- if (length(loeuf_vals) > 0) mean(loeuf_vals, na.rm = TRUE) else NA_real_
  loeuf_var  <- if (length(loeuf_vals) > 1) var(loeuf_vals, na.rm = TRUE) else NA_real_

  top_id <- get_top_ancestor(pw$id)
  top_name <- pw_names[id == top_id]$name
  if (length(top_name) == 0) top_name <- "Other"
  root_cat <- top_name[1]
  super <- if (root_cat %in% names(SUPER_MAP)) unname(SUPER_MAP[root_cat]) else "Other"

  data.table(
    pathway = pw$name,
    reactome_id = pw$id,
    n_genes = n,
    n_in_network = n_net,
    n_edges = n_edges,
    density = density,
    mean_loeuf = mean_loeuf,
    loeuf_variance = loeuf_var,
    insularity = insularity,
    root_category = root_cat,
    super_pathway = super
  )
}))

# Log-linear regression for predicted LOEUF
pw_regime_metrics[, log_density := log10(pmax(density, 1e-6))]
loeuf_fit <- lm(mean_loeuf ~ log_density, data = pw_regime_metrics[!is.na(mean_loeuf)])
pw_regime_metrics[, predicted_loeuf := predict(loeuf_fit, newdata = .SD)]
pw_regime_metrics[, residual := mean_loeuf - predicted_loeuf]

# Merge regime labels
pw_regime_metrics <- merge(pw_regime_metrics,
                            kappa_complete[, .(pathway, regime)],
                            by = "pathway", all.x = TRUE)
pw_regime_metrics[is.na(regime), regime := "Equilibrium"]

fwrite(pw_regime_metrics, file.path(OUTPUT, "pathway_regime_metrics.csv"))
cat(sprintf("  Saved pathway_regime_metrics.csv (%d rows)\n", nrow(pw_regime_metrics)))


# ==========================================================================
# SECTION 2: pathway_dfe_validation.csv
# ==========================================================================
cat("\n=== [2/9] Building pathway_dfe_validation.csv ===\n")

pw_dfe <- rbindlist(lapply(pw_genes_list, function(pw) {
  genes <- pw$genes

  shet_vals <- shet[gene %in% genes]$s_het
  mean_shet <- if (length(shet_vals) > 0) mean(shet_vals, na.rm = TRUE) else NA_real_
  n_shet <- sum(!is.na(shet_vals))

  misz_vals <- gnomad[gene %in% genes]$mis_z
  mean_misz <- if (length(misz_vals) > 0) mean(misz_vals, na.rm = TRUE) else NA_real_
  n_misz <- sum(!is.na(misz_vals))

  data.table(
    pathway = pw$name,
    mean_shet = mean_shet,
    n_shet = n_shet,
    mean_missense_z = mean_misz,
    n_misz = n_misz
  )
}))

# Merge with pw_classified columns
pw_dfe <- merge(pw_dfe,
                pw_classified[, .(pathway, regime, density, mean_loeuf,
                                  mean_paralog, mean_cds_length, clustering_score)],
                by = "pathway", all.x = TRUE)
pw_dfe[is.na(regime), regime := "Equilibrium"]

fwrite(pw_dfe, file.path(OUTPUT, "pathway_dfe_validation.csv"))
cat(sprintf("  Saved pathway_dfe_validation.csv (%d rows)\n", nrow(pw_dfe)))


# ==========================================================================
# SECTION 3: supp_mahalanobis_sensitivity.csv (rectangular threshold sweep)
# ==========================================================================
cat("\n=== [3/9] Building supp_mahalanobis_sensitivity.csv ===\n")

percentiles <- seq(0.75, 0.95, by = 0.01)

sens_results <- rbindlist(lapply(percentiles, function(p) {
  ru <- quantile(kappa_complete$mean_residual, p)
  rl <- quantile(kappa_complete$mean_residual, 1 - p)
  ku <- quantile(kappa_complete$deviation_kappa, p)
  kl <- quantile(kappa_complete$deviation_kappa, 1 - p)

  regimes <- classify_regimes_rect(kappa_complete$mean_residual,
                                    kappa_complete$deviation_kappa,
                                    ru, rl, ku, kl)
  counts <- table(factor(regimes, levels = REGIME_ORDER))

  data.table(
    percentile = p,
    Buffered = as.integer(counts["Buffered"]),
    Equilibrium = as.integer(counts["Equilibrium"]),
    Decoupled = as.integer(counts["Decoupled"]),
    `Topology-Amplified` = as.integer(counts["Topology-Amplified"]),
    total_outlier = as.integer(sum(counts) - counts["Equilibrium"])
  )
}))

fwrite(sens_results, file.path(OUTPUT, "supp_mahalanobis_sensitivity.csv"))
cat(sprintf("  Saved supp_mahalanobis_sensitivity.csv (%d rows)\n", nrow(sens_results)))


# ==========================================================================
# SECTION 4: supp_cohens_d_degree_bins.csv
# ==========================================================================
cat("\n=== [4/9] Building supp_cohens_d_degree_bins.csv ===\n")

# Assign each gene to regime of its most constrained pathway (lowest mean_loeuf)
gene_pw_regime <- rbindlist(lapply(pw_genes_list, function(pw) {
  data.table(gene = pw$genes, pathway = pw$name)
}))
gene_pw_regime <- merge(gene_pw_regime,
                         pw_classified[, .(pathway, regime, mean_loeuf)],
                         by = "pathway")
gene_pw_regime <- gene_pw_regime[!is.na(regime)]
gene_regime <- gene_pw_regime[order(mean_loeuf)][, .SD[1], by = gene][, .(gene, regime)]

# Merge with degree and s_het
gene_cd <- merge(gene_degree[degree > 0], shet, by = "gene")
gene_cd <- merge(gene_cd, gene_regime, by = "gene")
gene_cd[, log_degree := log10(degree)]

# 5 equal-width bins on log10(degree)
breaks <- seq(min(gene_cd$log_degree), max(gene_cd$log_degree), length.out = 6)
gene_cd[, bin := cut(log_degree, breaks = breaks, include.lowest = TRUE)]
gene_cd[, bin_idx := as.integer(bin)]

# Extract bin boundaries
bin_levels <- levels(gene_cd$bin)

cohens_d_fn <- function(x, y) {
  n1 <- length(x); n2 <- length(y)
  if (n1 < 2 || n2 < 2) return(NA_real_)
  s <- sqrt(((n1 - 1) * var(x) + (n2 - 1) * var(y)) / (n1 + n2 - 2))
  if (s == 0) return(NA_real_)
  (mean(x) - mean(y)) / s
}

cd_results <- rbindlist(lapply(seq_along(bin_levels), function(i) {
  bl <- bin_levels[i]
  sub <- gene_cd[bin == bl]
  dec_vals <- sub[regime == "Decoupled"]$s_het
  eq_vals  <- sub[regime == "Equilibrium"]$s_het

  n_dec <- length(dec_vals)
  n_eq  <- length(eq_vals)

  if (n_dec < 10 || n_eq < 10) return(NULL)

  d_obs <- cohens_d_fn(dec_vals, eq_vals)

  # Bootstrap CI
  set.seed(RNG_SEED + i)
  d_boot <- replicate(N_BOOTSTRAP, {
    x <- sample(dec_vals, n_dec, replace = TRUE)
    y <- sample(eq_vals, n_eq, replace = TRUE)
    cohens_d_fn(x, y)
  })
  ci <- quantile(d_boot, c(0.025, 0.975), na.rm = TRUE)

  # Mann-Whitney
  wt <- wilcox.test(dec_vals, eq_vals, alternative = "greater")

  # Parse bin boundaries
  nums <- as.numeric(regmatches(bl, gregexpr("[0-9eE.+-]+", bl))[[1]])
  bin_left <- nums[1]
  bin_right <- nums[2]
  deg_label <- sprintf("%d-%d", round(10^bin_left), round(10^bin_right))

  data.table(
    degree_bin = deg_label,
    bin_left = bin_left,
    bin_right = bin_right,
    n_decoupled = n_dec,
    n_equilibrium = n_eq,
    median_shet_dec = median(dec_vals),
    median_shet_eq = median(eq_vals),
    cohens_d = d_obs,
    d_ci_lo = ci[1],
    d_ci_hi = ci[2],
    mannwhitney_p = wt$p.value
  )
}))

fwrite(cd_results, file.path(OUTPUT, "supp_cohens_d_degree_bins.csv"))
cat(sprintf("  Saved supp_cohens_d_degree_bins.csv (%d rows)\n", nrow(cd_results)))


# ==========================================================================
# SECTION 5: supp_pathway_independence.csv
# ==========================================================================
cat("\n=== [5/9] Building supp_pathway_independence.csv ===\n")

# Build gene set lookup from GMT
pw_gene_sets <- lapply(pw_genes_list, function(pw) pw$genes)
names(pw_gene_sets) <- sapply(pw_genes_list, function(pw) pw$name)
pw_id_map <- setNames(
  sapply(pw_genes_list, function(pw) pw$id),
  sapply(pw_genes_list, function(pw) pw$name)
)

# Get regime for each pathway
pw_regime_dt <- kappa_complete[, .(pathway, regime)]

indep_results <- rbindlist(lapply(REGIME_ORDER, function(reg) {
  pws <- pw_regime_dt[regime == reg]$pathway
  pws <- pws[pws %in% names(pw_gene_sets)]
  n_pw <- length(pws)

  if (n_pw <= 1) {
    return(data.table(
      pathway = pws,
      reactome_id = pw_id_map[pws],
      regime = reg,
      n_genes = sapply(pws, function(p) length(pw_gene_sets[[p]])),
      cluster_id = paste0(reg, "_1"),
      cluster_size = 1L,
      is_cluster_representative = TRUE
    ))
  }

  # Pairwise Jaccard distance matrix
  dist_mat <- matrix(0, n_pw, n_pw)
  for (i in seq_len(n_pw - 1)) {
    for (j in (i + 1):n_pw) {
      a <- pw_gene_sets[[pws[i]]]
      b <- pw_gene_sets[[pws[j]]]
      jaccard <- length(intersect(a, b)) / length(union(a, b))
      dist_mat[i, j] <- 1 - jaccard
      dist_mat[j, i] <- 1 - jaccard
    }
  }
  rownames(dist_mat) <- pws
  colnames(dist_mat) <- pws

  # Single-linkage clustering
  hc <- hclust(as.dist(dist_mat), method = "single")
  clusters <- cutree(hc, h = 1 - JACCARD_THR)

  # Build result
  pw_dt <- data.table(
    pathway = pws,
    reactome_id = pw_id_map[pws],
    regime = reg,
    n_genes = sapply(pws, function(p) length(pw_gene_sets[[p]])),
    cluster_raw = clusters
  )
  pw_dt[, cluster_id := paste0(reg, "_", cluster_raw)]
  pw_dt[, cluster_size := .N, by = cluster_id]
  pw_dt[, is_cluster_representative := n_genes == max(n_genes), by = cluster_id]
  # Break ties: keep first
  pw_dt[, is_cluster_representative := seq_len(.N) == which.max(n_genes), by = cluster_id]
  pw_dt[, cluster_raw := NULL]
  pw_dt
}))

fwrite(indep_results, file.path(OUTPUT, "supp_pathway_independence.csv"))
cat(sprintf("  Saved supp_pathway_independence.csv (%d rows)\n", nrow(indep_results)))

# Print summary
indep_summary <- indep_results[, .(raw = .N, effective = uniqueN(cluster_id)), by = regime]
cat("  Regime independence summary:\n")
print(indep_summary)


# ==========================================================================
# SECTION 6: supp_jackknife_stability.csv
# ==========================================================================
cat("\n=== [6/9] Building supp_jackknife_stability.csv ===\n")

# For each pathway, leave-one-out and reclassify
jk_pathways <- kappa_complete[!is.na(mean_residual) & !is.na(deviation_kappa) &
                                n_genes >= MIN_GENES_KAPPA]

jk_results <- rbindlist(lapply(seq_len(nrow(jk_pathways)), function(idx) {
  pw_name <- jk_pathways$pathway[idx]
  orig_regime <- jk_pathways$regime[idx]
  orig_mr <- jk_pathways$mean_residual[idx]
  orig_dk <- jk_pathways$deviation_kappa[idx]

  # Get gene-level data for this pathway
  pw_info <- pw_genes_list[[pw_name]]
  if (is.null(pw_info)) return(NULL)

  gd <- gene_for_tcm[gene %in% pw_info$genes]
  gd <- gd[!is.na(residual) & !is.na(log_degree)]
  n <- nrow(gd)
  if (n < MIN_GENES_KAPPA) return(NULL)

  # Leave-one-out
  jk_mr <- numeric(n)
  jk_dk <- numeric(n)
  jk_regime <- character(n)

  for (i in seq_len(n)) {
    gd_loo <- gd[-i]
    jk_mr[i] <- mean(gd_loo$residual)
    if (nrow(gd_loo) >= 3 && var(gd_loo$log_degree) > 0) {
      X <- cbind(1, gd_loo$log_degree)
      fit <- .lm.fit(X, gd_loo$residual)
      jk_dk[i] <- fit$coefficients[2]
    } else {
      jk_dk[i] <- NA_real_
    }
    jk_regime[i] <- classify_regimes_rect(jk_mr[i], jk_dk[i])
  }

  stability <- mean(jk_regime == orig_regime, na.rm = TRUE)
  influential_idx <- which(jk_regime != orig_regime)
  n_influential <- length(influential_idx)
  influential_genes <- if (n_influential > 0) {
    paste(head(gd$gene[influential_idx], 10), collapse = ";")
  } else {
    ""
  }

  data.table(
    pathway = pw_name,
    original_regime = orig_regime,
    n_genes = n,
    mean_residual = orig_mr,
    deviation_kappa = orig_dk,
    stability = stability,
    n_influential = n_influential,
    influential_genes = influential_genes,
    mr_range = max(jk_mr) - min(jk_mr),
    dk_range = max(jk_dk, na.rm = TRUE) - min(jk_dk, na.rm = TRUE)
  )
}))

fwrite(jk_results, file.path(OUTPUT, "supp_jackknife_stability.csv"))
cat(sprintf("  Saved supp_jackknife_stability.csv (%d rows)\n", nrow(jk_results)))


# ==========================================================================
# SECTION 7: STRING replication
# ==========================================================================
cat("\n=== [7/9] Building supp_string_replication.csv + supp_string_concordance.csv ===\n")

# 7a. Load STRING network (experimental channel only, physical interactions)
cat("  Loading STRING v12.0 detailed links...\n")
string_info <- fread(file.path(DATA, "9606.protein.info.v12.0.txt.gz"))
id_to_gene <- setNames(string_info$preferred_name, string_info[[1]])

string_links <- fread(file.path(DATA, "9606.protein.links.detailed.v12.0.txt.gz"))
string_links <- string_links[experimental >= STRING_EXP_THR]
string_links[, gene1 := id_to_gene[protein1]]
string_links[, gene2 := id_to_gene[protein2]]
string_links <- string_links[!is.na(gene1) & !is.na(gene2) & gene1 != gene2]
cat(sprintf("  STRING: %d interactions (experimental >= %d)\n", nrow(string_links), STRING_EXP_THR))

# 7b. STRING degree per gene
string_deg <- rbind(
  string_links[, .(gene = gene1)],
  string_links[, .(gene = gene2)]
)[, .(degree = .N), by = gene]

# 7c. Genome-wide baseline on STRING
string_gene <- merge(string_deg[degree > 0], shet, by = "gene")
string_gene[, log_degree := log10(degree)]
string_baseline <- lm(s_het ~ log_degree, data = string_gene)
string_gene[, residual := residuals(string_baseline)]
cat(sprintf("  STRING baseline: %d genes, R² = %.4f\n",
            nrow(string_gene), summary(string_baseline)$r.squared))

# 7d. Per-pathway TCM on STRING
string_pw <- rbindlist(lapply(pw_genes_list, function(pw) {
  tcm <- compute_pathway_tcm(pw$genes, string_gene)
  if (is.na(tcm$mr)) return(NULL)
  data.table(pathway = pw$name, mean_residual = tcm$mr,
             deviation_kappa = tcm$dk, n_genes = tcm$n)
}))

# 7e. Rectangular thresholds on STRING data
string_complete <- string_pw[!is.na(mean_residual) & !is.na(deviation_kappa)]
s_ru <- quantile(string_complete$mean_residual, OUTLIER_PROB)
s_rl <- quantile(string_complete$mean_residual, 1 - OUTLIER_PROB)
s_ku <- quantile(string_complete$deviation_kappa, OUTLIER_PROB)
s_kl <- quantile(string_complete$deviation_kappa, 1 - OUTLIER_PROB)
string_complete[, regime := classify_regimes_rect(mean_residual, deviation_kappa,
                                                   s_ru, s_rl, s_ku, s_kl)]

fwrite(string_complete, file.path(OUTPUT, "supp_string_replication.csv"))
cat(sprintf("  Saved supp_string_replication.csv (%d rows)\n", nrow(string_complete)))

# 7f. Concordance
conc <- merge(
  string_complete[, .(pathway, mean_residual_string = mean_residual,
                       deviation_kappa_string = deviation_kappa,
                       regime_string = regime)],
  kappa_complete[, .(pathway, mean_residual_biogrid = mean_residual,
                      deviation_kappa_biogrid = deviation_kappa,
                      regime_biogrid = regime)],
  by = "pathway"
)

agreement <- mean(conc$regime_string == conc$regime_biogrid)
r_mr <- cor(conc$mean_residual_string, conc$mean_residual_biogrid)
r_dk <- cor(conc$deviation_kappa_string, conc$deviation_kappa_biogrid)
cat(sprintf("  Concordance: %.1f%% agreement, r(mr)=%.2f, r(dk)=%.2f\n",
            agreement * 100, r_mr, r_dk))

fwrite(conc, file.path(OUTPUT, "supp_string_concordance.csv"))
cat(sprintf("  Saved supp_string_concordance.csv (%d rows)\n", nrow(conc)))


# ==========================================================================
# SECTION 8: supp_permutation_null.csv (500 permutations, parallelised)
# ==========================================================================
cat("\n=== [8/9] Building supp_permutation_null.csv (500 permutations) ===\n")

# Pathway sizes (only pathways with >= MIN_GENES_KAPPA genes in gene_for_tcm)
pw_sizes <- sapply(pw_genes_list, function(pw) {
  sum(pw$genes %in% gene_for_tcm$gene)
})
pw_sizes <- pw_sizes[pw_sizes >= MIN_GENES_KAPPA]
cat(sprintf("  %d pathways with >= %d genes, total pool = %d genes\n",
            length(pw_sizes), MIN_GENES_KAPPA, nrow(gene_for_tcm)))

gene_pool <- gene_for_tcm$gene

# Set up parallel
n_workers <- min(parallel::detectCores() - 1, 8)
n_workers <- max(n_workers, 1)
plan(multisession, workers = n_workers)
cat(sprintf("  Running %d permutations on %d workers...\n", N_PERMUTATIONS, n_workers))

perm_results <- future_lapply(seq_len(N_PERMUTATIONS), function(perm_i) {
  set.seed(RNG_SEED + perm_i)

  # Shuffle gene pool
  shuffled <- sample(gene_pool)

  # Assign to pathway-sized groups
  sizes <- unname(pw_sizes)
  cum_sizes <- cumsum(sizes)
  starts <- c(1L, cum_sizes[-length(cum_sizes)] + 1L)

  # Compute TCM for each permuted group
  mrs <- numeric(length(sizes))
  dks <- numeric(length(sizes))
  n_fitted <- 0L

  for (j in seq_along(sizes)) {
    idx <- starts[j]:cum_sizes[j]
    genes_j <- shuffled[idx]
    gd <- gene_for_tcm[gene %in% genes_j]
    gd <- gd[!is.na(residual) & !is.na(log_degree)]

    if (nrow(gd) < MIN_GENES_KAPPA) {
      mrs[j] <- NA_real_
      dks[j] <- NA_real_
      next
    }

    n_fitted <- n_fitted + 1L
    mrs[j] <- mean(gd$residual)
    if (nrow(gd) >= 3 && var(gd$log_degree) > 0) {
      X <- cbind(1, gd$log_degree)
      fit <- .lm.fit(X, gd$residual)
      dks[j] <- fit$coefficients[2]
    } else {
      dks[j] <- NA_real_
    }
  }

  # Remove NAs
  valid <- !is.na(mrs) & !is.na(dks)
  mrs_v <- mrs[valid]
  dks_v <- dks[valid]

  if (length(mrs_v) < 10) {
    return(data.table(
      permutation = perm_i, n_fitted = n_fitted,
      n_outlier = 0L, n_TA = 0L, n_Dec = 0L, n_Buf = 0L,
      spread_mr = NA_real_, spread_dk = NA_real_
    ))
  }

  # Rectangular thresholds on PERMUTED data
  p_ru <- quantile(mrs_v, OUTLIER_PROB)
  p_rl <- quantile(mrs_v, 1 - OUTLIER_PROB)
  p_ku <- quantile(dks_v, OUTLIER_PROB)
  p_kl <- quantile(dks_v, 1 - OUTLIER_PROB)

  regimes <- fifelse(
    mrs_v > p_ru & dks_v > p_ku, "Topology-Amplified",
    fifelse(mrs_v > p_ru & dks_v < p_kl, "Decoupled",
    fifelse(mrs_v < p_rl & dks_v < p_kl, "Buffered",
    "Equilibrium")))

  counts <- table(factor(regimes, levels = REGIME_ORDER))

  data.table(
    permutation = perm_i,
    n_fitted = n_fitted,
    n_outlier = as.integer(sum(counts) - counts["Equilibrium"]),
    n_TA = as.integer(counts["Topology-Amplified"]),
    n_Dec = as.integer(counts["Decoupled"]),
    n_Buf = as.integer(counts["Buffered"]),
    spread_mr = sd(mrs_v),
    spread_dk = sd(dks_v)
  )
}, future.seed = TRUE)

plan(sequential)

perm_dt <- rbindlist(perm_results)
fwrite(perm_dt, file.path(OUTPUT, "supp_permutation_null.csv"))
cat(sprintf("  Saved supp_permutation_null.csv (%d rows)\n", nrow(perm_dt)))
