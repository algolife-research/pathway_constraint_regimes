suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(readxl)
  library(fgsea)
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

REGIME_ORDER <- c("Buffered", "Equilibrium", "Decoupled", "Topology-Amplified")


# ==========================================================================
# 1. Gene-level topology metrics
# ==========================================================================
cat("=== [1/7] Building PPI network & computing gene topology ===\n")

# Load BioGRID physical interactions (human only, no self-loops)
biogrid_zip <- file.path(DATA, "biogrid_latest.tab3.zip")
biogrid_file <- unzip(biogrid_zip, list = TRUE)$Name[1]
bg <- fread(cmd = sprintf('unzip -p "%s" "%s"', biogrid_zip, biogrid_file),
            sep = "\t", header = TRUE, select = c(8, 9, 13, 16, 17),
            col.names = c("geneA", "geneB", "system_type", "taxA", "taxB"))

bg <- bg[system_type == "physical" & taxA == 9606 & taxB == 9606 & geneA != geneB]
bg <- unique(bg[, .(geneA, geneB)])
cat(sprintf("  BioGRID: %d physical interactions\n", nrow(bg)))

# Build undirected graph
g <- graph_from_data_frame(bg, directed = FALSE)
g <- simplify(g)  # remove multi-edges
comps <- components(g)
giant_nodes <- names(comps$membership[comps$membership == which.max(comps$csize)])
g <- induced_subgraph(g, giant_nodes)
cat(sprintf("  Giant component: %d nodes, %d edges\n", vcount(g), ecount(g)))

# Compute topology metrics
cat("  Computing centrality metrics...\n")
n <- vcount(g)
deg_cent <- degree(g) / (n - 1)
clust_coeff <- transitivity(g, type = "local", isolates = "zero")
kcore_idx <- coreness(g)

# Betweenness with k=500 approximation: use cutoff
cat("  Computing betweenness...\n")
btwn <- estimate_betweenness(g, cutoff = 3, directed = FALSE)
btwn_norm <- btwn / ((n - 1) * (n - 2) / 2)

# Load gnomAD LOEUF
gnomad <- fread(file.path(DATA, "gnomad_constraint.txt.bgz"),
                select = c("gene", "transcript", "oe_lof_upper", "pLI", "mis_z"))
# Keep canonical transcript (lowest LOEUF per gene)
gnomad <- gnomad[!is.na(oe_lof_upper)][order(oe_lof_upper)][
  , .SD[1], by = gene]
setnames(gnomad, "oe_lof_upper", "loeuf")

# PubMed counts
pubmed <- fread(file.path(DATA, "pubmed_counts_per_gene.tsv"))

# Assemble gene-level table
gene_topo <- data.table(
  gene = V(g)$name,
  degree_centrality = deg_cent,
  clustering_coeff = clust_coeff,
  betweenness = btwn_norm,
  kcore = kcore_idx
)
gene_topo <- merge(gene_topo, gnomad[, .(gene, loeuf)], by = "gene", all.x = TRUE)
gene_topo <- merge(gene_topo, pubmed, by = "gene", all.x = TRUE)

fwrite(gene_topo, file.path(OUTPUT, "gene_topology_metrics.csv"))
cat(sprintf("  Saved gene_topology_metrics.csv (%d genes)\n", nrow(gene_topo)))


# ==========================================================================
# 2. Pathway-level metrics & TCM
# ==========================================================================
cat("\n=== [2/7] Computing pathway metrics & Topological Coupling Model ===\n")

# Load Reactome GMT
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

# Filter to 10-300 genes
pw_genes_list <- Filter(function(x) length(x$genes) >= 10 & length(x$genes) <= 300,
                         pw_genes_list)
cat(sprintf("  %d Reactome pathways (10-300 genes)\n", length(pw_genes_list)))

# Load Reactome hierarchy for super-pathway classification
hierarchy <- fread(file.path(DATA, "ReactomePathwaysRelation.txt"),
                   header = FALSE, col.names = c("parent", "child"))
pw_names <- fread(file.path(DATA, "ReactomePathways.txt"),
                  header = FALSE, col.names = c("id", "name", "species"))
pw_names <- pw_names[species == "Homo sapiens"]

# Find top-level ancestors
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

# Super-pathway mapping
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

# Load auxiliary data
shet <- as.data.table(read_excel(
  file.path(DATA, "Cassa_2017_41588_2017_BFng3831_MOESM71_ESM.xlsx")))
shet <- shet[!is.na(s_het), .(gene = gene_symbol, s_het)]

paralogs <- fread(file.path(DATA, "paralog_counts.tsv"))
cds <- fread(file.path(DATA, "ensembl_cds_length.tsv"))
if ("cds_length" %in% names(cds)) {
  cds_max <- cds[, .(cds_length = max(cds_length, na.rm = TRUE)), by = gene]
} else {
  # Try gencode version
  cds_alt <- fread(file.path(DATA, "gencode_cds_lengths.tsv"))
  cds_max <- cds_alt[, .(cds_length = max(cds_length, na.rm = TRUE)), by = gene]
}
coords <- fread(file.path(DATA, "gene_coordinates.tsv"))

# Network genes set for fast lookup
net_genes <- V(g)$name

# Compute pathway-level metrics
cat("  Computing per-pathway metrics...\n")
pw_metrics <- rbindlist(lapply(pw_genes_list, function(pw) {
  genes <- pw$genes
  genes_in_net <- intersect(genes, net_genes)
  n <- length(genes)
  n_net <- length(genes_in_net)

  # Density: edges among pathway genes / possible edges
  if (n_net >= 2) {
    sg <- induced_subgraph(g, genes_in_net)
    e_int <- ecount(sg)
    density <- 2 * e_int / (n_net * (n_net - 1))
    # Insularity
    total_deg <- sum(degree(g, genes_in_net))
    insularity <- if (total_deg > 0) 2 * e_int / total_deg else 0
  } else {
    density <- 0
    insularity <- 0
  }

  # Mean LOEUF
  loeuf_vals <- gnomad[gene %in% genes]$loeuf
  mean_loeuf <- if (length(loeuf_vals) > 0) mean(loeuf_vals, na.rm = TRUE) else NA_real_
  loeuf_var <- if (length(loeuf_vals) > 1) var(loeuf_vals, na.rm = TRUE) else NA_real_

  # Mean paralog
  para_vals <- paralogs[gene %in% genes]$paralog_count
  mean_paralog <- if (length(para_vals) > 0) mean(para_vals, na.rm = TRUE) else NA_real_

  # Mean CDS length
  cds_vals <- cds_max[gene %in% genes]$cds_length
  mean_cds <- if (length(cds_vals) > 0) mean(cds_vals, na.rm = TRUE) else NA_real_

  # Genomic clustering score
  gene_coords <- coords[gene %in% genes]
  if (nrow(gene_coords) >= 2) {
    pairs <- combn(seq_len(nrow(gene_coords)), 2)
    same_chr <- gene_coords$chromosome[pairs[1, ]] == gene_coords$chromosome[pairs[2, ]]
    close <- abs(gene_coords$start_position[pairs[1, ]] -
                   gene_coords$start_position[pairs[2, ]]) < 1e6
    clust_score <- sum(same_chr & close) / ncol(pairs)
  } else {
    clust_score <- NA_real_
  }

  # Mean PubMed
  pub_vals <- pubmed[gene %in% genes]$pubmed_count
  mean_pub <- if (length(pub_vals) > 0) mean(pub_vals, na.rm = TRUE) else NA_real_

  # Super-pathway
  top_id <- get_top_ancestor(pw$id)
  top_name <- pw_names[id == top_id]$name
  if (length(top_name) == 0) top_name <- "Other"
  super <- if (top_name %in% names(SUPER_MAP)) SUPER_MAP[top_name] else "Other"

  data.table(
    pathway = pw$name,
    reactome_id = pw$id,
    n_genes = n,
    density = density,
    mean_loeuf = mean_loeuf,
    loeuf_variance = loeuf_var,
    insularity = insularity,
    super_pathway = unname(super),
    mean_paralog = mean_paralog,
    mean_cds_length = mean_cds,
    clustering_score = clust_score,
    mean_pubmed = mean_pub
  )
}))

cat(sprintf("  %d pathways with metrics\n", nrow(pw_metrics)))


# ---- TCM: genome-wide baseline -------------------------------------------
cat("  Fitting genome-wide baseline s_het ~ log10(degree)...\n")

# Gene-level: degree + s_het
gene_degree <- data.table(
  gene = V(g)$name,
  degree = degree(g)
)
gene_for_tcm <- merge(gene_degree, shet, by = "gene")
gene_for_tcm <- gene_for_tcm[degree > 0]
gene_for_tcm[, log_degree := log10(degree)]

baseline <- lm(s_het ~ log_degree, data = gene_for_tcm)
cat(sprintf("  Baseline: s_het = %.4f + %.4f * log10(degree), R² = %.4f\n",
            coef(baseline)[1], coef(baseline)[2], summary(baseline)$r.squared))

gene_for_tcm[, residual := residuals(baseline)]

# ---- TCM: per-pathway kappa deviations -----------------------------------
cat("  Computing per-pathway kappa deviations...\n")
MIN_GENES_KAPPA <- 15

pw_kappa_list <- lapply(pw_genes_list, function(pw) {
  genes_with_data <- gene_for_tcm[gene %in% pw$genes]
  if (nrow(genes_with_data) < MIN_GENES_KAPPA) return(NULL)

  # Mean residual
  mean_res <- mean(genes_with_data$residual)

  # Deviation kappa: slope of residual ~ log10(degree) within pathway
  if (nrow(genes_with_data) >= 3 && var(genes_with_data$log_degree) > 0) {
    kappa_fit <- lm(residual ~ log_degree, data = genes_with_data)
    dev_kappa <- coef(kappa_fit)[2]
    dev_kappa_p <- summary(kappa_fit)$coefficients[2, 4]
    r2 <- summary(kappa_fit)$r.squared
  } else {
    dev_kappa <- NA_real_
    dev_kappa_p <- NA_real_
    r2 <- NA_real_
  }

  data.table(
    pathway = pw$name,
    mean_residual = mean_res,
    deviation_kappa = unname(dev_kappa),
    deviation_kappa_p = unname(dev_kappa_p),
    r2 = unname(r2),
    n_genes = nrow(genes_with_data),
    mean_shet = mean(genes_with_data$s_het),
    mean_degree = mean(genes_with_data$degree)
  )
})
pw_kappa <- rbindlist(Filter(Negate(is.null), pw_kappa_list))

# Add density, mean_loeuf, super_pathway from pw_metrics
pw_kappa <- merge(pw_kappa,
                  pw_metrics[, .(pathway, density, mean_loeuf, super_pathway)],
                  by = "pathway", all.x = TRUE)

cat(sprintf("  %d pathways with TCM fit\n", nrow(pw_kappa)))

# ---- Rectangular double-outlier regime assignment -------------------------
cat("  Classifying regimes via rectangular marginal thresholds...\n")
kappa_complete <- pw_kappa[!is.na(mean_residual) & !is.na(deviation_kappa)]

# Marginal percentile thresholds: a pathway must be extreme on BOTH axes
OUTLIER_PROB <- 0.90
resid_upper <- quantile(kappa_complete$mean_residual, OUTLIER_PROB)
resid_lower <- quantile(kappa_complete$mean_residual, 1 - OUTLIER_PROB)
kappa_upper <- quantile(kappa_complete$deviation_kappa, OUTLIER_PROB)
kappa_lower <- quantile(kappa_complete$deviation_kappa, 1 - OUTLIER_PROB)

cat(sprintf("  Thresholds: resid [%.4f, %.4f], kappa [%.4f, %.4f]\n",
            resid_lower, resid_upper, kappa_lower, kappa_upper))

kappa_complete[, regime := fifelse(
  mean_residual > resid_upper & deviation_kappa > kappa_upper, "Topology-Amplified",
  fifelse(
    mean_residual > resid_upper & deviation_kappa < kappa_lower, "Decoupled",
    fifelse(
      mean_residual < resid_lower & deviation_kappa < kappa_lower, "Buffered",
      "Equilibrium"
    )
  )
)]

cat("  Regime counts:\n")
print(kappa_complete[, .N, by = regime][order(-N)])

pw_kappa_out <- copy(kappa_complete)
fwrite(pw_kappa_out, file.path(OUTPUT, "pathway_kappa_model.csv"))
cat(sprintf("  Saved pathway_kappa_model.csv (%d rows)\n", nrow(pw_kappa_out)))


# ==========================================================================
# 3. Full pathway_regimes_classified.csv (with tau, pLI, complex_fraction)
# ==========================================================================
cat("\n=== [3/7] Building pathway_regimes_classified.csv ===\n")

# Merge regime labels onto full pathway metrics
regime_labels <- kappa_complete[, .(pathway, regime)]
pw_classified <- merge(pw_metrics, regime_labels, by = "pathway", all.x = TRUE)
pw_classified[is.na(regime), regime := "Equilibrium"]

# Compute residual from LOEUF for the full set
loeuf_resid <- merge(pw_metrics[, .(pathway, mean_loeuf)],
                     pw_kappa_out[, .(pathway, mean_residual)],
                     by = "pathway", all.x = TRUE)
pw_classified <- merge(pw_classified,
                       loeuf_resid[, .(pathway, residual = mean_residual)],
                       by = "pathway", all.x = TRUE)

# Add tissue specificity (Tau)
tau_df <- fread(file.path(DATA, "gtex_tau_scores.tsv"))

# Add pLI (from gnomAD)
pli_map <- gnomad[!is.na(pLI), .(gene, pLI)]

# Add CORUM complex genes
corum_genes <- fread(file.path(DATA, "corum_complex_genes.tsv"))$gene

# Compute per-pathway: mean_tau, mean_pli, complex_fraction
cat("  Adding mean_tau, mean_pli, complex_fraction...\n")
tau_lookup <- setNames(tau_df$tau, tau_df$gene)
pli_lookup <- setNames(pli_map$pLI, pli_map$gene)

extra_metrics <- rbindlist(lapply(pw_genes_list, function(pw) {
  genes <- pw$genes

  # Mean tau
  tau_vals <- tau_lookup[intersect(genes, names(tau_lookup))]
  tau_vals <- tau_vals[!is.na(tau_vals)]
  mt <- if (length(tau_vals) > 0) mean(tau_vals) else NA_real_

  # Mean pLI
  pli_vals <- pli_lookup[intersect(genes, names(pli_lookup))]
  pli_vals <- pli_vals[!is.na(pli_vals)]
  mp <- if (length(pli_vals) > 0) mean(pli_vals) else NA_real_

  # Complex fraction
  n <- length(genes)
  cf <- if (n > 0) sum(genes %in% corum_genes) / n else NA_real_

  data.table(pathway = pw$name, mean_tau = mt, mean_pli = mp, complex_fraction = cf)
}))

pw_classified <- merge(pw_classified, extra_metrics, by = "pathway", all.x = TRUE)
pw_classified <- pw_classified[order(-residual, na.last = TRUE)]

fwrite(pw_classified, file.path(OUTPUT, "pathway_regimes_classified.csv"))
cat(sprintf("  Saved pathway_regimes_classified.csv (%d rows, %d columns)\n",
            nrow(pw_classified), ncol(pw_classified)))


# ==========================================================================
# 4-5. GSEA on residual and kappa axes
# ==========================================================================
cat("\n=== [4/7] Running GSEA (residual axis) ===\n")

# Build gene-level rankings
# Each gene → mean pathway-level mean_residual (across its pathways)
gene_pw_map <- rbindlist(lapply(pw_genes_list, function(pw) {
  data.table(gene = pw$genes, pathway = pw$name)
}))
gene_pw_map <- merge(gene_pw_map,
                     pw_kappa_out[, .(pathway, mean_residual, deviation_kappa)],
                     by = "pathway")

# Residual axis: mean of mean_residual across gene's pathways
gene_alpha_rank <- gene_pw_map[, .(score = mean(mean_residual, na.rm = TRUE)), by = gene]
gene_alpha_rank <- gene_alpha_rank[!is.na(score)][order(-score)]
alpha_ranks <- setNames(gene_alpha_rank$score, gene_alpha_rank$gene)

# Kappa axis: mean of deviation_kappa across gene's pathways
gene_kappa_rank <- gene_pw_map[, .(score = mean(deviation_kappa, na.rm = TRUE)), by = gene]
gene_kappa_rank <- gene_kappa_rank[!is.na(score)][order(-score)]
kappa_ranks <- setNames(gene_kappa_rank$score, gene_kappa_rank$gene)

# Fetch gene sets from Enrichr (matching the Python pipeline)
fetch_enrichr_gmt <- function(lib_name) {
  url <- sprintf("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=%s",
                 lib_name)
  cache <- file.path(DATA, paste0(lib_name, ".gmt"))
  if (!file.exists(cache)) {
    cat(sprintf("  Downloading %s from Enrichr...\n", lib_name))
    download.file(url, cache, quiet = TRUE)
  }
  lines <- readLines(cache)
  sets <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    parts[parts != "" & seq_along(parts) > 2]
  })
  names(sets) <- sapply(lines, function(x) strsplit(x, "\t")[[1]][1])
  Filter(function(x) length(x) > 0, sets)
}

go_sets <- fetch_enrichr_gmt("GO_Biological_Process_2023")
kegg_sets <- fetch_enrichr_gmt("KEGG_2021_Human")
all_sets <- c(go_sets, kegg_sets)

# Run fgsea on alpha ranks
cat(sprintf("  %d genes ranked, %d gene sets\n", length(alpha_ranks), length(all_sets)))
set.seed(42)
gsea_alpha <- fgsea(
  pathways = all_sets,
  stats = alpha_ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 1000
)

gsea_alpha_out <- gsea_alpha[, .(
  Name = pathway,
  Term = pathway,
  ES = ES,
  NES = NES,
  `NOM p-val` = pval,
  `FDR q-val` = padj,
  `FWER p-val` = pval,
  `Tag %` = NA_real_,
  `Gene %` = NA_real_,
  Lead_genes = sapply(leadingEdge, function(x) paste(x, collapse = ";")),
  library = fifelse(grepl("^KEGG_", pathway), "KEGG_2021_Human",
                    "GO_Biological_Process_2023")
)]

fwrite(gsea_alpha_out, file.path(OUTPUT, "gsea_alpha_rank_results.csv"))
cat(sprintf("  Saved gsea_alpha_rank_results.csv (%d terms, %d sig at FDR<0.05)\n",
            nrow(gsea_alpha_out), sum(gsea_alpha_out$`FDR q-val` < 0.05, na.rm = TRUE)))

cat("\n=== [5/7] Running GSEA (kappa axis) ===\n")
set.seed(42)
gsea_kappa <- fgsea(
  pathways = all_sets,
  stats = kappa_ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 1000
)

gsea_kappa_out <- gsea_kappa[, .(
  Name = pathway,
  Term = pathway,
  ES = ES,
  NES = NES,
  `NOM p-val` = pval,
  `FDR q-val` = padj,
  `FWER p-val` = pval,
  `Tag %` = NA_real_,
  `Gene %` = NA_real_,
  Lead_genes = sapply(leadingEdge, function(x) paste(x, collapse = ";")),
  library = fifelse(grepl("^KEGG_", pathway), "KEGG_2021_Human",
                    "GO_Biological_Process_2023")
)]

fwrite(gsea_kappa_out, file.path(OUTPUT, "gsea_kappa_rank_results.csv"))
cat(sprintf("  Saved gsea_kappa_rank_results.csv (%d terms, %d sig at FDR<0.05)\n",
            nrow(gsea_kappa_out), sum(gsea_kappa_out$`FDR q-val` < 0.05, na.rm = TRUE)))


# ==========================================================================
# 6. HPO enrichment
# ==========================================================================
cat("\n=== [6/7] HPO phenotype enrichment ===\n")

# Gene-to-regime assignment (each gene → regime of its most constrained pathway)
gene_regime <- merge(
  gene_pw_map[, .(gene, pathway)],
  pw_classified[, .(pathway, regime, mean_loeuf)],
  by = "pathway"
)[order(mean_loeuf)][, .SD[1], by = gene][, .(gene, regime)]

# Load HPO
hpo <- fread(file.path(DATA, "genes_to_phenotype.txt"))
setnames(hpo, c("ncbi_gene_id", "gene_symbol", "hpo_id", "hpo_name",
                "frequency", "disease_id"))

# Filter HPO terms with >= 20 annotated genes
hpo_terms <- hpo[, .(n_genes = uniqueN(gene_symbol)), by = .(hpo_id, hpo_name)]
hpo_terms <- hpo_terms[n_genes >= 20]

# Merge HPO with regime
hpo_gene <- merge(hpo[, .(gene_symbol, hpo_id, hpo_name)],
                  gene_regime, by.x = "gene_symbol", by.y = "gene")
hpo_gene <- unique(hpo_gene)

# Enrichment: outlier regime vs Equilibrium (Fisher's exact test)
outlier_regimes <- c("Buffered", "Decoupled", "Topology-Amplified")

equil_genes <- gene_regime[regime == "Equilibrium"]$gene
n_equil <- length(equil_genes)

hpo_results <- rbindlist(lapply(outlier_regimes, function(reg) {
  reg_genes <- gene_regime[regime == reg]$gene
  n_reg <- length(reg_genes)

  rbindlist(lapply(seq_len(nrow(hpo_terms)), function(i) {
    hid <- hpo_terms$hpo_id[i]
    hname <- hpo_terms$hpo_name[i]
    hpo_genes_all <- unique(hpo[hpo_id == hid]$gene_symbol)

    reg_hits <- sum(reg_genes %in% hpo_genes_all)
    equil_hits <- sum(equil_genes %in% hpo_genes_all)

    # 2x2 table: regime hit/miss vs equilibrium hit/miss
    mat <- matrix(c(reg_hits, n_reg - reg_hits,
                     equil_hits, n_equil - equil_hits), nrow = 2)
    ft <- fisher.test(mat, alternative = "two.sided")

    data.table(
      hpo_id = hid,
      hpo_name = hname,
      regime = reg,
      regime_hits = reg_hits,
      regime_total = n_reg,
      equil_hits = equil_hits,
      equil_total = n_equil,
      regime_pct = 100 * reg_hits / n_reg,
      equil_pct = 100 * equil_hits / n_equil,
      odds_ratio = ft$estimate,
      p_value = ft$p.value
    )
  }))
}))

# FDR correction within each regime
hpo_results[, fdr := p.adjust(p_value, method = "BH"), by = regime]

fwrite(hpo_results, file.path(OUTPUT, "hpo_enrichment_results.csv"))
cat(sprintf("  Saved hpo_enrichment_results.csv (%d rows, %d sig at FDR<0.1)\n",
            nrow(hpo_results), sum(hpo_results$fdr < 0.1)))


# ==========================================================================
# 7. ClinGen haploinsufficiency overlap
# ==========================================================================
cat("\n=== [7/7] ClinGen haploinsufficiency overlap ===\n")

clingen <- fread(file.path(DATA, "ClinGen_gene_curation_list_GRCh38.tsv"),
                 skip = 5, sep = "\t")
# Score = 3 means sufficient evidence
hi_genes <- clingen[`Haploinsufficiency Score` == "3"]$`#Gene Symbol`
hi_genes <- intersect(hi_genes, gene_regime$gene)
cat(sprintf("  %d ClinGen HI (score=3) genes in network\n", length(hi_genes)))

clingen_results <- gene_regime[, .(
  hi_count = sum(gene %in% hi_genes),
  total = .N
), by = regime]
clingen_results[, pct := 100 * hi_count / total]

# Odds ratio vs Equilibrium
equil_hi <- clingen_results[regime == "Equilibrium"]$hi_count
equil_non <- clingen_results[regime == "Equilibrium"]$total - equil_hi

clingen_results[, `:=`(
  odds_ratio = mapply(function(hi, tot) {
    non <- tot - hi
    if (equil_hi == 0 || non == 0) return(NA_real_)
    (hi / non) / (equil_hi / equil_non)
  }, hi_count, total),
  p_value = mapply(function(hi, tot) {
    non <- tot - hi
    mat <- matrix(c(hi, non, equil_hi, equil_non), nrow = 2)
    fisher.test(mat)$p.value
  }, hi_count, total)
)]

fwrite(clingen_results, file.path(OUTPUT, "clingen_overlap_results.csv"))
cat(sprintf("  Saved clingen_overlap_results.csv\n"))
print(clingen_results)

# ==========================================================================
# 8. Multivariate OLS: mean s_het ~ density + confounders
# ==========================================================================
cat("\n=== [8/8] Multivariate OLS regression (s_het) ===\n")

# Merge mean_shet from pw_kappa into pw_metrics
ols_df <- merge(pw_metrics, pw_kappa[, .(pathway, mean_shet)],
                by = "pathway", all.x = TRUE)
ols_df <- ols_df[!is.na(mean_shet) & !is.na(density) &
                   !is.na(mean_paralog) & !is.na(mean_cds_length) &
                   !is.na(clustering_score)]

cat(sprintf("  %d pathways in primary model\n", nrow(ols_df)))

# Standardise predictors
predictors_primary <- c("density", "mean_paralog", "mean_cds_length", "clustering_score")
predictors_sens    <- c(predictors_primary, "mean_pubmed")

standardise <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Primary model
ols_primary_df <- copy(ols_df)
for (v in predictors_primary) ols_primary_df[[v]] <- standardise(ols_primary_df[[v]])
fit_primary <- lm(mean_shet ~ density + mean_paralog + mean_cds_length + clustering_score,
                  data = ols_primary_df)
s_primary <- summary(fit_primary)
cat(sprintf("  Primary: R² = %.4f, Adj R² = %.4f\n",
            s_primary$r.squared, s_primary$adj.r.squared))

# Sensitivity model (+ PubMed)
ols_sens_df <- ols_df[!is.na(mean_pubmed)]
for (v in predictors_sens) ols_sens_df[[v]] <- standardise(ols_sens_df[[v]])
fit_sens <- lm(mean_shet ~ density + mean_paralog + mean_cds_length +
                 clustering_score + mean_pubmed,
               data = ols_sens_df)
s_sens <- summary(fit_sens)
cat(sprintf("  Sensitivity: R² = %.4f, Adj R² = %.4f\n",
            s_sens$r.squared, s_sens$adj.r.squared))

# Build output table
coef_primary <- s_primary$coefficients
coef_sens    <- s_sens$coefficients

table_rows <- lapply(predictors_sens, function(v) {
  p_beta <- if (v %in% rownames(coef_primary)) coef_primary[v, "Estimate"] else NA_real_
  p_se   <- if (v %in% rownames(coef_primary)) coef_primary[v, "Std. Error"] else NA_real_
  p_p    <- if (v %in% rownames(coef_primary)) coef_primary[v, "Pr(>|t|)"] else NA_real_
  s_beta <- coef_sens[v, "Estimate"]
  s_se   <- coef_sens[v, "Std. Error"]
  s_p    <- coef_sens[v, "Pr(>|t|)"]
  data.table(predictor = v,
             beta_primary = p_beta, se_primary = p_se, p_primary = p_p,
             beta_sensitivity = s_beta, se_sensitivity = s_se, p_sensitivity = s_p)
})
table_ols <- rbindlist(table_rows)

# Add model-level stats as attributes in header comment
fwrite(table_ols, file.path(OUTPUT, "table1_ols_shet.csv"))
cat(sprintf("  Saved table1_ols_shet.csv\n"))

# Also save model summary stats
model_stats <- data.table(
  model = c("primary", "sensitivity"),
  n = c(nrow(ols_primary_df), nrow(ols_sens_df)),
  r_squared = c(s_primary$r.squared, s_sens$r.squared),
  adj_r_squared = c(s_primary$adj.r.squared, s_sens$adj.r.squared),
  f_statistic = c(s_primary$fstatistic[1], s_sens$fstatistic[1]),
  f_p_value = c(
    pf(s_primary$fstatistic[1], s_primary$fstatistic[2], s_primary$fstatistic[3], lower.tail = FALSE),
    pf(s_sens$fstatistic[1], s_sens$fstatistic[2], s_sens$fstatistic[3], lower.tail = FALSE)
  )
)
fwrite(model_stats, file.path(OUTPUT, "table1_ols_shet_summary.csv"))
cat(sprintf("  Saved table1_ols_shet_summary.csv\n"))
print(table_ols)
print(model_stats)

cat("\nDONE.\n")
