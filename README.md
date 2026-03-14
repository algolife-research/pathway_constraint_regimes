# pathway_constraint_regimes

Code and data for:

> Gouy (2026). Beyond the Cost of Complexity: Gene-Intrinsic versus
> Topological Constraint in Regulatory Networks.

## Repository structure

```
scripts/          R scripts
data/             Input data (see data/README.md for sources)
output/           Intermediate CSV tables
figures/          Publication figures
```

## Requirements

**R >= 4.3** with the following packages:

```r
install.packages(c(
  "data.table", "igraph", "readxl", "fgsea",
  "tidyverse", "patchwork", "ggrepel", "ggsignif",
  "scales", "ggcorrplot", "future.apply"
))
```

## Rerunning the pipeline

### 1. Obtain input data

Download all primary data files listed in [`data/README.md`](data/README.md)
into the `data/` directory.

### 2. Run the analysis pipeline

All scripts auto-resolve paths relative to their own location, so they
can be run from anywhere. Execute them in this order:

```bash
# Step 1: Core analysis — builds PPI network, computes topology metrics,
#         fits the Topological Coupling Model, classifies regimes, runs
#         GSEA, HPO enrichment, ClinGen overlap, and OLS regression.
#         Outputs 9 CSVs to output/.
Rscript scripts/01_main_analysis.R

# Step 2: Supplementary analysis — permutation null model, jackknife
#         stability, STRING replication, Jaccard independence clustering,
#         DFE validation, threshold sensitivity.
#         Outputs 9 CSVs to output/.
Rscript scripts/02_prepare_supp_tables.R

# Step 3 (optional): LaTeX tables for manuscript.
#         Outputs .tex files to manuscript/.
Rscript scripts/03_prepare_supp_latex_tables.R

# Step 4: Main figures 1-4.
Rscript scripts/04_figures_main.R

# Step 5: Supplementary figures S1-S7.
Rscript scripts/05_figures_supp.R
```

## Output files

### Intermediate tables (`output/`)

| File | Description |
|------|-------------|
| `gene_topology_metrics.csv` | Per-gene centrality metrics in BioGRID PPI |
| `pathway_kappa_model.csv` | Per-pathway TCM fit (mean residual, deviation kappa, regime) |
| `pathway_regimes_classified.csv` | Full pathway metadata with regime labels |
| `gsea_alpha_rank_results.csv` | GSEA on residual axis (GO + KEGG) |
| `gsea_kappa_rank_results.csv` | GSEA on kappa axis (GO + KEGG) |
| `hpo_enrichment_results.csv` | HPO phenotype enrichment per regime |
| `clingen_overlap_results.csv` | ClinGen haploinsufficiency overlap per regime |
| `table1_ols_shet.csv` | OLS regression coefficients |
| `table1_ols_shet_summary.csv` | OLS model-level statistics |
| `pathway_regime_metrics.csv` | Regime-level summary statistics |
| `pathway_dfe_validation.csv` | DFE parameter validation |
| `supp_mahalanobis_sensitivity.csv` | Rectangular threshold sensitivity |
| `supp_cohens_d_degree_bins.csv` | Cohen's d by degree bins |
| `supp_pathway_independence.csv` | Jaccard independence clustering |
| `supp_jackknife_stability.csv` | Leave-one-out jackknife stability |
| `supp_string_replication.csv` | STRING network replication |
| `supp_string_concordance.csv` | STRING-BioGRID regime concordance |
| `supp_permutation_null.csv` | Permutation null distribution |

### Figures (`figures/`)

| File | Description |
|------|-------------|
| `figure1_kappa_model.png` | TCM deviation space and regime distributions |
| `figure2_structural_context.png` | Dosage sensitivity and structural context by regime |
| `figure3_hub_fallacy.png` | Topology-fitness decoupling |
| `figure4_functional_validation.png` | HPO and ClinGen functional validation |
| `supp_s1_topological_landscape.png` | Gene-level topological landscape |
| `supp_s2_pathway_landscape.png` | Pathway landscape and bias exclusion |
| `supp_s3_dfe_validation.png` | DFE parameter validation |
| `supp_s4_permutation_null.png` | Permutation null model |
| `supp_s5_jackknife_stability.png` | Jackknife stability |
| `supp_s6_string_replication.png` | STRING network replication |
| `supp_s7_pathway_independence.png` | Pathway independence clustering |

## License

See [LICENSE](LICENSE).
