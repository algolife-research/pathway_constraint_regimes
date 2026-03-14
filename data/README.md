# Data directory

All data files in this directory are gitignored. This README documents each file, its source, and how to obtain it.

## Primary data (manual download)

These files must be downloaded manually before running any scripts.

### `Cassa_2017_41588_2017_BFng3831_MOESM71_ESM.xlsx`

Heterozygous selection coefficients (s_het) from Cassa et al. 2017,
*Nature Genetics*.
- **Source:** Supplementary Data from
  [doi:10.1038/ng.3831](https://doi.org/10.1038/ng.3831)
- **Download:** Manual, from the paper's supplementary materials on
  Nature's website.
- **Used by:** `prepare_tables.R`, `prepare_supp_tables.R`,
  `figures_main.R`, `figures_supp.R`, `string_threshold_sensitivity.R`,
  `string_channel_sensitivity.R`

### `biogrid_latest.tab3.zip`

BioGRID latest release, all organisms (filtered to human physical
interactions at runtime).
- **Source:** [BioGRID](https://thebiogrid.org/)
- **URL:**
  `https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ALL-LATEST.tab3.zip`
- **Download:** Manual or `wget`.
- **Used by:** `prepare_tables.R`, `prepare_supp_tables.R`

### `gnomad_constraint.txt.bgz`

gnomAD v2.1.1 gene constraint metrics (LOEUF, pLI, missense Z-scores).
- **Source:** [gnomAD](https://gnomad.broadinstitute.org/)
- **URL:**
  `https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz`
- **Download:** Manual or `wget`.
- **Used by:** `prepare_tables.R`, `prepare_supp_tables.R`

### `9606.protein.links.detailed.v12.0.txt.gz`

STRING v12.0 human protein-protein interaction links with per-channel
evidence scores (experimental, database, textmining, coexpression, etc.).
- **Source:** [STRING database](https://string-db.org/cgi/download?sessionId=&species_text=Homo+sapiens)
- **Download:** Manual, from the STRING download page for species 9606.
- **Used by:** `prepare_supp_tables.R`, `string_channel_sensitivity.R`

### `9606.protein.info.v12.0.txt.gz`

STRING v12.0 human protein info (Ensembl protein ID to gene name mapping).
- **Source:** [STRING database](https://string-db.org/cgi/download?sessionId=&species_text=Homo+sapiens)
- **Download:** Manual, from the STRING download page for species 9606.
- **Used by:** `prepare_supp_tables.R`, `string_threshold_sensitivity.R`,
  `string_channel_sensitivity.R`

### `ensembl_cds_length.tsv`

CDS lengths per gene from Ensembl BioMart.
- **Source:** [Ensembl BioMart](https://www.ensembl.org/biomart)
- **Download:** Manual export from BioMart (Gene stable ID, Gene name,
  CDS length).
- **Used by:** `prepare_tables.R`, `prepare_supp_tables.R`

### `genes_to_phenotype.txt`

HPO gene-to-phenotype annotations.
- **Source:** [Human Phenotype Ontology](https://hpo.jax.org/)
- **URL:**
  `https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt`
- **Download:** Manual or `wget`.
- **Used by:** `prepare_tables.R`

### `ClinGen_gene_curation_list_GRCh38.tsv`

ClinGen gene curation list with haploinsufficiency and triplosensitivity
scores.
- **Source:** [ClinGen](https://clinicalgenome.org/)
- **URL:**
  `https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv`
- **Download:** Manual or `wget`.
- **Used by:** `prepare_tables.R`

### Reactome files

Three files from Reactome are required:

- **`ReactomePathways.gmt.zip`** — Pathway gene sets in GMT format.
  URL: `https://reactome.org/download/current/ReactomePathways.gmt.zip`
- **`ReactomePathwaysRelation.txt`** — Pathway hierarchy
  (parent-child relationships).
  URL: `https://reactome.org/download/current/ReactomePathwaysRelation.txt`
- **`ReactomePathways.txt`** — Pathway names and species.
  URL: `https://reactome.org/download/current/ReactomePathways.txt`
