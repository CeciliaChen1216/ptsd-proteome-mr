# Proteome-wide Mendelian randomization of plasma proteins for PTSD

Code and analysis pipeline for a proteome-wide two-sample Mendelian randomization
(MR) study of plasma proteins and post-traumatic stress disorder (PTSD), extended
with cross-disorder Genomic SEM (a common internalizing factor over PTSD, major
depressive disorder [MDD], and anxiety disorders), colocalization, eQTL-SMR,
brain-pQTL, and childhood-maltreatment secondary-phenotype analyses.

---

## Overview

Plasma cis-pQTL instruments from the UK Biobank Pharma Proteomics Project
(UKB-PPP) are tested against the PGC-PTSD Freeze 3 GWAS to prioritize candidate
proteins, which are then followed up with:

- proteome-wide FDR/Bonferroni screening against PTSD, MDD, and anxiety disorders,
  with a per-outcome instrument-harmonization flow (rsID match + allele flip only;
  no proxy/position/palindromic recovery in the primary screen);
- Genomic SEM: a common internalizing factor (PTSD + MDD + anxiety), protein-to-
  common-factor MR, and Q_SNP heterogeneity testing;
- Bayesian colocalization (coloc) with a 1000 Genomes European LD reference;
- expression-QTL summary-data-based MR (eQTL-SMR + HEIDI) in brain (BrainMeta)
  and blood (Westra);
- cross-platform pQTL concordance with deCODE;
- brain pQTL replication (Wingo 2025) and observational/cell-type annotation
  (Human Protein Atlas; Daskalakis et al. 2024);
- a childhood-maltreatment secondary-phenotype (sensitivity) analysis;
- conditional/COJO analysis of the chromosome 15 FES/FURIN locus;
- drug-repurposing and target-tractability annotation (exploratory).

Cross-disorder effects are summarized with the common-factor model rather than
with categorical labels; earlier "transdiagnostic / PTSD-predominant"
classifications were removed and are not used in the manuscript or supplement.

---

## Data provenance

The manuscript's numeric results depend on the exact release/version of each
input GWAS. The versions used are:

| Dataset | File used | N (cases / controls) | Source |
|---|---|---|---|
| UKB-PPP cis-pQTLs | `ukbppp_cis_pqtl.rds` (54,219 participants; 1,860 proteins with genome-wide sentinel cis-pQTL) | — | Synapse syn51364943 (Sun et al. 2023, Nature) |
| PGC-PTSD Freeze 3 (EUR case-control) | `eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz` (VCF file date 2024-12-04) | 137,136 / 1,085,746 | <https://www.med.unc.edu/pgc/download-results/> (Nievergelt et al. 2024, Nat Genet) |
| PGC MDD 2025 (EUR, no-UKBB, no-23andMe) | `daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz` | 357,636 / 1,281,936 | Psychiatric Genomics Consortium (Dowsett et al. 2025, Cell) |
| PGC Anxiety 2026 (fullANX, wo-UTAH) | `ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz` | 122,083 / 729,602 | Psychiatric Genomics Consortium (Strom et al. 2026, Nat Genet) |

For PTSD, per-SNP `NCASE`, `NCON`, `NEFF`, and `NTOT` values are read directly
from the VCF INFO fields; the study-level totals shown above match the source
VCF `##nCase` and `##nControl` metadata, which correspond to the headline
European-ancestry meta-analysis reported in Nievergelt et al. 2024. Sample overlap
between UKB-PPP (a UK Biobank sub-cohort, n = 54,219) and the two UK Biobank
cohorts in PGC-PTSD Freeze 3 (`ukbb`: 10,913 cases + 124,888 controls; `ukb2`:
9,882 cases + 249,993 controls) is bounded at ~4.4% of the outcome sample and
does not exceed the strong-instrument regime for MR under Burgess et al. (2016);
this is discussed as a limitation in the manuscript.

---

## Repository layout

Active scripts (run in numeric order; numbering is not fully contiguous):

```
00_config.R                 Local paths (NOT committed; see template below)
00_config_template.R        Path template — copy to 00_config.R and edit
01_MR_screening_LOGGED.R    Proteome-wide MR + FDR + per-outcome harmonization-flow log
02_drug_repurposing.R       Known-drug / tractability annotation
03_deep_analyses.R          Additional downstream analyses
04_coloc_multiMR.R          Colocalization + multi-instrument MR
05_sensitivity_analyses.R   MR sensitivity analyses
06_biological_annotation.R  Biological / pathway annotation
07_mediation_network_MR.R   Mediation and protein-protein network MR
08_eQTL_SMR.R               Brain & blood eQTL-SMR (+ HEIDI)
09_figures_CLEAN.R          Main + supplementary figures (Fig 3A/3B volcano+forest,
                            Fig 4A/4B; Supp Fig 1, 2, 6, 7). Declassified
                            Fig 2 & Supp Fig 1; Fig 4B network MR de-emphasized
                            with Bonferroni note; Supp Fig 2 = four-state
                            evidence matrix (not summed).
10_tables.R                 Main-text tables
11_external_replication.R   deCODE cross-platform concordance
12_power_aware_classification.R  Legacy cross-disorder classification (retained as an
                            upstream candidate-list source only; its classification
                            output is no longer used in any figure or table)
13_CD40_target_safety.R     CD40 target-safety annotation
14_observational_validation.R    Observational / tissue validation
15_Figure3_celltype.R       Supplementary Figure S8 (cell-type expression &
                            PTSD dysregulation; formerly Figure 3)
16_celltype_brainregion_enrichment.R  Cell-type & brain-region enrichment;
                            writes Supp Fig 4 and region/cell-type summary tables
17_conditional_chr15_COJO.R Conditional / COJO analysis of the chr15 FES/FURIN locus
18_trauma_negative_control_MR.R        Childhood-maltreatment secondary-phenotype MR
18b_finalize_trauma_negative_control.R Finalize the secondary-phenotype analysis
19_brain_pqtl_replication_wingo2025.R  Brain pQTL replication (Wingo 2025)
20_generate_supplementary_tables.R     Supplementary tables (adds S20, S21;
                            other supplementary tables are curated separately
                            during revision — see Data and code availability)
21_append_conditional_chr15_table.R    Append chr15 conditional table
22_figure5.R                Figure 5 (KHK cross-tissue + childhood-maltreatment
                            secondary phenotype; outputs Figure_5_KHK_cross_tissue.*)
24_phewas_classify_supp.R   Supplementary pheWAS scan (Table S23)
25_genomicsem_stage0_feasibility.R     Genomic SEM: LDSC / feasibility. Uses the
                            PTSD case proportion from the source VCF metadata
                            (sample.prev = 137,136 / 1,222,882 = 0.1121).
26_genomicsem_stage0b_modelfix.R       Genomic SEM: model-fit (Heywood handling).
                            Retained for reproducibility of the constrained-model
                            decision; not required to rerun for the submission.
27_genomicsem_stage1_commonfactorGWAS.R Genomic SEM Stage 1: sumstats harmonization
                            (Block 1). Optional genome-wide factor GWAS (Block 2)
                            is OFF by default and NOT used for the manuscript;
                            if enabled it uses userGWAS with the explicit
                            constrained model (MDD residual fixed 0). Sample.prev
                            is not used in Block 1, so this script was not
                            rerun for the sample.prev correction.
28_genomicsem_stage12_factorMR.R       Genomic SEM Stage 1+2 (instrument-restricted):
                            userGWAS with the explicit constrained model
                            (MDD residual fixed 0) + Q_SNP=TRUE; protein-to-
                            common-factor Wald-ratio MR. Direct instruments only
                            (no proxy substitution).
make_figure2_SEM.R          MAIN Figure 2 (common internalizing factor: A
                            loadings, B protein-to-common-factor MR volcano,
                            C instrument-level Q_SNP distribution). R/ggplot;
                            reads stage12_factorMR_results.csv.
```

---

## Setup

1. Copy `00_config_template.R` to `00_config.R` and edit the paths to point at
   your local data directories. `00_config.R` is git-ignored because it contains
   machine-specific absolute paths.

2. Required external datasets (cited in the manuscript Data availability
   statement) and where the config expects them:

   | Variable          | Dataset |
   |-------------------|---------|
   | `pqtl_path`       | UKB-PPP cis-pQTL instruments (`ukbppp_cis_pqtl.rds`) |
   | `ptsd_path` / `mdd_path` / `anx_path` | PGC-PTSD Freeze 3, PGC-MDD 2025 (EUR, no-UKBB, no-23andMe), PGC-anxiety 2026 GWAS summary statistics (see Data provenance) |
   | `ld_ref_path`     | 1000 Genomes Phase 3 EUR (plink bfile prefix; SMR uses n = 503) |
   | `brain_eqtl_dir`  | BrainMeta cis-eQTL summary (SMR `.besd/.epi/.esi`) |
   | `blood_eqtl_dir`  | Westra blood cis-eQTL (SMR format) |
   | `hpa_dir`         | Human Protein Atlas TSVs (see below) |
   | `decode_supp`     | deCODE variant-annotation file |
   | `result_dir`      | Output directory for CSV/RDS results |
   | `fig_dir`         | Output directory for figures |
   | `ldsc_dir`        | LDSC `eur_w_ld_chr/` directory (trailing slash) |
   | `hm3_path`        | LDSC HapMap3 SNP list (`w_hm3.snplist`) |
   | `ptsd_smr_ma`     | PTSD `.ma` file derived from the VCF (SMR-formatted) |

   Genomic SEM (scripts 25-29) additionally reads the munged sumstats /
   LDSC output written under `result_dir/genomicsem/`.

3. Human Protein Atlas files used by `16_celltype_brainregion_enrichment.R`
   (place the unzipped `.tsv` in `hpa_dir`):

   - `rna_single_cell_type.tsv` — single-cell-type expression (nCPM);
     <https://www.proteinatlas.org/download/tsv/rna_single_cell_type.tsv.zip>
   - `rna_brain_gtex.tsv` — brain-region expression (nTPM);
     <https://www.proteinatlas.org/download/tsv/rna_brain_gtex.tsv.zip>

   HPA download URLs change periodically; the script falls back gracefully if a
   file is missing and prints the expected location.

---

## Software

- R (>= 4.5); analyses were run under R 4.5.2
- Key packages: `TwoSampleMR` (0.7.0), `GenomicSEM` (0.0.5), `coloc`,
  `data.table`, `dplyr`, `tidyr`, `ggplot2`, `patchwork`, `ggrepel`, `ggtext`,
  `readr`, `openxlsx`, `readxl`
- SMR (v1.3.1) for eQTL-SMR (external binary; path set in config)
- plink (for the LD reference) where required
- All figure scripts are R (including `make_figure2_SEM.R` for the main Figure 2).

---

## Figures

| Manuscript figure | Produced by |
|-------------------|-------------|
| Figure 2 (A loadings, B factor-MR volcano, C Q_SNP distribution) | `make_figure2_SEM.R` |
| Figure 3 (A volcano, B per-disorder forest) | `09_figures_CLEAN.R` |
| Figure 4 (A eQTL-SMR, B network MR) | `09_figures_CLEAN.R` |
| Figure 5 (KHK cross-tissue + childhood-maltreatment secondary phenotype) | `22_figure5.R` |
| Supplementary Figure 1 (cross-disorder pathway [NES] similarity) | `09_figures_CLEAN.R` |
| Supplementary Figure 2 (four-state evidence matrix) | `09_figures_CLEAN.R` |
| Supplementary Figure 4 (baseline brain cell-type expression) | `16_celltype_brainregion_enrichment.R` |
| Supplementary Figure 6 (CRP mediation) | `09_figures_CLEAN.R` |
| Supplementary Figure 7 (drug repurposing / tractability, exploratory) | `09_figures_CLEAN.R` |
| Supplementary Figure 8 (cell-type expression & PTSD dysregulation; formerly Figure 3) | `15_Figure3_celltype.R` |

**Supplementary Figures 3 and 5.** Supplementary Figure 3 (within-gene regional
expression of CD40, FURIN, and SIRPA) was generated from Human Protein Atlas
brain-region expression data (values shown as within-gene Z-scores), and
Supplementary Figure 5 (cell-type differential expression in PTSD) from the
single-nucleus RNA-seq data of Daskalakis et al. (2024). The underlying source
data are cited in the manuscript Data availability statement;
`16_celltype_brainregion_enrichment.R` writes the corresponding region- and
cell-type-level summary tables. Plotting code for these two descriptive figures
is available from the corresponding author on request.

---

## Notes on the common-factor (Genomic SEM) analyses

- The common internalizing factor is estimated by diagonally weighted least
  squares (DWLS). In the unconstrained model the MDD residual variance was a
  non-significant negative (Heywood) estimate and was fixed to zero in the
  reported model.
- The common-factor SNP effects and Q_SNP were estimated with `userGWAS()` using
  an EXPLICIT constrained model (`F1 =~ NA*PTSD + MDD + ANX`; `F1 ~~ 1*F1`;
  `MDD ~~ 0*MDD`; `F1 ~ SNP`) that matches the trait-level report, run on the
  instrument SNP set only (`Q_SNP = TRUE`; genome-wide estimation not performed).
  The default `commonfactorGWAS()` was NOT used because it silently free-estimates
  the MDD residual (a different model).
- Standardized loadings from the constrained model (PTSD 0.921, MDD 1.018 at
  the Heywood boundary, anxiety 0.893) are the primary SEM output reported in
  the manuscript. The factor is MDD-weighted rather than symmetrically
  transdiagnostic: PTSD's standardized loading is lower than MDD's, and PTSD
  retains 15.1% residual genetic variance not captured by the common factor
  (P = 1.2 × 10⁻⁴). Anxiety retains 20.2% residual (P = 2.7 × 10⁻⁶).
- Protein-to-factor MR is a per-protein Wald ratio (factor-SNP effect / cis-pQTL
  effect). A valid instrument therefore requires a single variant carrying BOTH a
  cis-pQTL and a common-factor estimate. Four PTSD-prioritized proteins are not
  testable in this framework: SIRPA, FURIN, and UBE2L6 (their sentinel cis-pQTL is
  absent from the common-factor sumstats while the high-LD proxy lacks cis-pQTL
  data, so no single variant carries both effects), and CD101 (no adequate high-LD
  proxy). This leaves 391 testable proteins; **24** reach factor-MR BH-FDR < 0.05,
  and all six testable PTSD-prioritized candidates (AKT3, CD40, CGREF1, FES, KHK,
  SNX18) are FDR-significant. **Thirty-two** proteins (32/391 = 8.2%, close to
  the ~5% expected under the null) show nominal Q_SNP heterogeneity but none
  survive Q_SNP FDR correction. Among the six PTSD-prioritized candidates only
  SNX18 (Q_SNP P = 0.034) reaches nominal significance; the other five have
  Q_SNP P ∈ [0.11, 0.30]. Instrument allele harmonization was verified (all 391
  matched on the full effect/other-allele pair; no palindromic instruments).

---

## Notes on reproducibility

- `00_config.R` holds machine-specific paths and is not committed; use
  `00_config_template.R` as the starting point.
- Some external datasets (UKB-PPP, deCODE, BrainMeta, Westra, HPA, PGC GWAS)
  are controlled-access or large and are not redistributed here; obtain them
  from their original sources as cited in the manuscript.
- The COJO LD reference at the chr15 locus comprised 633 1000 Genomes European
  individuals (the standard Phase 3 European set, n = 503, plus additional
  1000 Genomes European [IBS/CEU] individuals); the eQTL-SMR analyses used the
  standard Phase 3 European unrelated set (n = 503).
- **LDSC / Genomic SEM sample-size and prevalence conventions (script 25).**
  The N supplied to `munge()` / `ldsc()` differs by trait and the sample
  prevalence is set to match it:
    * **PTSD** — total sample size (per-SNP N in the `.ma` file, up to 1,222,882);
      sample.prev = case proportion = 137,136 / 1,222,882 = **0.1121**
      (from the source VCF `##nCase` / `##nControl` metadata).
    * **MDD** — sum-of-cohort effective N (829,250; daner `.neff` column);
      sample.prev = 0.5 (case/control imbalance already absorbed into N_eff).
    * **ANX** — sum-of-cohort effective N (390,020 = 2 × Neff_half);
      sample.prev = 0.5.

  Population (lifetime, EUR) prevalences are approximate (PTSD 0.07, MDD 0.15,
  ANX 0.10). Standardized loadings and the genetic-correlation matrix are robust
  to these liability-scale choices; only non-standardized quantities (e.g. the
  liability-scale PTSD residual variance ≈ 0.007) depend on them.

- HPA expression values can differ slightly between HPA data releases; cell-type
  and brain-region figures reflect the release downloaded at analysis time.

---

## Update log

### 2026-07-18 — PTSD sample.prev correction and downstream reruns

**Cause.** A prior version of `00_config.R` contained an incorrect PTSD case
count in a comment (118,203 / 1,072,525). That value was propagated as
`sample.prev = 0.0993` in `25_genomicsem_stage0_feasibility.R`. Inspection of
the source VCF (`eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz`) confirmed
`##nCase = 137,136` and `##nControl = 1,085,746`, matching the Nievergelt et al.
2024 Nature Genetics headline. The 118,203 number was the mode of the per-SNP
`NCASE` column for a subset of SNPs, not the study-level total.

**Correction.** `sample.prev` was set to 137,136 / 1,222,882 = 0.1121, and
`n_ptsd` in `00_config.R` was updated to 1,222,882.

**Reruns performed.**

1. **Script 25** (Stage 0 LDSC + common-factor model) — `LDSCoutput_stage0.rds`
   regenerated (~25 min: munging is the slow step).
2. **Script 27** (Stage 1) — not rerun: `sumstats()` harmonization does not use
   `sample.prev`, and `stage1_sumstats.rds` remained valid.
3. **Script 28** (Stage 2 factor MR + Q_SNP on 387 instrument SNPs) — 
   `stage12_factorMR_results.csv` regenerated (~2.5 min serial).
4. **`make_figure2_SEM.R`** — Figure 2 regenerated with corrected numbers.
5. Supplementary tables S24 / S25 / S26 updated in `Additional_file_2-supplementary_tables.xlsx`.

**Impact on manuscript numbers.**

| Quantity | Before | After |
|---|---|---|
| PTSD sample size in Methods | 118,203 / 1,072,525 | 137,136 / 1,085,746 |
| PTSD standardized loading | 0.926 | 0.921 |
| PTSD residual variance proportion | 14.2% (P = 5 × 10⁻⁴) | 15.1% (P = 1.2 × 10⁻⁴) |
| PTSD unstandardized residual | 0.003 | 0.007 |
| MDD standardized loading | 1.004 | 1.018 |
| Anxiety standardized loading | 0.897 | 0.893 |
| Anxiety residual variance proportion | 0.195 | 0.202 |
| Genetic correlations (PTSD-MDD / MDD-ANX / PTSD-ANX) | 0.94 / 0.91 / 0.82 | 0.938 / 0.909 / 0.823 |
| Factor-MR FDR<0.05 proteins | 28 | 24 |
| Nominal Q_SNP proteins | 38 | 32 |
| Nominal Q_SNP among 6 candidates | 6/6 | 1/6 (SNX18 only, P = 0.034) |

**Impact on conclusions.** All six PTSD-prioritized candidates remain factor-MR
FDR<0.05. Common-factor structure, KHK/CGREF1 and CD40 prioritization, and
FES/FURIN conditional interpretation are unchanged. The interpretation of
Q_SNP among the six candidates was strengthened, not weakened: what had appeared
as a 6/6 nominal-deviation pattern (previously discussed as consistent with
PTSD-selection bias) resolved into 1/6, which matches the ~5% null expectation
across the 391 testable proteins.

**Descriptive Z-score concordance metric removed.** A prior version of the
manuscript reported instrument-restricted SNP-effect Z-score correlations
between an optional genome-wide factor GWAS and each trait (0.61 / 0.91 / 0.82
for PTSD / MDD / anxiety). That analysis relied on Block 2 of Script 27
(genome-wide factor GWAS), which is disabled by default and was not rerun after
the sample.prev correction. The manuscript now describes factor–trait
asymmetry using the standardized loadings themselves, and Supplementary Table
S24 lists loadings and residual variance proportions rather than Z-score
concordance.

---

## Data and code availability

Summary statistics and external datasets are described in the manuscript Data
availability statement. UKB-PPP pQTL data are available via the UKB-PPP
consortium (controlled access). PGC GWAS summary statistics (PTSD Freeze 3,
MDD 2025, anxiety 2026) are available from the Psychiatric Genomics Consortium
(see Data provenance for exact file names and versions).

**Reproducibility note.** Several upstream inputs (UKB-PPP pQTL, PGC GWAS,
deCODE, BrainMeta, Westra, HPA) are controlled-access or too large to
redistribute, so the pipeline is not runnable end-to-end from raw data by third
parties without those access approvals. The primary proteome-wide MR screen is
run by `01_MR_screening_LOGGED.R`, and the derived per-outcome results are
written to `mr_all_outcomes.rds` / the harmonization-flow table. The
common-factor results are written to `stage12_factorMR_results.csv` by
Script 28. The final supplementary workbook (Additional File 2) was curated
from these derived outputs during revision and is provided as a frozen
deliverable rather than regenerated by a single end-to-end script. Accordingly,
the code-availability statement should read: "Scripts for all analyses are
provided, together with the derived summary-level inputs required to reproduce
the reported tables and figures where redistribution is permitted." The SMR
outputs are `eqtl_smr_brain_blood.csv` (per-tissue SMR/HEIDI) and
`eqtl_smr_brain_blood_fdr.csv` (adds per-tissue BH-FDR and graded verdicts).

The key frozen intermediates that reproduce the reported tables/figures are:

- `mr_all_outcomes.rds` — proteome-wide MR for PTSD/MDD/anxiety
- `stage12_factorMR_results.csv` — common-factor MR + Q_SNP for 391 proteins
- `LDSCoutput_stage0.rds` — Genomic SEM covariance (regenerated 2026-07-18 with
  corrected sample.prev)
- `stage1_sumstats.rds` — harmonized multivariate sumstats used by Script 28
- Harmonization-flow table for the proteome-wide screen
