###############################################################################
# 00_config_template.R — path & parameter configuration TEMPLATE
#
# Usage:
#   1. Copy this file to 00_config.R  (which is git-ignored).
#   2. Replace every /path/to/... below with an actual path on your machine.
#   3. Windows users: use forward slashes ("D:/PTSD/...") not backslashes.
#   4. Directory paths ending in / must keep the trailing slash where noted.
#
# All analysis scripts source 00_config.R at the top; no absolute path
# should appear anywhere else in the codebase.
###############################################################################

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Input data paths
# ═══════════════════════════════════════════════════════════════════════════════

## ---- Exposure: UKB-PPP cis-pQTL ------------------------------------------------
pqtl_path      <- "/path/to/ukbppp_cis_pqtl.rds"
pqtl_base      <- "/path/to/protein_full_gwas/"        # per-chr protein GWAS folder

## ---- Outcomes: three psychiatric-disorder GWAS ---------------------------------
# PTSD: PGC-PTSD Freeze 3 EUR case-control (Nievergelt 2024 Nat Genet)
# Confirm VCF metadata declares ##nCase = 137,136 / ##nControl = 1,085,746.
ptsd_path      <- "/path/to/eur_ptsdcasecontrol_pcs_v4_aug3_2021.vcf.gz"

# MDD: PGC MDD 2025 EUR, no-23andMe & no-UKBB (Dowsett 2025 Cell)
# Column names in .neff daner should contain FRQ_A_357636 / FRQ_U_1281936.
mdd_path       <- "/path/to/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz"

# Anxiety: PGC Anxiety 2026 fullANX, wo-UTAH (Strom 2026 Nat Genet)
# Column names in daner should contain FRQ_A_122083 / FRQ_U_729602.
anx_path       <- "/path/to/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz"

## ---- Secondary / sensitivity phenotypes ---------------------------------------
pcl_path       <- "/path/to/pcl_quantitative.gz"          # PCL continuous phenotype
aam_path       <- "/path/to/aam_ptsd.gz"                   # PGC-PTSD African-American
hna_path       <- "/path/to/hna_ptsd.gz"                   # PGC-PTSD Hispanic/Native
crp_path       <- "/path/to/crp_gwas/"                     # CRP GWAS folder

## ---- LD reference (1000G Phase 3 EUR plink prefix; SMR uses n = 503) -----------
ld_ref_path    <- "/path/to/1000G_v3/EUR"                  # plink prefix (no extension)

## ---- Output directories -------------------------------------------------------
cs_dir         <- "/path/to/susie_credible_sets"
result_dir     <- "/path/to/results/"                      # trailing slash required
fig_dir        <- "/path/to/Figures/"
table_dir      <- "/path/to/results/tables/"
validation_dir <- "/path/to/validation_data/"
hpa_dir        <- "/path/to/HPA/"

# ═══════════════════════════════════════════════════════════════════════════════
# 2. External software / reference files (LDSC / SMR / GenomicSEM)
# ═══════════════════════════════════════════════════════════════════════════════

plink_bin      <- "plink"                                  # or absolute path to plink.exe
smr_bin        <- "/path/to/smr-1.3.1/smr"                 # SMR v1.3.1 executable

# LDSC references (used by scripts 25, 27)
ldsc_dir       <- "/path/to/ldsc/eur_w_ld_chr/"            # trailing slash required
hm3_path       <- "/path/to/ldsc/w_hm3.snplist"

# GenomicSEM sumstats() reference (used by script 27 Block 1)
ref_sumstats   <- "/path/to/reference.custom.maf.0.005.txt.gz"

# PTSD .ma file derived from the source VCF (SMR-format; used by scripts 25/27/28)
ptsd_smr_ma    <- "/path/to/results/ptsd_freeze3_smr.ma"

# deCODE (Ferkingstad 2021) variant annotation (used by script 11)
decode_supp    <- "/path/to/deCODE/assocvariants.annotated.txt.gz"

# ═══════════════════════════════════════════════════════════════════════════════
# 3. eQTL references
# ═══════════════════════════════════════════════════════════════════════════════

brain_eqtl_dir <- "/path/to/BrainMeta_cis_eqtl_summary/"
blood_eqtl_dir <- "/path/to/blood_eqtl/"                   # westra_eqtl_hg19.{besd,epi,esi}

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Sample sizes (metadata constants; not used to compute results directly.
#    Verify these match the source GWAS files before running.)
# ═══════════════════════════════════════════════════════════════════════════════

n_pqtl   <- 54219        # UKB-PPP (Sun 2023, Nature)
n_ptsd   <- 1222882      # 137,136 cases + 1,085,746 controls
                         # (Nievergelt 2024 Nat Genet; VCF ##nCase / ##nControl)
n_mdd    <- 1639572      # 357,636 + 1,281,936  (Dowsett 2025 Cell)
n_anx    <- 851685       # 122,083 + 729,602    (Strom 2026 Nat Genet)
n_decode <- 35559        # deCODE Ferkingstad 2021

# NOTE on PTSD sample.prev used inside script 25:
# script 25 uses PTSD_SAMPPREV <- 137136 / (137136 + 1085746) = 0.1121
# (derived at runtime; do NOT hardcode it here.)

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Candidate protein list
# ═══════════════════════════════════════════════════════════════════════════════

candidates <- c("AKT3", "CD40", "CGREF1", "FES", "FURIN",
                "SIRPA", "CD101", "KHK", "SNX18", "UBE2L6")

# Sub-folder names under `pqtl_base` for the ten prioritized proteins.
# Adjust to match your local naming if downloaded from UKB-PPP directly.
protein_folders <- c(
  AKT3   = "AKT3_Q9Y243_OID21197_v1_Oncology",
  CD40   = "CD40_P25942_OID20724_v1_Inflammation",
  CGREF1 = "CGREF1_Q99674_OID20152_v1_Cardiometabolic",
  FES    = "FES_P07332_OID21207_v1_Oncology",
  FURIN  = "FURIN_P09958_OID21514_v1_Oncology",
  SIRPA  = "SIRPA_P78324_OID20304_v1_Cardiometabolic",
  CD101  = "CD101_Q93033_OID31480_v1_Oncology_II",
  KHK    = "KHK_P50053_OID30241_v1_Cardiometabolic_II",
  SNX18  = "SNX18_Q96RF0_OID31476_v1_Oncology_II",
  UBE2L6 = "UBE2L6_O14933_OID30321_v1_Cardiometabolic_II"
)

# ═══════════════════════════════════════════════════════════════════════════════
# 6. Analysis parameters
# ═══════════════════════════════════════════════════════════════════════════════

fdr_threshold   <- 0.05
nominal_p       <- 0.05
f_stat_min      <- 10
coloc_pph4_min  <- 0.8
heidi_p_min     <- 0.05
cis_window      <- 500000   # ±500kb cis definition
presso_n_perm   <- 5000
steiger_alpha   <- 0.05

cat("✓ 00_config.R loaded (template — replace /path/to/... entries before use)\n")
