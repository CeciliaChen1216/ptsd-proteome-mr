###############################################################################
# 25_genomicsem_stage0_feasibility.R
#
# Genomic SEM — STAGE 0 feasibility probe (cheap, no commonfactorGWAS).
#
# Question: does a common internalizing factor (PTSD + MDD + anxiety) hold, and
# does PTSD load lower than MDD/anxiety (i.e. is there a meaningful PTSD-specific
# component worth modelling)? If yes -> green light for Stage 1.
#
# Steps:
#   1. munge() the three GWAS to the LDSC .sumstats format
#   2. ldsc() to estimate the genetic covariance / correlation matrix
#   3. fit a one-factor model (commonfactor) and inspect loadings
#   4. diagnostics: rg matrix, heritability Z, loadings, PTSD-specific variance
#
# NOTE on sample size: MDD/anxiety are case-control daner files; effective N is
# passed explicitly (Neff), which is the correct LDSC scaling under unknown
# sample overlap. PTSD uses the per-SNP N in its .ma file. Stage 0 conclusions
# (factor structure, relative loadings) rest on the genetic-correlation matrix,
# which is robust to N scaling; revisit N carefully before any Stage 1 GWAS.
###############################################################################

suppressMessages({
  library(GenomicSEM)
  library(data.table)
})
source("00_config.R")

## ---- paths -----------------------------------------------------------------
hm3     <- hm3_path
out_dir <- file.path(result_dir, "genomicsem")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ptsd_ma <- ptsd_smr_ma
mdd_f   <- mdd_path
anx_f   <- anx_path

## ---- effective sample sizes (case-control) ---------------------------------
# MDD:  Neff column = 829249.6
# ANX:  Neff = 2 * Neff_half = 2 * 195010.1 = 390020.2
N_MDD <- 829250
N_ANX <- 390020

## ---- PTSD sample prevalence -----------------------------------------------
# PTSD uses per-SNP TOTAL N from the .ma file (converted from the source VCF's
# NCASE + NCON INFO fields), so LDSC receives a variable N per SNP for PTSD.
# sample.prev must still be a scalar; use the study-level case proportion
# reported in the VCF metadata (##nCase / ##nControl):
#   Ncase   = 137,136
#   Ncontrol = 1,085,746
#   Ntotal  = 1,222,882   (Nievergelt et al. 2024, Nat Genet)
#   sample.prev = 137136 / 1222882 = 0.1121
N_PTSD_CASE   <- 137136
N_PTSD_CTRL   <- 1085746
PTSD_SAMPPREV <- N_PTSD_CASE / (N_PTSD_CASE + N_PTSD_CTRL)   # 0.11215

## ---- 0. prep PTSD .ma into a munge-friendly file ---------------------------
# .ma columns: SNP A1 A2 freq b se p N  -> rename for munge clarity
ptsd <- fread(ptsd_ma)
setnames(ptsd,
         c("SNP","A1","A2","freq","b","se","p","N"),
         c("SNP","A1","A2","MAF","BETA","SE","P","N"))
ptsd_prep <- file.path(out_dir, "ptsd_freeze3_for_munge.txt")
fwrite(ptsd, ptsd_prep, sep = "\t")
cat("PTSD prepped:", nrow(ptsd), "SNPs ->", ptsd_prep, "\n")

## ---- 1. munge --------------------------------------------------------------
# munge writes <trait>.sumstats.gz into the working directory
setwd(out_dir)
files  <- c(ptsd_prep, mdd_f, anx_f)
trait_names <- c("PTSD", "MDD", "ANX")
N_vec  <- c(NA, N_MDD, N_ANX)   # NA = use per-SNP N column (PTSD .ma)

cat("\n=== munge() ===\n")
munge(files        = files,
      hm3          = hm3,
      trait.names  = trait_names,
      N            = N_vec,
      info.filter  = 0.9,
      maf.filter   = 0.01)

## ---- 2. LDSC genetic covariance --------------------------------------------
cat("\n=== ldsc() ===\n")
sumstats_files <- paste0(trait_names, ".sumstats.gz")
# sample prevalence / population prevalence for the case-control liability scale.
# All three are case-control; PTSD uses TOTAL per-SNP N (.ma), so its sample.prev
# is the actual case proportion, NOT NA. MDD/ANX use effective N (daner).
# Sample prevalence MUST match the type of N passed to munge()/ldsc():
#   PTSD  -> TOTAL sample size (per-SNP N in the .ma file, ~1,222,882 max).
#            sample.prev = real case proportion = 137,136 / 1,222,882 = 0.1121.
#            (Numbers from the source VCF metadata ##nCase / ##nControl,
#             consistent with Nievergelt et al. 2024 Nat Genet.)
#   MDD   -> SUM-OF-COHORT EFFECTIVE N (829,250; daner .neff). When effective N is
#            supplied, GenomicSEM expects sample.prev = 0.5 (the case/control
#            imbalance is already absorbed into Neff). Do NOT use the raw 0.2181.
#   ANX   -> SUM-OF-COHORT EFFECTIVE N (390,020 = 2 x Neff_half). sample.prev = 0.5.
# Confirmed in this script: munge() receives Neff for MDD/ANX (N_MDD = 829250,
#   N_ANX = 390020) and the total per-SNP N for PTSD, so the sample.prev values
#   below match the N type passed to each trait.
# Population (lifetime, EUR) prevalences are approximate (PTSD ~0.07, MDD ~0.15,
# ANX ~0.10); vary them in a sensitivity analysis. Standardized loadings and the
# genetic-correlation matrix are robust to these liability-scale choices;
# unstandardized (liability-scale) variances will scale with sample.prev.
LDSCoutput <- ldsc(
  traits          = sumstats_files,
  sample.prev     = c(PTSD_SAMPPREV, 0.5, 0.5),
  population.prev = c(0.07,          0.15, 0.10),
  ld              = ldsc_dir,
  wld             = ldsc_dir,
  trait.names     = trait_names
)
saveRDS(LDSCoutput, file.path(out_dir, "LDSCoutput_stage0.rds"))

## ---- 3. genetic correlation matrix + heritability --------------------------
S  <- LDSCoutput$S          # genetic covariance
V  <- LDSCoutput$V          # sampling covariance
# standardize S -> genetic correlation
D  <- sqrt(diag(diag(S)))
rg <- solve(D) %*% S %*% solve(D)
dimnames(rg) <- list(trait_names, trait_names)

cat("\n========== GENETIC CORRELATION MATRIX (rg) ==========\n")
print(round(rg, 3))

cat("\n========== HERITABILITY (diagonal of S, observed/liability) ==========\n")
h2 <- diag(S)
# heritability Z from V (diagonal sampling variances of the h2 estimates)
h2_se <- sqrt(diag(V)[1:length(trait_names)])
h2_Z  <- h2 / h2_se
for (i in seq_along(trait_names))
  cat(sprintf("  %-5s  h2 = %.4f  SE = %.4f  Z = %.2f\n",
              trait_names[i], h2[i], h2_se[i], h2_Z[i]))
cat("  (h2 Z > ~4 is desirable for a stable indicator)\n")

## ---- 4. one-factor (common internalizing) model ----------------------------
cat("\n========== COMMON FACTOR MODEL ==========\n")
cfm <- commonfactor(covstruc = LDSCoutput, estimation = "DWLS")
print(cfm$results)

# Extract standardized loadings for the diagnostic readout
load_tbl <- cfm$results[cfm$results$op == "=~", ]
cat("\n========== STANDARDIZED LOADINGS (key diagnostic) ==========\n")
print(load_tbl[, c("lhs","op","rhs","STD_Genotype","STD_Genotype_SE")])

cat("\n========== STAGE 0 VERDICT GUIDE ==========\n")
cat("GREEN light to proceed to Stage 1 if:\n")
cat("  - model converged and all three loadings are positive;\n")
cat("  - loadings are sizeable (>~0.4) for MDD and ANX;\n")
cat("  - PTSD loading is NOTABLY LOWER than MDD/ANX (=> PTSD-specific variance),\n")
cat("    OR rg(PTSD, others) is clearly < rg(MDD, ANX);\n")
cat("  - h2 Z for all three is healthy (>~4), PTSD especially.\n")
cat("RED / caution if: model fails, a loading is near zero or negative,\n")
cat("  PTSD rg with others is ~1 (no specificity), or PTSD h2 Z is weak.\n")
cat("\nIf borderline: add neuroticism / depressive-symptom GWAS as extra\n")
cat("indicators to stabilise the factor before Stage 1.\n")

cat("\nSaved: LDSCoutput_stage0.rds in", out_dir, "\n")
cat("Stage 0 complete.\n")
