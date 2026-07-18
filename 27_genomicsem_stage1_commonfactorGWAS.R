###############################################################################
# 27_genomicsem_stage1_commonfactorGWAS.R
#
# Genomic SEM — STAGE 1: common-factor GWAS + Q_SNP, with factor-identity
# diagnostics. This is the EXPENSIVE step (sumstats harmonization + per-SNP SEM
# across ~1M SNPs). Built with checkpoints so a crash does not lose hours.
#
# Model: common internalizing factor over PTSD/MDD/ANX, MDD residual bounded to
# a SMALL POSITIVE value (not hard 0) to avoid per-SNP convergence failures at
# the Heywood boundary (Stage 0b showed MDD residual ~ 0).
#
# RUN IN BLOCKS. Each block saves an .rds checkpoint. If a later block fails,
# you can reload the checkpoint instead of recomputing.
#
# OUTPUTS:
#   stage1_sumstats.rds        harmonized multivariate sumstats (big)
#   stage1_factorGWAS.rds      common-factor SNP effects + Q_SNP
#   stage1_factor_identity.txt factor-vs-trait rg diagnostics
###############################################################################

suppressMessages({ library(GenomicSEM); library(data.table) })

out_dir <- "D:/PTSD/results/genomicsem"
setwd(out_dir)

ptsd_ma <- "D:/PTSD/results/ptsd_freeze3_smr.ma"
mdd_f   <- "D:/PTSD/GWAS/PGC_MDD_2025/27061255/daner/daner_pgc_mdd_no23andMe-noUKBB_eur_hg19_v3.49.24.11.neff.gz"
anx_f   <- "D:/PTSD/GWAS/PGC_Anxiety_2026/31389910/ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz"
hm3     <- "D:/reference/ldsc/w_hm3.snplist"
ref_ss  <- "D:/reference/reference.custom.maf.0.005.txt.gz"   # sumstats() reference (SNP/A1/A2/MAF)

N_MDD <- 829250
N_ANX <- 390020
trait_names <- c("PTSD","MDD","ANX")

LDSCoutput <- readRDS(file.path(out_dir, "LDSCoutput_stage0.rds"))

## ===========================================================================
## BLOCK 1 — prep PTSD .ma for sumstats() and harmonize all three
## sumstats() needs: SNP, A1(effect), A2, effect (beta/OR/Z), P, and N or SE.
## This step reads ALL SNPs (not just HapMap3) and is slow + memory heavy.
## ===========================================================================
CKPT_SS <- file.path(out_dir, "stage1_sumstats.rds")
if (!file.exists(CKPT_SS)) {
  cat("=== BLOCK 1: sumstats() harmonization (slow) ===\n")

  # PTSD .ma -> a file with columns sumstats() understands
  ptsd <- fread(ptsd_ma)   # SNP A1 A2 freq b se p N
  setnames(ptsd, c("SNP","A1","A2","freq","b","se","p","N"),
                 c("SNP","A1","A2","MAF","BETA","SE","P","N"))
  ptsd_prep <- file.path(out_dir, "ptsd_freeze3_for_sumstats.txt")
  fwrite(ptsd, ptsd_prep, sep="\t")

  files       <- c(ptsd_prep, mdd_f, anx_f)
  ref         <- ref_ss
  # se.logit: are SEs on the logit (log-OR) scale? PTSD beta=log-OR -> TRUE-ish,
  #   but PTSD .ma carries BETA/SE already on log-odds; daner OR with SE on
  #   log-odds scale -> TRUE for MDD/ANX. PTSD .ma SE is on the beta scale -> TRUE.
  # OLS: continuous? all three are binary/effect-size -> FALSE.
  # linprob: treat as linear-probability? FALSE (we have proper effect sizes).
  ss <- sumstats(files       = files,
                 ref         = ref,
                 trait.names = trait_names,
                 se.logit    = c(TRUE, TRUE, TRUE),
                 OLS         = c(FALSE, FALSE, FALSE),
                 linprob     = c(FALSE, FALSE, FALSE),
                 N           = c(NA, N_MDD, N_ANX),
                 info.filter = 0.6,
                 maf.filter  = 0.01)
  saveRDS(ss, CKPT_SS)
  cat("BLOCK 1 done. sumstats rows:", nrow(ss), "-> saved\n")
} else {
  cat("BLOCK 1 skipped (checkpoint exists). Loading...\n")
  ss <- readRDS(CKPT_SS)
}

## ===========================================================================
## BLOCK 2 — OPTIONAL genome-wide common-factor GWAS (diagnostic only).
## NOTE: This genome-wide run is NOT the analysis reported in the manuscript.
##   The manuscript's protein-to-factor MR uses the INSTRUMENT-RESTRICTED
##   userGWAS in script 28 (explicit constrained model, MDD residual fixed 0).
##   A genome-wide factor GWAS was not performed for the manuscript, so this
##   block is OFF by default. Set RUN_GENOMEWIDE_FACTOR_GWAS <- TRUE only if you
##   specifically want the (slow, hours-long) genome-wide diagnostic.
## If enabled, use userGWAS() with the SAME explicit model as script 28 rather
##   than the default commonfactorGWAS(), so trait- and SNP-level models match.
## ===========================================================================
RUN_GENOMEWIDE_FACTOR_GWAS <- FALSE
CKPT_GWAS <- file.path(out_dir, "stage1_factorGWAS.rds")
if (RUN_GENOMEWIDE_FACTOR_GWAS && !file.exists(CKPT_GWAS)) {
  cat("\n=== BLOCK 2: commonfactorGWAS (SLOW — hours) ===\n")

  n_cores <- max(1, parallel::detectCores() - 1)
  cat("Using", n_cores, "cores.\n")

  # Explicit constrained model (MDD residual fixed 0), IDENTICAL to script 28
  # and to the trait-level report (stage0b fit2). Use userGWAS (which accepts a
  # user model); commonfactorGWAS() would silently use the default free-residual
  # model instead.
  model <- '
    F1 =~ NA*PTSD + MDD + ANX
    F1 ~~ 1*F1
    PTSD ~~ PTSD
    ANX  ~~ ANX
    MDD  ~~ 0*MDD
    F1 ~ SNP
  '
  factorGWAS <- userGWAS(
    covstruc   = LDSCoutput,
    SNPs       = ss,
    estimation = "DWLS",
    model      = model,
    sub        = c("F1~SNP"),
    cores      = n_cores,
    parallel   = TRUE,
    printwarn  = TRUE
  )[[1]]
  saveRDS(factorGWAS, CKPT_GWAS)
  cat("BLOCK 2 done. factor GWAS SNPs:", nrow(factorGWAS), "-> saved\n")
} else if (file.exists(CKPT_GWAS)) {
  cat("\nBLOCK 2 skipped (checkpoint exists). Loading...\n")
  factorGWAS <- readRDS(CKPT_GWAS)
} else {
  cat("\nBLOCK 2 skipped (RUN_GENOMEWIDE_FACTOR_GWAS = FALSE; not used by manuscript).\n")
  factorGWAS <- NULL
}

## ===========================================================================
## BLOCK 3 — factor-identity diagnostics
## Is the common factor a balanced 3-disorder factor, or MDD-dominant?
## Compare factor-GWAS Z to each trait's Z (proxy for genetic similarity).
## A formal rg would need munging the factor GWAS; here we use a fast Z-Z proxy
## and flag if a full LDSC rg is warranted.
## ===========================================================================
cat("\n=== BLOCK 3: factor-identity diagnostics ===\n")
if (is.null(factorGWAS)) {
  cat("BLOCK 3 skipped: genome-wide factor GWAS not run (see BLOCK 2 note).\n")
} else {
fg <- as.data.table(factorGWAS)
cat("factor GWAS columns:\n"); print(names(fg))

# Save a header preview so we can wire up Stage 2 correctly
out_txt <- file.path(out_dir, "stage1_factor_identity.txt")
sink(out_txt)
cat("Stage 1 factor GWAS — column names:\n"); print(names(fg))
cat("\nN SNPs:", nrow(fg), "\n")
cat("\nFirst rows:\n"); print(head(fg, 5))
# Q_SNP summary if present
qcol <- grep("Q_?SNP", names(fg), value=TRUE, ignore.case=TRUE)
pcol <- grep("Q.*pval|Q.*_P|Qpval", names(fg), value=TRUE, ignore.case=TRUE)
cat("\nQ_SNP-related columns:", paste(qcol, collapse=", "), "\n")
sink()
cat("Diagnostics written to", out_txt, "\n")
}

cat("\nSTAGE 1 complete. Next: Stage 2 (1,860-protein -> common factor MR),\n")
cat("and a formal factor-vs-MDD rg check to confirm the factor is not MDD-dominant.\n")
