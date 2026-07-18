###############################################################################
# 28_genomicsem_stage12_factorMR_userGWAS.R   (corrected)
#
# Fixes vs the previous 28:
#  (1) MODEL CONSISTENCY: uses userGWAS() with the SAME explicit constrained
#      model as the trait-level report (MDD residual fixed to 0; stage0b fit2),
#      instead of the default commonfactorGWAS() which free-estimates the MDD
#      residual. commonfactorGWAS() cannot take a user model, so the previous
#      run silently used a different (default) model than the manuscript reports.
#  (2) PROXY REMOVED: the earlier proxy swap paired a proxy's factor effect with
#      the sentinel's pQTL beta (an invalid Wald ratio), because the proxies are
#      absent from the cis-pQTL table. Proxies are therefore NOT used; SIRPA,
#      FURIN, UBE2L6 are not testable (no single SNP has BOTH a cis-pQTL and a
#      common-factor estimate), and CD101 has no adequate proxy. Only direct
#      instruments (variant present in BOTH pQTL and the factor sumstats) are used.
#  (3) SE COLUMN FIX: factor-effect SE is taken from userGWAS 'SE' (the previous
#      code read a non-existent 'se' column).
#
# Outputs:
#   stage12_factor_snpeffects_userGWAS.rds   common-factor effect + Q_SNP
#   stage12_factorMR_results.csv             protein -> common factor MR (Wald)
#   stage12_not_testable.csv                 candidates without a valid instrument
###############################################################################

suppressMessages({ library(GenomicSEM); library(data.table) })

out_dir <- "D:/PTSD/results/genomicsem"
setwd(out_dir)

LDSCoutput <- readRDS(file.path(out_dir, "LDSCoutput_stage0.rds"))
ss   <- readRDS(file.path(out_dir, "stage1_sumstats.rds"))         # from 27 BLOCK 1
pqtl <- as.data.table(readRDS("D:/PTSD/pQTL/ukbppp_cis_pqtl.rds"))

## ---- instrument set: DIRECT instruments only (present in both pQTL and ss) ---
# No proxy substitution: a valid Wald ratio needs ONE variant carrying both a
# cis-pQTL effect and a common-factor effect. The 3 candidate proxies are absent
# from the pQTL table, so they cannot form a valid instrument.
inst_direct <- intersect(unique(pqtl$SNP), ss$SNP)
ss_inst     <- ss[ss$SNP %in% inst_direct, ]
cat("Direct instrument SNPs to run:", nrow(ss_inst), "\n")

## ---- explicit constrained model (== trait-level stage0b fit2) ---------------
# MDD residual FIXED to 0 (Heywood boundary in the unconstrained model).
model <- '
  F1 =~ NA*PTSD + MDD + ANX
  F1 ~~ 1*F1
  PTSD ~~ PTSD
  ANX  ~~ ANX
  MDD  ~~ 0*MDD
  F1 ~ SNP
'

## ---- STAGE 1 (restricted): userGWAS with the explicit model ----------------
CKPT <- file.path(out_dir, "stage12_factor_snpeffects_userGWAS.rds")
if (!file.exists(CKPT)) {
  cat("\n=== userGWAS on instrument SNPs (explicit model; serial) ===\n")
  t0 <- Sys.time()
  gwas_list <- userGWAS(
    covstruc   = LDSCoutput,
    SNPs       = ss_inst,
    estimation = "DWLS",
    model      = model,
    sub        = c("F1~SNP"),   # return only the SNP -> factor effect
    Q_SNP      = TRUE,          # compute Q_SNP heterogeneity under THIS model
    cores      = 1,
    parallel   = FALSE,
    printwarn  = TRUE
  )
  cat("Elapsed:", round(as.numeric(Sys.time()-t0, units="mins"),1), "min\n")
  saveRDS(gwas_list, CKPT)
} else {
  cat("userGWAS checkpoint exists; loading.\n")
  gwas_list <- readRDS(CKPT)
}

# userGWAS(sub=, Q_SNP=TRUE) returns a LIST; take [[1]]. With Q_SNP=TRUE the
# output carries the model-consistent Q_SNP heterogeneity index and its P value.
fsnp <- as.data.table(gwas_list[[1]])
cat("\nFactor-SNP-effect columns:\n"); print(names(fsnp))

# est/SE: the SNP->factor effect and its SE (from THIS constrained-model run).
est_col <- if ("est" %in% names(fsnp)) "est" else grep("^est",   names(fsnp), value=TRUE, ignore.case=TRUE)[1]
se_col  <- if ("SE"  %in% names(fsnp)) "SE"  else grep("^SE$|^se$",names(fsnp), value=TRUE, ignore.case=TRUE)[1]
# Q_SNP + P: auto-detect the column names produced by this GenomicSEM version.
q_col   <- grep("^Q_?SNP$|^Q$",               names(fsnp), value=TRUE, ignore.case=TRUE)[1]
qp_col  <- grep("Q_?SNP_?pval|Q_?SNP_?P$|^Q_pval$|Qpval", names(fsnp), value=TRUE, ignore.case=TRUE)[1]
stopifnot(!is.na(est_col), !is.na(se_col))
if (is.na(q_col) || is.na(qp_col)) {
  cat("WARNING: Q_SNP columns not auto-detected. Columns present:\n")
  print(names(fsnp))
  stop("Set q_col/qp_col manually from the printed names above (need Q_SNP + its P).")
}
cat(sprintf("Using columns: est='%s' SE='%s' Q_SNP='%s' Q_SNP_pval='%s'\n",
            est_col, se_col, q_col, qp_col))

fkey <- fsnp[, .(SNP,
                 beta_factor = get(est_col),
                 se_factor   = get(se_col),
                 Q      = get(q_col),
                 Q_pval = get(qp_col))]

# clean character alleles from the harmonized sumstats (A1 = effect allele in ss)
ss_al <- as.data.table(ss)[, .(SNP, ssA1 = toupper(A1), ssA2 = toupper(A2))]
fkey  <- merge(fkey, ss_al, by = "SNP")

## ---- STAGE 2: protein -> common factor MR (Wald ratio) ---------------------
stopifnot(uniqueN(pqtl$protein) == nrow(pqtl))        # exactly one sentinel per protein
pqtl_use <- copy(pqtl)
pqtl_use[, outcome_snp := SNP]                        # NO proxy swap

res <- merge(pqtl_use, fkey, by.x = "outcome_snp", by.y = "SNP",
             all.x = FALSE, suffixes = c("_pqtl","_factor"))

# allele harmonization: align factor effect to the pQTL effect allele,
# using the CLEAN alleles (ssA1 = effect allele in the factor sumstats).
res[, aligned := fifelse(toupper(effect_allele) == ssA1,  1,
                  fifelse(toupper(effect_allele) == ssA2, -1, NA_real_))]

## --- diagnostics (non-destructive): full allele-pair check + palindromic flag -
res[, `:=`(EA = toupper(effect_allele), OA = toupper(other_allele))]
res[, direct  := (EA == ssA1 & OA == ssA2)]
res[, reverse := (EA == ssA2 & OA == ssA1)]
res[, pair_ok := direct | reverse]
res[, palindromic := (EA=="A"&OA=="T")|(EA=="T"&OA=="A")|(EA=="C"&OA=="G")|(EA=="G"&OA=="C")]
cat(sprintf("Allele check: %d aligned by EA; of these %d fail the full-pair (EA+OA) test; %d are palindromic (A/T or C/G).\n",
            sum(!is.na(res$aligned)),
            sum(!is.na(res$aligned) & !res$pair_ok, na.rm=TRUE),
            sum(!is.na(res$aligned) & res$palindromic, na.rm=TRUE)))
cat("  (Alignment below uses EA-identity, which reproduced the validated set;\n")
cat("   inspect the counts above to decide whether stricter full-pair/EAF handling changes anything.)\n")

res <- res[!is.na(aligned)]
res[, beta_factor_aligned := beta_factor * aligned]

# Wald ratio MR
res[, mr_beta := beta_factor_aligned / beta]
res[, mr_se   := abs(se_factor / beta)]
res[, mr_z    := mr_beta / mr_se]
res[, mr_p    := 2 * pnorm(-abs(mr_z))]

out <- res[, .(protein, cis_gene, outcome_snp,
               A1 = ssA1, A2 = ssA2,
               pqtl_beta = beta, pqtl_se = se,
               beta_factor = beta_factor_aligned, se_factor,
               mr_beta, mr_se, mr_z, mr_p, Q, Q_pval)]
setorder(out, mr_p)
fwrite(out, file.path(out_dir, "stage12_factorMR_results.csv"))

## ---- not-testable candidates (documented explicitly) -----------------------
not_testable <- data.table(
  protein = c("SIRPA","FURIN","UBE2L6","CD101"),
  sentinel = c("rs6136377","rs2071410","rs28362950","rs12130298"),
  proxy    = c("rs17775933","rs1573643","rs11603020", NA),
  proxy_r2 = c(0.988, 0.951, 0.979, 0.51),
  reason = c(
    rep("Sentinel cis-pQTL absent from common-factor sumstats and proxy absent from cis-pQTL data; no single variant carries both effects.", 3),
    "No adequate high-LD proxy (best r2 = 0.51); no valid Wald-ratio instrument."))
fwrite(not_testable, file.path(out_dir, "stage12_not_testable.csv"))

## ---- summary ----------------------------------------------------------------
cat("\nProtein -> common factor MR done. Testable proteins:", nrow(out), "\n")
bh <- p.adjust(out$mr_p, "BH")
cat("Factor-MR BH-FDR < 0.05:", sum(bh < 0.05), "\n")
cat("Nominal Q_SNP < 0.05:", sum(out$Q_pval < 0.05, na.rm=TRUE), "\n")
cat("Q_SNP BH-FDR < 0.05:", sum(p.adjust(out$Q_pval,"BH") < 0.05, na.rm=TRUE), "\n")

cands <- c("AKT3","CD40","CGREF1","FES","KHK","SNX18")   # the 6 testable candidates
out[, bh := bh]
cat("\n========== SIX TESTABLE CANDIDATES ==========\n")
print(out[cis_gene %in% cands | protein %in% cands,
          .(cis_gene, outcome_snp, mr_beta, mr_p, bh, Q_pval)])
cat("\nNot testable (see stage12_not_testable.csv): SIRPA, FURIN, UBE2L6, CD101\n")
