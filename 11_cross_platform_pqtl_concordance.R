###############################################################################
# 11_cross_platform_pqtl_concordance.R — Cross-platform cis-pQTL concordance (SIRPA)
#
# Background: Of the 10 candidate proteins, only SIRPA has a cis-pQTL in the
#   deCODE SomaScan plasma proteomics resource. CD40 / AKT3 / FES / KHK /
#   SNX18 / CD101 are not on the SomaScan panel; FURIN / CGREF1 / UBE2L6 have
#   SomaScan assays but only trans-pQTLs.
#
# Reported as CROSS-PLATFORM CONCORDANCE OF THE cis-pQTL EFFECT (not independent MR replication)
# (UKB-PPP Olink vs deCODE SomaScan), NOT as an independent deCODE MR estimate.
#
# Why no separate deCODE MR / meta-analysis:
#   A deCODE-specific Wald-ratio MR cannot be computed here without a deCODE
#   instrument run against the SAME PTSD outcome GWAS. Reconstructing a deCODE
#   MR by rescaling the UKB-PPP MR estimate (beta_outcome = ukb_beta_mr *
#   ukb_beta_pqtl, then dividing by decode_beta_pqtl) is algebraically circular:
#   the pQTL effect cancels in the Wald-ratio test statistic, so the deCODE
#   z-statistic (and P) is identical to UKB-PPP and carries no independent
#   information. Meta-analysing the two as if independent would understate the
#   standard error. Both were therefore removed; only the UKB-PPP MR estimate
#   is reported as primary, alongside cross-platform pQTL concordance.
#
# Data: Ferkingstad et al. 2021 Nat Genet, Supplementary Table S2 (sheet ST02)
#       SIRPA cis-pQTL: rs6136377, beta = -1.26, -log10(P) = 6920.5
#       VERIFIED: allele orientation checked against the deCODE variant-annotation
#       file (assocvariants.annotated.txt.gz). rs6136377 chr20:1915642,
#       effectAllele = G, otherAllele = A, EAF = 0.360 - identical orientation to
#       UKB-PPP (G/A, EAF 0.382). The two platforms are directly comparable and
#       concordant in sign and magnitude (UKB-PPP beta = -1.234; deCODE beta = -1.26).
#
# Output: decode_cross_platform_concordance.csv
###############################################################################
source("00_config.R")

library(dplyr)
library(readr)

mr_rds_path <- file.path(result_dir, "mr_all_outcomes.rds")

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 11: Cross-platform cis-pQTL concordance (SIRPA) ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# UKB-PPP SIRPA cis-pQTL instrument (Olink)
pqtl_inst  <- readRDS(pqtl_path)
sirpa_inst <- pqtl_inst[pqtl_inst$SNP == "rs6136377", ]
ukb_beta_pqtl <- sirpa_inst$beta[1]   # -1.234
ukb_se_pqtl   <- sirpa_inst$se[1]

# UKB-PPP primary MR result (SIRPA -> PTSD), reported as-is (primary estimate)
mr_orig  <- readRDS(mr_rds_path)
sirpa_mr <- mr_orig[grepl("SIRPA", mr_orig$protein, ignore.case = TRUE) &
                    grepl("PTSD|ptsd", mr_orig$outcome, ignore.case = TRUE), ]
ukb_beta_mr <- sirpa_mr$mr_beta[1]
ukb_se_mr   <- sirpa_mr$mr_se[1]
ukb_p_mr    <- 2 * pnorm(-abs(ukb_beta_mr / ukb_se_mr))

# deCODE SIRPA cis-pQTL (SomaScan; Ferkingstad et al. 2021, Supp Table S2).
# Magnitude from the published deCODE summary statistics; allele orientation and
# EAF verified against the deCODE annotation file to match UKB-PPP (see header).
decode_beta_pqtl <- -1.26

# Cross-platform concordance of the cis-pQTL effect (direction + magnitude)
direction_concordant <- sign(ukb_beta_pqtl) == sign(decode_beta_pqtl)

cat(sprintf("  cis-pQTL rs6136377  UKB-PPP beta = %.3f | deCODE beta = %.3f | direction concordant = %s\n",
            ukb_beta_pqtl, decode_beta_pqtl, direction_concordant))
cat(sprintf("  UKB-PPP MR (primary):  OR = %.4f,  P = %s\n",
            exp(ukb_beta_mr), signif(ukb_p_mr, 3)))

res <- data.frame(
  protein                   = "SIRPA",
  snp                       = "rs6136377",
  pqtl_beta_ukb             = ukb_beta_pqtl,
  pqtl_beta_decode          = decode_beta_pqtl,
  pqtl_direction_concordant = direction_concordant,
  MR_OR_ukb                 = round(exp(ukb_beta_mr), 4),
  MR_P_ukb                  = signif(ukb_p_mr, 4),
  note = paste0("Cross-platform cis-pQTL concordance only. deCODE-specific MR ",
                "and fixed-effect meta-analysis intentionally not computed: a ",
                "reconstructed deCODE Wald ratio is algebraically non-independent ",
                "of the UKB-PPP estimate (pQTL effect cancels) and provides no ",
                "additional evidence.")
)
write_csv(res, file.path(result_dir, "decode_cross_platform_concordance.csv"))
cat("\n✓ decode_cross_platform_concordance.csv (pQTL concordance; deCODE MR/meta removed)\n")
