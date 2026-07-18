###############################################################################
# 08b_add_smr_fdr.R
# Adds per-tissue BH-FDR and a graded verdict to the existing SMR results,
# WITHOUT rerunning smr.exe. Reads eqtl_smr_brain_blood.csv (which
# already contains p_SMR / p_HEIDI) and rewrites it with two new columns:
#   smr_fdr       BH-FDR of p_SMR computed SEPARATELY within brain and blood
#   verdict_fdr   graded verdict: FDR-supported / nominal only / HEIDI rejected /
#                 not significant  (reviewer request #8)
###############################################################################
source("00_config.R")

f_in  <- file.path(result_dir, "eqtl_smr_brain_blood.csv")
f_out <- file.path(result_dir, "eqtl_smr_brain_blood_fdr.csv")
es <- read.csv(f_in, stringsAsFactors = FALSE)
es$p_SMR   <- suppressWarnings(as.numeric(es$p_SMR))
es$p_HEIDI <- suppressWarnings(as.numeric(es$p_HEIDI))

# --- per-tissue BH-FDR over genes that returned an SMR result ----------------
es$smr_fdr <- NA_real_
for (tis in unique(es$tissue)) {
  ix <- which(es$tissue == tis & is.finite(es$p_SMR))
  if (length(ix) > 0) es$smr_fdr[ix] <- p.adjust(es$p_SMR[ix], method = "BH")
}

# --- graded verdict ----------------------------------------------------------
grade_fn <- function(p_smr, fdr, p_heidi, nsnp) {
  if (is.na(p_smr) || p_smr >= 0.05) return("Not significant")
  heidi_rej <- !is.na(p_heidi) && !is.na(nsnp) && nsnp >= 3 && p_heidi < 0.05
  if (heidi_rej) return("SMR significant, HEIDI rejected (possible linkage)")
  heidi_ok <- !is.na(p_heidi) && !is.na(nsnp) && nsnp >= 3 && p_heidi >= 0.05
  tier <- if (!is.na(fdr) && fdr < 0.05) "FDR-supported" else "nominal only"
  if (!heidi_ok) return(paste0("SMR significant (", tier, "), HEIDI untestable"))
  paste0("Strong evidence (SMR + HEIDI, ", tier, ")")
}
es$verdict_fdr <- mapply(grade_fn, es$p_SMR, es$smr_fdr, es$p_HEIDI, es$n_heidi)

# --- summary + write ---------------------------------------------------------
cat("\n== SMR per-tissue FDR summary ==\n")
for (tis in unique(es$tissue)) {
  sub <- es[es$tissue == tis & is.finite(es$p_SMR), ]
  cat(sprintf("  %-6s: tested=%d  nominal p<0.05=%d  FDR<0.05=%d\n",
              tis, nrow(sub), sum(sub$p_SMR < 0.05, na.rm=TRUE),
              sum(sub$smr_fdr < 0.05, na.rm=TRUE)))
}
cat("\nPer-gene graded verdicts:\n")
print(es[, c("gene","tissue","p_SMR","smr_fdr","p_HEIDI","verdict_fdr")], row.names = FALSE)

write.csv(es, f_out, row.names = FALSE)
cat("\n[OK] wrote", basename(f_out), "with smr_fdr + verdict_fdr\n")
