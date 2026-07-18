###############################################################################
# 05b_forward_multiIV_MR.R
# Forward multi-instrument MR (protein -> PTSD) using the SuSiE credible-set
# lead SNPs as instruments. Produces IVW, weighted median, and MR-Egger
# (with Egger intercept) plus Cochran's Q, via the MendelianRandomization
# package (standard implementations, not hand-coded).
#
# WHY THIS SCRIPT EXISTS
#   - Reviewer request: report full IVW / weighted-median / MR-Egger effects,
#     SE, P, Cochran Q and the Egger intercept for the multi-instrument MR.
#   - It also fills a pipeline gap: 06_biological_annotation.R reads
#     `multi_instrument_mr.csv`, which no other script produced.
#
# INSTRUMENTS
#   get_lead_snps() takes the highest-PIP SNP of EACH SuSiE credible set
#   (script 04/05 pipeline; L=10). Different credible sets are, by SuSiE
#   construction, approximately independent signals, so the number of
#   instruments equals the number of credible sets for that protein (it is
#   NOT fixed at five). Residual between-set LD is not modelled with a full
#   SNP correlation matrix; IVW/Egger SEs are therefore approximate and the
#   multi-instrument MR is reported as a supportive, not a primary, analysis.
#
# OUTPUT
#   multi_instrument_mr.csv  (long format: one row per protein x method)
#     columns: protein, n_iv, method, beta, se, pval, Q, Q_pval,
#              egger_intercept, egger_intercept_se, egger_intercept_p
###############################################################################

suppressMessages({
  library(dplyr); library(readr); library(data.table)
  library(MendelianRandomization)
})

source("00_config.R")

# --- shared helper: lead SNP (highest PIP) of each SuSiE credible set --------
get_lead_snps <- function(protein, cs_dir) {
  cs_files <- list.files(cs_dir, pattern = paste0("^", protein, "_pqtl_CS"),
                         full.names = TRUE)
  if (length(cs_files) == 0) return(data.table())
  rbindlist(lapply(cs_files, function(f) { cs <- fread(f); cs[order(-pip)][1] }))
}

coloc_full <- readRDS(file.path(result_dir, "susie_coloc_full.rds"))

run_forward_multiIV <- function() {
  cat("\n== Forward multi-instrument MR (protein -> PTSD) ==\n")
  out <- list()

  for (protein in candidates) {
    cat(sprintf("== %s ==\n", protein))

    cf <- coloc_full[[protein]]
    if (is.null(cf) || is.null(cf$merged_common)) {
      cat("  no merged_common; skipped\n")
      out[[protein]] <- tibble(protein = protein, n_iv = 0, method = NA,
                               beta = NA, se = NA, pval = NA, Q = NA, Q_pval = NA,
                               egger_intercept = NA, egger_intercept_se = NA,
                               egger_intercept_p = NA, note = "no region data")
      next
    }

    leads  <- get_lead_snps(protein, cs_dir)
    merged <- as.data.frame(cf$merged_common)
    iv <- merged %>% filter(rsid %in% leads$rsid) %>%
      select(rsid, beta_pqtl, se_pqtl, beta_gwas, se_gwas) %>%
      filter(is.finite(beta_pqtl), is.finite(se_pqtl),
             is.finite(beta_gwas), is.finite(se_gwas), se_pqtl > 0, se_gwas > 0)

    n_iv <- nrow(iv)
    cat(sprintf("  instruments (credible-set leads): %d\n", n_iv))

    # MR-Egger needs >= 3 instruments; IVW/WM need >= 2.
    if (n_iv < 3) {
      cat("  <3 instruments; Egger not estimable\n")
      out[[protein]] <- tibble(protein = protein, n_iv = n_iv, method = NA,
                               beta = NA, se = NA, pval = NA, Q = NA, Q_pval = NA,
                               egger_intercept = NA, egger_intercept_se = NA,
                               egger_intercept_p = NA,
                               note = "fewer than 3 credible-set instruments")
      next
    }

    mri <- mr_input(bx = iv$beta_pqtl, bxse = iv$se_pqtl,
                    by = iv$beta_gwas, byse = iv$se_gwas, snps = iv$rsid)

    ivw   <- MendelianRandomization::mr_ivw(mri)
    wm    <- MendelianRandomization::mr_median(mri, weighting = "weighted")
    egg   <- MendelianRandomization::mr_egger(mri)

    rows <- bind_rows(
      tibble(protein = protein, n_iv = n_iv, method = "IVW",
             beta = ivw@Estimate, se = ivw@StdError, pval = ivw@Pvalue,
             Q = ivw@Heter.Stat[1], Q_pval = ivw@Heter.Stat[2],
             egger_intercept = NA, egger_intercept_se = NA, egger_intercept_p = NA,
             note = ""),
      tibble(protein = protein, n_iv = n_iv, method = "Weighted median",
             beta = wm@Estimate, se = wm@StdError, pval = wm@Pvalue,
             Q = NA, Q_pval = NA,
             egger_intercept = NA, egger_intercept_se = NA, egger_intercept_p = NA,
             note = ""),
      tibble(protein = protein, n_iv = n_iv, method = "MR-Egger",
             beta = egg@Estimate, se = egg@StdError.Est, pval = egg@Pvalue.Est,
             Q = egg@Heter.Stat[1], Q_pval = egg@Heter.Stat[2],
             egger_intercept = egg@Intercept,
             egger_intercept_se = egg@StdError.Int,
             egger_intercept_p = egg@Pvalue.Int, note = "")
    )
    out[[protein]] <- rows

    cat(sprintf("  IVW P=%.3g | WM P=%.3g | Egger P=%.3g | Egger intercept P=%.3g\n\n",
                ivw@Pvalue, wm@Pvalue, egg@Pvalue.Est, egg@Pvalue.Int))
  }

  res <- bind_rows(out)
  write_csv(res, file.path(result_dir, "multi_instrument_mr.csv"))
  cat("[OK] multi_instrument_mr.csv (IVW / Weighted median / MR-Egger + intercept + Q)\n")
  invisible(res)
}

run_forward_multiIV()
