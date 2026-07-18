###############################################################################
# 01_MR_screening.R — Proteome-wide MR + FDR + Cross-disease Classification
#
# 1,860 UKB-PPP proteins × 3 outcome GWAS (PTSD, MDD, Anxiety)
# Wald ratio MR with allele harmonization
# FDR correction per outcome (Benjamini-Hochberg)
# Cross-disease classification of PTSD candidates
#
# Outputs:
#   mr_all_outcomes.rds         — full MR result set (3 outcomes × ~1,300 proteins)
#   ptsd_candidates_primary.csv — 10 FDR-significant PTSD candidates (primary screen)
#   mr_summary.csv              — summary statistics
###############################################################################

source("00_config.R")
library(data.table)
library(dplyr)
library(readr)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 01: Proteome-wide MR Screening                ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# ── 1. Read pQTL instruments ──
cat("Step 1: Loading pQTL instruments...\n")
pqtl <- readRDS(pqtl_path)
cat("  Proteins:", length(unique(pqtl$protein)), "\n")
cat("  Total instruments:", nrow(pqtl), "\n\n")

# ── 2. Read outcome GWAS ──
cat("Step 2: Loading outcome GWAS...\n")
outcomes <- list(
  PTSD_freeze3  = list(path = ptsd_path, name = "PTSD_freeze3"),
  MDD           = list(path = mdd_path,  name = "MDD"),
  Anxiety       = list(path = anx_path,  name = "Anxiety")
)

# ── 3. Wald ratio MR for each protein × outcome ──
cat("Step 3: Running Wald ratio MR...\n\n")

# Returns list(status, res). Math is IDENTICAL to the original screen;
# only a harmonization status is recorded so the 1,860 -> analysed flow can be
# tabulated. status in {analyzed_concordant, analyzed_flipped,
# dropped_snp_absent, dropped_allele_incompatible}.
run_wald_mr <- function(pqtl_row, outcome_gwas) {
  snp <- pqtl_row$SNP

  # 在outcome GWAS中查找SNP
  outcome_snp <- outcome_gwas[outcome_gwas$SNP == snp | outcome_gwas$rsid == snp, ]
  if (nrow(outcome_snp) == 0) return(list(status = "dropped_snp_absent", res = NULL))

  beta_exp <- pqtl_row$beta
  se_exp   <- pqtl_row$se
  beta_out <- outcome_snp$beta[1]
  se_out   <- outcome_snp$se[1]

  # Allele harmonization
  status <- "analyzed_concordant"
  if (toupper(pqtl_row$effect_allele) != toupper(outcome_snp$effect_allele[1])) {
    if (toupper(pqtl_row$effect_allele) == toupper(outcome_snp$other_allele[1])) {
      beta_out <- -beta_out
      status   <- "analyzed_flipped"
    } else {
      return(list(status = "dropped_allele_incompatible", res = NULL))  # Incompatible alleles
    }
  }

  # Wald ratio
  mr_beta <- beta_out / beta_exp
  mr_se   <- abs(se_out / beta_exp)
  mr_pval <- 2 * pnorm(-abs(mr_beta / mr_se))

  res <- data.frame(
    protein = pqtl_row$protein,
    SNP = snp,
    mr_beta = mr_beta, mr_se = mr_se, mr_pval = mr_pval,
    OR = exp(mr_beta),
    OR_lower = exp(mr_beta - 1.96 * mr_se),
    OR_upper = exp(mr_beta + 1.96 * mr_se),
    stringsAsFactors = FALSE
  )
  list(status = status, res = res)
}

# 对每个outcome执行MR (此处为框架代码, 实际运行需读取GWAS文件)
all_results <- list()
harm_list   <- list()   # harmonization log: protein x outcome x status
for (out_name in names(outcomes)) {
  cat("  Processing:", out_name, "...\n")
  out_path <- outcomes[[out_name]]$path

  # 读取outcome GWAS
  outcome_gwas <- fread(out_path)

  # 对每个蛋白执行Wald ratio MR
  proteins <- unique(pqtl$protein)
  out_each <- lapply(proteins, function(prot) {
    inst <- pqtl[pqtl$protein == prot, ][1, ]  # lead instrument
    o <- tryCatch(run_wald_mr(inst, outcome_gwas),
                  error = function(e) list(status = "dropped_error", res = NULL))
    res <- o$res
    if (!is.null(res)) res$outcome <- out_name
    list(res = res, status = o$status, protein = prot)
  })

  res_list <- lapply(out_each, `[[`, "res")
  all_results[[out_name]] <- bind_rows(res_list)
  harm_list[[out_name]] <- data.frame(
    protein = vapply(out_each, `[[`, character(1), "protein"),
    outcome = out_name,
    status  = vapply(out_each, `[[`, character(1), "status"),
    stringsAsFactors = FALSE
  )
  cat("    Matched:", nrow(all_results[[out_name]]), "proteins\n")
}

# 合并所有结果
mr_all <- bind_rows(all_results)

# ── 3b. Harmonization flow table (1,860 -> analysed) ───────────────────────
# Faithful to the screen's actual logic: direct rsID lookup + allele flip,
# with no proxy / position-matching / palindromic-EAF recovery. Losses are
# therefore either (i) SNP absent from the outcome GWAS or (ii) allele
# incompatible after attempting an effect/other-allele flip.
harm_log <- bind_rows(harm_list)
harm_flow <- harm_log %>%
  group_by(outcome) %>%
  summarise(
    eligible_instruments      = dplyr::n(),
    dropped_snp_absent        = sum(status == "dropped_snp_absent"),
    dropped_allele_incompat   = sum(status == "dropped_allele_incompatible"),
    analysed_total            = sum(grepl("^analyzed", status)),
    of_which_allele_concordant= sum(status == "analyzed_concordant"),
    of_which_allele_flipped   = sum(status == "analyzed_flipped"),
    .groups = "drop"
  )
write_csv(harm_log,  file.path(result_dir, "harmonization_log.csv"))
write_csv(harm_flow, file.path(result_dir, "harmonization_flow.csv"))
cat("\n  Harmonization flow (per outcome):\n")
print(as.data.frame(harm_flow))
cat("  ✓ harmonization_flow.csv / harmonization_log.csv\n\n")

# ── 4. FDR correction per outcome ──
cat("\nStep 4: FDR correction...\n")
mr_all <- mr_all %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(mr_pval, method = "BH")) %>%
  ungroup()

# 统计FDR显著数
for (out in unique(mr_all$outcome)) {
  n_sig <- sum(mr_all$fdr[mr_all$outcome == out] < fdr_threshold)
  cat("  ", out, ": ", n_sig, " FDR < 0.05\n")
}

# ── 5. Cross-disease classification ──
cat("\nStep 5: Cross-disease classification...\n")
ptsd_sig <- mr_all %>%
  filter(outcome == "PTSD_freeze3", fdr < fdr_threshold)

if (nrow(ptsd_sig) > 0) {
  ptsd_sig <- ptsd_sig %>%
    left_join(
      mr_all %>% filter(outcome == "MDD") %>% select(protein, pval_mdd = mr_pval),
      by = "protein"
    ) %>%
    left_join(
      mr_all %>% filter(outcome == "Anxiety") %>% select(protein, pval_anx = mr_pval),
      by = "protein"
    ) %>%
    mutate(
      # Descriptive nominal-overlap flags only (not a formal classification;
      # cross-disorder inference is based on Genomic SEM, not these flags).
      nominal_mdd = pval_mdd < nominal_p,
      nominal_anx = pval_anx < nominal_p
    )
}

# ── 6. 保存 ──
cat("\nStep 6: Saving results...\n")
saveRDS(mr_all, file.path(result_dir, "mr_all_outcomes.rds"))
write_csv(ptsd_sig, file.path(result_dir, "ptsd_candidates_primary.csv"))
write_csv(
  mr_all %>% group_by(outcome) %>%
    summarise(n_tested = n(), n_fdr = sum(fdr < 0.05), .groups = "drop"),
  file.path(result_dir, "mr_summary.csv")
)

cat("  ✓ mr_all_outcomes.rds\n")
cat("  ✓ ptsd_candidates_primary.csv\n")
cat("  ✓ mr_summary.csv\n")
cat("\n  PTSD FDR candidates:", nrow(ptsd_sig), "\n")
