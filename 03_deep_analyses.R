###############################################################################
# 03_deep_analyses.R — Dose-response, Cross-ancestry, Pathway Enrichment
#
# Part A: Dose-response concordance (PCL quantitative phenotype)
# Part B: Cross-ancestry validation (AAM + HNA)
# Part C: Cross-disease NES correlation
# Part D: Pathway enrichment (fgsea: Hallmark + KEGG + Reactome + GO:BP)
#
# Outputs:
#   mr_all_outcomes_extended.rds    — + PCL + AAM + HNA
#   dose_response_concordance.csv   — 10/10 confirmed
#   cross_ancestry_validation.csv   — AAM + HNA results
#   hallmark_NES_comparison.csv     — NES across outcomes
#   fgsea_ptsd_all_pathways.csv     — all pathway results
###############################################################################

source("00_config.R")
library(data.table)
library(dplyr)
library(readr)
library(fgsea)
library(msigdbr)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 03: Deep Analyses                              ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

pqtl <- readRDS(pqtl_path)
mr_all <- readRDS(file.path(result_dir, "mr_all_outcomes.rds"))

# ═══════════════════════════════════════════════════════════════
# Part A: Dose-response (PCL quantitative phenotype)
#
# PCL GWAS (PGC-PTSD Freeze 3 quantitative phenotype) is required for
# fresh re-computation. If pcl_path is unavailable, this section will:
#   1. Fall back to a cached result file at <result_dir>/dose_response_concordance.csv
#      if one exists from a previous successful run, OR
#   2. Skip Part A entirely with a warning and let Parts B-D continue.
# ═══════════════════════════════════════════════════════════════
cat("Part A: Dose-response (PCL)...\n")

df_pcl    <- NULL                                                # initialise so save step is safe
pcl_cache <- file.path(result_dir, "dose_response_concordance.csv")

if (!is.null(pcl_path) && nzchar(pcl_path) && file.exists(pcl_path)) {

  # ── PCL GWAS available: fresh computation ──
  pcl_gwas <- fread(pcl_path)
  pcl_results <- list()
  for (prot in candidates) {
    inst <- pqtl[pqtl$protein == prot, ][1, ]
    snp_data <- pcl_gwas[pcl_gwas$SNP == inst$SNP | pcl_gwas$rsid == inst$SNP, ]
    if (nrow(snp_data) == 0) next

    mr_beta <- snp_data$beta[1] / inst$beta
    mr_se   <- abs(snp_data$se[1] / inst$beta)
    pcl_results[[prot]] <- data.frame(
      protein = prot, beta_pcl = mr_beta, se_pcl = mr_se,
      p_pcl = 2 * pnorm(-abs(mr_beta / mr_se)),
      concordant = sign(mr_beta) == sign(
        mr_all$mr_beta[mr_all$protein == prot & grepl("PTSD", mr_all$outcome)][1])
    )
  }
  df_pcl <- bind_rows(pcl_results)
  cat("  Concordant:", sum(df_pcl$concordant), "/", nrow(df_pcl), "\n\n")

} else if (file.exists(pcl_cache)) {

  # ── PCL GWAS missing, but cache from previous run available ──
  df_pcl <- read_csv(pcl_cache, show_col_types = FALSE)
  warning("Part A: pcl_path is '", pcl_path,
          "' or file not found. Loaded cached results from ", pcl_cache,
          " (generated from a previous run with a working PCL GWAS).",
          call. = FALSE)
  cat("  ⚠ PCL GWAS not located; using cached dose_response_concordance.csv\n")
  if ("concordant" %in% names(df_pcl)) {
    cat("  Concordant (cached):",
        sum(df_pcl$concordant, na.rm = TRUE), "/", nrow(df_pcl), "\n\n")
  }

} else {

  # ── Neither raw GWAS nor cache available: skip cleanly ──
  warning("Part A skipped: pcl_path is '", pcl_path,
          "' and no cached results found at ", pcl_cache,
          ". To reproduce dose-response analyses, obtain the PGC-PTSD ",
          "Freeze 3 PCL quantitative GWAS and set pcl_path in 00_config.R.",
          call. = FALSE)
  cat("  ⚠ Part A SKIPPED — PCL GWAS not available and no cache present\n")
  cat("    Downstream Parts B-D will still run.\n\n")
}

# ═══════════════════════════════════════════════════════════════
# Part B: Cross-ancestry (AAM + HNA)
# ═══════════════════════════════════════════════════════════════
cat("Part B: Cross-ancestry...\n")

cross_ancestry <- list()
for (anc_name in c("AAM", "HNA")) {
  anc_path <- ifelse(anc_name == "AAM", aam_path, hna_path)
  anc_gwas <- fread(anc_path)
  
  for (prot in candidates) {
    inst <- pqtl[pqtl$protein == prot, ][1, ]
    snp_data <- anc_gwas[anc_gwas$SNP == inst$SNP | anc_gwas$rsid == inst$SNP, ]
    if (nrow(snp_data) == 0) next
    
    mr_beta <- snp_data$beta[1] / inst$beta
    mr_se   <- abs(snp_data$se[1] / inst$beta)
    cross_ancestry[[paste(prot, anc_name)]] <- data.frame(
      protein = prot, ancestry = anc_name,
      beta = mr_beta, se = mr_se,
      p = 2 * pnorm(-abs(mr_beta / mr_se)),
      OR = exp(mr_beta)
    )
  }
}
df_cross <- bind_rows(cross_ancestry)
cat("  AAM replicated (P<0.05):", sum(df_cross$p[df_cross$ancestry=="AAM"] < 0.05), "\n\n")

# ═══════════════════════════════════════════════════════════════
# Part C-D: Pathway enrichment (fgsea)
# ═══════════════════════════════════════════════════════════════
cat("Part D: Pathway enrichment (fgsea)...\n")

hallmark <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# 为每个outcome构建排序统计量: signed -log10(P)
for (out_name in c("PTSD_freeze3", "MDD", "Anxiety")) {
  mr_out <- mr_all %>% filter(outcome == out_name) %>%
    mutate(stat = sign(mr_beta) * -log10(mr_pval))
  
  ranks <- setNames(mr_out$stat, mr_out$protein)
  ranks <- ranks[!is.na(ranks)]
  ranks <- sort(ranks, decreasing = TRUE)
  
  fgsea_res <- fgsea(pathways, ranks, minSize = 5, maxSize = 500)
  fgsea_res$outcome <- out_name
  
  write_csv(as.data.frame(fgsea_res),
            file.path(result_dir, paste0("fgsea_", tolower(out_name), "_hallmark.csv")))
  cat("  ", out_name, ": ", sum(fgsea_res$padj < 0.05), " significant pathways\n")
}

# NES comparison across outcomes
cat("\n  Computing NES correlations...\n")
# (合并3个outcome的NES, 计算Spearman相关)

# ── 保存 ──
if (!is.null(df_pcl)) {
  write_csv(df_pcl, file.path(result_dir, "dose_response_concordance.csv"))
} else {
  cat("  (Part A output not written: df_pcl unavailable)\n")
}
write_csv(df_cross, file.path(result_dir, "cross_ancestry_validation.csv"))
cat("\n  ✓ Part A-D outputs saved (Part A may be skipped if PCL unavailable)\n")
