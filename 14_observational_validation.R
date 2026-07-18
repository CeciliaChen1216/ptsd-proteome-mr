###############################################################################
# 14_observational_validation.R — Observational Protein-Level Validation
#
# Two independent observational sources for the candidate proteins:
#   (1) UK Biobank Olink plasma × PTSD (Daskalakis et al. 2024 Science)
#       > 16,000 participants, PTSD case-control + symptom score
#       Data source: Supplementary Table S10A (seq10)
#
#   (2) Post-mortem brain tissue mass-spectrometry proteomics
#       (Wang et al. 2025 Genome Medicine)
#       N = 66, DLPFC + sgPFC, tandem mass spectrometry
#       Data source: Additional File 2
#
# Outputs:
#   observational_validation_ukb.csv
#   observational_validation_brain.csv
#   observational_validation_summary.csv
###############################################################################
source("00_config.R")

library(openxlsx)
library(dplyr)
library(readr)

candidates  <- c("AKT3","CD40","CGREF1","FES","FURIN",
                 "SIRPA","CD101","KHK","SNX18","UBE2L6")

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 14: Observational Protein-Level Validation     ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# MR方向 (用于一致性判断)
mr_direction <- data.frame(
  gene = candidates,
  mr_or = c(0.754, 0.953, 0.970, 1.200, 0.913,
            0.987, 1.023, 1.038, 0.816, 1.109),
  mr_dir = c("protective","protective","protective","risk","protective",
             "protective","risk","risk","protective","risk"),
  stringsAsFactors = FALSE
)

# ═══════════════════════════════════════════════════════════════
# 1. UKB Olink × PTSD (Daskalakis 2024)
# ═══════════════════════════════════════════════════════════════
# 需要先用Excel将 adh3707_Suppl. Excel_seq10_v6.xlsb 另存为 .xlsx

cat("═══ Part 1: UKB Plasma Protein × PTSD ═══\n\n")

f10 <- file.path(validation_dir, "Daskalakis_2024_Science/adh3707_Suppl. Excel_seq10_v6.xlsx")
genes_search <- c(candidates, "IGSF2")  # IGSF2 = CD101

if (file.exists(f10)) {
  # Case-control
  s_caco <- read.xlsx(f10, sheet = "S10A-1 (PTSD CACO)", startRow = 1)
  hit_caco <- s_caco[s_caco$gene %in% genes_search, ]
  
  cat("PTSD Case-Control:\n")
  if (nrow(hit_caco) > 0) {
    print(hit_caco[, c("gene","BETA","SE","P","p.adj","N_PTSDCase","N_PTSDControl")])
  } else {
    cat("  No candidates in significant set\n")
  }
  
  # Continuous score
  s_score <- read.xlsx(f10, sheet = "S10A-2 (PTSD Score)", startRow = 1)
  hit_score <- s_score[s_score$gene %in% genes_search, ]
  
  cat("\nPTSD Symptom Score:\n")
  if (nrow(hit_score) > 0) {
    print(hit_score[, c("gene","BETA","SE","P","p.adj")])
  } else {
    cat("  No candidates in significant set\n")
  }
  
  # 合并UKB结果
  ukb_results <- bind_rows(
    if (nrow(hit_caco) > 0) hit_caco %>% 
      select(gene, BETA, SE, P, p.adj) %>% mutate(phenotype = "PTSD_CaseControl"),
    if (nrow(hit_score) > 0) hit_score %>% 
      select(gene, BETA, SE, P, p.adj) %>% mutate(phenotype = "PTSD_Score")
  )
  
  if (nrow(ukb_results) > 0) {
    # 方向一致性判断
    ukb_results <- ukb_results %>%
      left_join(mr_direction, by = "gene") %>%
      mutate(
        obs_direction = ifelse(BETA > 0, "elevated_in_PTSD", "reduced_in_PTSD"),
        # MR protective + obs elevated = reverse causation (PTSD → protein ↑)
        # MR protective + obs reduced = concordant (protein ↓ → PTSD ↑)
        interpretation = case_when(
          mr_dir == "protective" & BETA > 0 ~ "Reverse causation (compensatory upregulation)",
          mr_dir == "protective" & BETA < 0 ~ "Concordant (lower protein → higher PTSD risk)",
          mr_dir == "risk" & BETA > 0       ~ "Concordant (higher protein → higher PTSD risk)",
          mr_dir == "risk" & BETA < 0       ~ "Reverse causation",
          TRUE ~ "Unclear"
        )
      )
    
    cat("\nDirectional interpretation:\n")
    for (i in 1:nrow(ukb_results)) {
      cat("  ", ukb_results$gene[i], " (", ukb_results$phenotype[i], "):",
          " beta=", round(ukb_results$BETA[i], 4),
          " P=", signif(ukb_results$P[i], 3),
          " | MR:", ukb_results$mr_dir[i],
          " | ", ukb_results$interpretation[i], "\n")
    }
    
    write_csv(ukb_results, file.path(result_dir, "observational_validation_ukb.csv"))
    cat("\n✓ observational_validation_ukb.csv\n")
  }
} else {
  cat("⚠ File not found:", f10, "\n")
  cat("  Open adh3707_Suppl. Excel_seq10_v6.xlsb in Excel and Save As .xlsx first\n")
}

# ═══════════════════════════════════════════════════════════════
# 2. 死后脑组织蛋白质谱 (Wang et al. 2025)
# ═══════════════════════════════════════════════════════════════

cat("\n═══ Part 2: Brain Protein MS (Wang 2025) ═══\n\n")

f_wang <- file.path(validation_dir, "Wang_2025_GenomeMed/13073_2025_1473_MOESM2_ESM.xlsx")

if (file.exists(f_wang)) {
  brain_results <- list()
  
  for (region in c("dlPFC", "sgPFC")) {
    s <- read.xlsx(f_wang, sheet = region, startRow = 1)
    hit <- s[s$Gene.Name %in% genes_search, ]
    
    if (nrow(hit) > 0) {
      for (i in 1:nrow(hit)) {
        brain_results[[paste(hit$Gene.Name[i], region)]] <- data.frame(
          gene = hit$Gene.Name[i],
          brain_region = region,
          PTSD_logFC = hit$PTSD.logFC[i],
          PTSD_P = hit$PTSD.P.Value[i],
          PTSD_FDR = hit$PTSD.adj.P.Val[i],
          MDD_logFC = hit$MDD.logFC[i],
          MDD_P = hit$MDD.P.Value[i],
          stringsAsFactors = FALSE
        )
        cat("  ", hit$Gene.Name[i], " (", region, "):",
            " logFC=", round(hit$PTSD.logFC[i], 4),
            " P=", signif(hit$PTSD.P.Value[i], 3), "\n")
      }
    }
  }
  
  df_brain <- bind_rows(brain_results)
  if (nrow(df_brain) > 0) {
    df_brain <- df_brain %>%
      left_join(mr_direction, by = "gene") %>%
      mutate(
        direction_vs_MR = case_when(
          mr_dir == "protective" & PTSD_logFC < 0 ~ "Concordant (reduced in PTSD)",
          mr_dir == "protective" & PTSD_logFC > 0 ~ "Discordant",
          mr_dir == "risk" & PTSD_logFC > 0       ~ "Concordant (elevated in PTSD)",
          mr_dir == "risk" & PTSD_logFC < 0       ~ "Discordant",
          TRUE ~ "Unclear"
        )
      )
    write_csv(df_brain, file.path(result_dir, "observational_validation_brain.csv"))
    cat("\n✓ observational_validation_brain.csv\n")
  }
} else {
  cat("⚠ File not found:", f_wang, "\n")
}

# ═══════════════════════════════════════════════════════════════
# 3. 小胶质细胞 snRNA-seq (Daskalakis 2024)
# ═══════════════════════════════════════════════════════════════

cat("\n═══ Part 3: Microglia snRNA-seq (Daskalakis 2024) ═══\n\n")

f8 <- file.path(validation_dir, "Daskalakis_2024_Science/adh3707_Suppl. Excel_seq8_v5.xlsx")

if (file.exists(f8)) {
  micro <- read.xlsx(f8, sheet = "S8A-5 (Microglia_PTSD)", startRow = 1)
  micro_hit <- micro[micro$genes %in% genes_search,
                     c("genes","beta","pval","FDR")]
  cat("Microglia PTSD DEGs:\n")
  print(micro_hit)
} else {
  cat("⚠ seq8 file not found\n")
}

# ═══════════════════════════════════════════════════════════════
# 4. 汇总
# ═══════════════════════════════════════════════════════════════

cat("\n═══ Summary: observational protein validation ═══\n\n")

summary_obs <- data.frame(
  Protein = candidates,
  UKB_PTSD_CACO = c("—","—","—","—",
                    "beta=0.131, P=1.49e-6",
                    "—","—","—","—","—"),
  UKB_PTSD_Score = c("—","—",
                     "beta=0.008, P=1.78e-4",
                     "—",
                     "beta=0.014, P=3.75e-12",
                     "—","—","—","—","—"),
  Brain_DLPFC = c("logFC=-0.009, P=0.74","Not detected","Not detected",
                  "Not detected","Not detected",
                  "logFC=-0.044, P=0.14","Not detected",
                  "Not detected","Not detected","Not detected"),
  stringsAsFactors = FALSE
)

print(summary_obs, row.names = FALSE)
write_csv(summary_obs, file.path(result_dir, "observational_validation_summary.csv"))
cat("\n✓ observational_validation_summary.csv\n")
cat("\nDone!\n")
