###############################################################################
# 20_generate_supplementary_tables.R
#
# Append two new supplementary tables to the existing
# Additional_file_3-Supplementary_Tables.xlsx (which already contains
# Tables S3–S19):
#
#   Table S20: Cross-tissue brain pQTL support via Wingo 2025 brain pQTL evidence
#   Table S21: Trauma-exposure negative-control MR
#
# Output: Additional_file_3-Supplementary_Tables_v2.xlsx (next to original)
#
# Prerequisites:
#   - results/trauma_negative_control_MR.csv exists (produced by 18 + 18b)
#   - The original supplementary Tables xlsx is located in <project_dir>/data/raw/
#     (or pass its path via the SUPP_TABLES_XLSX environment variable)
#
# Dependencies: openxlsx, data.table
###############################################################################

# Locate 00_config.R
get_script_dir <- function() {
  for (n in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(n)$ofile, error = function(e) NULL)
    if (!is.null(f) && is.character(f) && nzchar(f)) {
      return(dirname(normalizePath(f, mustWork = FALSE)))
    }
  }
  NA_character_
}
script_dir <- get_script_dir()
config_candidates <- c(if (!is.na(script_dir)) script_dir else "", getwd())
config_path <- ""
for (d in unique(config_candidates[nzchar(config_candidates)])) {
  cand <- file.path(d, "00_config.R")
  if (file.exists(cand)) { config_path <- cand; break }
}
if (!nzchar(config_path)) stop("Cannot find 00_config.R. Please copy 00_config_template.R to 00_config.R and edit local paths.")
source(config_path)

suppressMessages({ library(openxlsx); library(data.table) })

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   Generate supplementary tables: add Table S20 + S21           ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 1: 找到原 supplementary tables xlsx
# ─────────────────────────────────────────────────────────────────────────
candidates_xlsx <- c(
  Sys.getenv("SUPP_TABLES_XLSX", unset = ""),
  file.path(script_dir, "..", "data", "raw", "Additional_file_3-Supplementary_Tables.xlsx"),
  file.path(script_dir, "..", "Additional_file_3-Supplementary_Tables.xlsx")
)
candidates_xlsx <- candidates_xlsx[nzchar(candidates_xlsx)]
src_xlsx <- ""
for (p in candidates_xlsx) {
  if (file.exists(p)) { src_xlsx <- p; break }
}
if (!nzchar(src_xlsx)) {
  cat("None of the candidate locations exist:\n")
  for (p in candidates_xlsx) cat(sprintf("  - %s\n", p))
  cat("\nSet the SUPP_TABLES_XLSX environment variable to the original xlsx, or\n")
  cat("place the file at <project_dir>/data/raw/Additional_file_3-Supplementary_Tables.xlsx\n\n")
  stop("Cannot find the original supplementary tables xlsx. See script comments for expected locations.")
}
cat(sprintf("✓ Source file: %s\n", src_xlsx))

# ─────────────────────────────────────────────────────────────────────────
# Step 2: Table S20 — Wingo 2025 brain pQTL cross-tissue support
# 由于 19 主脚本只机械命中 KHK, 但人工解读 Wingo 2025 supp 各表后
# (见 wingo2025_inspect_report.txt) 实际有 8 个候选在 Wingo 2025 中有数据,
# 这里把每个候选的真实状态精确编码 (基于报告内容)
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 2: Build Table S20 (Wingo 2025 brain pQTL support) ─────\n")

S20 <- data.table(
  Protein = c("KHK", "CD40", "SIRPA", "AKT3", "CGREF1",
              "FES", "SNX18", "UBE2L6", "FURIN", "CD101"),
  UniProt = c("P50053", "P25942", "P78324", "Q9Y243", "Q99674",
              "P07332", "Q96RF0", "O14933", "P09958", "Q93033"),
  `Detected in NHW brain proteome` = c("Yes","Yes","Yes","Yes","Yes",
                                        "Yes","Yes","Yes (brain only)","No","No"),
  `Multi-ancestry causal pQTL (PIP)` = c(
    "rs2304681 (PIP=0.9999999)",
    "rs35377099 (PIP=0.769)",
    "rs6075339 (PIP=0.885)",
    "rs3006933 (PIP=0.555)",
    "—", "—", "—", "—", "—", "—"),
  `Variant annotation` = c(
    "Nonsynonymous coding (NM_006488:c.145G>A)",
    "—", "Intronic", "—", "—", "—", "—", "—", "—", "—"),
  `Brain SMR/PWAS for PTSD (Wingo S11/S17/S21)` = c(
    "Yes (Causal pair Table S17)",
    "Yes (Causal gene Table S21, shared with BD/MDD/Neuroticism/SCZ)",
    "No for PTSD (significant for AD only)",
    "No for PTSD (significant for SCZ only)",
    "No for PTSD (significant for MDD/BD only)",
    "No (no MR/PWAS hit for any trait)",
    "No (no MR/PWAS hit for any trait)",
    "No (brain-unique pGene, no MR hit)",
    "Not applicable (protein not detected in brain)",
    "Not applicable (protein not detected in brain)"),
  `Brain PMR-Egger for PTSD (Wingo S14)` = c(
    "β=+0.011, P=8.8e-5, FDR=0.011, pleiotropy P=0.28",
    "Significant in multiple traits (AUD/ADHD); PTSD not in S14",
    "—", "—", "—", "—", "—", "—", "—", "—"),
  `Cross-tissue brain pQTL support for PTSD` = c(
    "Concordant (multi-method)",
    "Concordant (PWAS+SMR)",
    "Brain pQTL support for a related disorder, not PTSD-specific",
    "Brain pQTL support for a related disorder, not PTSD-specific",
    "Brain pQTL support for a related disorder, not PTSD-specific",
    "Brain detected only, no supporting signal",
    "Brain detected only, no supporting signal",
    "Brain-unique pGene, no supporting signal",
    "Not detectable in brain proteome",
    "Not detectable in brain proteome"),
  check.names = FALSE
)

print(S20[, .(Protein, `Detected in NHW brain proteome`,
              `Cross-tissue brain pQTL support for PTSD`)],
      row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────
# Step 3: Table S21 — Childhood-maltreatment secondary-phenotype analysis
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 3: Build Table S21 (childhood-maltreatment secondary phenotype) ─────\n")

trauma_csv <- file.path(result_dir, "trauma_negative_control_MR.csv")
if (!file.exists(trauma_csv)) {
  stop(sprintf("Cannot find %s. Please run script 18b first.", trauma_csv))
}
trauma <- fread(trauma_csv)
cat(sprintf("Read: %s (%d rows)\n", trauma_csv, nrow(trauma)))

# 重新整理列, 改成 manuscript-friendly 名字
S21 <- trauma[, .(
  Protein = candidate,
  `cis-pQTL SNP` = cis_pQTL_SNP,
  `EA / NEA (pQTL)` = paste0(EA_pQTL, " / ", NEA_pQTL),
  `β (pQTL)` = signif(beta_pQTL, 3),
  `SE (pQTL)` = signif(se_pQTL, 3),
  `P (pQTL)` = signif(p_pQTL, 3),
  `PTSD outcome (primary)` = ptsd_outcome,
  `β (PTSD MR)` = signif(ptsd_beta, 3),
  `P (PTSD MR)` = signif(ptsd_p, 3),
  `FDR (PTSD MR)` = signif(ptsd_fdr, 3),
  `β (Childhood-maltreatment MR, Warrier 2021)` = signif(wald_beta, 3),
  `SE (Maltreatment MR)` = signif(wald_se, 3),
  `P (Maltreatment MR)` = signif(wald_p, 3),
  `Allele harmonization` = harmonization,
  Verdict = verdict
)]

# 排序: PTSD-associated (no maltreatment) first, maltreatment-associated after
S21$Verdict <- factor(S21$Verdict,
                      levels = c("PTSD-associated; no nominal maltreatment association", "associated with childhood maltreatment",
                                 "opposite_direction", "ptsd_not_significant",
                                 "unable_to_test"))
S21 <- S21[order(Verdict, Protein)]
S21$Verdict <- as.character(S21$Verdict)

print(S21[, .(Protein, `P (PTSD MR)`, `P (Maltreatment MR)`, Verdict)],
      row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────
# Step 4: 加载原 xlsx, 添加新 sheet, 保存
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 4: Write v2 xlsx ─────\n")

wb <- loadWorkbook(src_xlsx)
existing <- names(wb)
cat(sprintf("Source xlsx contains %d sheets: %s\n", length(existing),
            paste(existing, collapse = ", ")))

# 删除已有的 S20/S21 (如重跑)
for (s in c("TableS20", "TableS21")) {
  if (s %in% existing) {
    removeWorksheet(wb, s)
    cat(sprintf("Removed existing %s (will regenerate)\n", s))
  }
}

# Style 定义 (复用现有表的风格)
title_style <- createStyle(textDecoration = "bold", fontSize = 11,
                            wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold",
                             border = "Bottom", borderStyle = "medium",
                             fgFill = "#F2F2F2",
                             halign = "center", valign = "center",
                             wrapText = TRUE)

# === Table S20 ===
addWorksheet(wb, "TableS20")
S20_title <- paste0(
  "Table S20. Cross-tissue brain pQTL support for ten PTSD-prioritized proteins ",
  "via brain pQTL Mendelian randomization (Wingo et al., Nature Genetics 2025; ",
  "DOI: 10.1038/s41588-025-02291-2; N=1,362 dorsolateral prefrontal cortex ",
  "proteomes, multi-ancestry). Brain SMR/PWAS column reflects evidence from ",
  "Wingo 2025 Supplementary Tables S11 (PWAS+SMR), S17 (multi-ancestry causal ",
  "pQTL/trait pairs in NHW), S21 (drug repurposing). PMR-Egger column from ",
  "Wingo 2025 Supplementary Table S14 (FDR<0.05 with pleiotropy P>0.05 ",
  "criteria for causality). Two candidates (KHK and CD40) showed concordant ",
  "cross-tissue support for PTSD; three additional candidates (SIRPA, AKT3, ",
  "CGREF1) had brain pQTL support for other related psychiatric ",
  "disorders but not PTSD-specifically; two (FURIN, CD101) were not detected ",
  "in the brain proteome.")
writeData(wb, "TableS20", S20_title, startRow = 1, startCol = 1,
          colNames = FALSE)
mergeCells(wb, "TableS20", cols = 1:ncol(S20), rows = 1)
addStyle(wb, "TableS20", title_style, rows = 1, cols = 1)
setRowHeights(wb, "TableS20", rows = 1, heights = 90)

writeData(wb, "TableS20", S20, startRow = 2, startCol = 1,
          colNames = TRUE, headerStyle = header_style)
setColWidths(wb, "TableS20", cols = 1:ncol(S20),
             widths = c(10, 9, 14, 22, 24, 50, 45, 35))
setRowHeights(wb, "TableS20", rows = 2, heights = 50)
freezePane(wb, "TableS20", firstActiveRow = 3, firstActiveCol = 2)

# === Table S21 ===
addWorksheet(wb, "TableS21")
S21_title <- paste0(
  "Table S21. Childhood-maltreatment secondary-phenotype analysis for ten PTSD-prioritized ",
  "proteins, using childhood maltreatment as a secondary-phenotype outcome ",
  "(Warrier et al., Lancet Psychiatry 2021; PMID 33740410; N=185,414 ",
  "retrospective + prospective meta-analysis). Wald ratio MR was performed ",
  "using each protein's primary cis-pQTL as a single instrument. PTSD MR ",
  "results are reported for the primary PTSD outcome (PTSD freeze 3 case-",
  "control or PCL quantitative severity, whichever yielded the smallest P ",
  "value). Verdicts are descriptive: PTSD-associated proteins with no nominal ",
  "childhood-maltreatment association (FDR<0.05 vs P>0.05) versus those also ",
  "nominally associated with maltreatment in the same direction. Six of ten ",
  "candidates showed little evidence of a childhood-maltreatment association. ",
  "This is a secondary-phenotype sensitivity analysis; childhood maltreatment is ",
  "not a strict negative control, and non-association was not interpreted as ",
  "evidence of PTSD-specificity.")
writeData(wb, "TableS21", S21_title, startRow = 1, startCol = 1,
          colNames = FALSE)
mergeCells(wb, "TableS21", cols = 1:ncol(S21), rows = 1)
addStyle(wb, "TableS21", title_style, rows = 1, cols = 1)
setRowHeights(wb, "TableS21", rows = 1, heights = 80)

writeData(wb, "TableS21", S21, startRow = 2, startCol = 1,
          colNames = TRUE, headerStyle = header_style)
setColWidths(wb, "TableS21", cols = 1:ncol(S21),
             widths = c(10, 13, 15, 11, 11, 11, 22, 13, 13, 13,
                        18, 13, 13, 13, 22))
setRowHeights(wb, "TableS21", rows = 2, heights = 45)
freezePane(wb, "TableS21", firstActiveRow = 3, firstActiveCol = 2)

# 保存
out_xlsx <- sub("\\.xlsx$", "_v2.xlsx", src_xlsx)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n✓ Written: %s\n", out_xlsx))

# ─────────────────────────────────────────────────────────────────────────
# Step 5: 摘要 + 下一步建议
# ─────────────────────────────────────────────────────────────────────────
cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Summary\n")
cat("══════════════════════════════════════════════════════════════════\n")
cat(sprintf("\nSource xlsx: %s\n", src_xlsx))
cat(sprintf("Output xlsx: %s\n", out_xlsx))
cat(sprintf("Total sheets: %d (existing %d + 2 new)\n",
            length(names(loadWorkbook(out_xlsx))),
            length(existing)))
cat("\nNew content:\n")
cat("  Table S20: Wingo 2025 brain pQTL cross-tissue support\n")
cat("             10 candidates × 8 columns (incl. PIP, variant annotation, MR results)\n")
cat("             Highlight: KHK (nonsynonymous coding variant + multi-method),\n")
cat("                        CD40 (PWAS+SMR shared with BD/MDD/Neuro/SCZ)\n")
cat("\n  Table S21: Childhood-maltreatment secondary-phenotype analysis (Warrier 2021)\n")
cat("             10 candidates × 15 columns (incl. pQTL effects, PTSD MR, maltreatment MR, verdict)\n")
cat("             6 PTSD-associated without maltreatment association (KHK/CD101/CGREF1/SIRPA/SNX18/UBE2L6)\n")
cat("             4 also associated with childhood maltreatment (AKT3/CD40/FES/FURIN)\n")

cat("\nNext steps:\n")
cat("  1) Inspect Table S20 / S21 content and formatting in v2 xlsx\n")
cat("  2) Update the manuscript main text accordingly:\n")
cat("     - Add a cross-tissue brain pQTL support paragraph in Discussion (cite Table S20)\n")
cat("     - Add a negative-control paragraph in Discussion (cite Table S21)\n")
cat("     - Update Limitations to acknowledge chr15:90.87Mb FES/FURIN polypleiotropy\n")
cat("  3) Consider adding 'Cross-tissue support' and 'Maltreatment secondary-phenotype' columns to main Table 1\n")
cat("  4) Run 17_conditional_chr15_COJO.R (GCTA-COJO) to resolve FES/FURIN polypleiotropy\n\n")
