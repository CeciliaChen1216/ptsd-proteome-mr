###############################################################################
# 21_append_conditional_chr15_table.R
#
# Reads results/conditional_chr15_MR.csv (output of script 17) and appends
# Table S22 (chromosome 15 conditional MR) to the supplementary tables xlsx
# produced by script 20.
#
# Input source xlsx is located via, in order:
#   1. SUPP_TABLES_XLSX environment variable (if set);
#   2. <project_dir>/data/raw/;
#   3. the repository root.
# The script prefers files that already contain TableS20 and TableS21
# (i.e. the v2 output from script 20) and writes a *_v3.xlsx alongside.
###############################################################################

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

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Updating supplementary tables: append Table S22 (chr15 conditional MR)\n")
cat("══════════════════════════════════════════════════════════════════\n\n")

# Step 1: Locate any candidate supplementary tables xlsx
search_dirs <- c(
  dirname(Sys.getenv("SUPP_TABLES_XLSX", unset = "")),
  file.path(script_dir, "..", "data", "raw"),
  file.path(script_dir, "..")
)
search_dirs <- search_dirs[nzchar(search_dirs) & dir.exists(search_dirs)]

patterns <- c(
  "(?i)^additional[ _-]+file[ _-]*3.*supp.*tables?.*\\.xlsx$",
  "(?i).*supplementary[ _-]*tables?.*\\.xlsx$"
)

all_hits <- character(0)
for (d in search_dirs) {
  files <- list.files(d, pattern = "\\.xlsx$",
                      ignore.case = TRUE, full.names = TRUE)
  for (pat in patterns) {
    hits <- files[grepl(pat, basename(files), perl = TRUE)]
    # 排除已输出的 _v3 输出, 避免循环
    hits <- hits[!grepl("_v3\\.xlsx$", basename(hits), ignore.case = TRUE)]
    all_hits <- c(all_hits, hits)
  }
}
all_hits <- unique(all_hits)

if (length(all_hits) == 0) {
  cat("No supplementary tables xlsx found. Searched directories:\n")
  for (d in search_dirs) cat(sprintf("  %s\n", d))
  stop("Cannot find the source supplementary tables xlsx (expected output of script 20).")
}

# 优先选含 TableS20/S21 的文件 (说明已加 cross-tissue + trauma 子表)
src_xlsx <- ""
for (f in all_hits) {
  test_wb <- tryCatch(loadWorkbook(f), error = function(e) NULL)
  if (is.null(test_wb)) next
  sheets <- names(test_wb)
  if ("TableS20" %in% sheets && "TableS21" %in% sheets) {
    src_xlsx <- f
    break
  }
}
if (!nzchar(src_xlsx)) {
  cat("Found these xlsx files but none contain both TableS20 and TableS21:\n")
  for (f in all_hits) cat(sprintf("  %s\n", f))
  cat("\nPlease confirm script 20 has been run first.\n")
  cat("Falling back to the first hit (only if you are sure):\n  ", all_hits[1], "\n")
  src_xlsx <- all_hits[1]
}
cat(sprintf("✓ Source: %s\n", src_xlsx))

# Step 2: 读 conditional 结果
cond_csv <- file.path(result_dir, "conditional_chr15_MR.csv")
if (!file.exists(cond_csv)) stop("Cannot find conditional_chr15_MR.csv. Please run script 17 first.")
cond <- fread(cond_csv)
cat(sprintf("✓ Read conditional results: %d rows\n", nrow(cond)))

# Step 3: 整理 Table S22
S22 <- data.table(
  `Test` = cond$label,
  `Query protein` = cond$query_protein,
  `Conditioned on` = cond$cond_on_protein,
  `Query cis-pQTL` = cond$query_SNP,
  `pQTL β` = signif(cond$pQTL_b, 3),
  `pQTL SE` = signif(cond$pQTL_se, 3),
  `pQTL P` = signif(cond$pQTL_p, 3),
  `Unconditional outcome β` = signif(cond$uncond_outcome_b, 3),
  `Unconditional outcome P` = signif(cond$uncond_outcome_p, 3),
  `Conditional outcome β` = signif(cond$cond_outcome_b, 3),
  `Conditional outcome P` = signif(cond$cond_outcome_p, 3),
  `Wald ratio MR β (unconditional)` = signif(cond$mr_uncond_beta, 3),
  `Wald ratio MR P (unconditional)` = signif(cond$mr_uncond_p, 3),
  `Wald ratio MR β (conditional)` = signif(cond$mr_cond_beta, 3),
  `Wald ratio MR P (conditional)` = signif(cond$mr_cond_p, 3),
  `Allele flipped` = cond$flipped,
  `Verdict` = cond$verdict
)

# Step 4: 写入 xlsx
wb <- loadWorkbook(src_xlsx)
existing <- names(wb)
cat(sprintf("\nExisting sheets (%d): %s\n", length(existing),
            paste(existing, collapse = ", ")))
if ("TableS22" %in% existing) {
  removeWorksheet(wb, "TableS22")
  cat("(Removing existing TableS22 to regenerate)\n")
}

title_style <- createStyle(textDecoration = "bold", fontSize = 11,
                            wrapText = TRUE)
header_style <- createStyle(textDecoration = "bold",
                             border = "Bottom", borderStyle = "medium",
                             fgFill = "#F2F2F2",
                             halign = "center", valign = "center",
                             wrapText = TRUE)

addWorksheet(wb, "TableS22")
S22_title <- paste0(
  "Table S22. Conditional Mendelian randomization at the chromosome 15:90.87 Mb ",
  "FES/FURIN locus. GCTA-COJO conditional analysis was performed using the ",
  "1000 Genomes Phase 3 EUR reference panel (N = 633 individuals) restricted to ",
  "chromosome 15:89,877,710-91,885,812 (±1 Mb around the FES and FURIN cis-pQTLs). ",
  "For each of two outcomes (PTSD freeze 3 and Warrier et al. 2021 childhood ",
  "maltreatment), the effect at one protein's cis-pQTL was tested after ",
  "conditioning on the other protein's cis-pQTL. Wald ratio MR estimates use the ",
  "protein's primary cis-pQTL as a single instrument. The FES effect on PTSD ",
  "remained significant after conditioning on FURIN (conditional P = 0.002), ",
  "whereas the FURIN effect attenuated to non-significance after conditioning ",
  "on FES (conditional P = 0.78), indicating that any contribution to PTSD risk ",
  "at this locus is more parsimoniously attributable to FES than to FURIN. ",
  "The FES instrument also remained nominally associated with childhood ",
  "maltreatment after conditioning on FURIN (conditional P = 0.03), supporting ",
  "the locus's polypleiotropic interpretation. ",
  "Verdicts: robust_to_conditioning = significant before and after conditioning; ",
  "attenuated_by_LD = significant only before conditioning."
)
writeData(wb, "TableS22", S22_title, startRow = 1, startCol = 1,
          colNames = FALSE)
mergeCells(wb, "TableS22", cols = 1:ncol(S22), rows = 1)
addStyle(wb, "TableS22", title_style, rows = 1, cols = 1)
setRowHeights(wb, "TableS22", rows = 1, heights = 110)

writeData(wb, "TableS22", S22, startRow = 2, startCol = 1,
          colNames = TRUE, headerStyle = header_style)
setColWidths(wb, "TableS22", cols = 1:ncol(S22),
             widths = c(22, 12, 14, 14, 10, 10, 11, 18, 18, 17, 17,
                        20, 20, 19, 19, 11, 22))
setRowHeights(wb, "TableS22", rows = 2, heights = 60)
freezePane(wb, "TableS22", firstActiveRow = 3, firstActiveCol = 2)

# 输出文件名: 加 _v3 后缀 (即使原文件没有版本号)
out_xlsx <- sub("\\.xlsx$", "_v3.xlsx", src_xlsx)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat(sprintf("\n✓ Written: %s\n", out_xlsx))
cat(sprintf("Sheets: existing %d + new Table S22 = %d\n",
            length(existing), length(names(loadWorkbook(out_xlsx)))))

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Table S22 preview:\n")
cat("══════════════════════════════════════════════════════════════════\n")
print(S22[, .(Test, Verdict,
              `Conditional outcome P`,
              `Wald ratio MR P (conditional)`)],
      row.names = FALSE)
