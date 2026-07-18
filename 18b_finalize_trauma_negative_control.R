###############################################################################
# 18b_finalize_trauma_negative_control.R
#
# Finalization step for script 18. Reads the trauma-MR results object produced
# by script 18, merges them with the PTSD primary-MR results using an extended
# column-name parser (recognising the `mr_beta` / `mr_pval` prefixes used in
# `mr_all_outcomes_extended.rds`), and assigns the final verdict
# (post_pathology / behavioral_confound / opposite_direction).
#
# Run order: 18 first, then 18b. The output overwrites
#   results/trauma_negative_control_MR.csv
#   results/trauma_negative_control_MR.rds
# with the verdict-annotated table used by Table S21.
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

suppressMessages({ library(data.table); library(dplyr) })

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   18b. Finalize trauma NC verdict                              ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 1: 加载已生成的 trauma MR 主表
# ─────────────────────────────────────────────────────────────────────────
prev_path <- file.path(result_dir, "trauma_negative_control_MR.rds")
if (!file.exists(prev_path)) {
  stop(sprintf("Cannot find %s. Please run script 18 first.", prev_path))
}
prev <- readRDS(prev_path)
results <- as.data.table(prev$results)
cat(sprintf("✓ Loaded: %s (%d rows)\n", prev_path, nrow(results)))

# 清掉旧的 PTSD 列 (要重填)
# Clear ALL previously merged PTSD columns (including any .x/.y suffixes left by
# a prior merge) so re-running this script does not accumulate duplicate columns.
for (col in c("ptsd_beta", "ptsd_p", "ptsd_fdr", "ptsd_outcome", "verdict",
              "ptsd_outcome.x", "ptsd_fdr.x", "ptsd_beta.x", "ptsd_p.x",
              "ptsd_outcome.y", "ptsd_fdr.y", "ptsd_beta.y", "ptsd_p.y")) {
  if (col %in% colnames(results)) results[[col]] <- NULL
}

# ─────────────────────────────────────────────────────────────────────────
# Step 2: 读 PTSD MR, 用扩展的列名匹配
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 2: Load PTSD MR (mr_all_outcomes_extended.rds) ─────\n")
mr_path <- file.path(result_dir, "mr_all_outcomes_extended.rds")
mr <- as.data.table(readRDS(mr_path))
cat(sprintf("Shape: %d rows × %d columns\n", nrow(mr), ncol(mr)))
cat(sprintf("Column names: %s\n", paste(colnames(mr), collapse = ", ")))

# outcome 列的 unique values (确认 PTSD 编码)
cat("\nUnique values in the outcome column:\n")
print(unique(mr$outcome))

# 扩展的列识别
cn <- colnames(mr)
out_col <- grep("^outcome$|trait|disease", cn, value = TRUE, ignore.case = TRUE)[1]
pro_col <- grep("^protein$|^exposure$|^gene$", cn, value = TRUE, ignore.case = TRUE)[1]
b_col   <- grep("^mr_beta$|^beta$|^b$|mr.?beta", cn, value = TRUE, ignore.case = TRUE)[1]
p_col   <- grep("^mr_pval$|^pval$|^p$|mr.?pval|mr.?p$", cn, value = TRUE, ignore.case = TRUE)[1]
fdr_col <- grep("^fdr$|^p_adj|^padj", cn, value = TRUE, ignore.case = TRUE)[1]

cat(sprintf("\nIdentified columns: outcome=%s, protein=%s, beta=%s, p=%s, fdr=%s\n",
            out_col, pro_col, b_col, p_col, fdr_col))

# ─────────────────────────────────────────────────────────────────────────
# Step 3: 提取 PTSD × candidates
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 3: Extract PTSD MR for the 10 candidates ─────\n")

ptsd_mask <- grepl("ptsd|stress", mr[[out_col]], ignore.case = TRUE)
sub <- mr[ptsd_mask & mr[[pro_col]] %in% candidates]
cat(sprintf("Rows for PTSD outcomes × candidate proteins: %d\n", nrow(sub)))

if (nrow(sub) == 0) {
  cat("\n⚠ No PTSD-related rows matched the candidate proteins; check outcome column matching.\n")
  cat("Trying case-insensitive match on 'ptsd':\n")
  print(unique(mr[[out_col]])[grep("ptsd", unique(mr[[out_col]]),
                                    ignore.case = TRUE)])
}

# 每蛋白每 outcome 可能有多行 (multi-instrument), 但同时只一个 SNP/protein/outcome,
# 应该只 1 行。如多行, 取 P 最小的那行 (对于 PTSD primary outcome)
# 先按 P 排序, 再按 (protein, outcome) 取首行
sub <- sub[order(get(p_col))]
sub_uniq <- sub[, .SD[1], by = .(get(pro_col), get(out_col))]
setnames(sub_uniq, c("get", "get.1"), c(pro_col, out_col))

cat("\nCandidate-protein hits across PTSD outcomes:\n")
print(sub_uniq[, .(protein = get(pro_col), outcome = get(out_col),
                   beta = signif(get(b_col), 3),
                   p = signif(get(p_col), 3))],
      row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────
# Step 4: Merge the PRE-SPECIFIED primary PTSD outcome (freeze-3 case-control)
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 4: Use pre-specified primary PTSD outcome (freeze-3) and merge into trauma table ─────\n")

# Pre-specify PTSD freeze-3 case-control as the single primary outcome for every
# candidate. This is the primary PTSD GWAS used throughout the study. Previously
# this step selected, per protein, whichever PTSD outcome gave the smallest P
# (freeze-3 vs PCL quantitative severity), which is an outcome-selection step
# that can bias the reported PTSD association toward significance and inflate the
# contrast underlying the negative-control verdict. All ten candidates are
# FDR-significant under freeze-3, so fixing the outcome removes the selection step
# without changing any verdict. PCL quantitative severity is retained only as a
# secondary outcome (see secondary columns / supplementary reporting).
PRIMARY_OUTCOME <- "PTSD_freeze3"
primary <- sub_uniq[get(out_col) == PRIMARY_OUTCOME]

# Safety fallback: if a candidate lacks a freeze-3 row, fall back to smallest-P
# for that candidate only, with a warning (should not trigger for the 10 candidates).
missing_primary <- setdiff(unique(sub_uniq[[pro_col]]), primary[[pro_col]])
if (length(missing_primary) > 0) {
  warning("No freeze-3 row for: ", paste(missing_primary, collapse = ", "),
          " - falling back to smallest-P for these candidates.")
  fb <- sub_uniq[get(pro_col) %in% missing_primary][, .SD[which.min(get(p_col))], by = pro_col]
  primary <- rbind(primary, fb)
}
# 标准化列
ptsd_summary <- data.table(
  candidate = primary[[pro_col]],
  ptsd_outcome = primary[[out_col]],
  ptsd_beta = as.numeric(primary[[b_col]]),
  ptsd_p = as.numeric(primary[[p_col]])
)
if (!is.na(fdr_col)) {
  ptsd_summary$ptsd_fdr <- as.numeric(primary[[fdr_col]])
}

cat("Final PTSD MR (pre-specified freeze-3 primary outcome):\n")
print(ptsd_summary, row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────
# Step 5: Verdict
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 5: Compute final verdicts ─────\n")

# 用 fdr 而不是裸 p 来定 PTSD 显著性 (如果有 fdr 列的话)
ptsd_sig_col <- if ("ptsd_fdr" %in% colnames(ptsd_summary)) "ptsd_fdr" else "ptsd_p"
cat(sprintf("PTSD significance criterion: %s < 0.05\n", ptsd_sig_col))

results <- merge(results, ptsd_summary, by = "candidate",
                 all.x = TRUE, sort = FALSE)

classify <- function(ptsd_b, ptsd_sig, trauma_b, trauma_p,
                     alpha_ptsd = 0.05, alpha_trauma = 0.05) {
  if (is.na(trauma_b) || is.na(trauma_p)) return("unable_to_test")
  if (is.na(ptsd_sig) || ptsd_sig > alpha_ptsd) return("ptsd_not_significant")
  if (trauma_p > alpha_trauma) return("post_pathology")        # ✅
  if (!is.na(ptsd_b)) {
    if (sign(ptsd_b) == sign(trauma_b)) return("behavioral_confound")  # ⚠
    else return("opposite_direction")
  }
  "trauma_significant_unclassified"
}

ptsd_sig_vec <- results[[ptsd_sig_col]]
results[, verdict := mapply(classify, ptsd_beta, ptsd_sig_vec,
                            wald_beta, wald_p)]

# ─────────────────────────────────────────────────────────────────────────
# Step 6: 输出
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 6: Write outputs ─────\n")

out_csv <- file.path(result_dir, "trauma_negative_control_MR.csv")
fwrite(results, out_csv)
cat(sprintf("✓ Overwritten: %s\n", out_csv))

out_rds <- file.path(result_dir, "trauma_negative_control_MR.rds")
saveRDS(list(results = results,
             instruments = prev$instruments,
             trauma_gwas_files = prev$trauma_gwas_files,
             ptsd_summary = ptsd_summary,
             ptsd_sig_col_used = ptsd_sig_col),
        out_rds)
cat(sprintf("✓ Overwritten: %s\n", out_rds))

# ─────────────────────────────────────────────────────────────────────────
# Step 7: 解读报告
# ─────────────────────────────────────────────────────────────────────────
report <- character(0)
addln <- function(...) report <<- c(report, sprintf(...))

addln("══════════════════════════════════════════════════════════════════")
addln("Trauma Exposure Negative Control MR — 解读 (v2)")
addln("══════════════════════════════════════════════════════════════════")
addln("")
addln("Trauma exposure GWAS: Warrier 2021 Lancet Psy childhood maltreatment")
addln("                       (N=185,414, retrospective+prospective meta)")
addln("PTSD MR 来源: %s", file.path(result_dir, "mr_all_outcomes_extended.rds"))
addln("PTSD 显著性判据: %s < 0.05", ptsd_sig_col)
addln("")
addln("Verdict 逻辑:")
addln("  post_pathology      = PTSD 显著 + trauma 不显著  ✅ 真创伤后病理")
addln("  behavioral_confound = PTSD 显著 + trauma 同向显著 ⚠ 行为学混杂")
addln("  opposite_direction  = PTSD 显著 + trauma 反向显著  复杂 pleiotropy")
addln("  ptsd_not_significant = PTSD MR 主分析未达显著 (不在测试范围)")
addln("  unable_to_test      = trauma SNP 缺失 / harmonization 失败")
addln("")

for (g_name in unique(results$trauma_gwas)) {
  addln("──────────────────────────────────────────────────────────────────")
  addln("Trauma GWAS: %s", g_name)
  addln("──────────────────────────────────────────────────────────────────")
  sub_r <- results[trauma_gwas == g_name][order(verdict, candidate)]
  for (i in seq_len(nrow(sub_r))) {
    r <- sub_r[i]
    addln("")
    addln("  %s (cis-pQTL %s, %s>%s)",
          r$candidate, r$cis_pQTL_SNP, r$EA_pQTL, r$NEA_pQTL)
    if (is.na(r$wald_p)) {
      addln("    [无法 MR: %s]", r$harmonization)
      next
    }
    fdr_str <- if ("ptsd_fdr" %in% colnames(r))
                 sprintf(", FDR=%s",
                         if (is.na(r$ptsd_fdr)) "—" else signif(r$ptsd_fdr, 3))
               else ""
    addln("    PTSD MR (%s):%s β=%s, P=%s%s",
          if (!is.na(r$ptsd_outcome)) r$ptsd_outcome else "—",
          strrep(" ", max(0, 18 - nchar(if (!is.na(r$ptsd_outcome)) r$ptsd_outcome else "—"))),
          if (is.na(r$ptsd_beta)) "—" else signif(r$ptsd_beta, 3),
          if (is.na(r$ptsd_p)) "—" else signif(r$ptsd_p, 3),
          fdr_str)
    addln("    Trauma exposure MR:    β=%s, SE=%s, P=%s",
          signif(r$wald_beta, 3),
          signif(r$wald_se, 3),
          signif(r$wald_p, 3))
    addln("    Verdict:               %s", r$verdict)
  }
}

addln("")
addln("══════════════════════════════════════════════════════════════════")
addln("Summary 统计")
addln("══════════════════════════════════════════════════════════════════")
verdict_tbl <- as.data.frame(table(results$verdict))
for (i in seq_len(nrow(verdict_tbl))) {
  addln("  %-32s %d", verdict_tbl[i, 1], verdict_tbl[i, 2])
}

addln("")
addln("──────────────────────────────────────────────────────────────────")
addln("关键蛋白汇总 (Tier 1 候选, 按 verdict 分类):")
addln("──────────────────────────────────────────────────────────────────")
for (v in c("post_pathology", "behavioral_confound", "opposite_direction",
            "ptsd_not_significant", "unable_to_test")) {
  matched <- results[verdict == v]
  if (nrow(matched) > 0) {
    addln("")
    addln("  [%s]", v)
    for (i in seq_len(nrow(matched))) {
      m <- matched[i]
      addln("    %s   PTSD β=%s P=%s | Trauma β=%s P=%s",
            sprintf("%-7s", m$candidate),
            if (is.na(m$ptsd_beta)) "—" else signif(m$ptsd_beta, 3),
            if (is.na(m$ptsd_p)) "—" else signif(m$ptsd_p, 3),
            signif(m$wald_beta, 3),
            signif(m$wald_p, 3))
    }
  }
}

report_path <- file.path(result_dir, "trauma_negative_control_summary.txt")
writeLines(report, report_path)
cat(sprintf("✓ Overwritten: %s\n\n", report_path))

# 控制台打印
cat(paste(report, collapse = "\n"), "\n\n")
