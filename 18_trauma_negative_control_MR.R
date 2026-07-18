###############################################################################
# 18_trauma_negative_control_MR.R   (v1)
#
# Purpose:
#   Trauma-exposure negative-control MR (additional validation analysis).
#
#   PTSD = traumatic exposure + post-trauma pathology.
#   A candidate protein that is significant for both PTSD AND trauma exposure
#   may reflect a behavioural confound (genetically influencing exposure
#   likelihood); a protein significant for PTSD only points to genuine
#   post-trauma pathology.
#
#   Primary trauma GWAS: Warrier 2021 Lancet Psy childhood maltreatment
#                         (N = 185,414, retrospective + prospective meta-analysis)
#   Optional sensitivity: any additional trauma-exposure GWAS in the trauma_dir.
#
# Outputs:
#   results/trauma_negative_control_MR.csv     main table (10 candidates × ≥1 trauma GWAS)
#   results/trauma_negative_control_MR.rds     full R object
#   results/trauma_negative_control_summary.txt interpretation report
#
# R package dependencies: data.table, dplyr (tidyverse subset)
###############################################################################

# ─────────────────────────────────────────────────────────────────────────
# Step 0: Locate 00_config.R
# ─────────────────────────────────────────────────────────────────────────
get_script_dir <- function() {
  for (n in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(n)$ofile, error = function(e) NULL)
    if (!is.null(f) && is.character(f) && nzchar(f)) {
      return(dirname(normalizePath(f, mustWork = FALSE)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                  error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grepl("^--file=", args)]
  if (length(fa) > 0) {
    p <- sub("^--file=", "", fa[1])
    return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  NA_character_
}

script_dir <- get_script_dir()
config_candidates <- c(
  if (!is.na(script_dir)) script_dir else "",
  getwd()
)
config_candidates <- unique(config_candidates[nzchar(config_candidates)])
config_path <- ""
for (d in config_candidates) {
  cand <- file.path(d, "00_config.R")
  if (file.exists(cand)) { config_path <- cand; break }
}
if (!nzchar(config_path)) stop("Cannot find 00_config.R. Please copy 00_config_template.R to 00_config.R and edit local paths.")
cat(sprintf("✓ Loaded config: %s\n", config_path))
source(config_path)

suppressMessages({
  needed <- c("data.table", "dplyr")
  miss <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
  if (length(miss) > 0) {
    stop(sprintf("Missing R packages. Please install: install.packages(c(%s))",
                 paste(sprintf('"%s"', miss), collapse = ", ")))
  }
  library(data.table)
  library(dplyr)
})

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   18. Trauma Exposure — Negative Control MR                     ║\n")
cat("║       (Warrier 2021 childhood maltreatment + sensitivity)       ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 1: 准备目录, 检查 trauma exposure GWAS 数据
# ─────────────────────────────────────────────────────────────────────────
trauma_dir <- file.path(validation_dir, "trauma_exposure")
if (!dir.exists(trauma_dir)) {
  dir.create(trauma_dir, recursive = TRUE, showWarnings = FALSE)
}
cat(sprintf("Trauma exposure data directory: %s\n\n", trauma_dir))

# 自动列出该目录下所有可能的 GWAS 文件
existing_files <- list.files(trauma_dir, pattern = "\\.(gz|tsv|txt)$",
                             ignore.case = TRUE, full.names = TRUE)
if (length(existing_files) > 0) {
  cat("Files detected:\n")
  for (f in existing_files) {
    cat(sprintf("  %s  (%.1f MB)\n",
                basename(f), file.info(f)$size / 1024^2))
  }
  cat("\n")
}

# 智能识别 Warrier 2021 文件
warrier_pattern <- "warrier|childhood.*maltreatment|maltreatment.*childhood|child.*maltr"
warrier_candidates <- existing_files[grepl(warrier_pattern, basename(existing_files),
                                            ignore.case = TRUE)]
warrier_path <- if (length(warrier_candidates) > 0) warrier_candidates[1] else ""

# Coleman / SLE / lifetime trauma 等可选
sensitivity_pattern <- "coleman|stressful|sle|life.event|trauma.exposure|lifetime.trauma"
sensitivity_candidates <- existing_files[
  grepl(sensitivity_pattern, basename(existing_files), ignore.case = TRUE)]
# 排除已识别为 Warrier 的
sensitivity_candidates <- setdiff(sensitivity_candidates, warrier_path)

# 如果 Warrier 没找到, 给用户清晰指引并停止
if (!nzchar(warrier_path)) {
  cat("❌ Warrier 2021 childhood maltreatment GWAS not found.\n\n")
  cat("Please download manually (Cambridge Apollo does not allow programmatic downloads):\n")
  cat("  1) Open URL: https://www.repository.cam.ac.uk/handle/1810/318326\n")
  cat("     Title: 'Data for: Gene-environment correlations and causal effects\n")
  cat("            of childhood maltreatment...'\n")
  cat("     Authors: Warrier V et al. 2021 Lancet Psychiatry\n")
  cat("     N = 185,414\n\n")
  cat("  2) Download the sumstats file from the Files section (.txt.gz, ~100-300 MB)\n\n")
  cat(sprintf("  3) Place it under: %s\n", trauma_dir))
  cat("     Suggested filename: warrier_2021_childhood_maltreatment.txt.gz\n\n")
  cat("  4) Re-run this script\n\n")
  stop("Warrier 2021 childhood maltreatment GWAS not found. Please download it and place it under data/raw/trauma_exposure/.")
}
cat(sprintf("✓ Primary analysis input: %s (%.1f MB)\n",
            basename(warrier_path), file.info(warrier_path)$size / 1024^2))

if (length(sensitivity_candidates) > 0) {
  cat(sprintf("✓ Also found %d sensitivity GWAS files:\n",
              length(sensitivity_candidates)))
  for (f in sensitivity_candidates) {
    cat(sprintf("    %s\n", basename(f)))
  }
} else {
  cat("(No sensitivity GWAS files found; running primary Warrier analysis only.\n")
  cat(" To add sensitivity GWAS, place files in this directory.)\n")
}
cat("\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 2: 加载 cis-pQTL instruments (10 候选)
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 2: Load candidate-protein cis-pQTLs ─────\n")

pqtl <- readRDS(pqtl_path)
pqtl_dt <- as.data.table(pqtl)

# 取 10 候选
inst <- pqtl_dt[protein %in% candidates]
cat(sprintf("Extracting instruments for %d candidate proteins from %d cis-pQTLs\n",
            nrow(pqtl_dt), nrow(inst)))

# 每个候选可能有多个 cis-pQTL (主+次), 这里取每个蛋白 P 值最小的 (primary instrument)
inst <- inst[order(protein, pval)][, .SD[1], by = protein]
cat(sprintf("One smallest-P SNP per protein; total %d instruments:\n\n",
            nrow(inst)))
print(inst[, .(protein, SNP, chr, pos_hg38,
               effect_allele, other_allele, beta, se, pval, F_stat)],
      row.names = FALSE)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 3: 加载并标准化 trauma exposure GWAS (智能列名匹配)
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 3: Load trauma GWAS and standardize columns ─────\n")

# 把任意 GWAS 标准化到 (SNP, chr, pos, EA, NEA, EAF, beta, se, p, n)
standardize_gwas <- function(path, label = "GWAS") {
  cat(sprintf("\n[%s] Reading: %s\n", label, basename(path)))
  df <- tryCatch(fread(path, showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("  ⚠ Failed to read or empty table\n"))
    return(NULL)
  }
  cat(sprintf("  Raw shape: %d rows × %d columns\n", nrow(df), ncol(df)))
  cat(sprintf("  Raw column names: %s\n", paste(colnames(df), collapse = ", ")))

  cn <- tolower(colnames(df))
  pick <- function(patterns, must = TRUE) {
    for (pat in patterns) {
      hit <- which(grepl(pat, cn))
      if (length(hit) > 0) return(colnames(df)[hit[1]])
    }
    if (must) NA_character_ else NA_character_
  }

  col_snp  <- pick(c("^snp$", "^rsid$", "^rs_id$", "^markername$",
                     "^variant_id$", "^variantid$", "^id$"))
  col_chr  <- pick(c("^chr$", "^chrom$", "^chromosome$", "^#chrom$"))
  col_pos  <- pick(c("^pos$", "^bp$", "^position$", "^base_pair_location$",
                     "^pos_b37$", "^pos_b38$"))
  col_ea   <- pick(c("^a1$", "^ea$", "^effect_allele$", "^allele1$",
                     "^alt$", "^allele_b$", "^a_1$"))
  col_nea  <- pick(c("^a2$", "^nea$", "^other_allele$", "^non_effect_allele$",
                     "^allele0$", "^allele2$", "^ref$", "^allele_a$",
                     "^a_0$"))
  col_beta <- pick(c("^beta$", "^b$", "^effect$", "^estimate$",
                     "^beta_meta$", "^or$"))   # OR 后面对数处理
  col_se   <- pick(c("^se$", "^standard_error$", "^stderr$", "^sebeta$"))
  col_p    <- pick(c("^p$", "^pval$", "^p_value$", "^p\\.value$",
                     "^p_bolt_lmm$", "^p_meta$"))
  col_eaf  <- pick(c("^eaf$", "^a1freq$", "^maf$", "^freq$",
                     "^effect_allele_frequency$", "^a1_freq$"))
  col_n    <- pick(c("^n$", "^n_eff$", "^neff$", "^samplesize$",
                     "^sample_size$", "^n_total$"))

  cat("  Identified columns:\n")
  cat(sprintf("    SNP=%s, chr=%s, pos=%s, EA=%s, NEA=%s\n",
              col_snp, col_chr, col_pos, col_ea, col_nea))
  cat(sprintf("    beta=%s, se=%s, p=%s, eaf=%s, n=%s\n",
              col_beta, col_se, col_p, col_eaf, col_n))

  if (is.na(col_snp) || is.na(col_ea) || is.na(col_nea) ||
      is.na(col_beta) || is.na(col_se) || is.na(col_p)) {
    cat("  ⚠ Required columns missing; cannot proceed\n")
    return(NULL)
  }

  out <- data.table(
    SNP = as.character(df[[col_snp]]),
    chr = if (!is.na(col_chr)) df[[col_chr]] else NA,
    pos = if (!is.na(col_pos)) df[[col_pos]] else NA,
    EA  = toupper(as.character(df[[col_ea]])),
    NEA = toupper(as.character(df[[col_nea]])),
    EAF = if (!is.na(col_eaf)) as.numeric(df[[col_eaf]]) else NA_real_,
    beta_raw = as.numeric(df[[col_beta]]),
    se = as.numeric(df[[col_se]]),
    p = as.numeric(df[[col_p]]),
    n_obs = if (!is.na(col_n)) as.numeric(df[[col_n]]) else NA_real_
  )

  # OR → log(OR) 检测: 如果 col_beta 列名是 OR 或 odds_ratio, 转 log
  if (grepl("^or$|odds.?ratio", col_beta, ignore.case = TRUE)) {
    cat("  Note: beta column appears to be OR; taking log()\n")
    out$beta <- log(out$beta_raw)
  } else {
    out$beta <- out$beta_raw
  }
  out$beta_raw <- NULL

  # 过滤无效行
  before <- nrow(out)
  out <- out[!is.na(SNP) & !is.na(beta) & !is.na(se) & !is.na(p) &
             nchar(EA) >= 1 & nchar(NEA) >= 1]
  cat(sprintf("  After cleaning: %d rows (filtered %d missing/abnormal rows)\n",
              nrow(out), before - nrow(out)))
  out
}

trauma_gwas <- list()
trauma_gwas$Warrier_2021_childhood_maltreatment <- standardize_gwas(
  warrier_path, "Warrier 2021")
if (is.null(trauma_gwas$Warrier_2021_childhood_maltreatment)) {
  stop("Failed to parse the Warrier 2021 GWAS file.")
}

# Sensitivity GWAS (可选)
for (f in sensitivity_candidates) {
  bn <- basename(f)
  label <- gsub("\\.(gz|tsv|txt|tsv.gz|txt.gz)$", "", bn, ignore.case = TRUE)
  res <- standardize_gwas(f, label)
  if (!is.null(res)) trauma_gwas[[label]] <- res
}

cat(sprintf("\nLoaded %d trauma-exposure GWAS in total\n\n", length(trauma_gwas)))

# ─────────────────────────────────────────────────────────────────────────
# Step 4: SNP harmonization + Wald ratio MR
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 4: Extract instrument effects in trauma GWAS; run Wald ratios ─────\n")

# Wald ratio: β_outcome / β_exposure
# SE: 用 delta method 一阶近似 SE = sqrt((se_out^2 + (β_out/β_exp)^2 * se_exp^2) / β_exp^2)
# 简化版 SE ≈ se_out / |β_exp| (假设 β_exp 测量误差小, 适合 strong instrument)
wald_ratio <- function(b_exp, se_exp, b_out, se_out) {
  ratio <- b_out / b_exp
  # 一阶 delta
  ratio_se <- sqrt((se_out^2 / b_exp^2) +
                   (b_out^2 * se_exp^2 / b_exp^4))
  z <- ratio / ratio_se
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  list(beta = ratio, se = ratio_se, z = z, p = p)
}

results_all <- list()
for (g_name in names(trauma_gwas)) {
  cat(sprintf("\n[%s]\n", g_name))
  trauma <- trauma_gwas[[g_name]]
  out_rows <- list()
  for (i in seq_len(nrow(inst))) {
    prot <- inst$protein[i]
    snp  <- inst$SNP[i]
    ea_e <- toupper(inst$effect_allele[i])
    nea_e <- toupper(inst$other_allele[i])
    b_e  <- inst$beta[i]; se_e <- inst$se[i]; p_e <- inst$pval[i]

    hit <- trauma[SNP == snp]
    if (nrow(hit) == 0) {
      out_rows[[i]] <- data.table(
        candidate = prot, cis_pQTL_SNP = snp,
        EA_pQTL = ea_e, NEA_pQTL = nea_e,
        beta_pQTL = b_e, se_pQTL = se_e, p_pQTL = p_e,
        EA_trauma = NA_character_, NEA_trauma = NA_character_,
        beta_trauma_raw = NA_real_, se_trauma = NA_real_, p_trauma = NA_real_,
        flipped = NA, harmonization = "SNP_not_in_trauma_GWAS",
        wald_beta = NA_real_, wald_se = NA_real_, wald_p = NA_real_
      )
      next
    }
    if (nrow(hit) > 1) hit <- hit[1]

    ea_t <- hit$EA; nea_t <- hit$NEA
    b_t <- hit$beta; se_t <- hit$se; p_t <- hit$p

    # Harmonization
    flip <- NA; status <- ""
    if (ea_e == ea_t && nea_e == nea_t) {
      flip <- FALSE
      status <- "aligned"
    } else if (ea_e == nea_t && nea_e == ea_t) {
      flip <- TRUE
      b_t <- -b_t
      status <- "flipped"
    } else {
      # alleles 不一致 (比如 ambiguous palindromic)
      flip <- NA
      status <- sprintf("allele_mismatch(pQTL %s/%s vs trauma %s/%s)",
                        ea_e, nea_e, ea_t, nea_t)
    }

    if (!is.na(flip)) {
      mr <- wald_ratio(b_e, se_e, b_t, se_t)
      out_rows[[i]] <- data.table(
        candidate = prot, cis_pQTL_SNP = snp,
        EA_pQTL = ea_e, NEA_pQTL = nea_e,
        beta_pQTL = b_e, se_pQTL = se_e, p_pQTL = p_e,
        EA_trauma = ea_t, NEA_trauma = nea_t,
        beta_trauma_raw = if (flip) -b_t else b_t,
        se_trauma = se_t, p_trauma = p_t,
        flipped = flip, harmonization = status,
        wald_beta = mr$beta, wald_se = mr$se, wald_p = mr$p
      )
    } else {
      out_rows[[i]] <- data.table(
        candidate = prot, cis_pQTL_SNP = snp,
        EA_pQTL = ea_e, NEA_pQTL = nea_e,
        beta_pQTL = b_e, se_pQTL = se_e, p_pQTL = p_e,
        EA_trauma = ea_t, NEA_trauma = nea_t,
        beta_trauma_raw = b_t, se_trauma = se_t, p_trauma = p_t,
        flipped = NA, harmonization = status,
        wald_beta = NA_real_, wald_se = NA_real_, wald_p = NA_real_
      )
    }
  }
  res_dt <- rbindlist(out_rows, fill = TRUE)
  res_dt[, trauma_gwas := g_name]
  results_all[[g_name]] <- res_dt

  # 简洁打印
  cat(sprintf("  %d/%d instruments successfully harmonized\n",
              sum(!is.na(res_dt$wald_p)), nrow(res_dt)))
  print(res_dt[, .(candidate, cis_pQTL_SNP,
                   wald_beta = signif(wald_beta, 3),
                   wald_se = signif(wald_se, 3),
                   wald_p = signif(wald_p, 3),
                   harmonization)],
        row.names = FALSE)
}

results_combined <- rbindlist(results_all, fill = TRUE)

# ─────────────────────────────────────────────────────────────────────────
# Step 5: 跟 PTSD 主分析比较, 给 verdict
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 5: Compare with primary PTSD MR, assign verdict ─────\n")

# 试图从已生成结果读 PTSD MR
ptsd_mr_files <- c(
  file.path(result_dir, "mr_all_outcomes_extended.rds"),
  file.path(result_dir, "mr_main.rds"),
  file.path(result_dir, "mr_results.rds")
)
ptsd_mr_path <- ptsd_mr_files[file.exists(ptsd_mr_files)][1]
ptsd_mr <- if (!is.na(ptsd_mr_path) && nzchar(ptsd_mr_path)) {
  cat(sprintf("Reading PTSD MR results from %s\n", ptsd_mr_path))
  tryCatch(readRDS(ptsd_mr_path), error = function(e) NULL)
} else {
  NULL
}

# 尝试提取 PTSD-specific MR 结果 (10 candidates)
ptsd_summary <- data.table(candidate = candidates,
                            ptsd_beta = NA_real_, ptsd_p = NA_real_)

if (!is.null(ptsd_mr)) {
  ptsd_dt <- if (is.data.frame(ptsd_mr)) as.data.table(ptsd_mr) else NULL
  if (!is.null(ptsd_dt)) {
    cn <- colnames(ptsd_dt)
    out_col <- grep("outcome|trait|disease", cn, value = TRUE, ignore.case = TRUE)[1]
    pro_col <- grep("^protein$|^exposure$|^gene$", cn, value = TRUE,
                    ignore.case = TRUE)[1]
    b_col   <- grep("^beta$|^b$|^effect$", cn, value = TRUE, ignore.case = TRUE)[1]
    p_col   <- grep("^p$|^pval$|^p_value$|^p\\.value$", cn, value = TRUE,
                    ignore.case = TRUE)[1]
    if (!is.na(out_col) && !is.na(pro_col)) {
      ptsd_mask <- grepl("ptsd|stress.*disorder", as.character(ptsd_dt[[out_col]]),
                         ignore.case = TRUE)
      sub <- ptsd_dt[ptsd_mask & ptsd_dt[[pro_col]] %in% candidates]
      if (nrow(sub) > 0 && !is.na(b_col)) {
        # 每蛋白取一行 (P 最小)
        sub <- sub[order(get(p_col))][, .SD[1], by = pro_col]
        for (g in candidates) {
          row <- sub[get(pro_col) == g]
          if (nrow(row) > 0) {
            ptsd_summary[candidate == g, ptsd_beta := as.numeric(row[[b_col]])]
            if (!is.na(p_col)) {
              ptsd_summary[candidate == g, ptsd_p := as.numeric(row[[p_col]])]
            }
          }
        }
      }
    }
  }
}
cat("PTSD MR (from primary analysis):\n")
print(ptsd_summary, row.names = FALSE)
cat("\n")

# Verdict logic
classify <- function(ptsd_b, ptsd_p, trauma_b, trauma_p,
                     alpha_ptsd = 0.05, alpha_trauma = 0.05) {
  if (is.na(trauma_b) || is.na(trauma_p)) return("unable_to_test")
  if (is.na(ptsd_p) || ptsd_p > alpha_ptsd) return("ptsd_not_significant")
  # PTSD 显著 → 看 trauma
  if (trauma_p > alpha_trauma) {
    return("post_pathology")  # ✅ 期望: 蛋白参与创伤后病理, 而非创伤暴露倾向
  }
  # trauma 也显著
  if (!is.na(ptsd_b)) {
    if (sign(ptsd_b) == sign(trauma_b)) {
      return("behavioral_confound")  # ⚠ 同向: 蛋白可能影响创伤暴露倾向
    } else {
      return("opposite_direction")    # 异向: 复杂 pleiotropy
    }
  }
  "trauma_significant_unclassified"
}

# 把 PTSD 信息合并到 results_combined
results_combined <- merge(results_combined, ptsd_summary,
                          by = "candidate", all.x = TRUE)
results_combined[, verdict := mapply(classify,
                                      ptsd_beta, ptsd_p,
                                      wald_beta, wald_p)]

# ─────────────────────────────────────────────────────────────────────────
# Step 6: 输出
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 6: Write outputs ─────\n")

out_csv <- file.path(result_dir, "trauma_negative_control_MR.csv")
fwrite(results_combined, out_csv)
cat(sprintf("✓ Written: %s\n", out_csv))

out_rds <- file.path(result_dir, "trauma_negative_control_MR.rds")
saveRDS(list(results = results_combined,
             instruments = inst,
             trauma_gwas_files = sapply(names(trauma_gwas),
                                         function(x) basename(switch(
                                           x,
                                           "Warrier_2021_childhood_maltreatment" = warrier_path,
                                           ""))),
             ptsd_summary = ptsd_summary),
        out_rds)
cat(sprintf("✓ Written: %s\n", out_rds))

# ─────────────────────────────────────────────────────────────────────────
# Step 7: 解读报告
# ─────────────────────────────────────────────────────────────────────────
report <- character(0)
addln <- function(...) report <<- c(report, sprintf(...))

addln("══════════════════════════════════════════════════════════════════")
addln("Trauma Exposure Negative Control MR — 解读")
addln("══════════════════════════════════════════════════════════════════")
addln("")
addln("假说测试逻辑:")
addln("  PTSD = 经历创伤 + 创伤后产生病理反应")
addln("  - 蛋白对 PTSD 显著 + 对 trauma 不显著  →  post-trauma pathology ✅")
addln("    (蛋白参与创伤后应激病理过程, 不是行为学风险)")
addln("  - 蛋白对 PTSD 显著 + 对 trauma 同向显著  →  behavioral confound ⚠")
addln("    (蛋白可能影响创伤暴露倾向, 是行为学混杂因素)")
addln("  - 蛋白对 PTSD 显著 + 对 trauma 反向显著  →  opposite direction")
addln("    (复杂 pleiotropy, 需进一步审视)")
addln("")

for (g_name in unique(results_combined$trauma_gwas)) {
  addln("──────────────────────────────────────────────────────────────────")
  addln("Trauma GWAS: %s", g_name)
  addln("──────────────────────────────────────────────────────────────────")

  sub <- results_combined[trauma_gwas == g_name]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    addln("")
    addln("  %s (cis-pQTL %s, %s>%s)",
          r$candidate, r$cis_pQTL_SNP, r$EA_pQTL, r$NEA_pQTL)
    if (is.na(r$wald_p)) {
      addln("    [无法 MR: %s]", r$harmonization)
      next
    }
    addln("    PTSD MR (主分析):     β=%s, P=%s",
          if (is.na(r$ptsd_beta)) "—" else signif(r$ptsd_beta, 3),
          if (is.na(r$ptsd_p)) "—" else signif(r$ptsd_p, 3))
    addln("    Trauma exposure MR:   β=%s, SE=%s, P=%s",
          signif(r$wald_beta, 3),
          signif(r$wald_se, 3),
          signif(r$wald_p, 3))
    addln("    Verdict:              %s", r$verdict)
  }
}

addln("")
addln("══════════════════════════════════════════════════════════════════")
addln("Summary 统计")
addln("══════════════════════════════════════════════════════════════════")
verdict_tbl <- as.data.frame(table(results_combined$verdict))
for (i in seq_len(nrow(verdict_tbl))) {
  addln("  %-32s %d", verdict_tbl[i, 1], verdict_tbl[i, 2])
}

report_path <- file.path(result_dir, "trauma_negative_control_summary.txt")
writeLines(report, report_path)
cat(sprintf("✓ Written: %s\n\n", report_path))

# 同时控制台打印
cat(paste(report, collapse = "\n"), "\n\n")

cat("Next steps:\n")
cat("  1) Inspect verdicts for PTSD-significant proteins in trauma_negative_control_MR.csv\n")
cat("  2) Proteins flagged 'post_pathology' can be highlighted in the manuscript:\n")
cat("     \"Importantly, [protein] showed no significant effect on childhood\n")
cat("      maltreatment (Warrier 2021, P=...), supporting its role in post-\n")
cat("      trauma stress pathology rather than behavioral exposure risk.\"\n")
cat("  3) Proteins flagged 'behavioral_confound' should be discussed in Limitations\n")
