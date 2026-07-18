###############################################################################
# 17_conditional_chr15_COJO.R
#
# GCTA-COJO conditional Mendelian randomization at the chromosome 15:90.87 Mb
# FES/FURIN locus. For each of two outcomes (PTSD freeze 3 case-control and the
# Warrier 2021 childhood maltreatment GWAS), the effect at one protein's
# cis-pQTL is tested after conditioning on the other protein's cis-pQTL.
#
# Implementation note: the chromosome 15 region is extracted from the 1000G EUR
# reference panel in base R via readBin/writeBin, so the 11 GB bfile is never
# fully loaded into memory. GCTA is then run on the resulting ~10 MB sub-bfile.
#
# Inputs (from 00_config.R):
#   pqtl_path, ld_ref_prefix, ptsd_ma_path, validation_dir, gcta_path
#
# Outputs:
#   results/conditional_chr15_MR.csv
#   results/conditional_chr15_MR.rds
#   results/cojo_chr15/  (intermediate GCTA outputs)
#
# PLINK BED v1 layout (SNP-major):
#   3 magic bytes (0x6c 0x1b 0x01), then ceil(n_ind/4) bytes per SNP.
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
suppressMessages({ library(data.table) })

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   17. chr15 GCTA-COJO Conditional Analysis                    ║\n")
cat("║       (region extracted via base R; GCTA runs on small bfile)   ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# Step 0
cat("───── Step 0: Check tools and data ─────\n")
required <- list(
  GCTA = gcta_path,
  `LD ref bed` = paste0(ld_ref_prefix, ".bed"),
  `LD ref bim` = paste0(ld_ref_prefix, ".bim"),
  `LD ref fam` = paste0(ld_ref_prefix, ".fam"),
  `PTSD .ma` = ptsd_ma_path
)
for (k in names(required)) {
  ok <- file.exists(required[[k]])
  cat(sprintf("  %s %-15s %s\n", if (ok) "✓" else "❌", k, required[[k]]))
  if (!ok) stop(sprintf("Missing required config variable: %s", required[[k]]))
}
trauma_raw <- file.path(validation_dir, "trauma_exposure",
                        "Retro_prospective_meta_childhoodmaltreatment.txt")
if (!file.exists(trauma_raw)) {
  d <- file.path(validation_dir, "trauma_exposure")
  hits <- list.files(d, pattern = "\\.txt$", ignore.case = TRUE,
                     full.names = TRUE, recursive = TRUE)
  if (length(hits) > 0) trauma_raw <- hits[1]
}
if (!file.exists(trauma_raw)) stop("Cannot locate the trauma exposure GWAS file (Warrier 2021).")
cat(sprintf("  ✓ Trauma raw GWAS: %s (%.0f MB)\n",
            basename(trauma_raw), file.info(trauma_raw)$size / 1024^2))
cojo_dir <- file.path(result_dir, "cojo_chr15")
dir.create(cojo_dir, showWarnings = FALSE, recursive = TRUE)

# Step 1
cat("\n───── Step 1: Locate FES + FURIN cis-pQTLs ─────\n")
pqtl <- as.data.table(readRDS(pqtl_path))
fes_inst   <- pqtl[protein == "FES"  ][order(pval)][1]
furin_inst <- pqtl[protein == "FURIN"][order(pval)][1]
cond_snps <- list(
  FES   = list(snp = fes_inst$SNP, chr = 15, pos = fes_inst$pos_hg38,
               EA = fes_inst$effect_allele,
               beta = fes_inst$beta, se = fes_inst$se, pval = fes_inst$pval),
  FURIN = list(snp = furin_inst$SNP, chr = 15, pos = furin_inst$pos_hg38,
               EA = furin_inst$effect_allele,
               beta = furin_inst$beta, se = furin_inst$se, pval = furin_inst$pval)
)
cat(sprintf("FES top: %s @ %d, EA=%s\n",
            cond_snps$FES$snp, cond_snps$FES$pos, cond_snps$FES$EA))
cat(sprintf("FURIN top: %s @ %d, EA=%s\n",
            cond_snps$FURIN$snp, cond_snps$FURIN$pos, cond_snps$FURIN$EA))
region_lo <- min(cond_snps$FES$pos, cond_snps$FURIN$pos) - 1000000
region_hi <- max(cond_snps$FES$pos, cond_snps$FURIN$pos) + 1000000
cat(sprintf("\nchr15 region: %d-%d (hg38, ±1Mb)\n", region_lo, region_hi))

# Step 2: 流式抽 chr15 子 bfile
cat("\n───── Step 2: Stream-extract chr15 region sub-bfile ─────\n")
small_prefix <- file.path(cojo_dir, "1000G_EUR_chr15_region")
small_bed <- paste0(small_prefix, ".bed")
small_bim <- paste0(small_prefix, ".bim")
small_fam <- paste0(small_prefix, ".fam")

if (file.exists(small_bed) && file.exists(small_bim) &&
    file.exists(small_fam) && file.info(small_bed)$size > 1000) {
  cat(sprintf("  Cache hit: %s (%.1f MB)\n",
              basename(small_bed), file.info(small_bed)$size / 1024^2))
} else {
  cat("  Reading .bim to locate chr15 region indices...\n")
  t0 <- Sys.time()
  bim_full <- fread(paste0(ld_ref_prefix, ".bim"),
                    col.names = c("chr", "id", "cm", "bp", "a1", "a2"),
                    showProgress = FALSE)
  cat(sprintf("    .bim full: %d rows, %.1f s\n", nrow(bim_full),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  bim_full[, row_idx := .I]
  mask <- (bim_full$chr == "15" | bim_full$chr == 15) &
          bim_full$bp >= region_lo & bim_full$bp <= region_hi
  bim_keep <- bim_full[mask]
  cat(sprintf("    chr15 ±1Mb region: %d SNPs (indices %d - %d, contiguous: %s)\n",
              nrow(bim_keep), min(bim_keep$row_idx), max(bim_keep$row_idx),
              if (all(diff(bim_keep$row_idx) == 1)) "连续" else "不连续"))
  if (nrow(bim_keep) == 0) stop("No SNPs found in the chr15 region after LD reference filtering.")
  for (p in names(cond_snps)) {
    if (nrow(bim_keep[id == cond_snps[[p]]$snp]) == 0)
      stop(sprintf("Conditioning SNP %s is not within the chr15 region.", cond_snps[[p]]$snp))
  }

  fam_lines <- readLines(paste0(ld_ref_prefix, ".fam"))
  n_ind <- length(fam_lines)
  n_bytes_per_snp <- ceiling(n_ind / 4)
  cat(sprintf("    .fam: %d samples, %d bytes per SNP\n",
              n_ind, n_bytes_per_snp))

  rm(bim_full); gc()  # 释放 .bim 全表内存

  cat("\n  Stream-reading .bed file...\n")
  t0 <- Sys.time()
  con_in <- file(paste0(ld_ref_prefix, ".bed"), "rb")
  con_out <- file(small_bed, "wb")
  magic <- readBin(con_in, "raw", n = 3L)
  expected_magic <- as.raw(c(0x6c, 0x1b, 0x01))
  if (!identical(magic, expected_magic)) {
    close(con_in); close(con_out)
    stop(sprintf("PLINK BED magic bytes do not match: 0x%02x 0x%02x 0x%02x",
                 as.integer(magic[1]), as.integer(magic[2]),
                 as.integer(magic[3])))
  }
  writeBin(expected_magic, con_out)

  indices <- bim_keep$row_idx
  is_contiguous <- all(diff(indices) == 1)
  if (is_contiguous) {
    cat(sprintf("    Indices contiguous; block copy (%d SNP × %d bytes = %.1f MB)\n",
                length(indices), n_bytes_per_snp,
                length(indices) * n_bytes_per_snp / 1024^2))
    start_idx <- indices[1]
    file_offset <- 3 + (start_idx - 1) * n_bytes_per_snp
    seek(con_in, file_offset)
    chunk <- readBin(con_in, "raw",
                     n = length(indices) * n_bytes_per_snp)
    writeBin(chunk, con_out)
  } else {
    cat(sprintf("    Indices non-contiguous; per-SNP extraction (%d SNPs)\n", length(indices)))
    pb <- txtProgressBar(min = 0, max = length(indices), style = 3)
    for (i in seq_along(indices)) {
      seek(con_in, 3 + (indices[i] - 1) * n_bytes_per_snp)
      writeBin(readBin(con_in, "raw", n = n_bytes_per_snp), con_out)
      if (i %% 5000 == 0) setTxtProgressBar(pb, i)
    }
    setTxtProgressBar(pb, length(indices)); close(pb); cat("\n")
  }
  close(con_in); close(con_out)

  bim_out <- bim_keep[, .(chr, id, cm, bp, a1, a2)]
  fwrite(bim_out, small_bim, sep = "\t", col.names = FALSE, quote = FALSE)
  file.copy(paste0(ld_ref_prefix, ".fam"), small_fam, overwrite = TRUE)

  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("\n  ✓ Done (%.1f s)\n", dt))
  cat(sprintf("    %s (%.1f MB)\n", basename(small_bed),
              file.info(small_bed)$size / 1024^2))
}

bim_small <- fread(small_bim,
                   col.names = c("chr", "id", "cm", "bp", "a1", "a2"))
cat(sprintf("\n  Small bfile: %d SNPs, %d samples\n", nrow(bim_small),
            length(readLines(small_fam))))
for (p in names(cond_snps)) {
  hit <- bim_small[id == cond_snps[[p]]$snp]
  if (nrow(hit) == 0) stop(sprintf("SNP %s not found in the extracted chr15 sub-bfile.", cond_snps[[p]]$snp))
  cat(sprintf("  ✓ %s SNP %s (bp=%d, %s/%s)\n",
              p, cond_snps[[p]]$snp, hit$bp, hit$a1, hit$a2))
}

# Step 3: trauma .ma (跟 v5 一致)
cat("\n───── Step 3: Prepare trauma .ma file ─────\n")
trauma_ma_path <- file.path(cojo_dir, "trauma_warrier_chr15.ma")
if (file.exists(trauma_ma_path) && file.info(trauma_ma_path)$size > 1000) {
  cat(sprintf("  Cache hit: %s\n", trauma_ma_path))
} else {
  cat("  Reading PTSD .ma to borrow allele frequencies...\n")
  ptsd_ma <- fread(ptsd_ma_path, showProgress = FALSE)
  ptsd_region <- ptsd_ma[SNP %in% bim_small$id]
  freq_lookup <- ptsd_region[, .(SNP, ptsd_A1 = toupper(A1),
                                  ptsd_A2 = toupper(A2), ptsd_freq = freq)]
  rm(ptsd_ma); gc()
  cat("  Reading Warrier chr15 region...\n")
  trauma_full <- fread(trauma_raw, showProgress = FALSE)
  trauma_chr15 <- trauma_full[CHR == 15 & SNP %in% bim_small$id]
  rm(trauma_full); gc()
  m <- merge(trauma_chr15, freq_lookup, by = "SNP")
  m[, freq := fifelse(toupper(A1) == ptsd_A1, ptsd_freq,
              fifelse(toupper(A1) == ptsd_A2, 1 - ptsd_freq,
                      NA_real_))]
  m <- m[!is.na(freq) & freq > 0 & freq < 1]
  ma_out <- m[, .(SNP, A1 = toupper(A1), A2 = toupper(A2),
                   freq, b = BETA, se = SE, p = P, N = 185414)]
  fwrite(ma_out, trauma_ma_path, sep = " ", quote = FALSE)
  cat(sprintf("  ✓ Written: %d rows\n", nrow(ma_out)))
}

# Step 4: GCTA-COJO × 4
cat("\n───── Step 4: GCTA-COJO conditional × 4 (small bfile) ─────\n")
run_cojo <- function(ma_path, cond_snp, out_prefix, label) {
  cond_list <- file.path(cojo_dir,
                          sprintf("cond_%s.snplist", basename(out_prefix)))
  writeLines(cond_snp, cond_list)
  args <- c("--bfile", shQuote(small_prefix),
            "--chr", "15",
            "--cojo-file", shQuote(ma_path),
            "--cojo-cond", shQuote(cond_list),
            "--out", shQuote(out_prefix),
            "--thread-num", "4")
  cat(sprintf("\n  [%s]\n    cond on: %s\n", label, cond_snp))
  t0 <- Sys.time()
  log <- system2(gcta_path, args = args, stdout = TRUE, stderr = TRUE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out_file <- paste0(out_prefix, ".cma.cojo")
  if (file.exists(out_file)) {
    cat(sprintf("    ✓ Done (%.1f s)\n", dt))
    return(out_file)
  } else {
    cat(sprintf("    ❌ No .cma.cojo generated; tail of log:\n"))
    for (ln in tail(log, 20)) cat(sprintf("      %s\n", ln))
    return(NA_character_)
  }
}
cojo_results <- list()
cojo_results$PTSD_condFES <- run_cojo(ptsd_ma_path, cond_snps$FES$snp,
  file.path(cojo_dir, "ptsd_cond_FES"), "PTSD | cond FES")
cojo_results$PTSD_condFURIN <- run_cojo(ptsd_ma_path, cond_snps$FURIN$snp,
  file.path(cojo_dir, "ptsd_cond_FURIN"), "PTSD | cond FURIN")
cojo_results$Trauma_condFES <- run_cojo(trauma_ma_path, cond_snps$FES$snp,
  file.path(cojo_dir, "trauma_cond_FES"), "Trauma | cond FES")
cojo_results$Trauma_condFURIN <- run_cojo(trauma_ma_path, cond_snps$FURIN$snp,
  file.path(cojo_dir, "trauma_cond_FURIN"), "Trauma | cond FURIN")

# Step 5-7: 解析 + MR + verdict + 报告
cat("\n───── Step 5: Parse results and conditional MR ─────\n")
parse_cma <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  fread(path)
}
wald_ratio <- function(b_exp, se_exp, b_out, se_out) {
  ratio <- b_out / b_exp
  ratio_se <- sqrt((se_out^2/b_exp^2) + (b_out^2 * se_exp^2/b_exp^4))
  z <- ratio / ratio_se
  list(beta = ratio, se = ratio_se,
       p = 2 * pnorm(abs(z), lower.tail = FALSE))
}
placeholder_row <- function(label, qp, qi, cp, note) {
  data.table(label = label, query_protein = qp, cond_on_protein = cp,
             query_SNP = qi$snp,
             pQTL_b = qi$beta, pQTL_se = qi$se, pQTL_p = qi$pval,
             uncond_outcome_b = NA_real_, uncond_outcome_se = NA_real_,
             uncond_outcome_p = NA_real_,
             cond_outcome_b = NA_real_, cond_outcome_se = NA_real_,
             cond_outcome_p = NA_real_,
             mr_uncond_beta = NA_real_, mr_uncond_se = NA_real_,
             mr_uncond_p = NA_real_,
             mr_cond_beta = NA_real_, mr_cond_se = NA_real_,
             mr_cond_p = NA_real_,
             flipped = NA, note = note)
}
run_conditional_mr <- function(cma_file, qp, qi, cp, ci, label) {
  cma <- parse_cma(cma_file)
  if (is.null(cma)) {
    cat(sprintf("\n  [%s] ⚠ COJO failed\n", label))
    return(placeholder_row(label, qp, qi, cp, "COJO_failed"))
  }
  cat(sprintf("\n  [%s] .cma.cojo rows: %d\n", label, nrow(cma)))
  hit <- cma[SNP == qi$snp]
  if (nrow(hit) == 0) {
    cat(sprintf("    ⚠ Query SNP %s was dropped by COJO\n", qi$snp))
    return(placeholder_row(label, qp, qi, cp, "query_SNP_dropped"))
  }
  cat(sprintf("    refA=%s | uncond b=%.4g, p=%.3g | cond bC=%.4g, pC=%.3g\n",
              hit$refA, hit$b, hit$p, hit$bC, hit$pC))
  flip <- FALSE
  cond_b <- hit$bC; cond_se <- hit$bC_se
  uncond_b <- hit$b; uncond_se <- hit$se
  if (toupper(hit$refA) != toupper(qi$EA)) {
    flip <- TRUE
    cond_b <- -cond_b; uncond_b <- -uncond_b
    cat(sprintf("    [flipped]\n"))
  }
  mr_u <- wald_ratio(qi$beta, qi$se, uncond_b, uncond_se)
  mr_c <- wald_ratio(qi$beta, qi$se, cond_b, cond_se)
  data.table(label = label, query_protein = qp, cond_on_protein = cp,
             query_SNP = qi$snp,
             pQTL_b = qi$beta, pQTL_se = qi$se, pQTL_p = qi$pval,
             uncond_outcome_b = uncond_b, uncond_outcome_se = uncond_se,
             uncond_outcome_p = hit$p,
             cond_outcome_b = cond_b, cond_outcome_se = cond_se,
             cond_outcome_p = hit$pC,
             mr_uncond_beta = signif(mr_u$beta, 4),
             mr_uncond_se = signif(mr_u$se, 4),
             mr_uncond_p = signif(mr_u$p, 3),
             mr_cond_beta = signif(mr_c$beta, 4),
             mr_cond_se = signif(mr_c$se, 4),
             mr_cond_p = signif(mr_c$p, 3),
             flipped = flip, note = "")
}
mr_table <- list()
mr_table$A <- run_conditional_mr(cojo_results$PTSD_condFES,
  "FURIN", cond_snps$FURIN, "FES", cond_snps$FES, "FURIN→PTSD | adj FES")
mr_table$B <- run_conditional_mr(cojo_results$PTSD_condFURIN,
  "FES", cond_snps$FES, "FURIN", cond_snps$FURIN, "FES→PTSD | adj FURIN")
mr_table$C <- run_conditional_mr(cojo_results$Trauma_condFES,
  "FURIN", cond_snps$FURIN, "FES", cond_snps$FES, "FURIN→Trauma | adj FES")
mr_table$D <- run_conditional_mr(cojo_results$Trauma_condFURIN,
  "FES", cond_snps$FES, "FURIN", cond_snps$FURIN, "FES→Trauma | adj FURIN")
mr_results <- rbindlist(mr_table, fill = TRUE)

cat("\n───── Step 6: Verdict ─────\n\n")
verdict_for_row <- function(uncond_p, cond_p, alpha = 0.05) {
  if (is.na(cond_p)) return("query_SNP_dropped_or_unable")
  if (is.na(uncond_p)) return("uncond_p_NA")
  if (uncond_p < alpha && cond_p < alpha) return("robust_to_conditioning")
  if (uncond_p < alpha && cond_p >= alpha) return("attenuated_by_LD")
  if (uncond_p >= alpha && cond_p < alpha) return("emerged_after_conditioning")
  if (uncond_p >= alpha && cond_p >= alpha) return("non_significant_either")
  "unclassified"
}
mr_results[, verdict := mapply(verdict_for_row, mr_uncond_p, mr_cond_p)]
print(mr_results[, .(label, mr_uncond_p, mr_cond_p, verdict, note)],
      row.names = FALSE)

out_csv <- file.path(result_dir, "conditional_chr15_MR.csv")
fwrite(mr_results, out_csv)
cat(sprintf("\n✓ Written: %s\n", out_csv))
saveRDS(list(mr_results = mr_results, cond_snps = cond_snps,
             cojo_files = cojo_results,
             region = list(chr = 15, lo = region_lo, hi = region_hi),
             small_bfile = small_prefix),
        file.path(result_dir, "conditional_chr15_MR.rds"))

# Step 7: 报告
report <- character(0)
addln <- function(...) report <<- c(report, sprintf(...))
addln("══════════════════════════════════════════════════════════════════")
addln("chr15 GCTA-COJO Conditional Analysis — 解读 (v6)")
addln("══════════════════════════════════════════════════════════════════")
addln("")
addln("方法: base R readBin 流式抽 chr15 子 bfile (避开 OOM), GCTA 跑小文件")
addln("Region: chr15:%s-%s (hg38, ±1Mb)",
       format(region_lo, big.mark=","), format(region_hi, big.mark=","))
addln("FES top:   %s (β=%.3f, P=%.2e)",
       cond_snps$FES$snp, cond_snps$FES$beta, cond_snps$FES$pval)
addln("FURIN top: %s (β=%.3f, P=%.2e)",
       cond_snps$FURIN$snp, cond_snps$FURIN$beta, cond_snps$FURIN$pval)
addln("")
addln("Verdict 逻辑:")
addln("  robust_to_conditioning   = uncond P<0.05 + cond P<0.05  ✅ 真信号")
addln("  attenuated_by_LD         = uncond P<0.05 + cond P≥0.05  ⚠ LD hitchhiker")
addln("  non_significant_either   = 两个都不显著")
addln("")
for (i in seq_len(nrow(mr_results))) {
  r <- mr_results[i]
  addln("──────────────────────────────────────────────────────────────────")
  addln("  %s", r$label)
  addln("──────────────────────────────────────────────────────────────────")
  addln("    Query SNP: %s (pQTL β=%.4f, SE=%.4f)",
        r$query_SNP, r$pQTL_b, r$pQTL_se)
  if (!is.na(r$mr_uncond_p)) {
    addln("    Outcome (uncond): β=%.4g, SE=%.4g, P=%.3g",
          r$uncond_outcome_b, r$uncond_outcome_se, r$uncond_outcome_p)
    addln("    Outcome (cond):   β=%.4g, SE=%.4g, P=%.3g",
          r$cond_outcome_b, r$cond_outcome_se, r$cond_outcome_p)
    addln("    Wald MR (uncond): β=%.4g, P=%.3g",
          r$mr_uncond_beta, r$mr_uncond_p)
    addln("    Wald MR (cond):   β=%.4g, P=%.3g",
          r$mr_cond_beta, r$mr_cond_p)
  } else {
    addln("    (COJO 失败或 SNP 被剔除)")
  }
  addln("    Verdict: %s", r$verdict)
  if (nzchar(r$note)) addln("    Note: %s", r$note)
  addln("")
}
addln("══════════════════════════════════════════════════════════════════")
addln("Manuscript 决策树")
addln("══════════════════════════════════════════════════════════════════")
addln("")
v_furin_ptsd <- mr_results[label == "FURIN→PTSD | adj FES", verdict]
addln("FURIN→PTSD | adj FES: %s", v_furin_ptsd)
addln("FES→PTSD | adj FURIN: %s",
      mr_results[label == "FES→PTSD | adj FURIN", verdict])
addln("FURIN→Trauma | adj FES: %s",
      mr_results[label == "FURIN→Trauma | adj FES", verdict])
addln("FES→Trauma | adj FURIN: %s",
      mr_results[label == "FES→Trauma | adj FURIN", verdict])
addln("")
if (length(v_furin_ptsd) > 0 && v_furin_ptsd == "robust_to_conditioning") {
  addln("✅ FURIN→PTSD robust to conditioning FES — 路线 B (但 manuscript 已写路线 A)")
} else if (length(v_furin_ptsd) > 0 && v_furin_ptsd == "attenuated_by_LD") {
  addln("⚠ FURIN→PTSD attenuated after adj FES — 已经在 manuscript 路线 A 中体现")
}

writeLines(report, file.path(result_dir, "conditional_chr15_summary.txt"))
cat(sprintf("\n✓ Written: %s\n\n",
            file.path(result_dir, "conditional_chr15_summary.txt")))
cat(paste(report, collapse = "\n"), "\n\n")
