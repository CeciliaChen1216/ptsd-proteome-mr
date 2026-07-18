###############################################################################
# 05_sensitivity_analyses.R
#
# Six sensitivity-analysis modules:
#   Part A: Steiger directionality test
#   Part B: MR-PRESSO outlier detection
#   Part C: Reverse MR (PTSD → protein)
#   Part D: Multivariable MR (FES/FURIN, CGREF1/KHK)
#   Part E: SMR + HEIDI
#   Part F: Coloc prior sensitivity
#
# Dependencies: dplyr, readr, data.table, MRPRESSO, MVMR, coloc
# Inputs:  pQTL instruments, primary MR results, coloc results, full protein
#          GWAS, PTSD GWAS (paths set in 00_config.R)
# Outputs: six CSV result tables written to result_dir
###############################################################################

library(data.table)
library(dplyr)
library(readr)

source("00_config.R")
# ═══════════════════════════════════════════════════════════════════════════════
# Global paths and parameters
# ═══════════════════════════════════════════════════════════════════════════════

mr_rds_path <- file.path(result_dir, "mr_all_outcomes.rds")
coloc_path  <- file.path(result_dir, "susie_coloc_full.rds")
pqtl_base   <- pqtl_base

n_pqtl <- 54219       # UKB-PPP样本量
n_ptsd <- 1222882     # 137,136 cases + 1,085,746 controls (Nievergelt 2024 Nat Genet; VCF ##nCase / ##nControl)

candidates <- c("AKT3", "CD40", "CGREF1", "FES", "FURIN",
                "SIRPA", "CD101", "KHK", "SNX18", "UBE2L6")

# 蛋白全GWAS文件夹映射 (已核实)
protein_folders <- c(
  AKT3   = "AKT3_Q9Y243_OID21197_v1_Oncology",
  CD40   = "CD40_P25942_OID20724_v1_Inflammation",
  CGREF1 = "CGREF1_Q99674_OID20152_v1_Cardiometabolic",
  FES    = "FES_P07332_OID21207_v1_Oncology",
  FURIN  = "FURIN_P09958_OID21514_v1_Oncology",
  SIRPA  = "SIRPA_P78324_OID20304_v1_Cardiometabolic",
  CD101  = "CD101_Q93033_OID31480_v1_Oncology_II",
  KHK    = "KHK_P50053_OID30241_v1_Cardiometabolic_II",
  SNX18  = "SNX18_Q96RF0_OID31476_v1_Oncology_II",
  UBE2L6 = "UBE2L6_O14933_OID30321_v1_Cardiometabolic_II"
)

# 共享辅助函数: 从CS文件获取每个credible set的lead SNP
get_lead_snps <- function(protein, cs_dir) {
  cs_files <- list.files(cs_dir, pattern = paste0("^", protein, "_pqtl_CS"),
                         full.names = TRUE)
  rbindlist(lapply(cs_files, function(f) {
    cs <- fread(f)
    cs[order(-pip)][1]
  }))
}


###############################################################################
# Part A: Steiger 方向性检验
###############################################################################

run_steiger <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part A: Steiger Directionality Test                  ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  pqtl       <- readRDS(pqtl_path)
  mr_results <- readRDS(mr_rds_path)
  
  # 合并pQTL和MR结果
  mr_ptsd <- mr_results %>%
    filter(grepl("PTSD|ptsd", outcome, ignore.case = TRUE),
           protein %in% candidates)
  
  dat <- mr_ptsd %>%
    inner_join(pqtl %>% select(protein, SNP, 
                                beta_exp = beta, se_exp = se, F_stat),
               by = c("protein", "SNP"))
  
  # 反推SNP-outcome关联: beta_out = mr_beta * beta_exp
  dat <- dat %>%
    mutate(beta_out = mr_beta * beta_exp,
           se_out   = mr_se * abs(beta_exp))
  
  # Steiger检验
  steiger_results <- dat %>%
    mutate(
      r2_exp = beta_exp^2 / (beta_exp^2 + n_pqtl * se_exp^2),
      r2_out = beta_out^2 / (beta_out^2 + n_ptsd * se_out^2),
      r2_ratio = r2_exp / r2_out,
      correct_direction = r2_exp > r2_out,
      # Steiger Z (Fisher z-transform)
      steiger_z = {
        z_e <- 0.5 * log((1 + sqrt(r2_exp)) / (1 - sqrt(r2_exp)))
        z_o <- 0.5 * log((1 + sqrt(r2_out)) / (1 - sqrt(r2_out)))
        (z_e - z_o) / sqrt(1/(n_pqtl - 3) + 1/(n_ptsd - 3))
      },
      steiger_p = pnorm(steiger_z, lower.tail = FALSE),
      interpretation = case_when(
        correct_direction & steiger_p < 0.05 ~ "✅ 正确方向 (蛋白→PTSD)",
        correct_direction ~ "— 方向正确但P不显著",
        TRUE ~ "⚠ 可能反向因果"
      )
    ) %>%
    select(protein, SNP, F_stat, r2_exposure = r2_exp, r2_outcome = r2_out,
           r2_ratio, correct_direction, steiger_z, steiger_p, interpretation) %>%
    arrange(steiger_p)
  
  write_csv(steiger_results, file.path(result_dir, "steiger_directionality_results.csv"))
  
  cat("Results:\n")
  steiger_results %>%
    select(protein, r2_exposure, r2_outcome, r2_ratio, steiger_p, interpretation) %>%
    mutate(r2_exposure = formatC(r2_exposure, format = "e", digits = 2),
           r2_outcome  = formatC(r2_outcome, format = "e", digits = 2),
           r2_ratio    = round(r2_ratio, 1),
           steiger_p   = formatC(steiger_p, format = "e", digits = 2)) %>%
    as.data.frame() %>%
    print(right = FALSE)
  
  cat(sprintf("\nDirectionally correct and significant: %d / %d\n",
              sum(steiger_results$correct_direction & steiger_results$steiger_p < 0.05),
              nrow(steiger_results)))
  
  invisible(steiger_results)
}


###############################################################################
# Part B: MR-PRESSO
###############################################################################

run_mr_presso <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part B: MR-PRESSO Outlier Detection                 ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  library(MRPRESSO)
  coloc_full <- readRDS(coloc_path)
  
  set.seed(42)
  presso_results <- list()
  presso_details <- list()
  
  for (protein in candidates) {
    cat(sprintf("══ %s ══\n", protein))
    
    leads <- get_lead_snps(protein, cs_dir)
    merged <- as.data.frame(coloc_full[[protein]]$merged_common)
    iv_data <- merged %>% filter(rsid %in% leads$rsid) %>%
      select(rsid, beta_pqtl, se_pqtl, beta_gwas, se_gwas)
    
    cat(sprintf("  Matched %d instruments\n", nrow(iv_data)))
    
    if (nrow(iv_data) < 4) {
      cat("  ⚠ Too few instruments; skipped\n\n")
      presso_results[[protein]] <- tibble(
        protein = protein, n_iv = nrow(iv_data),
        global_test_p = NA, raw_beta = NA, raw_p = NA,
        n_outliers = NA, corrected_p = NA, note = "IV数量不足")
      next
    }
    
    presso_input <- data.frame(bx = iv_data$beta_pqtl, by = iv_data$beta_gwas,
                               bxse = iv_data$se_pqtl, byse = iv_data$se_gwas)
    
    tryCatch({
      res <- mr_presso(BetaOutcome = "by", BetaExposure = "bx",
                       SdOutcome = "byse", SdExposure = "bxse",
                       OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                       data = presso_input, NbDistribution = 5000,
                       SignifThreshold = 0.05)
      
      main_res <- res$`Main MR results`
      global_p <- res$`MR-PRESSO results`$`Global Test`$Pvalue
      raw_row  <- main_res[main_res$`MR Analysis` == "Raw", ]
      corr_row <- main_res[main_res$`MR Analysis` == "Outlier-corrected", ]
      
      outlier_test <- res$`MR-PRESSO results`$`Outlier Test`
      n_outliers <- 0; outlier_snps <- ""
      if (!is.null(outlier_test) && is.data.frame(outlier_test)) {
        outlier_idx <- which(outlier_test$Pvalue < 0.05)
        n_outliers <- length(outlier_idx)
        outlier_snps <- paste(iv_data$rsid[outlier_idx], collapse = "; ")
      }
      
      presso_results[[protein]] <- tibble(
        protein = protein, n_iv = nrow(iv_data),
        global_test_p = global_p,
        raw_beta = as.numeric(raw_row$`Causal Estimate`),
        raw_p = as.numeric(raw_row$`P-value`),
        n_outliers = n_outliers, outlier_snps = outlier_snps,
        corrected_p = ifelse(nrow(corr_row) > 0, as.numeric(corr_row$`P-value`), NA),
        note = ifelse(global_p < 0.05, "⚠ 存在水平多效性", "✅ 无显著多效性"))
      
      presso_details[[protein]] <- res
      
      cat(sprintf("  Global P = %s, outliers = %d → %s\n\n",
                  formatC(global_p, format = "e", digits = 2), n_outliers,
                  ifelse(global_p < 0.05, "⚠", "✅")))
      
    }, error = function(e) {
      cat(sprintf("  ❌ Error: %s\n\n", e$message))
      presso_results[[protein]] <<- tibble(
        protein = protein, n_iv = nrow(iv_data),
        global_test_p = NA, raw_beta = NA, raw_p = NA,
        n_outliers = NA, corrected_p = NA, note = paste("Error:", e$message))
    })
  }
  
  presso_summary <- bind_rows(presso_results)
  write_csv(presso_summary, file.path(result_dir, "mr_presso_results.csv"))
  saveRDS(presso_details, file.path(result_dir, "mr_presso_full.rds"))
  
  cat(sprintf("No pleiotropy: %d / %d\n",
              sum(presso_summary$global_test_p >= 0.05, na.rm = TRUE), nrow(presso_summary)))
  
  invisible(presso_summary)
}


###############################################################################
# Part C: Reverse MR (PTSD → 蛋白)
###############################################################################

run_reverse_mr <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part C: Reverse MR (PTSD → Protein)                 ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  # Step 1: 读取PTSD GWAS
  cat("Reading PTSD GWAS...\n")
  con <- gzfile(ptsd_path)
  header_lines <- readLines(con, n = 1000)
  close(con)
  skip_n <- max(grep("^##", header_lines))
  
  ptsd_gwas <- fread(ptsd_path, skip = skip_n, header = TRUE, sep = "\t")
  setnames(ptsd_gwas, c("CHR", "SNP", "POS", "A1", "A2", "FREQ",
                         "BETA", "SE", "P", "NCASE", "NCON", "NEFF", "NTOT"))
  
  ptsd_gws <- ptsd_gwas[P < 5e-8]
  cat(sprintf("  GWS SNPs: %d\n", nrow(ptsd_gws)))
  rm(ptsd_gwas); gc()
  
  # Step 2: Distance-based pruning
  cat("Distance pruning (500kb)...\n")
  ptsd_gws <- ptsd_gws[order(P)]
  kept <- character(0)
  used_regions <- data.table(CHR = integer(), POS = integer())
  
  for (i in seq_len(nrow(ptsd_gws))) {
    row <- ptsd_gws[i]
    if (nrow(used_regions) == 0 ||
        !any(used_regions$CHR == row$CHR & abs(used_regions$POS - row$POS) < 500000)) {
      kept <- c(kept, row$SNP)
      used_regions <- rbind(used_regions, data.table(CHR = row$CHR, POS = row$POS))
    }
  }
  
  ptsd_instruments <- ptsd_gws[SNP %in% kept]
  cat(sprintf("  Independent instruments: %d\n\n", nrow(ptsd_instruments)))
  
  # Step 3: Reverse MR for each protein
  reverse_mr_results <- list()
  
  for (protein in candidates) {
    cat(sprintf("══ %s ══\n", protein))
    
    folder <- protein_folders[protein]
    inner_dir <- file.path(pqtl_base, folder, folder)
    
    if (!dir.exists(inner_dir)) {
      cat(sprintf("  ❌ Directory not found: %s\n\n", inner_dir))
      reverse_mr_results[[protein]] <- tibble(
        protein = protein, n_iv = NA, method = "IVW",
        beta = NA, se = NA, pval = NA, note = "文件夹不存在")
      next
    }
    
    needed_chrs <- unique(ptsd_instruments$CHR)
    pqtl_data <- rbindlist(lapply(needed_chrs, function(chr) {
      chr_files <- list.files(inner_dir, pattern = sprintf("discovery_chr%d_", chr),
                              full.names = TRUE)
      if (length(chr_files) == 0) return(NULL)
      dt <- fread(chr_files[1], select = c("CHROM", "GENPOS", "ID", "ALLELE0",
                                            "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"))
      id_parts <- tstrsplit(dt$ID, ":", fixed = TRUE)
      dt[, pos_hg19 := as.integer(id_parts[[2]])]
      dt
    }))
    
    if (nrow(pqtl_data) == 0) { cat("  ⚠ No data\n\n"); next }
    
    merged <- merge(
      ptsd_instruments[, .(SNP, CHR, POS, A1, A2, FREQ,
                           beta_exp = BETA, se_exp = SE, p_exp = P)],
      pqtl_data[, .(CHR = CHROM, POS = pos_hg19,
                    allele0 = ALLELE0, allele1 = ALLELE1,
                    beta_out = BETA, se_out = SE)],
      by = c("CHR", "POS"), all = FALSE)
    
    # Harmonization
    merged[, harmonized := FALSE]
    merged[A1 == allele1 & A2 == allele0, harmonized := TRUE]
    merged[A1 == allele0 & A2 == allele1, `:=`(beta_out = -beta_out, harmonized = TRUE)]
    merged <- merged[harmonized == TRUE]
    
    cat(sprintf("  Harmonized: %d IVs\n", nrow(merged)))
    
    if (nrow(merged) < 3) {
      cat("  ⚠ Insufficient\n\n")
      reverse_mr_results[[protein]] <- tibble(
        protein = protein, n_iv = nrow(merged), method = "IVW",
        beta = NA, se = NA, pval = NA, note = "instruments不足")
      next
    }
    
    # IVW
    wald_ratios <- merged$beta_out / merged$beta_exp
    wald_ses    <- abs(merged$se_out / merged$beta_exp)
    weights     <- 1 / wald_ses^2
    ivw_beta    <- sum(weights * wald_ratios) / sum(weights)
    ivw_se      <- sqrt(1 / sum(weights))
    ivw_p       <- 2 * pnorm(abs(ivw_beta / ivw_se), lower.tail = FALSE)
    
    # Weighted median
    sorted_idx  <- order(wald_ratios)
    cum_weights <- cumsum(weights[sorted_idx]) / sum(weights)
    median_idx  <- min(which(cum_weights >= 0.5))
    wm_beta     <- wald_ratios[sorted_idx[median_idx]]
    set.seed(42)
    boot_est <- replicate(1000, {
      bi <- sample(seq_len(nrow(merged)), replace = TRUE)
      br <- wald_ratios[bi]; bw <- weights[bi]
      si <- order(br); cw <- cumsum(bw[si]) / sum(bw)
      br[si[min(which(cw >= 0.5))]]
    })
    wm_se <- sd(boot_est)
    wm_p  <- 2 * pnorm(abs(wm_beta / wm_se), lower.tail = FALSE)
    
    sig_label <- ifelse(ivw_p < 0.05, "⚠ 反向显著!", "✅ 反向不显著")
    cat(sprintf("  IVW: beta=%.4f, P=%s | WM: P=%s → %s\n\n",
                ivw_beta, formatC(ivw_p, format = "e", digits = 2),
                formatC(wm_p, format = "e", digits = 2), sig_label))
    
    reverse_mr_results[[protein]] <- bind_rows(
      tibble(protein = protein, n_iv = nrow(merged), method = "IVW",
             beta = ivw_beta, se = ivw_se, pval = ivw_p, note = sig_label),
      tibble(protein = protein, n_iv = nrow(merged), method = "Weighted Median",
             beta = wm_beta, se = wm_se, pval = wm_p, note = ""))
  }
  
  reverse_summary <- bind_rows(reverse_mr_results)
  write_csv(reverse_summary, file.path(result_dir, "reverse_mr_results.csv"))
  
  cat("Reverse MR summary (IVW):\n")
  reverse_summary %>% filter(method == "IVW", !is.na(pval)) %>%
    select(protein, n_iv, beta, pval, note) %>%
    mutate(beta = round(beta, 4),
           pval = formatC(pval, format = "e", digits = 2)) %>%
    as.data.frame() %>% print(right = FALSE)
  
  invisible(reverse_summary)
}


###############################################################################
# Part D: Multivariable MR
###############################################################################

run_mvmr <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part D: Multivariable MR (MVMR)                     ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  library(MVMR)
  coloc_full <- readRDS(coloc_path)
  
  # 辅助: 从蛋白全GWAS查找指定hg19位置的效应
  lookup_pqtl_effects <- function(protein, chr, positions_hg19) {
    folder <- protein_folders[protein]
    inner_dir <- file.path(pqtl_base, folder, folder)
    chr_files <- list.files(inner_dir, pattern = sprintf("discovery_chr%d_", chr),
                            full.names = TRUE)
    if (length(chr_files) == 0) return(NULL)
    
    dt <- fread(chr_files[1], select = c("CHROM", "GENPOS", "ID", "ALLELE0",
                                          "ALLELE1", "A1FREQ", "BETA", "SE", "LOG10P"))
    id_parts <- tstrsplit(dt$ID, ":", fixed = TRUE)
    dt[, pos_hg19 := as.integer(id_parts[[2]])]
    dt_matched <- dt[pos_hg19 %in% positions_hg19]
    
    setnames(dt_matched, c("BETA", "SE", "ALLELE0", "ALLELE1", "A1FREQ"),
             c(paste0("beta_", protein), paste0("se_", protein),
               paste0("a0_", protein), paste0("a1_", protein),
               paste0("freq_", protein)))
    
    dt_matched[, c("pos_hg19", paste0("beta_", protein), paste0("se_", protein),
                   paste0("a0_", protein), paste0("a1_", protein),
                   paste0("freq_", protein)), with = FALSE]
  }
  
  mvmr_pairs <- list(
    list(p1 = "FES",    p2 = "FURIN",  chr = 15),
    list(p1 = "CGREF1", p2 = "KHK",    chr = 2)
  )
  
  mvmr_results <- list()
  
  for (pair in mvmr_pairs) {
    p1 <- pair$p1; p2 <- pair$p2; chr <- pair$chr
    cat(sprintf("══ MVMR: %s + %s (chr%d) ══\n", p1, p2, chr))
    
    leads_p1 <- get_lead_snps(p1, cs_dir)
    leads_p2 <- get_lead_snps(p2, cs_dir)
    all_positions <- unique(c(leads_p1$pos_hg19, leads_p2$pos_hg19))
    
    eff_p1 <- lookup_pqtl_effects(p1, chr, all_positions)
    eff_p2 <- lookup_pqtl_effects(p2, chr, all_positions)
    if (is.null(eff_p1) || is.null(eff_p2)) { cat("  ❌ No data\n\n"); next }
    
    # PTSD效应
    merged_p1 <- as.data.table(coloc_full[[p1]]$merged_common)
    merged_p2 <- as.data.table(coloc_full[[p2]]$merged_common)
    ptsd_eff <- rbind(
      merged_p1[pos_hg19 %in% all_positions, .(pos_hg19, beta_ptsd = beta_gwas,
                                                se_ptsd = se_gwas, a1_ptsd = a1_gwas,
                                                a2_ptsd = a2_gwas)],
      merged_p2[pos_hg19 %in% all_positions, .(pos_hg19, beta_ptsd = beta_gwas,
                                                se_ptsd = se_gwas, a1_ptsd = a1_gwas,
                                                a2_ptsd = a2_gwas)]
    )[!duplicated(pos_hg19)]
    
    dat <- merge(merge(eff_p1, eff_p2, by = "pos_hg19"), ptsd_eff, by = "pos_hg19")
    cat(sprintf("  %d SNPs\n", nrow(dat)))
    if (nrow(dat) < 4) { cat("  ⚠ Insufficient\n\n"); next }
    
    # Harmonization: 翻转p2和PTSD使其与p1方向一致
    a1_p1 <- paste0("a1_", p1); a0_p1 <- paste0("a0_", p1)
    a1_p2 <- paste0("a1_", p2); a0_p2 <- paste0("a0_", p2)
    bp2    <- paste0("beta_", p2)
    
    flip_p2 <- dat[[a1_p1]] == dat[[a0_p2]] & dat[[a0_p1]] == dat[[a1_p2]]
    if (any(flip_p2)) dat[flip_p2, (bp2) := -get(bp2)]
    
    flip_ptsd <- dat[[a1_p1]] == dat$a2_ptsd & dat[[a0_p1]] == dat$a1_ptsd
    if (any(flip_ptsd)) dat[flip_ptsd, beta_ptsd := -beta_ptsd]
    
    # MVMR
    bp1 <- paste0("beta_", p1); sp1 <- paste0("se_", p1); sp2 <- paste0("se_", p2)
    betaX <- as.matrix(dat[, .(get(bp1), get(bp2))])
    seX   <- as.matrix(dat[, .(get(sp1), get(sp2))])
    colnames(betaX) <- colnames(seX) <- c(p1, p2)
    
    mvmr_input <- format_mvmr(BXGs = betaX, BYG = dat$beta_ptsd,
                               seBXGs = seX, seBYG = dat$se_ptsd)
    mvmr_res <- ivw_mvmr(mvmr_input)
    cond_f   <- strength_mvmr(mvmr_input, gencov = 0)
    het_test <- tryCatch(pleiotropy_mvmr(mvmr_input, gencov = 0), error = function(e) NULL)
    
    cat(sprintf("  Conditional F: %s=%.1f, %s=%.1f\n",
                p1, cond_f$exposure[1], p2, cond_f$exposure[2]))
    
    for (j in 1:2) {
      prot_name <- c(p1, p2)[j]
      sig <- ifelse(mvmr_res[j, 4] < 0.05, "✅ 独立显著", "— 不独立显著")
      cat(sprintf("  %s: beta=%.4f, P=%s → %s\n", prot_name, mvmr_res[j, 1],
                  formatC(mvmr_res[j, 4], format = "e", digits = 2), sig))
      
      mvmr_results[[paste(p1, p2, prot_name, sep = "_")]] <- tibble(
        pair = paste0(p1, " + ", p2), protein = prot_name,
        n_iv = nrow(dat), mvmr_beta = mvmr_res[j, 1], mvmr_se = mvmr_res[j, 2],
        mvmr_pval = mvmr_res[j, 4], mvmr_OR = exp(mvmr_res[j, 1]),
        cond_F = cond_f$exposure[j],
        Q_pval = ifelse(!is.null(het_test), het_test$Qpval, NA),
        status = sig)
    }
    cat("\n")
  }
  
  mvmr_summary <- bind_rows(mvmr_results)
  write_csv(mvmr_summary, file.path(result_dir, "mvmr_results.csv"))
  print(mvmr_summary)
  
  invisible(mvmr_summary)
}


###############################################################################
# Part E: SMR + HEIDI
###############################################################################

run_smr_heidi <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part E: SMR + HEIDI Test                             ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  coloc_full <- readRDS(coloc_path)
  results <- list()
  
  for (protein in candidates) {
    cat(sprintf("══ %s ══\n", protein))
    merged <- as.data.table(coloc_full[[protein]]$merged_common)
    
    # SMR: top pQTL
    merged[, z_pqtl := abs(beta_pqtl / se_pqtl)]
    top_idx <- which.max(merged$z_pqtl)
    top <- merged[top_idx]
    
    z_gwas <- top$beta_gwas / top$se_gwas
    z_pqtl <- top$beta_pqtl / top$se_pqtl
    chi2_smr <- z_gwas^2 * z_pqtl^2 / (z_gwas^2 + z_pqtl^2)
    p_smr <- pchisq(chi2_smr, df = 1, lower.tail = FALSE)
    b_smr <- top$beta_gwas / top$beta_pqtl
    se_smr <- sqrt(top$se_gwas^2 / top$beta_pqtl^2 +
                     top$beta_gwas^2 * top$se_pqtl^2 / top$beta_pqtl^4)
    
    # HEIDI: 区域内其他SNPs
    other <- merged[-top_idx][z_pqtl > 3]
    if (nrow(other) > 20) other <- other[order(-z_pqtl)][1:20]
    
    if (nrow(other) >= 3) {
      other[, d_i := beta_gwas / beta_pqtl]
      other[, d_diff := d_i - b_smr]
      other[, var_d_diff := se_gwas^2 / beta_pqtl^2 +
              top$se_gwas^2 / top$beta_pqtl^2]
      other[, heidi_chi2_i := d_diff^2 / var_d_diff]
      
      heidi_chi2 <- sum(other$heidi_chi2_i)
      heidi_p <- pchisq(heidi_chi2, df = nrow(other), lower.tail = FALSE)
      heidi_note <- ifelse(heidi_p > 0.05, "✅ shared causal variant",
                           "⚠ 可能linkage")
    } else {
      heidi_chi2 <- NA; heidi_p <- NA
      heidi_note <- "SNPs不足"
    }
    
    cat(sprintf("  SMR: P=%s | HEIDI: %d SNPs, P=%s → %s\n\n",
                formatC(p_smr, format = "e", digits = 2),
                nrow(other),
                ifelse(is.na(heidi_p), "NA", formatC(heidi_p, format = "e", digits = 2)),
                heidi_note))
    
    results[[protein]] <- tibble(
      protein = protein, top_snp = top$rsid, n_snps_region = nrow(merged),
      smr_beta = b_smr, smr_se = se_smr, smr_pval = p_smr, smr_OR = exp(b_smr),
      heidi_n_snps = nrow(other), heidi_pval = heidi_p, heidi_note = heidi_note,
      overall = case_when(
        p_smr < 0.05 & !is.na(heidi_p) & heidi_p > 0.05 ~ "✅ 强证据",
        p_smr < 0.05 & !is.na(heidi_p) & heidi_p <= 0.05 ~ "⚠ HEIDI拒绝",
        p_smr < 0.05 & is.na(heidi_p) ~ "SMR显著, HEIDI无法检验",
        TRUE ~ "— SMR不显著"))
  }
  
  smr_summary <- bind_rows(results)
  write_csv(smr_summary, file.path(result_dir, "smr_heidi_results.csv"))
  
  smr_summary %>%
    select(protein, smr_pval, smr_OR, heidi_pval, overall) %>%
    mutate(smr_OR = round(smr_OR, 3),
           smr_pval = formatC(smr_pval, format = "e", digits = 2),
           heidi_pval = ifelse(is.na(heidi_pval), "NA",
                               formatC(heidi_pval, format = "e", digits = 2))) %>%
    as.data.frame() %>% print(right = FALSE)
  
  invisible(smr_summary)
}


###############################################################################
# Part F: Coloc Prior Sensitivity
###############################################################################

run_coloc_sensitivity <- function() {
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  Part F: Coloc Prior Sensitivity Analysis             ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  
  library(coloc)
  coloc_full <- readRDS(coloc_path)
  
  p12_values <- c(1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4)
  p1 <- 1e-4; p2 <- 1e-4
  
  all_results <- list()
  
  for (protein in candidates) {
    cat(sprintf("══ %s ══\n", protein))
    merged <- as.data.table(coloc_full[[protein]]$merged_common)
    
    d1 <- list(beta = merged$beta_pqtl, varbeta = merged$se_pqtl^2,
               type = "quant", N = n_pqtl, snp = merged$rsid,
               MAF = pmin(merged$freq_pqtl, 1 - merged$freq_pqtl))
    d2 <- list(beta = merged$beta_gwas, varbeta = merged$se_gwas^2,
               type = "cc", N = n_ptsd, s = 137136 / n_ptsd,
               snp = merged$rsid,
               MAF = pmin(merged$freq_gwas, 1 - merged$freq_gwas))
    
    pph4_values <- sapply(p12_values, function(p12_val) {
      res <- tryCatch(coloc.abf(dataset1 = d1, dataset2 = d2,
                                 p1 = p1, p2 = p2, p12 = p12_val),
                      error = function(e) NULL)
      if (!is.null(res)) res$summary["PP.H4.abf"] else NA
    })
    
    for (i in seq_along(p12_values)) {
      label <- ifelse(p12_values[i] == 1e-5, " ← default", "")
      cat(sprintf("  p12=%s → PPH4=%.4f%s\n",
                  formatC(p12_values[i], format = "e", digits = 0),
                  pph4_values[i], label))
    }
    
    min_pph4 <- min(pph4_values, na.rm = TRUE)
    robust <- ifelse(min_pph4 > 0.8, "✅ 稳健",
                     ifelse(min_pph4 > 0.5, "— 中等", "⚠ prior敏感"))
    cat(sprintf("  → min=%.4f, %s\n\n", min_pph4, robust))
    
    result_row <- tibble(protein = protein, n_snps = nrow(merged))
    for (i in seq_along(p12_values)) {
      col_name <- sprintf("PPH4_p12_%s", gsub("[+]", "", formatC(p12_values[i], format = "e", digits = 0)))
      result_row[[col_name]] <- pph4_values[i]
    }
    result_row$min_PPH4 <- min_pph4
    result_row$robust <- robust
    all_results[[protein]] <- result_row
  }
  
  sensitivity <- bind_rows(all_results)
  write_csv(sensitivity, file.path(result_dir, "coloc_sensitivity_results.csv"))
  
  sensitivity %>%
    select(protein, n_snps, min_PPH4, robust) %>%
    mutate(min_PPH4 = round(min_PPH4, 4)) %>%
    as.data.frame() %>% print(right = FALSE)
  
  invisible(sensitivity)
}


###############################################################################
# 主函数: 按顺序运行所有分析
###############################################################################

run_all_sensitivity <- function() {
  cat("╔═══════════════════════════════════════════════════════╗\n")
  cat("║  PTSD Proteome MR — Sensitivity Analyses              ║\n")
  cat("║  10 FDR candidates × 6 validation methods             ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n")
  
  t0 <- Sys.time()
  
  steiger_res <- run_steiger()
  presso_res  <- run_mr_presso()
  reverse_res <- run_reverse_mr()
  mvmr_res    <- run_mvmr()
  smr_res     <- run_smr_heidi()
  coloc_res   <- run_coloc_sensitivity()
  
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  
  cat("\n╔═══════════════════════════════════════════════════════╗\n")
  cat("║  All done                                              ║\n")
  cat("╚═══════════════════════════════════════════════════════╝\n\n")
  cat(sprintf("Total time: %.1f minutes\n\n", as.numeric(elapsed)))
  cat("Output files:\n")
  cat("  steiger_directionality_results.csv\n")
  cat("  mr_presso_results.csv\n")
  cat("  reverse_mr_results.csv\n")
  cat("  mvmr_results.csv\n")
  cat("  smr_heidi_results.csv\n")
  cat("  coloc_sensitivity_results.csv\n")
  cat(sprintf("\nAll outputs saved to: %s\n", result_dir))
  
  invisible(list(steiger = steiger_res, presso = presso_res,
                 reverse = reverse_res, mvmr = mvmr_res,
                 smr = smr_res, coloc_sens = coloc_res))
}

# ── 运行 ──
# 可以整体运行:
#   run_all_sensitivity()
# 也可以单独运行某个模块:
#   run_steiger()
#   run_mr_presso()
#   run_reverse_mr()
#   run_mvmr()
#   run_smr_heidi()
#   run_coloc_sensitivity()
