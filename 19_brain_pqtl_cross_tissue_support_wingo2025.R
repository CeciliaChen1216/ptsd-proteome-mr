###############################################################################
# 19_brain_pqtl_cross_tissue_support_wingo2025.R
#
# cross-tissue support: scan the supplementary tables of Wingo et al.
# Nature Genetics 2025 (multi-ancestry brain pQTL atlas, N = 1,362 dlPFC) for
# the ten PTSD-prioritized candidate proteins, and align the brain-derived
# evidence with this study's plasma cis-pQTL MR results.
#
# This script identifies which Wingo supplementary sheets contain each
# candidate; the precise quantitative estimates (β, SE, P, FDR, PIP) used in
# manuscript Table S20 are then curated from the relevant Wingo sheets and
# encoded in script 20_generate_supplementary_tables.R.
#
# Inputs (from 00_config.R):
#   wingo2025_dir       — local cache for Wingo 2025 supplementary xlsx files
#   candidates,         — Tier-1 protein symbols
#   candidate_uniprot,
#   result_dir
#
# Outputs:
#   results/wingo2025_brain_pqtl_support.csv
#   results/wingo2025_brain_pqtl_support.rds
#   results/wingo2025_sheet_overview.csv
#   results/wingo2025_candidate_details/<protein>_wingo2025_rows.csv
###############################################################################

# ─────────────────────────────────────────────────────────────────────────
# Step 0: Locate 00_config.R in the same directory as this script
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
  Sys.getenv("PTSD_CODE_DIR", unset = ""),
  if (!is.na(script_dir)) script_dir else "",
  getwd()
)
config_candidates <- unique(config_candidates[nzchar(config_candidates)])

config_path <- ""
for (d in config_candidates) {
  cand <- file.path(d, "00_config.R")
  if (file.exists(cand)) { config_path <- cand; break }
}
if (!nzchar(config_path)) {
  cat("Searched locations:\n")
  for (d in config_candidates) cat(sprintf("  - %s\n", d))
  stop("Cannot find 00_config.R. Please copy 00_config_template.R to 00_config.R in the same directory as this script.")
}
cat(sprintf("✓ Loaded config: %s\n", config_path))
source(config_path)

needed_vars <- c("candidates", "candidate_uniprot", "wingo2025_dir", "result_dir")
missing_vars <- needed_vars[!sapply(needed_vars, exists)]
if (length(missing_vars) > 0) {
  stop(sprintf("Missing variables in 00_config.R: %s\nPlease use the current template.",
               paste(missing_vars, collapse = ", ")))
}

suppressMessages({
  needed <- c("openxlsx", "httr", "dplyr", "data.table")
  missing_pkg <- needed[!sapply(needed, requireNamespace, quietly = TRUE)]
  if (length(missing_pkg) > 0) {
    cat(sprintf("⚠ Please install required packages: install.packages(c(%s))\n",
                paste(sprintf('"%s"', missing_pkg), collapse = ", ")))
    stop("Missing required R packages.")
  }
  library(openxlsx)
  library(httr)
  library(dplyr)
  library(data.table)
})

cat("\n╔══════════════════════════════════════════════════════════════════╗\n")
cat("║   19. Wingo 2025 Brain pQTL — cross-tissue support Lookup ║\n")
cat("╚══════════════════════════════════════════════════════════════════╝\n\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 1: 准备目录
# ─────────────────────────────────────────────────────────────────────────
if (!dir.exists(wingo2025_dir)) {
  dir.create(wingo2025_dir, recursive = TRUE, showWarnings = FALSE)
}
cat(sprintf("Wingo 2025 data cache directory: %s\n\n", wingo2025_dir))

# ─────────────────────────────────────────────────────────────────────────
# Step 2: 下载所有 supplementary xlsx
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 2: Download Wingo 2025 supplementary tables ─────\n")

base_url <- paste0("https://static-content.springer.com/esm/",
                   "art%3A10.1038%2Fs41588-025-02291-2/MediaObjects/")

download_supp <- function(moesm_num, ext = "xlsx") {
  fname <- sprintf("41588_2025_2291_MOESM%d_ESM.%s", moesm_num, ext)
  local_path <- file.path(wingo2025_dir, fname)
  if (file.exists(local_path) && file.info(local_path)$size > 1000) {
    return(list(path = local_path, status = "cached"))
  }
  url <- paste0(base_url, fname)
  res <- tryCatch(
    GET(url, write_disk(local_path, overwrite = TRUE),
        timeout(120), user_agent("Mozilla/5.0 R-script")),
    error = function(e) NULL
  )
  if (is.null(res) || status_code(res) != 200) {
    if (file.exists(local_path)) unlink(local_path)
    return(list(path = NA_character_, status = "failed"))
  }
  if (file.info(local_path)$size < 1000) {
    unlink(local_path)
    return(list(path = NA_character_, status = "too_small"))
  }
  list(path = local_path, status = "downloaded")
}

downloaded <- list()
for (i in 3:30) {
  res <- download_supp(i, "xlsx")
  if (!is.na(res$path)) {
    downloaded[[as.character(i)]] <- res$path
    cat(sprintf("  ✓ MOESM%d: %s (%s)\n", i,
                basename(res$path), res$status))
  }
}

if (length(downloaded) == 0) {
  cat("\n❌ All automatic downloads failed. Please download manually and re-run:\n")
  cat("   1) Open https://www.nature.com/articles/s41588-025-02291-2#Sec22\n")
  cat("   2) Download all .xlsx files from the Supplementary Information section\n")
  cat(sprintf("   3) Place them under: %s\n\n", wingo2025_dir))
  stop("Failed to download the Wingo 2025 supplementary tables.")
}
cat(sprintf("\nDownloaded %d xlsx files (Wingo 2025 combines all supp tables in one file)\n\n",
            length(downloaded)))

# ─────────────────────────────────────────────────────────────────────────
# Step 3: 扫描所有 xlsx, 找含候选蛋白且与 PTSD 相关的内容
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 3: Scan supplementary content ─────\n")

# 评分函数: 返回 numeric (避免类型问题)
score_sheet_relevance <- function(df) {
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) return(0)
  blob <- tolower(paste(c(
    colnames(df),
    unlist(lapply(df, function(x) head(as.character(x), 50)))
  ), collapse = " "))
  ptsd_kw <- as.numeric(any(grepl("ptsd", blob))) +
             as.numeric(any(grepl("post.?traumatic", blob))) +
             as.numeric(any(grepl("stress.?disorder", blob)))
  mr_kw <- as.numeric(any(grepl("\\bmr\\b|mendelian|causal.+pair|causal.+protein", blob))) +
           as.numeric(any(grepl("\\bbeta\\b|\\bp.?value\\b|\\bfdr\\b|\\bse\\b", blob)))
  ptsd_kw * 10 + mr_kw
}

candidate_genes <- candidates
candidate_uniprots_lower <- tolower(candidate_uniprot)

# 列出每个 xlsx 的所有 sheet 名 (帮助人工查找 PTSD MR 表)
cat("\nSheet listing per xlsx:\n")
for (xlsx_path in downloaded) {
  fname <- basename(xlsx_path)
  sheet_names <- tryCatch(getSheetNames(xlsx_path),
                          error = function(e) character(0))
  cat(sprintf("\n  %s — %d sheets:\n", fname, length(sheet_names)))
  for (i in seq_along(sheet_names)) {
    cat(sprintf("    %2d. %s\n", i, sheet_names[i]))
  }
}
cat("\n")

# 实际扫描每个 sheet
all_sheets <- list()
for (xlsx_path in downloaded) {
  fname <- basename(xlsx_path)
  sheet_names <- tryCatch(getSheetNames(xlsx_path),
                          error = function(e) character(0))
  for (sn in sheet_names) {
    df <- tryCatch(
      read.xlsx(xlsx_path, sheet = sn, colNames = TRUE,
                detectDates = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(df)) next
    rel_score <- score_sheet_relevance(df)
    blob <- tolower(paste(unlist(df), collapse = " "))
    gene_hits <- vapply(candidate_genes, function(g)
      grepl(sprintf("\\b%s\\b", tolower(g)), blob), logical(1))
    uniprot_hits <- vapply(candidate_uniprots_lower, function(u)
      grepl(sprintf("\\b%s\\b", u), blob), logical(1))
    any_hit <- any(gene_hits) || any(uniprot_hits)

    all_sheets[[length(all_sheets) + 1]] <- list(
      file = fname, sheet = sn,
      n_row = as.integer(nrow(df)),
      n_col = as.integer(ncol(df)),
      cols = colnames(df),
      rel_score = as.numeric(rel_score),
      gene_hits = candidate_genes[gene_hits],
      uniprot_hits = names(candidate_uniprots_lower)[uniprot_hits],
      any_hit = any_hit, df = df
    )
  }
}

cat(sprintf("Scanned %d sheets in total.\n\n", length(all_sheets)))

# 全部用 numeric(1) / character(1) / logical(1) 严格匹配
sheet_overview <- data.frame(
  file      = vapply(all_sheets, function(s) as.character(s$file),    character(1)),
  sheet     = vapply(all_sheets, function(s) as.character(s$sheet),   character(1)),
  n_row     = vapply(all_sheets, function(s) as.numeric(s$n_row),     numeric(1)),
  n_col     = vapply(all_sheets, function(s) as.numeric(s$n_col),     numeric(1)),
  rel_score = vapply(all_sheets, function(s) as.numeric(s$rel_score), numeric(1)),
  any_hit   = vapply(all_sheets, function(s) as.logical(s$any_hit),   logical(1)),
  hit_genes = vapply(all_sheets, function(s)
                paste(unique(c(s$gene_hits, s$uniprot_hits)), collapse = "|"),
                character(1)),
  stringsAsFactors = FALSE
)
sheet_overview <- sheet_overview[order(-sheet_overview$rel_score,
                                        -as.integer(sheet_overview$any_hit)), ]

cat("Top sheets by relevance + candidate hits:\n")
print(head(sheet_overview, 20), row.names = FALSE)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 4: 提取 PTSD-related sheet 中的候选蛋白行
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 4: Extract candidate-protein records from Wingo 2025 ─────\n\n")

ptsd_all_sheets <- Filter(function(s) s$rel_score >= 10, all_sheets)
cat(sprintf("Number of PTSD-related sheets: %d\n\n", length(ptsd_all_sheets)))

extract_candidate_rows <- function(sheet_info, candidates_to_match) {
  df <- sheet_info$df
  row_blobs <- apply(df, 1, function(r)
    tolower(paste(as.character(r), collapse = " ")))
  out <- list()
  for (g in candidates_to_match) {
    hits <- which(grepl(sprintf("\\b%s\\b", tolower(g)), row_blobs))
    upr <- candidate_uniprot[g]
    if (!is.na(upr)) {
      hits_u <- which(grepl(sprintf("\\b%s\\b", tolower(upr)), row_blobs))
      hits <- unique(c(hits, hits_u))
    }
    if (length(hits) > 0) {
      sub <- df[hits, , drop = FALSE]
      sub$.candidate <- g
      sub$.source_file <- sheet_info$file
      sub$.source_sheet <- sheet_info$sheet
      out[[g]] <- sub
    }
  }
  out
}

candidate_rows <- list()
for (g in candidate_genes) candidate_rows[[g]] <- list()

for (s in ptsd_all_sheets) {
  rows <- extract_candidate_rows(s, candidate_genes)
  for (g in names(rows)) {
    candidate_rows[[g]][[length(candidate_rows[[g]]) + 1]] <- rows[[g]]
  }
}

candidate_summary <- data.frame(
  candidate = candidate_genes,
  uniprot = candidate_uniprot[candidate_genes],
  in_wingo2025_ptsd = FALSE,
  n_rows = 0L,
  source_files = "",
  raw_evidence = "",
  stringsAsFactors = FALSE
)

for (i in seq_along(candidate_genes)) {
  g <- candidate_genes[i]
  rows <- candidate_rows[[g]]
  if (length(rows) > 0) {
    n_total <- sum(vapply(rows, function(r) as.integer(nrow(r)), integer(1)))
    candidate_summary$in_wingo2025_ptsd[i] <- TRUE
    candidate_summary$n_rows[i] <- n_total
    files_used <- unique(unlist(lapply(rows, function(r) r$.source_file)))
    candidate_summary$source_files[i] <- paste(files_used, collapse = ";")
    first_row <- rows[[1]][1, ]
    keep <- which(!is.na(first_row) &
                  nchar(as.character(first_row)) > 0 &
                  !names(first_row) %in% c(".candidate", ".source_file",
                                           ".source_sheet"))
    keep <- keep[seq_len(min(8, length(keep)))]
    pairs <- paste(names(first_row)[keep],
                   as.character(first_row)[keep],
                   sep = "=", collapse = "; ")
    candidate_summary$raw_evidence[i] <- substr(pairs, 1, 250)
  }
}

cat("Cross-reference summary:\n")
print(candidate_summary[, c("candidate", "uniprot", "in_wingo2025_ptsd",
                            "n_rows", "source_files")],
      row.names = FALSE)
cat("\n")

# ─────────────────────────────────────────────────────────────────────────
# Step 5: 跟本研究外周 MR 结果对齐方向
# ─────────────────────────────────────────────────────────────────────────
cat("───── Step 5: Align with this study's plasma MR results ─────\n")

peripheral_mr_path <- file.path(result_dir, "mr_all_outcomes_extended.rds")
if (file.exists(peripheral_mr_path)) {
  per_mr <- tryCatch(readRDS(peripheral_mr_path), error = function(e) NULL)
} else {
  per_mr <- NULL
}

if (!is.null(per_mr)) {
  per_df <- if (is.data.frame(per_mr)) as.data.frame(per_mr) else NULL
  if (!is.null(per_df)) {
    outcome_col <- grep("outcome|trait|disease",
                        colnames(per_df), value = TRUE, ignore.case = TRUE)[1]
    protein_col <- grep("^protein$|^exposure$|^gene$",
                        colnames(per_df), value = TRUE, ignore.case = TRUE)[1]
    beta_col    <- grep("^beta$|^b$|^effect$|estimate",
                        colnames(per_df), value = TRUE, ignore.case = TRUE)[1]
    p_col       <- grep("^p$|^pval$|^p_value$|^p\\.value$",
                        colnames(per_df), value = TRUE, ignore.case = TRUE)[1]

    cat(sprintf("Plasma MR columns identified: outcome=%s, protein=%s, beta=%s, p=%s\n",
                outcome_col, protein_col, beta_col, p_col))

    if (!is.na(outcome_col) && !is.na(protein_col) && !is.na(beta_col)) {
      ptsd_mask <- grepl("ptsd|stress", per_df[[outcome_col]],
                         ignore.case = TRUE)
      per_ptsd <- per_df[ptsd_mask & per_df[[protein_col]] %in% candidate_genes, ]
      cat(sprintf("Plasma PTSD MR candidate-protein coverage: %d/%d\n",
                  length(unique(per_ptsd[[protein_col]])),
                  length(candidate_genes)))

      candidate_summary$peripheral_beta <- NA_real_
      candidate_summary$peripheral_p    <- NA_real_
      for (i in seq_len(nrow(candidate_summary))) {
        g <- candidate_summary$candidate[i]
        match_row <- per_ptsd[per_ptsd[[protein_col]] == g, , drop = FALSE]
        if (nrow(match_row) > 0) {
          candidate_summary$peripheral_beta[i] <- as.numeric(match_row[[beta_col]][1])
          if (!is.na(p_col))
            candidate_summary$peripheral_p[i] <- as.numeric(match_row[[p_col]][1])
        }
      }
    }
  }
} else {
  cat("⚠ mr_all_outcomes_extended.rds not found; skipping plasma MR alignment\n")
}

# ─────────────────────────────────────────────────────────────────────────
# Step 6: 输出
# ─────────────────────────────────────────────────────────────────────────
cat("\n───── Step 6: Write outputs ─────\n")

candidate_summary$cross_tissue_support <- vapply(seq_len(nrow(candidate_summary)),
  function(i) {
    if (!candidate_summary$in_wingo2025_ptsd[i]) "Not_in_Wingo2025_PTSD"
    else "Found_in_Wingo2025_PTSD_supp"
  }, character(1))

out_csv <- file.path(result_dir, "wingo2025_brain_pqtl_support.csv")
write.csv(candidate_summary, out_csv, row.names = FALSE)
cat(sprintf("✓ Written: %s\n", out_csv))

overview_csv <- file.path(result_dir, "wingo2025_sheet_overview.csv")
write.csv(sheet_overview, overview_csv, row.names = FALSE)
cat(sprintf("✓ Written: %s (relevance ranking across 29 sheets)\n", overview_csv))

out_rds <- file.path(result_dir, "wingo2025_brain_pqtl_support.rds")
saveRDS(list(summary = candidate_summary,
             raw_rows = candidate_rows,
             ptsd_sheets_overview = sheet_overview,
             downloaded_files = downloaded),
        out_rds)
cat(sprintf("✓ Written: %s (raw rows + sheet index)\n", out_rds))

detail_dir <- file.path(result_dir, "wingo2025_candidate_details")
dir.create(detail_dir, showWarnings = FALSE, recursive = TRUE)
for (g in candidate_genes) {
  rows <- candidate_rows[[g]]
  if (length(rows) > 0) {
    combined <- tryCatch(
      rbindlist(lapply(rows, as.data.table), fill = TRUE),
      error = function(e) NULL
    )
    if (!is.null(combined)) {
      out_path <- file.path(detail_dir, sprintf("%s_wingo2025_rows.csv", g))
      fwrite(combined, out_path)
    }
  }
}
cat(sprintf("✓ Written: %s/<protein>_wingo2025_rows.csv (per-protein details)\n",
            detail_dir))

cat("\n══════════════════════════════════════════════════════════════════\n")
cat("Summary — cross-tissue support via Wingo 2025 brain pQTL MR\n")
cat("══════════════════════════════════════════════════════════════════\n")
n_in_wingo <- sum(candidate_summary$in_wingo2025_ptsd)
cat(sprintf("\nOf the 10 candidate proteins, %d appear in Wingo 2025 PTSD-related supp tables:\n",
            n_in_wingo))
for (i in seq_len(nrow(candidate_summary))) {
  status_emoji <- if (candidate_summary$in_wingo2025_ptsd[i]) "✓" else "—"
  cat(sprintf("  %s %-8s (%s)\n",
              status_emoji,
              candidate_summary$candidate[i],
              candidate_summary$uniprot[i]))
}

cat("\nNext steps:\n")
cat("  1) Inspect results/wingo2025_candidate_details/<protein>_wingo2025_rows.csv\n")
cat("     Verify each row is a PTSD MR result (vs PWAS / coloc / SMR or other analyses)\n")
cat("  2) If a protein is significant in Wingo 2025 PTSD MR (FDR<0.05) with consistent direction:\n")
cat("     -> add a cross-tissue support paragraph in the manuscript discussion\n")
cat("  3) If non-significant but directionally consistent: interpret as power-limited (Wingo N=1362 vs UKB-PPP N=54k)\n")
cat("  4) If directions are opposite: discuss in Limitations as possible tissue-specific divergence\n\n")
