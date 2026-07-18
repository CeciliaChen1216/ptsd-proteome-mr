###############################################################################
# 08_eQTL_SMR.R — Brain (BrainMeta v2) + Blood (Westra) eQTL-SMR via SMR v1.3.1
# 对 10 个 PTSD-prioritized 基因跑标准 SMR+HEIDI（smr.exe），每基因取 cis-eQTL
# 最强(p_eQTL 最小)的 probe，按 FINAL verdict 规则判定。
# 依赖：smr.exe, 1000G EUR bfile, PTSD .ma, BrainMeta(per-chr)/Westra besd
###############################################################################
source("00_config.R")
result_dir <- "D:/PTSD/results"
library(data.table)

V     <- "D:/PTSD/validation_data/external_data"
smr   <- file.path(V,"smr","smr-1.3.1-win-x86_64","smr-1.3.1-win.exe")
bfile <- "D:/reference/1000G_v3/EUR"
ma    <- file.path(result_dir,"ptsd_freeze3_smr.ma")
brain_dir <- file.path(V,"brain_eqtl","BrainMeta_cis_eqtl_summary")
blood_besd<- file.path(V,"blood_eqtl","westra_eqtl_hg19")
tmp <- file.path(result_dir,"smr_runs"); dir.create(tmp, showWarnings=FALSE)

## 基因 → 染色体（用于选 brain per-chr besd）
gene_chr <- c(AKT3=1, CD40=20, CGREF1=2, FES=15, FURIN=15,
              SIRPA=20, CD101=1, KHK=2, SNX18=5, UBE2L6=11)

run_smr <- function(gene, besd, tag){
  out <- file.path(tmp, paste0(gene,"_",tag))
  args <- c("--bfile",shQuote(bfile),"--gwas-summary",shQuote(ma),
            "--beqtl-summary",shQuote(besd),"--gene",gene,
            "--peqtl-smr","5e-8","--out",shQuote(out),"--thread-num","4")
  system2(smr, args, stdout=file.path(tmp,paste0(gene,"_",tag,".log")), stderr=FALSE)
  sf <- paste0(out,".smr")
  if (!file.exists(sf)) return(NULL)
  d <- tryCatch(read.delim(sf), error=function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d[which.min(d$p_eQTL), ]          # 每基因取 cis-eQTL 最强的 probe
}

verdict_fn <- function(p_smr, p_heidi, nsnp){
  if (is.na(p_smr) || p_smr>=0.05) return("Not significant")
  if (is.na(p_heidi) || nsnp<3)    return("SMR significant, HEIDI untestable (n<3)")
  if (p_heidi<0.05)                return("SMR significant, HEIDI rejected (possible linkage)")
  "Strong evidence (SMR + HEIDI)"
}

rows <- list()
for (g in names(gene_chr)){
  cat("\n==", g, "==\n")
  ## Brain
  bb <- file.path(brain_dir, paste0("BrainMeta_cis_eQTL_chr", gene_chr[g]))
  rb <- run_smr(g, bb, "brain")
  ## Blood
  rl <- run_smr(g, blood_besd, "blood")
  for (tis in c("brain","blood")){
    r <- if (tis=="brain") rb else rl
    if (is.null(r)) { cat("  ",tis,": no result\n"); next }
    v <- verdict_fn(r$p_SMR, r$p_HEIDI, r$nsnp_HEIDI)
    cat(sprintf("  %-5s probe=%s top=%s b=%.3f pSMR=%.2e HEIDI=%.3f n=%d -> %s\n",
                tis, r$probeID, r$topSNP, r$b_SMR, r$p_SMR, r$p_HEIDI, r$nsnp_HEIDI, v))
    rows[[paste(g,tis)]] <- data.frame(
      gene=g, tissue=ifelse(tis=="brain","Brain (BrainMeta)","Blood (Westra)"),
      probe=r$probeID, top_snp=r$topSNP, b_SMR=r$b_SMR, p_SMR=r$p_SMR,
      n_heidi=r$nsnp_HEIDI, p_HEIDI=r$p_HEIDI, verdict=v)
  }
}
es <- do.call(rbind, rows); rownames(es)<-NULL
## ---- multiple-testing correction: FDR WITHIN each tissue -------------------
# Reviewer point: 10 genes x 2 tissues ~ up to 20 SMR tests; nominal p_SMR<0.05
# is not the same evidence tier as an FDR-supported result. Apply BH-FDR
# separately within brain and within blood (over genes that returned an SMR
# result), then grade the verdict into FDR-supported / nominal / not-significant.
es$p_SMR <- suppressWarnings(as.numeric(es$p_SMR))
es$smr_fdr <- NA_real_
for (tis in unique(es$tissue)) {
  ix <- which(es$tissue == tis & is.finite(es$p_SMR))
  if (length(ix) > 0) es$smr_fdr[ix] <- p.adjust(es$p_SMR[ix], method = "BH")
}

grade_fn <- function(p_smr, fdr, p_heidi, nsnp) {
  if (is.na(p_smr) || p_smr >= 0.05)                         return("Not significant")
  heidi_ok <- !(is.na(p_heidi) || nsnp < 3) && p_heidi >= 0.05
  heidi_rej <- !is.na(p_heidi) && nsnp >= 3 && p_heidi < 0.05
  if (heidi_rej)                                             return("SMR significant, HEIDI rejected (possible linkage)")
  tier <- if (!is.na(fdr) && fdr < 0.05) "FDR-supported" else "nominal only"
  if (!heidi_ok)                                            return(paste0("SMR significant (", tier, "), HEIDI untestable (n<3)"))
  paste0("Strong evidence (SMR + HEIDI, ", tier, ")")
}
es$verdict_fdr <- mapply(grade_fn, es$p_SMR, es$smr_fdr, es$p_HEIDI, es$n_heidi)

cat("\n== SMR per-tissue FDR summary ==\n")
for (tis in unique(es$tissue)) {
  sub <- es[es$tissue == tis & is.finite(es$p_SMR), ]
  cat(sprintf("  %s: %d tested, nominal p<0.05 = %d, FDR<0.05 = %d\n",
              tis, nrow(sub), sum(sub$p_SMR<0.05, na.rm=TRUE),
              sum(sub$smr_fdr<0.05, na.rm=TRUE)))
}

write.csv(es, file.path(result_dir,"eqtl_smr_brain_blood.csv"), row.names=FALSE)
cat("\n[OK] eqtl_smr_brain_blood.csv\n"); print(es)

## 计数
cat("\nStrong evidence 计数：\n")
print(table(es$tissue[grepl("Strong", es$verdict)]))