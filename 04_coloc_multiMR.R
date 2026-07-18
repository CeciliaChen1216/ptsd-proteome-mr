###############################################################################
# 04_coloc_multiMR.R — SuSiE fine-mapping + coloc.susie (matched 1000G EUR LD)
# 从协调后的区域级中间产物 region_merged/ 跑 coloc.susie。
# 上游 cis 区域提取依赖 UKB-PPP 受控数据 (Synapse syn51365303)，不随代码再分发。
###############################################################################
source("00_config.R")
result_dir <- "D:/PTSD/results"
region_dir <- file.path(result_dir, "region_merged")
suppressMessages({library(ieugwasr); library(susieR); library(coloc)
  library(data.table); library(dplyr)})
stopifnot(nchar(Sys.getenv("OPENGWAS_JWT")) > 20)

align_ld <- function(ld, mc){
  rn <- rownames(ld); base <- sub("_[ACGT]+_[ACGT]+$","",rn)
  a1 <- sub("^.*_([ACGT]+)_[ACGT]+$","\\1",rn)
  idx <- match(base, mc$rsid); keep <- !is.na(idx)
  ld <- ld[keep,keep,drop=FALSE]; a1 <- a1[keep]; idx <- idx[keep]
  mc2 <- mc[idx,,drop=FALSE]
  flip <- ifelse(a1==mc2$a1_pqtl, 1, ifelse(a1==mc2$a0_pqtl, -1, NA))
  ok <- !is.na(flip); ld <- ld[ok,ok,drop=FALSE]; mc2 <- mc2[ok,,drop=FALSE]; flip <- flip[ok]
  S <- diag(flip, length(flip)); Ra <- S %*% ld %*% S
  dimnames(Ra) <- list(mc2$rsid, mc2$rsid)
  list(R=Ra, mc=mc2)
}

cat("\n== Script 04 (revised): coloc.susie with matched LD ==\n")
full <- list(); rows <- list()
for (prot in candidates){
  cat("\n--", prot, "--\n")
  res <- tryCatch({
    mc <- as.data.frame(readRDS(file.path(region_dir, paste0(prot, "_merged.rds"))))
    ld <- ieugwasr::ld_matrix(mc$rsid, pop="EUR", with_alleles=TRUE)
    al <- align_ld(ld, mc); R <- al$R; m <- al$mc
    cat("  aligned SNP:", nrow(R), "/", nrow(mc), "\n")
    sp <- susie_rss(bhat=m$beta_pqtl, shat=m$se_pqtl, R=R, n=n_pqtl, L=10)
    sg <- susie_rss(bhat=m$beta_gwas, shat=m$se_gwas, R=R, n=n_ptsd, L=10)
    npq <- length(sp$sets$cs_index); ngw <- length(sg$sets$cs_index)
    cat("  CS  pQTL:", npq, " GWAS:", ngw, "\n")
    csm <- NULL; pph3 <- NA; pph4 <- NA
    minp <- min(m$pval_gwas, na.rm=TRUE)
    verdict <- paste0("Not assessable (no PTSD signal; min P=", signif(minp,2), ")")
    if (npq>=1 && ngw>=1){
      csm <- as.data.frame(coloc.susie(sp, sg)$summary)
      if (!is.null(csm) && nrow(csm)>0){
        sentinel <- m$rsid[which.max(abs(m$beta_pqtl/m$se_pqtl))]
        cand <- csm[csm$hit1==sentinel,,drop=FALSE]
        if (nrow(cand)==0) cand <- csm[which.max(csm$PP.H4.abf),,drop=FALSE]
        best <- cand[which.max(cand$PP.H4.abf),,drop=FALSE]
        pph3 <- best$PP.H3.abf; pph4 <- best$PP.H4.abf
        verdict <- if (pph4>=0.8) "Colocalized" else if (pph3>=0.8) "Distinct signals (no colocalization)" else "Inconclusive"
      }
    }
    full[[prot]] <- list(protein=prot, susie_pqtl=sp, susie_gwas=sg, R=R,
                         merged_common=m, coloc_summary=csm,
                         PP.H3=pph3, PP.H4=pph4, verdict=verdict)
    data.frame(protein=prot, n_snp=nrow(R), n_cs_pqtl=npq, n_cs_gwas=ngw,
               PP.H3=round(pph3,4), PP.H4=round(pph4,4), verdict=verdict)
  }, error=function(e){
    cat("  x", conditionMessage(e), "\n")
    data.frame(protein=prot, n_snp=NA, n_cs_pqtl=NA, n_cs_gwas=NA,
               PP.H3=NA, PP.H4=NA, verdict=paste("ERROR:", conditionMessage(e)))
  })
  rows[[prot]] <- res
}
saveRDS(full, file.path(result_dir, "susie_coloc_full.rds"))
res_tab <- do.call(rbind, rows); rownames(res_tab) <- NULL
write.csv(res_tab, file.path(result_dir, "coloc_susie_results.csv"), row.names=FALSE)
cat("\n[OK] susie_coloc_full.rds + coloc_susie_results.csv\n")
print(res_tab)