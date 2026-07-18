###############################################################################
# 26_genomicsem_stage0b_modelfix.R
#
# Genomic SEM — STAGE 0b: fix and stress-test the common-factor model BEFORE
# any commonfactorGWAS. Addresses the Heywood case (MDD residual < 0) and
# verifies the factor structure is stable under a non-negativity constraint.
#
# Decision rule to proceed to Stage 1:
#   - constrained model converges with no negative residual variances;
#   - PTSD loading remains ~0.9 and PTSD residual variance remains significant;
#   - no pathological standard errors;
#   - CFI/SRMR reasonable where estimable.
#
# IMPORTANT interpretive notes (kept conservative on purpose):
#   - The PTSD residual is "genetic variance NOT captured by the common factor",
#     NOT "PTSD-specific biology". It may reflect specific effects, trauma-exposure
#     or case-definition differences, sampling/diagnostic heterogeneity, or
#     measurement error. Do not label it "PTSD-specific" without locus-level support.
#   - ANX residual (0.20) is numerically larger than PTSD (0.15); do not claim PTSD
#     is "most specific" and do not dismiss ANX residual as noise without evidence.
###############################################################################

suppressMessages({ library(GenomicSEM); library(data.table) })

out_dir <- "D:/PTSD/results/genomicsem"
LDSCoutput <- readRDS(file.path(out_dir, "LDSCoutput_stage0.rds"))

trait_names <- c("PTSD","MDD","ANX")

## ---------- helper: print standardized solution cleanly --------------------
show_fit <- function(fit, tag){
  cat("\n========== MODEL:", tag, "==========\n")
  res <- fit$results
  print(res[, intersect(c("lhs","op","rhs",
                          "Unstandardized_Estimate","Unstandardized_SE",
                          "Standardized_Est","Standardized_SE","p_value"),
                        names(res))])
  cfi <- tryCatch(fit$CFI, error=function(e) NA)
  srmr<- tryCatch(fit$SRMR, error=function(e) NA)
  chisq<-tryCatch(fit$modelfit$chisq, error=function(e) NA)
  df  <- tryCatch(fit$modelfit$df, error=function(e) NA)
  cat(sprintf("\nFit: CFI=%s  SRMR=%s  chisq=%s  df=%s\n",
              format(cfi), format(srmr), format(chisq), format(df)))
}

## ---------- Model 1: unconstrained common factor (baseline) ----------------
m1 <- '
F1 =~ NA*PTSD + MDD + ANX
F1 ~~ 1*F1
'
fit1 <- usermodel(LDSCoutput, estimation="DWLS", model=m1,
                  CFIcalc=TRUE, std.lv=TRUE, imp_cov=FALSE)
show_fit(fit1, "1: unconstrained (baseline, expect MDD Heywood)")

## ---------- Model 2: MDD residual constrained to be >= 0 --------------------
# Fix MDD residual to 0 (boundary) — standard handling of the Heywood case.
m2 <- '
F1 =~ NA*PTSD + MDD + ANX
F1 ~~ 1*F1
MDD ~~ 0*MDD
'
fit2 <- usermodel(LDSCoutput, estimation="DWLS", model=m2,
                  CFIcalc=TRUE, std.lv=TRUE, imp_cov=FALSE)
show_fit(fit2, "2: MDD residual fixed to 0 (Heywood fix)")

## ---------- Model 3: all residuals constrained >= 0 (lower bound) -----------
# Soft non-negativity via lower bounds on residual variances.
m3 <- '
F1 =~ NA*PTSD + MDD + ANX
F1 ~~ 1*F1
PTSD ~~ a*PTSD
MDD  ~~ b*MDD
ANX  ~~ c*ANX
a > 0.001
b > 0.001
c > 0.001
'
fit3 <- tryCatch(
  usermodel(LDSCoutput, estimation="DWLS", model=m3,
            CFIcalc=TRUE, std.lv=TRUE, imp_cov=FALSE),
  error=function(e){ cat("\nModel 3 (bounded) error:", conditionMessage(e), "\n"); NULL })
if(!is.null(fit3)) show_fit(fit3, "3: residuals lower-bounded > 0.001")

## ---------- decision readout -----------------------------------------------
cat("\n========== STAGE 0b DECISION READOUT ==========\n")
get_load <- function(fit, trait){
  r<-fit$results; r<-r[r$op=="=~" & r$rhs==trait,]
  if(nrow(r)) r$Standardized_Est[1] else NA
}
get_resid <- function(fit, trait){
  r<-fit$results; r<-r[r$op=="~~" & r$lhs==trait & r$rhs==trait,]
  if(nrow(r)) r$Standardized_Est[1] else NA
}
get_resid_p <- function(fit, trait){
  r<-fit$results; r<-r[r$op=="~~" & r$lhs==trait & r$rhs==trait,]
  if(nrow(r)) r$p_value[1] else NA
}
cat(sprintf("Constrained model (MDD resid=0):\n"))
cat(sprintf("  PTSD loading = %.3f  (target ~0.9)\n", get_load(fit2,"PTSD")))
cat(sprintf("  PTSD residual = %.3f  (p = %s; should stay significant)\n",
            get_resid(fit2,"PTSD"), format(get_resid_p(fit2,"PTSD"))))
cat(sprintf("  ANX  loading = %.3f  residual = %.3f\n",
            get_load(fit2,"ANX"), get_resid(fit2,"ANX")))
cat("\nProceed to Stage 1 (commonfactorGWAS) only if:\n")
cat("  - constrained model converged, no negative residuals, no crazy SEs;\n")
cat("  - PTSD loading ~0.9 and PTSD residual still significant;\n")
cat("  - the SAME constrained model spec will be used in commonfactorGWAS.\n")

saveRDS(list(fit1=fit1, fit2=fit2, fit3=fit3),
        file.path(out_dir, "stage0b_models.rds"))
cat("\nSaved stage0b_models.rds. Stage 0b complete.\n")
