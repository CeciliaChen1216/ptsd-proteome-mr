###############################################################################
# 12_power_aware_classification.R — Power-aware cross-disorder classification
#
# Three-tier criteria:
#   (1) directional concordance across disorders
#   (2) Cochran's Q heterogeneity
#   (3) ratio of effect sizes between PTSD and the comparator disorder
# Output: power_aware_classification.csv
###############################################################################
source("00_config.R")

library(data.table)
library(dplyr)
library(readr)

mr_rds_path <- file.path(result_dir, "mr_all_outcomes.rds")
candidates  <- c("AKT3","CD40","CGREF1","FES","FURIN",
                 "SIRPA","CD101","KHK","SNX18","UBE2L6")

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 12: Power-Aware Classification                ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

mr_all <- readRDS(mr_rds_path)
ptsd <- mr_all %>% filter(protein %in% candidates, grepl("PTSD|ptsd", outcome, ignore.case=TRUE))
mdd  <- mr_all %>% filter(protein %in% candidates, grepl("MDD|mdd|depress", outcome, ignore.case=TRUE))
anx  <- mr_all %>% filter(protein %in% candidates, grepl("ANX|anx|anxiety", outcome, ignore.case=TRUE))

wide <- ptsd %>% select(protein, beta_ptsd=mr_beta, se_ptsd=mr_se, p_ptsd=mr_pval) %>%
  left_join(mdd %>% select(protein, beta_mdd=mr_beta, se_mdd=mr_se, p_mdd=mr_pval), by="protein") %>%
  left_join(anx %>% select(protein, beta_anx=mr_beta, se_anx=mr_se, p_anx=mr_pval), by="protein")

# Direction
wide <- wide %>% mutate(
  concordant_mdd = (sign(beta_ptsd)==sign(beta_mdd)),
  concordant_anx = (sign(beta_ptsd)==sign(beta_anx)))

# Cochran's Q
cochran_q <- function(b1,s1,b2,s2) {
  w1<-1/s1^2; w2<-1/s2^2; bf<-(b1*w1+b2*w2)/(w1+w2)
  Q<-w1*(b1-bf)^2+w2*(b2-bf)^2; pchisq(Q,1,lower.tail=FALSE)
}
wide$pQ_ptsd_mdd <- NA; wide$pQ_ptsd_anx <- NA; wide$I2_global <- NA
for(i in 1:nrow(wide)){
  if(!is.na(wide$beta_mdd[i])) wide$pQ_ptsd_mdd[i] <- cochran_q(wide$beta_ptsd[i],wide$se_ptsd[i],wide$beta_mdd[i],wide$se_mdd[i])
  if(!is.na(wide$beta_anx[i])) wide$pQ_ptsd_anx[i] <- cochran_q(wide$beta_ptsd[i],wide$se_ptsd[i],wide$beta_anx[i],wide$se_anx[i])
  if(!is.na(wide$beta_mdd[i])&!is.na(wide$beta_anx[i])){
    b<-c(wide$beta_ptsd[i],wide$beta_mdd[i],wide$beta_anx[i]); s<-c(wide$se_ptsd[i],wide$se_mdd[i],wide$se_anx[i])
    w<-1/s^2; bf<-sum(b*w)/sum(w); Q<-sum(w*(b-bf)^2); wide$I2_global[i]<-max(0,(Q-2)/Q)*100
  }
}

# Ratios
wide$ratio_mdd <- wide$beta_mdd/wide$beta_ptsd
wide$ratio_anx <- wide$beta_anx/wide$beta_ptsd

# Classification
wide <- wide %>% mutate(
  shared_mdd = concordant_mdd & (pQ_ptsd_mdd>=0.05|is.na(pQ_ptsd_mdd)),
  shared_anx = concordant_anx & (pQ_ptsd_anx>=0.05|is.na(pQ_ptsd_anx)),
  # A protein needs at least one comparator (MDD or anxiety) with usable data to be
  # classified. Previously, a protein lacking an MDD estimate (beta_mdd = NA) fell
  # through to "PTSD-predominant" by default, which mislabels missing data as PTSD
  # specificity (e.g. SIRPA, which has no MDD instrument but IS shared with anxiety).
  # We now require MDD data for a positive "PTSD-predominant" call and mark proteins
  # with an untestable MDD comparison as Unclassified.
  has_mdd = !is.na(beta_mdd),
  has_anx = !is.na(beta_anx),
  class_power_aware = case_when(
    !has_mdd & !has_anx                 ~ "Unclassified (no comparator data)",
    !has_mdd                            ~ "Unclassified (MDD untestable)",
    shared_mdd &  shared_anx            ~ "Transdiagnostic",
    shared_mdd & !shared_anx            ~ "PTSD-MDD shared",
    !shared_mdd & shared_anx            ~ "PTSD-Anxiety shared",
    TRUE                                ~ "PTSD-predominant"))

orig <- read_csv(file.path(result_dir,"ptsd_candidates_classified.csv"))
wide <- wide %>% left_join(orig %>% select(protein, class_original=category), by="protein")

for(i in 1:nrow(wide)) cat(sprintf("  %-8s | %-25s -> %-20s\n", wide$protein[i],
  ifelse(is.na(wide$class_original[i]),"—",wide$class_original[i]), wide$class_power_aware[i]))

# Final manuscript classification: collapse to 2 categories
# "Transdiagnostic" + shared subtypes → Transdiagnostic
# "PTSD-predominant" → Relative PTSD predominance
wide$class_manuscript <- dplyr::case_when(
  grepl("^Unclassified", wide$class_power_aware) ~ "Unclassified (insufficient comparator data)",
  wide$class_power_aware == "PTSD-predominant"   ~ "Relative PTSD predominance",
  TRUE                                            ~ "Transdiagnostic")
cat("\nManuscript classification:\n")
print(table(wide$class_manuscript))

write_csv(wide %>% select(protein,beta_ptsd,se_ptsd,p_ptsd,beta_mdd,se_mdd,p_mdd,
  beta_anx,se_anx,p_anx,concordant_mdd,concordant_anx,pQ_ptsd_mdd,pQ_ptsd_anx,
  I2_global,ratio_mdd,ratio_anx,class_original,class_power_aware,class_manuscript),
  file.path(result_dir,"power_aware_classification.csv"))
cat("\n\u2713 power_aware_classification.csv\n")
