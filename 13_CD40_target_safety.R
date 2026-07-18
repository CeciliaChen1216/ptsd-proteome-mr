###############################################################################
# 13_CD40_target_safety.R — CD40 on-target safety assessment (literature)
#
# Source: Zhao et al. 2023 Nat Immunol (SCALLOP consortium, rs1883832)
# Output: cd40_target_safety_literature.csv
###############################################################################
source("00_config.R")

library(dplyr)
library(readr)


cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 13: CD40 Target Safety (Literature)            ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

safety <- tribble(
  ~disease,                    ~domain,       ~direction_cd40_high, ~significance, ~source,
  "PTSD",                     "Psychiatric", "Protective (OR=0.953)","FDR<0.05",   "Present study",
  "MDD",                      "Psychiatric", "Protective (nominal)", "Nominal",    "Present study",
  "Anxiety disorders",        "Psychiatric", "Protective (nominal)", "Nominal",    "Present study",
  "Rheumatoid arthritis",     "Autoimmune",  "Protective",          "Significant","Zhao 2023 Nat Immunol",
  "Multiple sclerosis",       "Autoimmune",  "Risk-increasing",     "Significant","Zhao 2023 Nat Immunol",
  "Inflammatory bowel disease","Autoimmune", "Risk-increasing",     "Significant","Zhao 2023 Nat Immunol",
  "Ulcerative colitis",       "Autoimmune",  "Risk-increasing",     "Significant","Zhao 2023 Nat Immunol"
) %>% mutate(agonism_implication = case_when(
  direction_cd40_high=="Protective" ~ "Beneficial", TRUE ~ "Safety concern"))

print(safety[,c("disease","domain","direction_cd40_high","agonism_implication")])
write_csv(safety, file.path(result_dir,"cd40_target_safety_literature.csv"))
cat("\n✓ cd40_target_safety_literature.csv\n")
