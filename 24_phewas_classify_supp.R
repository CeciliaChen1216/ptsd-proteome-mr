###############################################################################
# 24_phewas_classify_supp.R
# Classify pheWAS associations and build Supplementary Table S23.
# Run AFTER 23_phewas_pleiotropy.R (reads phewas_full_associations.csv).
#
# Mechanism-aware classification:
#   cis_self     - eQTL/pQTL associations or "...levels" = the cis instrument
#                  acting on its own transcript/protein (NOT pleiotropy)
#   psychiatric  - PTSD / depression / anxiety / related internalizing &
#                  psychiatric traits = the substantive comparators
#   immune       - immune / inflammatory / blood-cell traits (mechanistically
#                  relevant to the neuroimmune hypothesis)
#   metabolic    - fructose / glucose / lipid / urate / adiposity (relevant to
#                  KHK fructose-metabolism biology)
#   other        - traits with no mechanistic link to the exposure = the
#                  category used to gauge potential HORIZONTAL pleiotropy
###############################################################################

suppressPackageStartupMessages({ library(dplyr); library(readr); library(tidyr) })

ROOT     <- "D:/PTSD"
IN_FULL  <- file.path(ROOT, "results", "phewas_full_associations.csv")
OUT_S23  <- file.path(ROOT, "results", "phewas_TableS23_classified.csv")
OUT_WIDE <- file.path(ROOT, "results", "phewas_TableS23_summary_wide.csv")

df <- readr::read_csv(IN_FULL, show_col_types = FALSE)

re_psych  <- "ptsd|post-trauma|stress disorder|depress|mdd|mood|anxiet|neurotic|internalizing|bipolar|schizophren|wellbeing|well-being|worry|feeling"
re_immune <- paste0("immun|inflamm|rheumat|arthrit|lupus|crohn|colitis|bowel|celiac|coeliac|psoriasis|thyro|graves|sclerosis|allerg|asthma|eczema|",
                    "diabetes type 1|t1d|autoimmun|leukocyte|lymphocyte|monocyte|neutrophil|eosinophil|basophil|white blood|c-reactive|cytokine|",
                    "interleukin|ig[gam]\\b|b cell|t cell|complement|blood cell|platelet|reticulocyte|corpuscular")
re_metab  <- "fructose|glucose|metaboli|lipid|cholesterol|triglycerid|hdl|ldl|urate|uric|bmi|body mass|adipos|weight|fat\\b|glycemic|hba1c|igf|diabetes type 2|t2d"

classify <- function(id, trait) {
  t <- tolower(trait)
  if (grepl("^eqtl-a|^prot-a|^prot-b", id) || grepl("levels$| levels", t)) return("cis_self")
  if (grepl(re_psych,  t)) return("psychiatric")
  if (grepl(re_immune, t)) return("immune")
  if (grepl(re_metab,  t)) return("metabolic")
  "other"
}

df <- df %>% mutate(category = mapply(classify, id, trait))

# long table (Supplementary S23) - sorted, with p formatted
s23 <- df %>%
  arrange(protein, category, p) %>%
  transmute(protein, rsid, category, trait,
            opengwas_id = id, beta, se, p, n)
readr::write_csv(s23, OUT_S23)

# wide per-candidate summary
wide <- df %>%
  count(protein, category) %>%
  pivot_wider(names_from = category, values_from = n, values_fill = 0)
for (cc in c("cis_self","psychiatric","immune","metabolic","other"))
  if (!cc %in% names(wide)) wide[[cc]] <- 0L
wide <- wide %>%
  transmute(protein,
            cis_self, psychiatric, immune, metabolic,
            other_unrelated = other,
            total = cis_self + psychiatric + immune + metabolic + other) %>%
  arrange(desc(other_unrelated))
readr::write_csv(wide, OUT_WIDE)

cat("=== Mechanism-aware pheWAS classification (per candidate) ===\n")
print(as.data.frame(wide), row.names = FALSE)
cat("\n'other_unrelated' = associations with NO mechanistic link to the exposure;",
    "\nthis is the column used to gauge potential horizontal pleiotropy.\n")
cat("\nWritten:\n  ", OUT_S23, "\n  ", OUT_WIDE, "\n")
###############################################################################
