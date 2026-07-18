###############################################################################
# 02_drug_repurposing.R — Open Targets Drug Annotation
#
# Query Open Targets GraphQL API for:
#   - Tractability (small molecule, antibody, PROTAC, other)
#   - Known drugs (clinical phase, MOA)
#   - Directional concordance (MR direction vs drug MOA)
#
# Outputs:
#   drug_target_summary.csv  — tractability per protein
#   drug_full_list.csv       — all drug-target pairs
#   Table3_drug_repurposing.csv — Phase II+ drugs
###############################################################################

source("00_config.R")
library(httr)
library(jsonlite)
library(dplyr)
library(readr)

cat("\n╔═══════════════════════════════════════════════════════╗\n")
cat("║  Script 02: Drug Repurposing (Open Targets)            ║\n")
cat("╚═══════════════════════════════════════════════════════╝\n\n")

# ── 1. Read candidate-protein MR directions ──
mr_all <- readRDS(file.path(result_dir, "mr_all_outcomes.rds"))
ptsd_mr <- mr_all %>%
  filter(grepl("PTSD", outcome), protein %in% candidates) %>%
  select(protein, mr_beta, OR) %>%
  mutate(mr_direction = ifelse(mr_beta > 0, "Risk", "Protective"))

# ── 2. Open Targets GraphQL查询 ──
query_open_targets <- function(ensembl_id) {
  query <- sprintf('{
    target(ensemblId: "%s") {
      id
      approvedSymbol
      tractability {
        label
        modality
        value
      }
      knownDrugs {
        rows {
          drug { name mechanismOfAction }
          phase
          status
          disease { name }
        }
      }
    }
  }', ensembl_id)
  
  res <- POST("https://api.platform.opentargets.org/api/v4/graphql",
              body = list(query = query), encode = "json",
              content_type_json())
  fromJSON(content(res, "text", encoding = "UTF-8"))
}

# Gene → Ensembl ID mapping
gene_ensembl <- c(
  AKT3   = "ENSG00000117020", CD40   = "ENSG00000101017",
  CGREF1 = "ENSG00000138028", FES    = "ENSG00000182511",
  FURIN  = "ENSG00000140564", SIRPA  = "ENSG00000198053",
  CD101  = "ENSG00000134256", KHK    = "ENSG00000138030",
  SNX18  = "ENSG00000178209", UBE2L6 = "ENSG00000156587"
)

# ── 3. 批量查询 ──
cat("Querying Open Targets API...\n\n")
drug_results <- list()
tractability_results <- list()

for (prot in candidates) {
  cat("  ", prot, "...")
  eid <- gene_ensembl[prot]
  
  ot <- tryCatch(query_open_targets(eid), error = function(e) NULL)
  if (is.null(ot) || is.null(ot$data$target)) { cat(" [not found]\n"); next }
  
  target <- ot$data$target
  
  # Tractability
  if (!is.null(target$tractability)) {
    tract <- as.data.frame(target$tractability)
    tract$protein <- prot
    tractability_results[[prot]] <- tract
  }
  
  # Known drugs
  if (!is.null(target$knownDrugs$rows) && length(target$knownDrugs$rows) > 0) {
    drugs <- target$knownDrugs$rows
    drugs_df <- data.frame(
      protein = prot,
      drug_name = sapply(drugs$drug, function(x) x$name),
      moa = sapply(drugs$drug, function(x) x$mechanismOfAction),
      phase = drugs$phase,
      stringsAsFactors = FALSE
    )
    drug_results[[prot]] <- drugs_df
    cat(" ", nrow(drugs_df), "drugs\n")
  } else {
    cat(" no drugs\n")
  }
  Sys.sleep(0.5)
}

# ── 4. Directional concordance ──
df_drugs <- bind_rows(drug_results)
df_tract <- bind_rows(tractability_results)

if (nrow(df_drugs) > 0) {
  df_drugs <- df_drugs %>%
    left_join(ptsd_mr %>% select(protein, mr_direction), by = "protein") %>%
    mutate(
      drug_direction = case_when(
        grepl("agonist|activator", moa, ignore.case = TRUE)   ~ "Activator",
        grepl("inhibitor|antagonist", moa, ignore.case = TRUE) ~ "Inhibitor",
        TRUE ~ "Other"
      ),
      concordant = case_when(
        mr_direction == "Protective" & drug_direction == "Activator" ~ TRUE,
        mr_direction == "Risk" & drug_direction == "Inhibitor"       ~ TRUE,
        TRUE ~ FALSE
      )
    )
}

# ── 5. 保存 ──
cat("\nSaving results...\n")
write_csv(df_tract, file.path(result_dir, "drug_target_summary.csv"))
write_csv(df_drugs, file.path(result_dir, "drug_full_list.csv"))
write_csv(df_drugs %>% filter(phase >= 2),
          file.path(result_dir, "Table3_drug_repurposing.csv"))

cat("  ✓ drug_target_summary.csv\n")
cat("  ✓ drug_full_list.csv\n")
cat("  ✓ Table3_drug_repurposing.csv\n")
