###############################################################################
# 09_figures.R — Publication-ready Figures
#
# Generates the panel content for several manuscript main figures. Output
# filenames preserve historical (draft) figure numbers; the mapping to the
# final manuscript numbering is documented in README.md.
#
# Design principles:
#   - Results read from upstream CSV/RDS where possible; curated fallbacks noted
#   - CD40 emphasis via ggtext::element_markdown(), not vectorized face
#   - Colors: blue #4575B4 (transdiagnostic), red #D73027 (predominance)
###############################################################################

source("00_config.R")
library(ggplot2)
library(ggrepel)
library(ggtext)
library(patchwork)
library(dplyr)
library(tidyr)

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Data ──
mr_all  <- readRDS(file.path(result_dir, "mr_all_outcomes.rds"))
ptsd_mr <- mr_all %>% filter(grepl("PTSD", outcome))
cand    <- read.csv(file.path(result_dir, "ptsd_candidates_classified.csv"))
ev      <- read.csv(file.path(result_dir, "evidence_summary_table.csv"))

col_t <- "#4575B4"
col_p <- "#D73027"

# ── Helper: markdown label for y-axis (bold CD40 only) ──
md_label <- function(x) ifelse(x == "CD40", "**CD40**", x)

# ════════════════════════════════════════
# Fig 3: PTSD proteome-wide MR volcano (3A) + per-disorder forest (3B)
# Declassified: no power-aware / transdiagnostic classification.
# ════════════════════════════════════════
cat("\nFig 3: Volcano (3A) + per-disorder forest (3B)\n")

col_fdr <- "#3b6ea5"   # blue = FDR-significant
col_ns  <- "#bdbdbd"   # grey = not significant

## ---- 2A: PTSD volcano ----
pd <- ptsd_mr %>%
  mutate(nlp = -log10(mr_pval),
         grp = ifelse(fdr < 0.05, "FDR<0.05", "NS"))
bonf      <- 0.05 / nrow(pd)
fdr_line  <- -log10(max(pd$mr_pval[pd$fdr < 0.05], na.rm = TRUE))
hits      <- pd %>% filter(fdr < 0.05)
bonf_hits <- pd %>% filter(mr_pval < bonf)

p2 <- ggplot(pd, aes(mr_beta, nlp)) +
  geom_point(aes(color = grp), size = 1.3, alpha = 0.6) +
  scale_color_manual(values = c("FDR<0.05" = col_fdr, "NS" = col_ns), name = NULL) +
  geom_point(data = bonf_hits, shape = 21, fill = NA, color = "black",
             size = 3.2, stroke = 1) +
  geom_hline(yintercept = fdr_line, linetype = "dashed", color = col_fdr, linewidth = 0.4) +
  geom_hline(yintercept = -log10(bonf), linetype = "dotted", color = "black", linewidth = 0.4) +
  geom_text_repel(data = hits, aes(label = protein), size = 3.4, fontface = "bold",
                  seed = 42, min.segment.length = 0, box.padding = 0.4) +
  annotate("text", x = max(pd$mr_beta, na.rm = TRUE), y = fdr_line,
           label = "FDR 0.05", hjust = 1, vjust = -0.4, size = 3, color = col_fdr) +
  annotate("text", x = max(pd$mr_beta, na.rm = TRUE), y = -log10(bonf),
           label = "Bonferroni", hjust = 1, vjust = -0.4, size = 3, color = "black") +
  labs(title = "PTSD proteome-wide MR (1,258 proteins)",
       x = "MR effect on PTSD liability (log-OR per SD pQTL)",
       y = expression(-log[10](italic(P)))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = c(0.12, 0.9))

ggsave(file.path(fig_dir, "Figure3A_volcano.pdf"), p2, width = 7, height = 5.5)
ggsave(file.path(fig_dir, "Figure3A_volcano.tiff"), p2, width = 7, height = 5.5,
       dpi = 300, compression = "lzw")

## ---- 2B: per-disorder MR estimates for the prioritized proteins ----
cat("Fig 3B: per-disorder forest\n")
cands <- unique(cand$protein)
lab   <- c(PTSD_freeze3 = "PTSD", MDD_adams2025 = "MDD", ANX_strom2026 = "Anxiety")
b2 <- mr_all %>%
  filter(protein %in% cands) %>%
  mutate(disorder = factor(lab[outcome], levels = c("PTSD", "MDD", "Anxiety")),
         lo = log(OR_lower), hi = log(OR_upper), sigflag = fdr < 0.05,
         protein = factor(protein, levels = rev(cands)))

p3 <- ggplot(b2, aes(mr_beta, protein, color = disorder, group = disorder)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                 position = position_dodge(width = 0.6), linewidth = 0.6) +
  geom_point(aes(fill = disorder),
             position = position_dodge(width = 0.6), size = 2.4, shape = 21,
             color = "white", stroke = 0.3) +
  geom_point(data = subset(b2, sigflag), aes(group = disorder),
             position = position_dodge(width = 0.6), size = 2.6, shape = 21,
             fill = NA, color = "black", stroke = 0.9) +
  scale_color_manual(values = c(PTSD = "#3b6ea5", MDD = "#e08214", Anxiety = "#4a9b6e"), name = NULL) +
  scale_fill_manual(values = c(PTSD = "#3b6ea5", MDD = "#e08214", Anxiety = "#4a9b6e"), guide = "none") +
  labs(title = "Per-disorder MR estimates for prioritized proteins",
       x = "MR effect (log-OR per SD increase in protein)", y = NULL,
       caption = "Dark-outlined point = FDR<0.05 in that disorder") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "Figure3B_forest.pdf"), p3, width = 7, height = 5)
ggsave(file.path(fig_dir, "Figure3B_forest.tiff"), p3, width = 7, height = 5,
       dpi = 300, compression = "lzw")

# -- Figure 3 combined (A volcano | B forest) --
Figure3 <- (p2 | p3) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(file.path(fig_dir, "Figure3_combined.pdf"),  Figure3, width = 13, height = 5.5)
ggsave(file.path(fig_dir, "Figure3_combined.tiff"), Figure3, width = 13, height = 5.5,
       dpi = 300, compression = "lzw")
cat("  -> Figure3_combined (A volcano | B forest)\n")

# ════════════════════════════════════════
# Fig 4 (draft) = Suppl Fig 2: Evidence matrix
# Revised: replaces the summed "pass/fail" tier check-mark figure with a
# four-state evidence matrix (Supportive / Cautionary / Inconclusive /
# Not testable). Cells are NOT summed into a composite score. The SMR+HEIDI
# column is filled from the real smr.exe re-analysis (eqtl_smr_brain_blood_fdr.csv);
# the other seven axes are taken from the verified supplementary tables.
# ════════════════════════════════════════
cat("Fig 4 / Suppl Fig 2: Evidence matrix\n")

# SMR+HEIDI column from real re-analysis: per gene take strongest verdict across tissues.
# FDR-aware (reviewer #8): brain and blood SMR were BH-FDR-corrected separately;
# a gene is scored Supportive only when it has an FDR-supported Strong (SMR+HEIDI)
# result in at least one tissue. HEIDI-rejected -> Cautionary; SMR-significant but
# nominal-only (not FDR-supported) -> Cautionary; SMR not significant -> Inconclusive.
smr_redone <- read.csv(file.path(result_dir, "eqtl_smr_brain_blood_fdr.csv"))
smr_col <- if ("verdict_fdr" %in% names(smr_redone)) "verdict_fdr" else "verdict"
smr_state <- function(g){
  v <- smr_redone[[smr_col]][smr_redone$gene == g]
  if (any(grepl("FDR-supported", v)))                    return("S")   # Strong + FDR
  if (any(grepl("rejected", v)))                         return("C")   # HEIDI rejected
  if (any(grepl("nominal only", v)))                     return("C")   # sig but not FDR
  if (any(grepl("^Strong", v)))                          return("S")   # back-compat (old verdict col)
  if (any(grepl("Not significant", v)))                  return("I")
  "N"
}
prot_order <- c("KHK","CGREF1","AKT3","FES","UBE2L6","CD40","FURIN","SNX18","SIRPA","CD101")
axes <- c("Primary MR","Steiger","Pleiotropy","Reverse MR",
          "Multi-IV IVW","SMR+HEIDI","Colocalization","Brain pQTL/PWAS")
# Seven hardcoded axes from verified tables; column 6 (SMR) overwritten by real data.
M <- rbind(
  KHK    = c("S","S","S","C","S",NA,"S","S"),
  CGREF1 = c("S","S","S","C","S",NA,"S","C"),
  AKT3   = c("S","S","S","S","S",NA,"S","C"),
  FES    = c("S","S","S","C","S",NA,"N","I"),
  UBE2L6 = c("S","S","S","C","S",NA,"C","I"),
  CD40   = c("S","S","S","C","S",NA,"N","S"),
  FURIN  = c("S","S","S","S","S",NA,"N","N"),
  SNX18  = c("S","S","S","S","S",NA,"N","I"),
  SIRPA  = c("S","S","S","C","S",NA,"N","C"),
  CD101  = c("S","S","C","S","S",NA,"N","N")
)
M[,6] <- sapply(rownames(M), smr_state)
colnames(M) <- axes

df4 <- as.data.frame(as.table(M)); names(df4) <- c("Protein","Axis","state")
df4$Protein <- factor(df4$Protein, levels = rev(prot_order))
df4$Axis    <- factor(df4$Axis, levels = axes)
state_lab <- c(S="Supportive", C="Cautionary", I="Inconclusive", N="Not testable")
df4$state <- factor(state_lab[as.character(df4$state)],
                    levels = c("Supportive","Cautionary","Inconclusive","Not testable"))
pal4 <- c("Supportive"="#2C7FB8","Cautionary"="#D95F0E",
          "Inconclusive"="#BDBDBD","Not testable"="#F0F0F0")

p4 <- ggplot(df4, aes(Axis, Protein, fill = state)) +
  geom_tile(color = "white", linewidth = 1.1) +
  scale_fill_manual(values = pal4, name = NULL, drop = FALSE) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  # coord_equal() removed: it squeezed the panel and pushed the rotated top
  # labels up into the title (the "串行" overlap). Aspect set via ggsave instead.
  labs(title = "Evidence matrix for ten PTSD-prioritized proteins",
       subtitle = "Eight complementary evidence domains; cells are not summed into a composite score",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x.top = element_text(angle = 0, hjust = 0.5, vjust = 0, face = "bold", size = 8.5),
    axis.text.y     = element_text(face = "bold.italic"),
    panel.grid      = element_blank(),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, color = "grey30"),
    plot.margin     = margin(t = 14, r = 20, b = 8, l = 8)   # horizontal top labels need little headroom
  )

ggsave(file.path(fig_dir, "SuppFig2_evidence_matrix.pdf"),  p4, width = 11.5, height = 8.2)
ggsave(file.path(fig_dir, "SuppFig2_evidence_matrix.tiff"), p4, width = 11.5, height = 8.2,
       dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "SuppFig2_evidence_matrix.png"),  p4, width = 11.5, height = 8.2, dpi = 200)
cat("  -> SuppFig2_evidence_matrix (.pdf/.tiff/.png)\n")

# ════════════════════════════════════════
# Fig 5: Cross-disorder
# ════════════════════════════════════════
cat("\nSuppFig1: cross-disorder pathway (NES) similarity\n")

nf <- file.path(result_dir, "hallmark_NES_comparison.csv")
if (file.exists(nf)) {
  ns <- read.csv(nf)
  pc <- grep("ptsd|PTSD", names(ns), ignore.case = TRUE, value = TRUE)[1]
  mc <- grep("mdd|MDD", names(ns), ignore.case = TRUE, value = TRUE)[1]
  ac <- grep("anx|ANX", names(ns), ignore.case = TRUE, value = TRUE)[1]
  rm <- cor(ns[[pc]], ns[[mc]], use = "complete.obs", method = "spearman")
  ra <- cor(ns[[pc]], ns[[ac]], use = "complete.obs", method = "spearman")

  p5B <- ggplot(ns, aes(.data[[pc]], .data[[mc]])) +
    geom_smooth(method = "lm", se = TRUE, color = col_t,
                fill = col_t, alpha = 0.08, linewidth = 0.6) +
    geom_point(size = 1.2, color = col_t, alpha = 0.4) +
    annotate("text", x = min(ns[[pc]], na.rm = TRUE), y = Inf,
      label = paste0("\u03C1 = ", round(rm, 3)),
      size = 3.8, hjust = 0, vjust = 1.5, fontface = "bold", color = "grey20") +
    labs(title = "A  PTSD vs MDD pathway similarity",
         x = "NES (PTSD)", y = "NES (MDD)") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

  p5C <- ggplot(ns, aes(.data[[pc]], .data[[ac]])) +
    geom_smooth(method = "lm", se = TRUE, color = col_p,
                fill = col_p, alpha = 0.08, linewidth = 0.6) +
    geom_point(size = 1.2, color = col_p, alpha = 0.4) +
    annotate("text", x = min(ns[[pc]], na.rm = TRUE), y = Inf,
      label = paste0("\u03C1 = ", round(ra, 3)),
      size = 3.8, hjust = 0, vjust = 1.5, fontface = "bold", color = "grey20") +
    labs(title = "B  PTSD vs anxiety pathway similarity",
         x = "NES (PTSD)", y = "NES (Anxiety)") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 11))

  f5 <- p5B | p5C
  ggsave(file.path(fig_dir, "SuppFig1_cross_disorder.pdf"), f5, width = 10, height = 5.5)
  ggsave(file.path(fig_dir, "SuppFig1_cross_disorder.tiff"), f5, width = 10, height = 5.5,
         dpi = 300, compression = "lzw")
}

# ════════════════════════════════════════
# Fig 6: Drug repurposing (data-driven)
# ════════════════════════════════════════
cat("Fig 6: Drug repurposing\n")

drug_file <- file.path(result_dir, "drug_target_summary.csv")
tract_file <- file.path(result_dir, "Table3_drug_repurposing.csv")

if (file.exists(drug_file)) {
  dd <- read.csv(drug_file)
  cat("Drug summary columns:", paste(names(dd), collapse = ", "), "\n")
  # Identify columns
  prot_d <- grep("protein|target|gene", names(dd), ignore.case = TRUE, value = TRUE)[1]
  ndrug_d <- grep("n_drug|count|total", names(dd), ignore.case = TRUE, value = TRUE)[1]
  dir_d <- grep("direct|protect|risk|mr_dir", names(dd), ignore.case = TRUE, value = TRUE)[1]
  
  if (!is.na(prot_d) && !is.na(ndrug_d)) {
    dd$protein <- dd[[prot_d]]
    dd$n_drugs <- as.numeric(dd[[ndrug_d]])
    if (!is.na(dir_d)) {
      dd$direction <- dd[[dir_d]]
    } else {
      # Derive from OR
      or_map <- setNames(cand$OR, cand$protein)
      dd$direction <- ifelse(or_map[dd$protein] < 1, "Protective", "Risk")
    }
  }
} else {
  # Fallback: build from known data
  or_map <- setNames(cand$OR, cand$protein)
  dd <- data.frame(
    protein = cand$protein,
    n_drugs = 0,
    direction = ifelse(or_map[cand$protein] < 1, "Protective", "Risk"))
  cat("NOTE: drug_target_summary.csv not found, using placeholder\n")
}

# Sort: CD40 bottom (visual endpoint), AKT3 next
pl6 <- c("CD40", "AKT3", setdiff(sort(dd$protein), c("CD40", "AKT3")))
dd <- dd %>% filter(protein %in% pl6)
dd$protein <- factor(dd$protein, levels = pl6)
dd$ba <- ifelse(dd$n_drugs == 0, 0.15, 0.85)
dd$label <- factor(as.character(dd$protein), levels = pl6)

p6A <- ggplot(dd, aes(n_drugs, label, fill = direction, alpha = ba)) +
  geom_col(width = 0.6, color = "grey30", linewidth = 0.3) +
  scale_alpha_identity() +
  geom_text(data = dd %>% filter(n_drugs > 0),
            aes(label = n_drugs), hjust = -0.3, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(Protective = col_t, Risk = col_p),
    name = NULL, labels = c("Protective association", "Risk association")) +
  scale_x_continuous(limits = c(0, 15), breaks = seq(0, 15, 5)) +
  labs(title = "A  Known drugs for prioritized targets",
       x = "Number of known drugs", y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom", legend.text = element_text(size = 8))

# Tractability panel — try to read from Table3, fallback to curated values
tract_built <- FALSE
if (file.exists(tract_file)) {
  tr <- read.csv(tract_file)
  cat("Tractability columns:", paste(names(tr), collapse = ", "), "\n")
  prot_t <- grep("protein|target|gene", names(tr), ignore.case = TRUE, value = TRUE)[1]
  sm_col <- grep("small.mol|SM|bucket_sm", names(tr), ignore.case = TRUE, value = TRUE)[1]
  ab_col <- grep("antibod|Ab|bucket_ab", names(tr), ignore.case = TRUE, value = TRUE)[1]
  pr_col <- grep("protac|PROTAC|bucket_pr", names(tr), ignore.case = TRUE, value = TRUE)[1]
  ot_col <- grep("other|bucket_ot", names(tr), ignore.case = TRUE, value = TRUE)[1]

  if (!is.na(prot_t) && !is.na(sm_col) && !is.na(ab_col)) {
    cat("  Building tractability from Table3\n")
    tm <- tr %>%
      filter(.data[[prot_t]] %in% pl6) %>%
      select(protein = all_of(prot_t),
             `Small molecule` = all_of(sm_col),
             Antibody = all_of(ab_col)) %>%
      mutate(PROTAC = if (!is.na(pr_col)) tr[[pr_col]][match(protein, tr[[prot_t]])] else FALSE,
             Other = if (!is.na(ot_col)) tr[[ot_col]][match(protein, tr[[prot_t]])] else FALSE) %>%
      pivot_longer(-protein, names_to = "modality", values_to = "tractable") %>%
      mutate(tractable = as.logical(tractable))
    tract_built <- TRUE
  }
}

if (!tract_built) {
  # FALLBACK: curated from Open Targets (accessed 2024-12).
  # If Table3 format changes, update this lookup or fix column matching above.
  cat("  NOTE: Table3 parsing failed; using curated fallback values (Open Targets 2024-12)\n")
  tract_known <- list(
    CD40  = c(FALSE, TRUE, FALSE, TRUE),
    AKT3  = c(TRUE, FALSE, TRUE, TRUE),
    CGREF1 = c(FALSE, FALSE, FALSE, FALSE),
    FES   = c(TRUE, FALSE, TRUE, FALSE),
    FURIN = c(TRUE, FALSE, FALSE, FALSE),
    KHK   = c(TRUE, FALSE, FALSE, FALSE),
    SNX18 = c(FALSE, FALSE, FALSE, FALSE),
    UBE2L6 = c(FALSE, FALSE, FALSE, FALSE),
    SIRPA = c(FALSE, TRUE, FALSE, FALSE),
    CD101 = c(FALSE, TRUE, FALSE, FALSE))
  tm <- data.frame(
    protein = rep(pl6, each = 4),
    modality = rep(c("Small molecule", "Antibody", "PROTAC", "Other"), length(pl6)),
    tractable = unlist(tract_known[pl6]))
}

tm$protein <- factor(tm$protein, levels = pl6)
tm$modality <- factor(tm$modality,
  levels = c("Small molecule", "Antibody", "PROTAC", "Other"))
# Plain protein names on both drug-panel axes; markdown/HTML emphasis did not
# render reliably here under patchwork (printed literal ** / <b> tags).
tm$label <- factor(as.character(tm$protein), levels = pl6)

p6B <- ggplot(tm, aes(modality, label, fill = tractable)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(tractable, "\u2713", "")),
            size = 3.5, color = "white", fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "#2E7D6F", "FALSE" = "#EBEBEB"),
    name = NULL, labels = c("Not tractable", "Tractable")) +
  labs(title = "B  Tractability across modalities", x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.text.x = element_text(size = 8, angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom", legend.text = element_text(size = 8))

f6 <- p6A | p6B
f6 <- f6 + plot_layout(widths = c(0.55, 0.45))
ggsave(file.path(fig_dir, "SuppFig7_drug_repurposing.pdf"), f6, width = 10, height = 5)
ggsave(file.path(fig_dir, "SuppFig7_drug_repurposing.tiff"), f6, width = 10, height = 5,
       dpi = 300, compression = "lzw")

# ════════════════════════════════════════
# Fig 7A (draft) = Fig 4A: eQTL-SMR
# Revised: real smr.exe results with FINAL verdicts. Star = SMR significant
# with HEIDI support; triangle = HEIDI rejected (possible linkage);
# blank = not significant or no cis-eQTL. Markers drawn with geom_point
# (shape 8 / 17) so rendering does not depend on system fonts.
# ════════════════════════════════════════
cat("Fig 7A / Fig 4A: eQTL-SMR\n")

sm <- read.csv(file.path(result_dir, "eqtl_smr_brain_blood_fdr.csv"))
go7 <- c("CD40","FURIN","SIRPA","CGREF1","FES","KHK","SNX18","UBE2L6","AKT3","CD101")
sm <- sm %>% filter(gene %in% go7)
sm$tissue_lab <- factor(ifelse(sm$tissue == "blood", "Blood (Westra)", "Brain (BrainMeta)"),
                        levels = c("Blood (Westra)", "Brain (BrainMeta)"))
sm$nlp <- ifelse(is.na(sm$p_SMR), NA, -log10(as.numeric(sm$p_SMR)))
sm$mark <- dplyr::case_when(
  grepl("^Strong", sm$verdict)  ~ "star",
  grepl("rejected", sm$verdict) ~ "tri",
  TRUE                          ~ "")
gene_lv <- rev(c("UBE2L6","SNX18","KHK","FES","CGREF1","SIRPA","FURIN","CD40","AKT3","CD101"))
sm$gene <- factor(sm$gene, levels = gene_lv)

p7A <- ggplot(sm, aes(tissue_lab, gene, fill = nlp)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_point(data = subset(sm, mark == "star"), shape = 8,  size = 3,   color = "white", stroke = 1.1) +
  geom_point(data = subset(sm, mark == "tri"),  shape = 17, size = 2.6, color = "white") +
  scale_fill_gradient(low = "#EDF8E9", high = "#1B7837",
                      name = expression(-log[10](P[SMR])), na.value = "grey90") +
  scale_x_discrete(position = "top") +
  labs(title = "Brain and blood eQTL-SMR", x = NULL, y = NULL,
       caption = "* SMR significant + HEIDI supported    \u25B3 HEIDI rejected") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 10)),
        axis.text.x.top = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 9, face = "bold"),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 7, hjust = 0, color = "grey30",
                                    margin = margin(t = 8)),
        legend.position = "bottom",
        legend.box.margin = margin(t = 4),
        legend.key.width = unit(1.2, "cm"),
        legend.key.height = unit(0.3, "cm"),
        plot.margin = margin(t = 5, r = 8, b = 12, l = 8))

ggsave(file.path(fig_dir, "Figure4A_eqtl_smr.pdf"), p7A, width = 5.5, height = 5.5)
ggsave(file.path(fig_dir, "Figure4A_eqtl_smr.tiff"), p7A, width = 5.5, height = 5.5,
       dpi = 300, compression = "lzw")

# ════════════════════════════════════════
# Fig 7B (draft) = Fig 4B: Network MR
# Revised: neutral title (the chr15 conditional analysis reassigns the signal
# to FES, so FURIN->CD40 is reported as nominal rather than "identified").
# ════════════════════════════════════════
cat("Fig 7B / Fig 4B: Network MR\n")

nt <- read.csv(file.path(result_dir, "network_mr_results.csv"))
nt <- nt %>%
  mutate(pair = paste0(upstream, " \u2192 ", downstream),
         bl = beta - 1.96 * se, bu = beta + 1.96 * se,
         sg = ifelse(significant, "Nominal (P<0.05, not Bonferroni)", "Not significant")) %>%
  arrange(desc(significant), desc(abs(beta)))
nt$pair <- factor(nt$pair, levels = rev(unique(nt$pair)))
nt$plb <- ifelse(nt$significant,
  paste0("P = ", format(round(nt$pval, 3), nsmall = 3)), "")

p7B <- ggplot(nt, aes(beta, pair, color = sg)) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_segment(aes(x = bl, xend = bu, y = pair, yend = pair), linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = plb), hjust = 1.15, vjust = -0.8,
            size = 3.2, show.legend = FALSE) +
  scale_color_manual(values = c("Nominal (P<0.05, not Bonferroni)" = "#737373",
    "Not significant" = "#BDBDBD"), name = NULL) +
  labs(title = "Hypothesis-driven network MR of prioritized protein pairs",
       x = expression(Causal~effect~(beta %+-% 95*"% CI")), y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        legend.position = "bottom",
        legend.box.margin = margin(t = 4),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "lines"),
        plot.margin = margin(t = 5, r = 8, b = 12, l = 8))

ggsave(file.path(fig_dir, "Figure4B_network_mr.pdf"), p7B, width = 6, height = 4)
ggsave(file.path(fig_dir, "Figure4B_network_mr.tiff"), p7B, width = 6, height = 4,
       dpi = 300, compression = "lzw")

# ── Figure 4 combined (A eQTL-SMR | B network MR, side by side) ──
Figure4 <- p7A + p7B +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))
ggsave(file.path(fig_dir, "Figure4_combined.pdf"),  Figure4, width = 11, height = 6.0)
ggsave(file.path(fig_dir, "Figure4_combined.tiff"), Figure4, width = 11, height = 6.0,
       dpi = 300, compression = "lzw")
cat("  -> Figure4_combined (A eQTL-SMR | B network MR)\n")

# ════════════════════════════════════════
# Fig 7C: CRP mediation
# ════════════════════════════════════════
cat("Fig 7C: CRP mediation\n")

md <- read.csv(file.path(result_dir, "mediation_crp_results.csv"))
md <- md %>%
  mutate(ap = abs(prop_mediated), ic = protein == "CD40") %>%
  arrange(ic, ap) %>%
  mutate(label = md_label(protein),
         label = factor(label, levels = label))

p7C <- ggplot(md, aes(prop_mediated * 100, label)) +
  geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_col(width = 0.5, fill = "#A0A0A0", color = "grey50", linewidth = 0.3) +
  labs(title = "No evidence of CRP-mediated effects",
       subtitle = "Global CRP \u2192 PTSD MR: P = 0.41",
       x = "% of total effect mediated via CRP", y = NULL) +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5,
                                     color = "grey40", face = "italic"),
        axis.text.y = element_markdown(size = 10))

ggsave(file.path(fig_dir, "SuppFig6_crp_mediation.pdf"), p7C, width = 5.5, height = 4)
ggsave(file.path(fig_dir, "SuppFig6_crp_mediation.tiff"), p7C, width = 5.5, height = 4,
       dpi = 300, compression = "lzw")

# ════════════════════════════════════════
# Supplementary Figure S8 (cell-type): see 15_Figure3_celltype.R
# ════════════════════════════════════════
cat("Supplementary Figure S8 (cell-type): see 15_Figure3_celltype.R\n")

cat("\n\u2713 All figures saved to", fig_dir, "\n")
