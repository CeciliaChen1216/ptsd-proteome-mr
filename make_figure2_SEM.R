###############################################################################
# make_figure2_SEM.R  —  MAIN Figure 2: common internalizing factor
#
# A  Standardized loadings of PTSD / MDD / anxiety on the common factor
#    (loadings as bars; residual variance proportions annotated as TEXT beside
#    the panel — NOT drawn below zero, to avoid implying negative variance).
# B  Protein -> common-factor MR (391 testable; 28 FDR<0.05; six candidates).
# C  Instrument-level deviation from the common-factor model (Q_SNP): ALL 391
#    proteins in grey; six PTSD-prioritized candidates highlighted; nominal and
#    Bonferroni reference lines.
#
# Reads stage12_factorMR_results.csv (userGWAS constrained model, real Q_SNP).
# Requires: ggplot2, patchwork, data.table, ggrepel.
###############################################################################

source("00_config.R")
suppressMessages({library(ggplot2); library(patchwork); library(data.table); library(ggrepel)})

df <- tryCatch(fread(file.path(result_dir, "genomicsem", "stage12_factorMR_results.csv")),
               error = function(e) data.table())
if (!nrow(df)) df <- fread(file.path(result_dir, "stage12_factorMR_results.csv"))

df[, mr_fdr := p.adjust(mr_p, "BH")]
df[, nlp    := -log10(mr_p)]
df[, q_nlp  := -log10(Q_pval)]
q_nominal <- sum(df$Q_pval < 0.05, na.rm = TRUE)
q_fdr     <- sum(p.adjust(df$Q_pval, "BH") < 0.05, na.rm = TRUE)
bonf      <- -log10(0.05 / nrow(df))          # Bonferroni line ~ 3.89

BLUE <- "#3b6ea5"; GREY <- "#c4c4c4"; RED <- "#c0504d"; ORANGE <- "#d98c3f"
cand <- c("AKT3","CD40","CGREF1","FES","KHK","SNX18")

## ---- Panel A: loadings as bars; residuals as TEXT (method B) ----------------
# Loadings and residuals from the updated Stage 0 constrained model
# (LDSCoutput_stage0.rds regenerated with corrected PTSD sample.prev = 0.1121).
la <- data.table(
  trait = factor(c("PTSD","MDD","Anxiety"), levels = c("PTSD","MDD","Anxiety")),
  load  = c(0.921, 1.018, 0.893))
# Residual note in the upper-right, using plotmath for a device-independent
# superscript in the P-value (unicode superscripts do not render on some R
# graphics devices and can appear as "10.4").
pA <- ggplot(la, aes(trait, load)) +
  geom_col(fill = BLUE, width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", load)), vjust = -0.6, size = 3.9) +
  # light background box (drawn first) then stacked text lines
  annotate("rect", xmin = 0.48, xmax = 2.52, ymin = 0.02, ymax = 0.265,
           fill = "#f7f7f7", colour = "#cccccc", linewidth = 0.3) +
  annotate("text", x = 0.58, y = 0.235, hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333", fontface = "bold", label = "Residual variance proportion") +
  annotate("text", x = 0.58, y = 0.175, hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333", parse = TRUE,
           label = "PTSD*':'~0.15*','~italic(P)==1.2%*%10^-4") +
  annotate("text", x = 0.58, y = 0.115, hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333", label = "MDD: fixed to 0") +
  annotate("text", x = 0.58, y = 0.055, hjust = 0, vjust = 1, size = 2.7,
           colour = "#333333", label = "Anxiety: 0.20") +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.5, 1.0),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A  Common internalizing factor",
       x = NULL, y = "Standardized loading") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        axis.text.x = element_text(face = "bold"))

## ---- Panel B: factor-MR volcano --------------------------------------------
df[, grp := ifelse(mr_fdr < 0.05, "FDR<0.05", "NS")]
lab_df <- df[protein %in% cand]

pB <- ggplot(df, aes(mr_beta, nlp)) +
  geom_point(data = df[grp == "NS"], colour = GREY, alpha = 0.6, size = 1.4) +
  geom_point(data = df[grp == "FDR<0.05"], colour = BLUE, alpha = 0.85, size = 1.7) +
  geom_point(data = lab_df, size = 3.1, shape = 21, fill = BLUE, stroke = 0.9, colour = "black") +
  geom_text_repel(data = lab_df, aes(label = protein), size = 3.6,
                  box.padding = 0.8, point.padding = 0.5, min.segment.length = 0,
                  max.overlaps = 20, seed = 1) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  labs(title = "B  Protein-to-common-factor MR (391 proteins)",
       subtitle = paste0(sum(df$mr_fdr < 0.05), " proteins FDR < 0.05"),
       x = "MR effect on common factor (Wald ratio)",
       y = expression(-log[10]~italic(P))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.subtitle = element_text(size = 9, colour = "grey30", hjust = 0.5))

## ---- Panel C: distribution of Q_SNP deviation across all 391 proteins ------
# Histogram of -log10(P_QSNP) for all testable proteins (y = protein count, a
# quantitative axis). The six PTSD-prioritized candidates are marked on the axis
# by a rug + labels. This shows directly that the great majority of proteins have
# small Q_SNP deviation and that none crosses the Bonferroni threshold.
#
# In the updated (corrected sample.prev) analysis only SNX18 crosses the nominal
# threshold; the remaining five candidates cluster together below nominal. To
# avoid overlapping labels we label SNX18 (crosses nominal) and FES (highest of
# the remaining five, useful anchor); the other four are shown as rug ticks only.
qc <- df[protein %in% cand][order(q_nlp)]
label_set <- c("SNX18","FES")
qc[, lab := ifelse(protein %in% label_set, protein, "")]
# staggered label heights (fraction of ymax); only the two labelled candidates
# need explicit positions.
lab_h_map <- c(SNX18 = 0.85, FES = 0.70)

ymax <- max(hist(df$q_nlp, breaks = 30, plot = FALSE)$counts)
qc[, lab_h := ymax * ifelse(protein %in% names(lab_h_map), lab_h_map[protein], 0.85)]

pC <- ggplot(df, aes(q_nlp)) +
  geom_histogram(bins = 30, fill = "#9fb8d4", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "#666666", linewidth = 0.7) +
  geom_vline(xintercept = bonf, linetype = "dotted", colour = "#444444", linewidth = 0.7) +
  # candidate rug ticks just above the x-axis
  geom_point(data = qc, aes(q_nlp, -ymax*0.045), colour = "black", fill = ORANGE,
             shape = 25, size = 2.6, stroke = 0.5) +
  # labels for four candidates, vertically staggered so KHK/CD40 (adjacent on x)
  # do not collide; short dotted leaders connect each label to its rug marker.
  geom_segment(data = qc[lab != ""],
               aes(x = q_nlp, xend = q_nlp, y = -ymax*0.02, yend = lab_h),
               colour = "#bbbbbb", linewidth = 0.25, linetype = "dotted") +
  geom_text(data = qc[lab != ""], aes(q_nlp, lab_h, label = lab),
            size = 3.1, fontface = "italic", vjust = 0,
            hjust = ifelse(qc[lab != ""]$q_nlp > 2.5, 1.05, -0.05)) +
  annotate("text", x = -log10(0.05), y = ymax*1.02, label = "nominal\nP = 0.05",
           size = 2.9, colour = "#666666", vjust = 1, hjust = 1.05) +
  annotate("text", x = bonf, y = ymax*1.02, label = "Bonferroni\nP = 0.05/391",
           size = 2.9, colour = "#444444", vjust = 1, hjust = 1.05) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(title = expression(bold("C  Instrument-level Q"[SNP]*" deviation")),
       x = expression(-log[10]~italic(P)[Q[SNP]]), y = "Number of proteins") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

## ---- assemble --------------------------------------------------------------
Figure2 <- pA + pB + pC + plot_layout(widths = c(0.9, 1.5, 1.25))

ggsave(file.path(fig_dir, "Figure2_SEM.pdf"),  Figure2, width = 17.5, height = 5.6)
ggsave(file.path(fig_dir, "Figure2_SEM.tiff"), Figure2, width = 17.5, height = 5.6,
       dpi = 300, compression = "lzw")
ggsave(file.path(fig_dir, "Figure2_SEM.png"),  Figure2, width = 17.5, height = 5.6, dpi = 150)
cat(sprintf("Figure 2 (SEM) written. nominal Q_SNP=%d, FDR Q_SNP=%d, Bonferroni line=%.2f\n",
            q_nominal, q_fdr, bonf))
