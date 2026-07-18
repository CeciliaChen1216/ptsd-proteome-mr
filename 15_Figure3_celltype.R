source("00_config.R")
###############################################################################
# 15_Figure3_celltype.R  (now Supplementary Figure S8)
#
# Cell-type expression heatmap (Supplementary Figure S8): baseline HPA single-
# nucleus expression for CD40, FURIN, SIRPA (Panel A) and PTSD-associated
# cell-type differential expression from Daskalakis et al. 2024 (Panel B).
###############################################################################

library(ggplot2)
library(patchwork)
library(dplyr)

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_fig8 <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 6)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "bottom",
    plot.margin = margin(5, 8, 5, 5)
  )

# 4 cell classes: Neuronal / Glial / Immune / Vascular
col_neuronal <- "#FFD54F"
col_glial    <- "#6D4C41"
col_immune   <- "#E65100"   # distinct orange for microglia
col_vascular <- "#7E57C2"

# ═══════════════════════════════════════════════════════
# A  Cell-type expression in human brain
# ═══════════════════════════════════════════════════════

ct_order <- c("Ex neurons", "In neurons", "Astrocyte",
              "Microglia", "Oligo", "OPC",
              "Endothelial", "Pericyte")

hpa <- data.frame(
  celltype = rep(ct_order, 2),
  gene = rep(c("CD40", "SIRPA"), each = 8),
  nCPM = c(
    1.5, 1.0, 6.0, 29.2, 0.5, 0.5, 22.0, 32.3,
    100, 80, 180, 120, 60, 35, 3, 3
  ),
  cell_class = rep(c("Neuronal","Neuronal","Glial","Immune",
                      "Glial","Glial","Vascular","Vascular"), 2)
)
hpa$celltype <- factor(hpa$celltype, levels = ct_order)
hpa$gene <- factor(hpa$gene, levels = c("CD40", "SIRPA"))
hpa$cell_class <- factor(hpa$cell_class,
                         levels = c("Neuronal","Glial","Immune","Vascular"))

pA <- ggplot(hpa, aes(x = celltype, y = nCPM, fill = cell_class)) +
  geom_col(width = 0.65, color = "grey30", linewidth = 0.3) +
  facet_wrap(~gene, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Neuronal" = col_neuronal,
                                "Glial"    = col_glial,
                                "Immune"   = col_immune,
                                "Vascular" = col_vascular),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "A  Cell-type expression in human brain",
       x = NULL, y = "nCPM") +
  theme_fig8 +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8.5),
    strip.text = element_text(face = "bold.italic", size = 11),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.background = element_blank()
  )

# ═══════════════════════════════════════════════════════
# B  PTSD-related cell-type dysregulation
# Unified y-axis (scales = "fixed") for fair comparison
# ═══════════════════════════════════════════════════════

ct7 <- c("Astro", "Endo", "Ex", "In", "Microglia", "Oligo", "OPC")

ptsd_de <- data.frame(
  celltype = rep(ct7, 2),
  gene = rep(c("CD40", "FURIN"), each = 7),
  beta = c(
    0.303, 0.054, NA, NA, 0.277, NA, NA,
    -0.613, 1.243, -0.202, 0.182, -0.057, -0.104, 0.183
  ),
  pval = c(
    0.463, 0.778, NA, NA, 0.347, NA, NA,
    0.057, 0.014, 0.182, 0.417, 0.860, 0.680, 0.482
  )
)
ptsd_de$celltype <- factor(ptsd_de$celltype, levels = ct7)
ptsd_de$gene <- factor(ptsd_de$gene, levels = c("CD40", "FURIN"))

det <- ptsd_de %>% filter(!is.na(beta))
det$sig <- case_when(
  det$pval < 0.05 ~ "P < 0.05",
  det$pval < 0.10 ~ "P < 0.10",
  TRUE ~ "n.s."
)

pB <- ggplot(det, aes(x = celltype, y = beta, fill = sig)) +
  geom_col(width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  facet_wrap(~gene, ncol = 2) +   # scales = "fixed" (default)
  scale_fill_manual(
    values = c("P < 0.05" = "#D73027",
               "P < 0.10" = "#FC8D59",
               "n.s."     = "#B0B0B0"),
    name = NULL
  ) +
  labs(title = "B  PTSD-related cell-type dysregulation",
       x = NULL,
       y = expression(beta~"(PTSD vs control)")) +
  theme_fig8 +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 8.5),
    strip.text = element_text(face = "bold.italic", size = 11),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.background = element_blank()
  )

# ═══════════════════════════════════════════════════════
# Assemble: A stacked on B (Panel C removed in current revision)
# ═══════════════════════════════════════════════════════

fig8 <- pA / pB

ggsave(file.path(fig_dir, "SuppFig8_celltype.pdf"),
       fig8, width = 7, height = 7.5)
ggsave(file.path(fig_dir, "SuppFig8_celltype.tiff"),
       fig8, width = 7, height = 7.5, dpi = 300, compression = "lzw")

cat("\u2713 Supplementary Figure S8 (cell-type) saved to", fig_dir, "\n")
