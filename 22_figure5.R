###############################################################################
# 22_figure5_KHK_and_negative_control.R
#
# Manuscript Figure 5:
#   Panel A: KHK cross-tissue MR forest plot
#   Panel B: PTSD symptom severity effect vs childhood maltreatment effect
#
# This cleaned version is suitable for Additional file / GitHub release.
# It removes local absolute paths and uses relative paths from 00_config.R.
#
# Expected input:
#   - 00_config.R defining at least result_dir; optionally fig_dir
#   - results/trauma_negative_control_MR.csv, if available
#
# Outputs:
#   - figures/Figure_5_KHK_cross_tissue.tiff
#   - figures/Figure_5_KHK_cross_tissue.pdf
#   - figures/Figure_5_KHK_cross_tissue.png
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 0. Setup
# ─────────────────────────────────────────────────────────────────────────────

get_script_dir <- function() {
  for (n in rev(seq_len(sys.nframe()))) {
    f <- tryCatch(sys.frame(n)$ofile, error = function(e) NULL)
    if (!is.null(f) && is.character(f) && nzchar(f)) {
      return(dirname(normalizePath(f, mustWork = FALSE)))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
}

script_dir <- get_script_dir()

config_path <- file.path(script_dir, "00_config.R")
if (!file.exists(config_path)) {
  config_path <- file.path(getwd(), "00_config.R")
}
if (!file.exists(config_path)) {
  stop("Cannot find 00_config.R. Copy 00_config_template.R to 00_config.R and edit local paths.")
}
source(config_path)

required_pkgs <- c("ggplot2", "data.table", "patchwork")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(patchwork)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) suppressPackageStartupMessages(library(ggrepel))

if (!exists("result_dir")) {
  result_dir <- file.path(dirname(script_dir), "results")
}
if (!exists("fig_dir")) {
  fig_dir <- file.path(dirname(script_dir), "figures")
}
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

format_p <- function(p) {
  ifelse(p < 0.001, sprintf("P = %.1e", p), sprintf("P = %.3f", p))
}

message("Generating Figure 5: KHK cross-tissue MR + childhood-maltreatment secondary-phenotype analysis")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Panel A: KHK cross-tissue MR forest plot
# ─────────────────────────────────────────────────────────────────────────────

khk <- data.table(
  Analysis = c(
    "Plasma -> PCL",
    "Plasma -> PTSD",
    "Brain -> PTSD",
    "Plasma -> CM"
  ),
  beta = c(0.0144, 0.0376, 0.0108, -0.00302),
  se   = c(0.00296, 0.01050, 0.00276, 0.00819),
  p    = c(1.17e-6, 3.47e-4, 8.83e-5, 0.713),
  Evidence = c(
    "Plasma",
    "Plasma",
    "Brain",
    "Childhood maltreatment (secondary phenotype)"
  )
)

khk[, ci_low  := beta - 1.96 * se]
khk[, ci_high := beta + 1.96 * se]
khk[, p_label := format_p(p)]
khk[, Analysis := factor(Analysis, levels = rev(Analysis))]

x_range <- max(abs(c(khk$ci_low, khk$ci_high)))

p_A <- ggplot(khk, aes(x = beta, y = Analysis, color = Evidence)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.5
  ) +
  geom_errorbarh(
    aes(xmin = ci_low, xmax = ci_high),
    height = 0.25,
    linewidth = 0.8
  ) +
  geom_point(size = 4) +
  geom_text(
    aes(label = p_label, x = ci_high + x_range * 0.08),
    hjust = 0,
    size = 3.4,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Plasma" = "#1f77b4",
      "Brain" = "#d62728",
      "Childhood maltreatment (secondary phenotype)" = "#7f7f7f"
    ),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.34))) +
  labs(
    x = "Effect estimate (beta, 95% CI)\nEstimates are on different scales; compare direction only, not magnitude.",
    y = NULL,
    subtitle = "A. Cross-tissue evidence for KHK",
    caption = "Instrument: rs2304681 (p.Val49Ile missense; MAF ~37%); brain pQTL posterior inclusion probability 0.99"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.6, "cm"),
    plot.subtitle = element_text(face = "bold", size = 12),
    plot.caption = element_text(
      hjust = 0,
      size = 8.5,
      face = "italic",
      color = "gray30",
      margin = margin(t = 8)
    ),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 11),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 2. Panel B: PTSD symptom severity effect vs childhood maltreatment effect
# ─────────────────────────────────────────────────────────────────────────────

trauma_csv <- file.path(result_dir, "trauma_negative_control_MR.csv")

if (file.exists(trauma_csv)) {
  tr <- fread(trauma_csv)
  required_cols <- c("candidate", "ptsd_beta", "wald_beta", "wald_se", "verdict")
  missing_cols <- setdiff(required_cols, names(tr))
  if (length(missing_cols) > 0) {
    stop("Input file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  scatter <- data.table(
    protein = tr$candidate,
    ptsd_beta = tr$ptsd_beta,
    ptsd_se = NA_real_,
    trauma_beta = tr$wald_beta,
    trauma_se = tr$wald_se,
    verdict = tr$verdict
  )

} else {
  warning("trauma_negative_control_MR.csv not found; using embedded values from the final analysis log.")

  scatter <- data.table(
    protein = c(
      "AKT3", "CD101", "CD40", "CGREF1", "FES",
      "FURIN", "KHK", "SIRPA", "SNX18", "UBE2L6"
    ),
    ptsd_beta = c(
      -0.283, 0.0228, -0.0478, -0.0303, 0.1827,
      -0.0911, 0.0376, -0.0133, -0.2033, 0.1036
    ),  # PTSD freeze-3 case-control betas (pre-specified primary outcome)
    ptsd_se = rep(NA_real_, 10),
    trauma_beta = c(
      -0.238, 0.000562, -0.0261, 0.00175, 0.180,
      -0.126, -0.00302, -0.00385, -0.0525, -0.00156
    ),
    trauma_se = c(
      0.0725, 0.00511, 0.00931, 0.00624, 0.0279,
      0.0195, 0.00819, 0.00286, 0.0452, 0.0153
    ),
    verdict = c(
      "behavioral_confound", "post_pathology", "behavioral_confound",
      "post_pathology", "behavioral_confound", "behavioral_confound",
      "post_pathology", "post_pathology", "post_pathology",
      "post_pathology"
    )
  )
}

scatter[verdict == "behavioral_confound", verdict_plot := "Childhood-maltreatment P < 0.05"]
scatter[verdict == "post_pathology", verdict_plot := "Childhood-maltreatment P > 0.05"]
scatter[verdict == "opposite_direction", verdict_plot := "Opposite direction"]
scatter[verdict == "ptsd_not_significant", verdict_plot := "PTSD not significant"]
scatter[is.na(verdict_plot), verdict_plot := verdict]

scatter[, is_highlight := protein %in% c("KHK", "CD40", "FES", "FURIN", "AKT3")]
scatter[, label_lab := ifelse(is_highlight, protein, "")]
scatter[, trauma_ci_low := trauma_beta - 1.96 * trauma_se]
scatter[, trauma_ci_high := trauma_beta + 1.96 * trauma_se]

verdict_colors <- c(
  "Childhood-maltreatment P > 0.05" = "#9E9E9E",
  "Childhood-maltreatment P < 0.05" = "#6a51a3",
  "Opposite direction" = "#F57C00",
  "PTSD not significant" = "#C7C7C7"
)

p_B_base <- ggplot(scatter, aes(x = ptsd_beta, y = trauma_beta)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.5
  ) +
  geom_errorbar(
    aes(
      ymin = trauma_ci_low,
      ymax = trauma_ci_high,
      color = verdict_plot
    ),
    width = 0,
    alpha = 0.45,
    linewidth = 0.7
  ) +
  geom_point(
    aes(color = verdict_plot, size = is_highlight),
    alpha = 0.95
  ) +
  scale_color_manual(
    values = verdict_colors,
    name = NULL,
    drop = FALSE
  ) +
  scale_size_manual(
    values = c("FALSE" = 3, "TRUE" = 5),
    guide = "none"
  ) +
  labs(
    x = "Effect on PTSD (freeze-3 case-control beta)",
    y = "Effect on childhood maltreatment (beta +/- 95% CI)",
    subtitle = "B. PTSD (freeze-3) versus childhood maltreatment effects"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.6, "cm"),
    plot.subtitle = element_text(face = "bold", size = 12),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11, margin = margin(r = 4)),
    axis.text = element_text(size = 10),
    plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
  )

label_data <- scatter[is_highlight == TRUE]

if (has_ggrepel) {
  p_B <- p_B_base +
    geom_text_repel(
      data = label_data,
      aes(label = label_lab),
      size = 4.2,
      fontface = "bold",
      color = "black",
      segment.size = 0.3,
      segment.color = "gray50",
      min.segment.length = 0.1,
      box.padding = 0.55,
      point.padding = 0.35,
      max.overlaps = Inf
    )
} else {
  p_B <- p_B_base +
    geom_text(
      data = label_data,
      aes(label = label_lab),
      size = 4.2,
      fontface = "bold",
      nudge_y = 0.015,
      color = "black"
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Combine and export
# ─────────────────────────────────────────────────────────────────────────────

fig5 <- p_A / p_B +
  plot_layout(heights = c(1.25, 1.25))

out_tiff <- file.path(fig_dir, "Figure_5_KHK_cross_tissue.tiff")
out_pdf  <- file.path(fig_dir, "Figure_5_KHK_cross_tissue.pdf")
out_png  <- file.path(fig_dir, "Figure_5_KHK_cross_tissue.png")

ggsave(
  filename = out_tiff,
  plot = fig5,
  width = 9.5,
  height = 11,
  dpi = 600,
  compression = "lzw",
  device = "tiff"
)

ggsave(
  filename = out_pdf,
  plot = fig5,
  width = 9.5,
  height = 11
)

ggsave(
  filename = out_png,
  plot = fig5,
  width = 9.5,
  height = 11,
  dpi = 250
)

message("Figure 5 generated:")
message("  TIFF: ", out_tiff)
message("  PDF:  ", out_pdf)
message("  PNG:  ", out_png)

###############################################################################
# Suggested legend:
#
# Figure 5. Cross-tissue and secondary-phenotype sensitivity analyses at the
# KHK/CGREF1 and FES/FURIN loci.
# (A) Forest plot of KHK Mendelian randomization estimates across plasma, brain,
# and a childhood-maltreatment secondary phenotype. Plasma cis-pQTL MR supported
# associations with PTSD symptom severity and PTSD case-control status. Independent
# brain pQTL-based MR showed a concordant association with PTSD, while the same
# plasma cis-pQTL showed no evidence of association with childhood maltreatment.
# The KHK instrument rs2304681 is a common missense variant (p.Val49Ile; minor
# allele frequency approximately 37%) with high posterior support in the brain
# pQTL atlas. Error bars indicate 95% confidence intervals.
# (B) Comparison of MR effect estimates for PTSD (freeze-3 case-control) and
# childhood maltreatment across the 10 PTSD-prioritized proteins. Points are
# coloured by whether the childhood-maltreatment association reached nominal
# P < 0.05; childhood maltreatment is treated as a secondary phenotype rather
# than a strict negative control, and non-significance is not interpreted as
# evidence of absence. KHK showed a PTSD association with a near-null childhood
# maltreatment estimate, whereas FES and FURIN showed same-direction effects on
# both outcomes, consistent with a polypleiotropic interpretation of the
# chromosome 15 FES/FURIN locus. Dashed lines indicate zero effect. CM, childhood
# maltreatment; PTSD, post-traumatic stress disorder; MR, Mendelian
# randomization; pQTL, protein quantitative trait locus; PCL, PTSD Checklist.
###############################################################################
