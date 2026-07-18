###############################################################################
# 06_biological_annotation.R
#
# Biological annotation of the prioritized candidate proteins.
#
# Modules:
#   Part A: Tissue / cell-type expression annotation (literature-based)
#   Part B: Biological pathway grouping
#   Part C: PTSD-relevance annotation
#   Part D: Composite evidence-strength score
#
# No external API dependencies; runs standalone.
# Outputs: protein_expression_annotation.csv, evidence_summary_table.csv
###############################################################################

library(dplyr)
library(readr)
source("00_config.R")


candidates <- c("AKT3", "CD40", "CGREF1", "FES", "FURIN",
                "SIRPA", "CD101", "KHK", "SNX18", "UBE2L6")


###############################################################################
# Part A: Tissue / cell-type expression
###############################################################################

cat("Part A: Tissue / cell-type expression annotation\n")
cat("════════════════════════════════════════\n\n")

expression_annotation <- tribble(
  ~protein, ~brain_expression, ~immune_expression, ~key_cell_types, ~primary_system,
  
  "AKT3",
  "高 (皮层、海马、杏仁核)",
  "低",
  "神经元 (兴奋性+抑制性), 星形胶质细胞",
  "CNS",
  
  "CD40",
  "低-中 (小胶质细胞)",
  "高 (B细胞、DC、单核细胞)",
  "B细胞, 树突状细胞, 小胶质细胞, 活化T细胞",
  "Immune",
  
  "CGREF1",
  "中 (广泛表达)",
  "低",
  "神经元, 内皮细胞",
  "CNS/Other",
  
  "FES",
  "低-中",
  "中 (髓系细胞)",
  "巨噬细胞, 肥大细胞, 中性粒细胞",
  "Immune",
  
  "FURIN",
  "中-高 (广泛表达)",
  "中 (广泛表达)",
  "神经元, 小胶质细胞, T细胞, 巨噬细胞",
  "CNS+Immune",
  
  "SIRPA",
  "高 (小胶质细胞特异)",
  "中 (巨噬细胞)",
  "小胶质细胞 (高), 巨噬细胞, 树突状细胞",
  "CNS (microglia)",
  
  "CD101",
  "极低",
  "中 (T细胞亚群)",
  "调节性T细胞 (Treg), 粒细胞",
  "Immune",
  
  "KHK",
  "低",
  "极低",
  "肝细胞 (极高), 肾小管",
  "Metabolic",
  
  "SNX18",
  "中 (广泛表达)",
  "低",
  "神经元, 广泛表达",
  "CNS/Other",
  
  "UBE2L6",
  "低-中",
  "中-高 (IFN诱导)",
  "单核细胞, 巨噬细胞, NK细胞",
  "Immune (IFN)"
)


###############################################################################
# Part B: 生物学通路
###############################################################################

cat("Part B: Biological pathway grouping\n")
cat("════════════════════════════════════════\n\n")

pathway_annotation <- tribble(
  ~protein, ~pathway_category, ~specific_pathway, ~ptsd_mechanism,
  
  "AKT3",
  "Neuroplasticity",
  "PI3K-AKT-mTOR",
  "Fear extinction和synaptic plasticity; AKT3 KO小鼠显示焦虑样行为; mTOR参与记忆巩固",
  
  "CD40",
  "Neuroinflammation",
  "CD40-CD40L → NF-κB",
  "调节神经炎症和小胶质细胞激活; PTSD患者外周CD40L水平升高; CD40 agonists可增强免疫应答",
  
  "CGREF1",
  "Cell growth",
  "Calcium signaling",
  "钙结合蛋白, 参与细胞生长调控; PTSD相关功能研究有限, 属novel finding",
  
  "FES",
  "Innate immunity",
  "TLR/cytokine signaling",
  "非受体酪氨酸激酶, 参与先天免疫信号和巨噬细胞激活",
  
  "FURIN",
  "Neuropeptide processing",
  "Proprotein convertase",
  "切割proNGF, proBDNF, TGF-β等底物; 参与HPA轴调控; 与SARS-CoV-2 spike蛋白切割相关",
  
  "SIRPA",
  "Neuroimmune interface",
  "CD47-SIRPα checkpoint",
  "调控小胶质细胞吞噬和突触修剪; 'don't eat me'信号; PTSD神经炎症假说的关键分子",
  
  "CD101",
  "Immune regulation",
  "Treg function",
  "免疫检查点分子; CD101+ Treg参与免疫抑制; PTSD免疫失调可能涉及Treg功能障碍",
  
  "KHK",
  "Metabolic",
  "Fructose metabolism",
  "果糖代谢关键酶; PTSD患者代谢综合征共病率高; 潜在的代谢共病桥梁",
  
  "SNX18",
  "Membrane trafficking",
  "Endocytosis/Autophagy",
  "参与受体内吞和膜运输; 调控自噬体形成; 神经元中参与突触受体转运",
  
  "UBE2L6",
  "IFN signaling",
  "ISG15 conjugation",
  "干扰素信号通路核心组件; ISG15 E2连接酶; PTSD与慢性IFN通路激活相关"
)


###############################################################################
# Part C: 综合证据表 (整合所有分析结果)
###############################################################################

cat("Part C: Composite evidence table\n")
cat("════════════════════════════════════════\n\n")

# 读取各分析结果
mr_primary   <- readRDS(file.path(result_dir, "mr_all_outcomes.rds"))
steiger      <- read_csv(file.path(result_dir, "steiger_directionality_results.csv"), show_col_types = FALSE)
presso       <- read_csv(file.path(result_dir, "mr_presso_results.csv"), show_col_types = FALSE)
reverse      <- read_csv(file.path(result_dir, "reverse_mr_results.csv"), show_col_types = FALSE)
mvmr         <- read_csv(file.path(result_dir, "mvmr_results.csv"), show_col_types = FALSE)
smr_heidi    <- read_csv(file.path(result_dir, "smr_heidi_results.csv"), show_col_types = FALSE)
multi_iv     <- read_csv(file.path(result_dir, "multi_instrument_mr.csv"), show_col_types = FALSE)
classified   <- read_csv(file.path(result_dir, "ptsd_candidates_classified.csv"), show_col_types = FALSE)

# 构建综合表
evidence_table <- tibble(protein = candidates) %>%
  # 基本MR结果
  left_join(
    mr_primary %>%
      filter(grepl("PTSD", outcome), protein %in% candidates) %>%
      select(protein, OR, mr_pval, fdr),
    by = "protein"
  ) %>%
  # 分类
  left_join(
    classified %>% select(protein, category = any_of(c("category", "classification"))),
    by = "protein"
  ) %>%
  # Steiger
  left_join(
    steiger %>% select(protein, steiger_correct = correct_direction, steiger_p = steiger_p),
    by = "protein"
  ) %>%
  # MR-PRESSO
  left_join(
    presso %>% select(protein, presso_global_p = global_test_p, presso_outliers = n_outliers),
    by = "protein"
  ) %>%
  # Reverse MR (IVW + WM)
  left_join(
    reverse %>% filter(method == "IVW") %>%
      select(protein, reverse_ivw_p = pval),
    by = "protein"
  ) %>%
  left_join(
    reverse %>% filter(method == "Weighted Median") %>%
      select(protein, reverse_wm_p = pval),
    by = "protein"
  ) %>%
  # Multi-instrument IVW
  left_join(
    multi_iv %>% filter(method == "IVW") %>%
      select(protein, multi_iv_p = pval, cochran_q_p),
    by = "protein"
  ) %>%
  # SMR + HEIDI
  left_join(
    smr_heidi %>% select(protein, smr_p = smr_pval, heidi_p = heidi_pval, smr_overall = overall),
    by = "protein"
  ) %>%
  # 通路注释
  left_join(
    pathway_annotation %>% select(protein, pathway_category, specific_pathway),
    by = "protein"
  ) %>%
  # 表达注释
  left_join(
    expression_annotation %>% select(protein, primary_system, key_cell_types),
    by = "protein"
  )

# 计算综合证据评分 (通过验证的分析数量)
evidence_table <- evidence_table %>%
  mutate(
    pass_steiger   = steiger_correct & steiger_p < 0.05,
    pass_presso    = presso_global_p >= 0.05,
    pass_reverse   = reverse_wm_p >= 0.05,    # 用WM判断
    pass_multi_iv  = multi_iv_p < 0.05,
    pass_smr       = grepl("强证据", smr_overall),
    pass_no_het    = cochran_q_p >= 0.05,
    
    n_passed = rowSums(across(starts_with("pass_"), ~replace_na(., FALSE))),
    
    evidence_tier = case_when(
      n_passed >= 5 ~ "Tier 1 (强)",
      n_passed >= 3 ~ "Tier 2 (中)",
      TRUE          ~ "Tier 3 (弱)"
    )
  )


###############################################################################
# Part D: 输出
###############################################################################

cat("Part D: Save results\n")
cat("════════════════════════════════════════\n\n")

# 保存
write_csv(expression_annotation, file.path(result_dir, "protein_expression_annotation.csv"))
write_csv(pathway_annotation, file.path(result_dir, "pathway_annotation.csv"))
write_csv(evidence_table, file.path(result_dir, "evidence_summary_table.csv"))

# 打印综合证据表
cat("═══════════════════════════════════════════════════════\n")
cat("Composite evidence score\n")
cat("═══════════════════════════════════════════════════════\n\n")

evidence_table %>%
  select(protein, OR, evidence_tier, n_passed,
         pass_steiger, pass_presso, pass_reverse,
         pass_multi_iv, pass_smr, pass_no_het) %>%
  mutate(OR = round(OR, 3)) %>%
  arrange(desc(n_passed)) %>%
  as.data.frame() %>%
  print(right = FALSE)

cat("\n")
cat("═══════════════════════════════════════════════════════\n")
cat("Biological narrative\n")
cat("═══════════════════════════════════════════════════════\n\n")
cat("1. Neuroinflammation / microglia axis: SIRPA + CD40 + FURIN\n")
cat("   -> microglial activation and aberrant synaptic pruning\n")
cat("2. Synaptic plasticity: AKT3 (PI3K-mTOR) + SNX18 (receptor trafficking)\n")
cat("   -> impaired fear extinction and memory consolidation\n")
cat("3. Interferon / immune dysregulation: UBE2L6 + FES + CD101\n")
cat("   -> chronic immune activation and Treg dysfunction\n")
cat("4. Metabolic comorbidity: KHK (fructose metabolism)\n")
cat("   -> bridge to PTSD-metabolic syndrome comorbidity\n")

cat(sprintf("\nAll outputs saved to: %s\n", result_dir))
cat("  protein_expression_annotation.csv\n")
cat("  pathway_annotation.csv\n")
cat("  evidence_summary_table.csv\n")
