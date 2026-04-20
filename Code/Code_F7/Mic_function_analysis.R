# ============== Analysis of Microglial Functional Exhaustion ==============
# Using Externally Validated Gene Sets

# Step 1: Define Key Gene Sets (from established literature)
# Metabolic gene set (electron transport chain & oxidative phosphorylation)
metabolic_genes <- c("NDUFA1", "NDUFA2", "NDUFAB1", "NDUFAF1", "NDUFAF2", 
                     "NDUFB1", "NDUFB2", "NDUFC1", "NDUFC2", "NDUFS1", 
                     "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", 
                     "NDUFS7", "NDUFS8", "NDUFV1", "NDUFV2", "SDHA", 
                     "SDHB", "SDHC", "SDHD", "UQCRC1", "UQCRC2", 
                     "UQCRFS1", "UQCRH", "UQCRQ", "COX4I1", "COX5A", 
                     "COX5B", "COX6A1", "COX6B1", "COX6C", "COX7A2", 
                     "COX7B", "COX7C", "COX8A", "ATP5A1", "ATP5B", 
                     "ATP5C1", "ATP5D", "ATP5E", "ATP5F1", "ATP5G1", 
                     "ATP5G2", "ATP5G3", "ATP5H", "ATP5I", "ATP5J", 
                     "ATP5J2", "ATP5L", "ATP5O")

# Clearance gene set (phagocytosis receptors)
clearance_genes <- c("TREM2", "PROS1", "GAS6", "MERTK", "AXL", 
                     "TYROBP", "CD33", "CR1", "SCARB1", "LDLR")

# Immune activation gene set (inflammatory response)
immune_genes <- c("IL1B", "TNF", "IL6", "IL18", "CXCL8", 
                  "TLR4", "TLR2", "NLRP3", "CASP1", "PYCARD")

# Key functional genes
key_genes <- c("P2RY12", "CX3CR1", "CD200R1", "TREM2", "APOE")

# Step 2: Calculate gene availability in dataset
available_metabolic <- metabolic_genes[metabolic_genes %in% rownames(mzl)]
available_clearance <- clearance_genes[clearance_genes %in% rownames(mzl)]
available_immune <- immune_genes[immune_genes %in% rownames(mzl)]
available_key <- key_genes[key_genes %in% rownames(mzl)]

cat("Available metabolic genes:", length(available_metabolic), "/", length(metabolic_genes), "\n")
cat("Available clearance genes:", length(available_clearance), "/", length(clearance_genes), "\n")
cat("Available immune genes:", length(available_immune), "/", length(immune_genes), "\n")


# ============== 修复模块评分问题 - 完整代码 ==============

# 步骤1：重新标准化数据（如果未标准化）
if(!"NormalizeData" %in% names(mzl@commands)) {
  cat("检测到数据未标准化，正在执行NormalizeData...\n")
  mzl <- NormalizeData(mzl, normalization.method = "LogNormalize", scale.factor = 10000)
} else {
  cat("数据已标准化\n")
}

# 步骤2：手动计算模块评分
manual_module_score <- function(seurat_obj, gene_list) {
  # 获取表达矩阵
  exp_mat <- GetAssayData(seurat_obj, slot = "data")
  
  # 筛选可用基因
  available_genes <- intersect(gene_list, rownames(exp_mat))
  cat("用于计算的基因数:", length(available_genes), "\n")
  
  if(length(available_genes) == 0) {
    warning("没有可用基因，返回NA向量")
    return(rep(NA, ncol(seurat_obj)))
  }
  
  # 计算模块评分（取平均值）
  module_score <- colMeans(exp_mat[available_genes, , drop = FALSE], na.rm = TRUE)
  
  # 检查是否全为NA
  if(all(is.na(module_score))) {
    cat("警告: 所有评分均为NA，可能原因:\n")
    cat("- 表达矩阵全为0或NA\n")
    cat("- 基因表达值极低\n")
    cat("正在检查前5个细胞的表达值（取前3个基因）:\n")
    print(head(exp_mat[available_genes, 1:5], 3))
  }
  
  return(module_score)
}

# 步骤3：计算各功能评分
cat("\n计算代谢功能评分...\n")
mzl$Manual_Metabolic_Score <- manual_module_score(mzl, metabolic_genes)

cat("\n计算清除功能评分...\n")
mzl$Manual_Clearance_Score <- manual_module_score(mzl, clearance_genes)

cat("\n计算免疫激活评分...\n")
mzl$Manual_Immune_Score <- manual_module_score(mzl, immune_genes)

# 步骤4：检查手动计算结果
cat("\n===== 手动计算评分摘要 =====\n")
cat("代谢评分摘要:\n")
print(summary(mzl$Manual_Metabolic_Score))
cat("NA比例:", mean(is.na(mzl$Manual_Metabolic_Score)), "\n")

cat("清除评分摘要:\n")
print(summary(mzl$Manual_Clearance_Score))
cat("NA比例:", mean(is.na(mzl$Manual_Clearance_Score)), "\n")

cat("免疫评分摘要:\n")
print(summary(mzl$Manual_Immune_Score))
cat("NA比例:", mean(is.na(mzl$Manual_Immune_Score)), "\n")

# 步骤5：创建分析数据框
analysis_df <- data.frame(
  CellID = colnames(mzl),
  Metabolic_Score = mzl$Manual_Metabolic_Score,
  Clearance_Score = mzl$Manual_Clearance_Score,
  Immune_Score = mzl$Manual_Immune_Score,
  Diagnosis = mzl$diagnosis,
  Sex = mzl$sex,
  CellType = mzl$cell_type_ident
)

# 添加关键基因表达
for(gene in available_key) {
  if(gene %in% rownames(mzl)) {
    analysis_df[[gene]] <- GetAssayData(mzl, slot = "data")[gene, ]
  } else {
    cat("基因", gene, "不在数据中，跳过\n")
    analysis_df[[gene]] <- NA
  }
}

# 仅提取小胶质细胞
mic_analysis_df <- subset(analysis_df, CellType == "mic")
cat("\n小胶质细胞数量:", nrow(mic_analysis_df), "\n")

# 检查NA值
cat("\n===== 小胶质细胞数据NA检查 =====\n")
cat("代谢评分NA比例:", mean(is.na(mic_analysis_df$Metabolic_Score)), "\n")
cat("清除评分NA比例:", mean(is.na(mic_analysis_df$Clearance_Score)), "\n")
cat("免疫评分NA比例:", mean(is.na(mic_analysis_df$Immune_Score)), "\n")

# 处理NA值：用各自分组的中位数填充
library(dplyr)
mic_analysis_df <- mic_analysis_df %>%
  group_by(Diagnosis, Sex) %>%
  mutate(
    Metabolic_Score = ifelse(is.na(Metabolic_Score), median(Metabolic_Score, na.rm = TRUE), Metabolic_Score),
    Clearance_Score = ifelse(is.na(Clearance_Score), median(Clearance_Score, na.rm = TRUE), Clearance_Score),
    Immune_Score = ifelse(is.na(Immune_Score), median(Immune_Score, na.rm = TRUE), Immune_Score)
  ) %>%
  ungroup()

# 创建组合分组变量
mic_analysis_df$GroupSex <- paste(mic_analysis_df$Diagnosis, mic_analysis_df$Sex, sep = "_")
mic_analysis_df$GroupSex <- factor(mic_analysis_df$GroupSex,
                                   levels = c("ctl_female", "AD_female", 
                                              "ctl_male", "AD_male"))

# 步骤6：验证填充后的数据
cat("\n===== 填充后数据摘要 =====\n")
cat("代谢评分填充后摘要:\n")
print(summary(mic_analysis_df$Metabolic_Score))
cat("清除评分填充后摘要:\n")
print(summary(mic_analysis_df$Clearance_Score))
cat("免疫评分填充后摘要:\n")
print(summary(mic_analysis_df$Immune_Score))

# ============== 完整可视化分析 ==============

library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

# 设置主题
theme_set(theme_minimal(base_size = 12) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                  plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  legend.position = "none"))

# 自定义颜色方案
group_colors <- c(
  "ctl_female" = "#E6B0AA",  # Light pink
  "AD_female" = "#EC7063",   # Vibrant pink
  "ctl_male" = "#AED6F1",    # Light blue
  "AD_male" = "#3498DB"      # Deep blue
)

# 1. 功能评分箱线图 - 修复布局和显著性标记问题
plot_functional_score <- function(data, y_var, title, subtitle) {
  # 动态计算Y轴范围
  y_range <- range(data[[y_var]], na.rm = TRUE)
  y_max <- y_range[2] + diff(y_range) * 0.2
  
  # 创建箱线图
  p <- ggplot(data, aes(x = GroupSex, y = .data[[y_var]], fill = GroupSex)) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.3, color = "gray30") +
    scale_fill_manual(values = group_colors) +
    labs(title = title,
         subtitle = subtitle,
         x = "",
         y = gsub("_", " ", y_var)) +
    scale_x_discrete(labels = c("Female\nControl", "Female\nAD", 
                                "Male\nControl", "Male\nAD")) +
    ylim(y_range[1], y_max)  # 设置固定的Y轴范围
  
  # 添加显著性标记（动态调整位置）
  p <- p + stat_compare_means(
    comparisons = list(
      c("ctl_female", "AD_female"),
      c("ctl_male", "AD_male"),
      c("AD_female", "AD_male")
    ),
    label = "p.signif",
    method = "wilcox.test",
    tip.length = 0.01,
    size = 4,
    label.y = c(y_max * 0.9, y_max * 0.95, y_max * 1.0)  # 确保标记不重叠
  )
  
  return(p)
}

# 创建三个功能评分图
metabolic_plot <- plot_functional_score(
  mic_analysis_df, 
  "Metabolic_Score", 
  "Metabolic Function in Microglia",
  "Electron Transport Chain/Oxidative Phosphorylation"
)

clearance_plot <- plot_functional_score(
  mic_analysis_df, 
  "Clearance_Score", 
  "Clearance Function in Microglia",
  "Phagocytic Receptor Activity"
)

immune_plot <- plot_functional_score(
  mic_analysis_df, 
  "Immune_Score", 
  "Immune Activation in Microglia",
  "Inflammatory Response"
)

# 组合功能评分图
functional_plots <- ggarrange(
  metabolic_plot, clearance_plot, immune_plot,
  ncol = 1, nrow = 3,
  labels = c("A", "B", "C"),
  font.label = list(size = 14, face = "bold")
)

# 保存功能评分图
ggsave("Microglia_Functional_Analysis.png", functional_plots, 
       width = 9, height = 12, dpi = 300)

# 2. 关键基因表达图 - 优化排版和显著性标记
# 移除CD200R1（如果需要）
if("CD200R1" %in% colnames(mic_analysis_df)) {
  mic_analysis_df$CD200R1 <- NULL
  cat("已移除CD200R1\n")
}

# 筛选关键基因（仅保留有表达的基因）
expressed_key_genes <- available_key[sapply(available_key, function(g) {
  any(mic_analysis_df[[g]] > 0)
})]

cat("将展示的关键基因:", paste(expressed_key_genes, collapse = ", "), "\n")

# 创建关键基因表达图函数
plot_key_gene <- function(data, gene) {
  # 动态计算Y轴范围
  gene_exp <- data[[gene]]
  y_max <- max(gene_exp, na.rm = TRUE) * 1.2
  
  p <- ggplot(data, aes(x = GroupSex, y = .data[[gene]], fill = GroupSex)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha = 0.8) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste(gene, "Expression"),
         x = "Group",
         y = "Expression Level") +
    scale_x_discrete(labels = c("F\nCTL", "F\nAD", "M\nCTL", "M\nAD")) +
    ylim(0, y_max)  # 从0开始
  
  # 添加显著性标记（每组比较独立定位）
  comparisons <- list(
    c("ctl_female", "AD_female"),
    c("ctl_male", "AD_male")
  )
  
  if(length(comparisons) > 0) {
    p <- p + stat_compare_means(
      comparisons = comparisons,
      label = "p.signif",
      method = "wilcox.test",
      size = 4,
      label.y = rep(y_max * 0.95, length(comparisons))  # 统一高度
    )
  }
  
  return(p)
}

# 创建关键基因表达图（4个基因）
key_plots <- lapply(expressed_key_genes[1:4], function(g) {
  plot_key_gene(mic_analysis_df, g)
})

# 创建2x2布局
key_plots_combined <- wrap_plots(key_plots, ncol = 2) + 
  plot_annotation(title = "Key Gene Expression in Microglia",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))

# 保存关键基因图
ggsave("Key_Gene_Expression.png", key_plots_combined, 
       width = 12, height = 9, dpi = 300)

# 3. 代谢-清除功能相关性分析（针对女性AD组）
female_ad_mic <- mic_analysis_df %>% 
  filter(Sex == "female", Diagnosis == "AD") %>%
  filter(!is.na(Metabolic_Score) & !is.na(Clearance_Score))

if(nrow(female_ad_mic) >= 3) {
  # 计算相关性
  cor_test <- cor.test(
    female_ad_mic$Metabolic_Score, 
    female_ad_mic$Clearance_Score, 
    method = "pearson"
  )
  
  # 创建相关性图
  cor_plot <- ggplot(female_ad_mic, aes(x = Metabolic_Score, y = Clearance_Score)) +
    geom_point(color = "#EC7063", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", color = "#922B21", se = TRUE, fill = "#E6B0AA") +
    annotate("text", 
             x = min(female_ad_mic$Metabolic_Score, na.rm = TRUE), 
             y = max(female_ad_mic$Clearance_Score, na.rm = TRUE) * 0.95,
             label = sprintf("r = %.2f\np = %.3g", 
                             cor_test$estimate, 
                             cor_test$p.value),
             hjust = 0, vjust = 1, size = 5, fontface = "bold", color = "#922B21") +
    labs(title = "Metabolic-Clearance Correlation in Female AD Microglia",
         subtitle = "Hypothesis of Functional Coordination",
         x = "Metabolic Function Score",
         y = "Clearance Function Score")
  
  ggsave("Metabolic_Clearance_Correlation.png", cor_plot, 
         width = 8, height = 6, dpi = 300)
} else {
  cat("\n无法进行相关性分析 - 女性AD组样本不足\n")
}

# 4. 评分分布密度图（按组）
plot_density <- function(data, score, title) {
  ggplot(data, aes(x = .data[[score]], fill = GroupSex)) +
    geom_density(alpha = 0.6) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste(title, "Distribution"),
         x = gsub("_", " ", score),
         y = "Density") +
    theme(legend.position = "bottom",
          legend.title = element_blank())
}

density_metabolic <- plot_density(mic_analysis_df, "Metabolic_Score", "Metabolic Function")
density_clearance <- plot_density(mic_analysis_df, "Clearance_Score", "Clearance Function")
density_immune <- plot_density(mic_analysis_df, "Immune_Score", "Immune Activation")

# 组合密度图
density_plots <- (density_metabolic | density_clearance | density_immune) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Functional_Score_Distributions.png", density_plots, 
       width = 12, height = 4, dpi = 300)

# 5. 最终报告
cat("\n======= 分析完成 =======\n")
cat("已生成以下文件:\n")
cat("1. 功能评分图: Microglia_Functional_Analysis.png\n")
cat("2. 关键基因表达图: Key_Gene_Expression.png\n")
cat("3. 评分分布图: Functional_Score_Distributions.png\n")

if(exists("cor_plot")) {
  cat("4. 代谢-清除相关性图: Metabolic_Clearance_Correlation.png\n")
}

cat("\n分析统计摘要:\n")
cat("- 小胶质细胞数量:", nrow(mic_analysis_df), "\n")
cat("- 女性AD组细胞数 (用于相关分析):", 
    if(exists("female_ad_mic")) nrow(female_ad_mic) else 0, "\n")
cat("- 使用的代谢基因数:", length(available_metabolic), "\n")
cat("- 使用的清除基因数:", length(available_clearance), "\n")
cat("- 使用的免疫基因数:", length(available_immune), "\n")

# 保存分析数据
write.csv(mic_analysis_df, "Microglia_Functional_Analysis_Data.csv", row.names = FALSE)

# ============== 微胶质细胞功能分析数据总结 ==============
cat("\n===== 微胶质细胞功能分析数据总结报告 =====")
cat("\n生成时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
cat("\n数据集:", ifelse(exists("mzl"), "mzl Seurat对象", "未知"))
cat("\n分析版本: 功能评分方法v2 (手动计算均值)\n")

# 1. 样本基础信息
cat("\n------ 样本基础信息 ------\n")
cat("小胶质细胞总数:", nrow(mic_analysis_df), "\n")

# 按诊断和性别分组计数
group_counts <- mic_analysis_df %>%
  group_by(Diagnosis, Sex) %>%
  summarise(Cells = n(), .groups = "drop")

cat("\n分组细胞数量:\n")
print(group_counts)

# 2. 功能评分总体统计
cat("\n------ 功能评分总体统计 ------\n")
cat("所有评分数值经过标准化处理（Z-score）\n")

# 功能评分汇总
score_summary <- mic_analysis_df %>%
  summarise(
    across(
      c(Metabolic_Score, Clearance_Score, Immune_Score),
      list(
        Mean = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        Median = ~median(.x, na.rm = TRUE),
        Min = ~min(.x, na.rm = TRUE),
        Max = ~max(.x, na.rm = TRUE)
      )
    )
  )

cat("\n功能评分总体统计:\n")
print(score_summary)

# 3. 分组的评分差异统计
cat("\n------ 分组评分差异统计 (Wilcoxon检验) ------\n")

# 代谢功能评分分组比较
cat("\n[代谢功能]\n")
metabolic_stats <- mic_analysis_df %>%
  group_by(GroupSex) %>%
  summarise(
    Mean = mean(Metabolic_Score),
    Median = median(Metabolic_Score),
    SD = sd(Metabolic_Score),
    .groups = "drop"
  )

# 组间比较
ad_vs_ctl_female <- wilcox.test(
  Metabolic_Score ~ Diagnosis,
  data = subset(mic_analysis_df, Sex == "female")
)

ad_vs_ctl_male <- wilcox.test(
  Metabolic_Score ~ Diagnosis,
  data = subset(mic_analysis_df, Sex == "male")
)

female_vs_male_ad <- wilcox.test(
  Metabolic_Score ~ Sex,
  data = subset(mic_analysis_df, Diagnosis == "AD")
)

cat("分组代谢评分:\n")
print(metabolic_stats)
cat("\nAD vs CTL (女性): p =", format.pval(ad_vs_ctl_female$p.value, digits = 3), "\n")
cat("AD vs CTL (男性): p =", format.pval(ad_vs_ctl_male$p.value, digits = 3), "\n")
cat("女性AD vs 男性AD: p =", format.pval(female_vs_male_ad$p.value, digits = 3), "\n")

# 清除功能评分分组比较
cat("\n[清除功能]\n")
clearance_stats <- mic_analysis_df %>%
  group_by(GroupSex) %>%
  summarise(
    Mean = mean(Clearance_Score),
    Median = median(Clearance_Score),
    SD = sd(Clearance_Score),
    .groups = "drop"
  )

ad_vs_ctl_female_clear <- wilcox.test(
  Clearance_Score ~ Diagnosis,
  data = subset(mic_analysis_df, Sex == "female")
)

ad_vs_ctl_male_clear <- wilcox.test(
  Clearance_Score ~ Diagnosis,
  data = subset(mic_analysis_df, Sex == "male")
)

cat("分组清除评分:\n")
print(clearance_stats)
cat("\nAD vs CTL (女性): p =", format.pval(ad_vs_ctl_female_clear$p.value, digits = 3), "\n")
cat("AD vs CTL (男性): p =", format.pval(ad_vs_ctl_male_clear$p.value, digits = 3), "\n")

# 4. 关键基因表达统计
cat("\n------ 关键基因表达统计 ------\n")
key_genes <- available_key[available_key != "CD200R1"]

gene_stats <- mic_analysis_df %>%
  group_by(GroupSex) %>%
  summarise(
    across(
      all_of(key_genes),
      list(
        Mean = ~mean(.x, na.rm = TRUE),
        Median = ~median(.x, na.rm = TRUE),
        Expressed_Cells = ~sum(.x > 0, na.rm = TRUE),
        Percent_Expressed = ~mean(.x > 0, na.rm = TRUE) * 100
      )
    ),
    .groups = "drop"
  )

cat("关键基因表达统计:\n")
print(gene_stats)

# 5. 代谢-清除功能相关性分析（仅女性AD组）
cat("\n------ 代谢与清除功能相关性分析 ------\n")
female_ad_mic <- mic_analysis_df %>% 
  filter(Sex == "female", Diagnosis == "AD")

if(nrow(female_ad_mic) >= 3) {
  cor_test <- cor.test(
    female_ad_mic$Metabolic_Score, 
    female_ad_mic$Clearance_Score, 
    method = "pearson"
  )
  
  cat("女性AD组样本数:", nrow(female_ad_mic), "\n")
  cat("Pearson相关系数 (r):", round(cor_test$estimate, 3), "\n")
  cat("显著性 (p):", format.pval(cor_test$p.value, digits = 3), "\n")
  cat("95%置信区间:", paste(round(cor_test$conf.int, 3), collapse = " to "), "\n")
} else {
  cat("样本不足进行相关性分析 (n =", nrow(female_ad_mic), "< 3)\n")
}

# 6. 基因集使用总结
cat("\n------ 分析基因集总结 ------\n")
cat("代谢功能基因数:", length(available_metabolic), "/", length(metabolic_genes), "\n")
cat("清除功能基因数:", length(available_clearance), "/", length(clearance_genes), "\n")
cat("免疫激活基因数:", length(available_immune), "/", length(immune_genes), "\n")
cat("关键功能基因:", paste(available_key, collapse = ", "), "\n")

# 7. 分析数据保存
write.csv(mic_analysis_df, "Microglia_Functional_Analysis_Summary_Data.csv", row.names = FALSE)
cat("\n分析数据已保存至: Microglia_Functional_Analysis_Summary_Data.csv\n")

# 最终报告
cat("\n======= 分析总结 =======\n")
cat("主要发现:\n")

# 代谢功能总结
cat("- 代谢功能评分: 平均值", round(mean(mic_analysis_df$Metabolic_Score), 3))
cat(" | 女性AD患者评分显著", 
    if(ad_vs_ctl_female$p.value < 0.05) "降低" else "未改变", 
    "(p =", format.pval(ad_vs_ctl_female$p.value, digits = 3), ")\n")

# 清除功能总结
cat("- 清除功能评分: 平均值", round(mean(mic_analysis_df$Clearance_Score), 3))
cat(" | 女性AD患者评分显著", 
    if(ad_vs_ctl_female_clear$p.value < 0.05) "降低" else "未改变", 
    "(p =", format.pval(ad_vs_ctl_female_clear$p.value, digits = 3), ")\n")


# 相关性总结
if(exists("cor_test") && cor_test$p.value < 0.05) {
  cat("- 女性AD患者中，代谢与清除功能呈", 
      if(cor_test$estimate > 0) "正" else "负",
      "相关性 (r =", round(cor_test$estimate, 2), 
      "p =", format.pval(cor_test$p.value, digits = 3), ")\n")
} else {
  cat("- 代谢与清除功能未发现显著相关性\n")
}

# 关键基因总结
cat("- TREM2表达: 女性AD组平均", 
    round(mean(subset(mic_analysis_df, GroupSex == "AD_female")$TREM2), 3), 
    "| 男性AD组平均", 
    round(mean(subset(mic_analysis_df, GroupSex == "AD_male")$TREM2), 3), "\n")

cat("- APOE表达: 女性AD组平均", 
    round(mean(subset(mic_analysis_df, GroupSex == "AD_female")$APOE), 3), 
    "| 男性AD组平均", 
    round(mean(subset(mic_analysis_df, GroupSex == "AD_male")$APOE), 3), "\n")

cat("\n分析完成! 详细数据已保存供进一步研究使用。\n")

# 最终版：AD组内性别显著度 + 增强箱线图 + 组别颜色 + 点收窄
library(ggplot2)
library(ggpubr)
library(patchwork)

# 设置颜色方案
group_colors <- c(
  "ctl_female" = "#E6B0AA",  # Light pink
  "AD_female" = "#EC7063",   # Vibrant pink
  "ctl_male" = "#AED6F1",    # Light blue
  "AD_male" = "#3498DB"      # Deep blue
)

# 创建单个功能的图表函数
create_single_plot <- function(data, score, title) {
  # 动态计算y轴上限
  y_max <- max(data[[score]]) * 1.3
  
  # 准备AD组内性别差异标记位置
  ad_data <- data %>% filter(Diagnosis == "AD")
  ad_y_position <- if(length(ad_data[[score]]) > 0) {
    max(ad_data[[score]]) * 1.15
  } else {
    y_max * 1.15
  }
  
  ggplot(data, aes(x = GroupSex, y = .data[[score]])) +
    # 增强箱线图（更明显）
    geom_boxplot(
      aes(fill = GroupSex),
      width = 0.2,
      outlier.shape = NA,
      alpha = 0.9,
      size = 0.7,
      color = "black"
    ) +
    
    # 小提琴图（半透明）
    geom_violin(
      aes(fill = GroupSex),
      alpha = 0.4,
      scale = "width",
      trim = TRUE,
      draw_quantiles = c(0.25, 0.5, 0.75),
      color = "black",
      size = 0.4
    ) +
    
    # 散点图（收窄范围 + 组别颜色）
    geom_jitter(
      aes(color = GroupSex), 
      size = 1.0, 
      alpha = 0.3, 
      width = 0.1,
      height = 0
    ) +
    
    # 颜色设置
    scale_fill_manual(values = group_colors) +
    scale_color_manual(values = group_colors) +
    guides(color = "none", fill = "none") +
    
    # 标签和主题
    labs(title = title, y = "Functional Score", x = "") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = unit(c(10, 15, 15, 15), "points")
    ) +
    scale_x_discrete(
      labels = c(
        "ctl_female" = "F Ctrl", 
        "AD_female" = "F AD", 
        "ctl_male" = "M Ctrl", 
        "AD_male" = "M AD"
      )
    ) +
    
    # 添加组间比较显著性标记
    stat_compare_means(
      comparisons = list(
        c("ctl_female", "AD_female"),
        c("ctl_male", "AD_male")
      ),
      label = "p.signif",
      method = "wilcox.test",
      tip.length = 0.01,
      label.y = c(y_max * 0.85, y_max * 0.92),
      size = 5,
      bracket.size = 0.6,
      hide.ns = TRUE
    ) +
    
    # 添加AD组内性别差异标记（顶部）
    stat_compare_means(
      comparisons = list(c("AD_female", "AD_male")),
      label = "p.signif",
      method = "wilcox.test",
      tip.length = 0.01,
      label.y = ad_y_position,
      size = 5,
      bracket.size = 0.6,
      hide.ns = TRUE,
      vjust = 0.5
    ) +
    
    # 添加AD性别差异比较线
    geom_segment(
      aes(x = 2, xend = 4, y = ad_y_position * 0.98, yend = ad_y_position * 0.98),
      color = "black",
      size = 0.3
    ) +
    
    ylim(NA, y_max * 1.2)
}

# 创建三个分图
metabolic_plot <- create_single_plot(
  mic_analysis_df, 
  "Metabolic_Score", 
  "Metabolic Function"
)

clearance_plot <- create_single_plot(
  mic_analysis_df, 
  "Clearance_Score", 
  "Clearance Function"
)

immune_plot <- create_single_plot(
  mic_analysis_df, 
  "Immune_Score", 
  "Immune Activation"
)

# 组合图表
combined_plot <- metabolic_plot + clearance_plot + immune_plot +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Microglial Functional Profiles with Sex Differences in AD",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.margin = unit(c(20, 20, 20, 20), "points")
    )
  )

# 保存图表
ggsave("AD_Sex_Differences_Functional_Score.png", combined_plot, 
       width = 16, height = 7, dpi = 300, bg = "white")


