# 加载必要的包
library(Seurat)
library(GSVA)
library(msigdbr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)

# 步骤1：从Seurat对象提取数据并确保分组正确 --------------------------------
# 提取元数据
metadata <- mzl@meta.data

# 使用diagnosis列作为分组依据
metadata$Group <- metadata$diagnosis

# 检查分组
cat("诊断分组:", unique(metadata$Group), "\n")
cat("各组样本量:\n")
print(table(metadata$Group))

# 提取表达矩阵（标准化数据）
expr_matrix <- GetAssayData(mzl, slot = "data")

# 确保细胞ID对齐
if(!identical(rownames(metadata), colnames(expr_matrix))) {
  cat("警告：细胞ID不匹配，正在修复...\n")
  rownames(metadata) <- colnames(expr_matrix)
}

# 步骤2：准备基因集 -----------------------------------------------------
# 获取HALLMARK雌激素反应基因集
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
estrogen_sets <- hallmark_sets[
  hallmark_sets$gs_name %in% c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
]

# 转换为GSVA需要的列表格式
gene_sets <- split(estrogen_sets$gene_symbol, estrogen_sets$gs_name)

# 步骤3：执行GSVA分析 ---------------------------------------------------
set.seed(123)
gsva_results <- gsva(
  as.matrix(expr_matrix),
  gene_sets,
  method = "ssgsea",
  kcdf = "Gaussian",
  min.sz = 1,
  parallel.sz = 4,
  verbose = TRUE
)

cat("GSVA结果维度:", dim(gsva_results), "\n")

# 步骤4：整合结果 ------------------------------------------------------
# 转置结果并添加元数据
gsva_df <- as.data.frame(t(gsva_results))
gsva_df$CellID <- rownames(gsva_df)
gsva_df <- merge(gsva_df, metadata, by.x = "CellID", by.y = "row.names")

# 步骤5：数据整理 ------------------------------------------------------
# 转换为长格式用于分析
plot_data <- reshape2::melt(
  gsva_df,
  id.vars = c("CellID", "Group", "cell_type_ident"),
  measure.vars = c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
  variable.name = "Pathway",
  value.name = "EnrichmentScore"
)

# 简化路径名称
plot_data$Pathway <- gsub("HALLMARK_", "", plot_data$Pathway)
plot_data$Pathway <- gsub("_RESPONSE", "", plot_data$Pathway)
plot_data$Pathway <- factor(
  plot_data$Pathway,
  levels = c("ESTROGEN_EARLY", "ESTROGEN_LATE"),
  labels = c("Early Response", "Late Response")
)

# 重命名列
plot_data <- plot_data %>%
  rename(
    CellType = cell_type_ident
  )

# 步骤6：计算效应量（均值差）及其置信区间 --------------------------------
effect_data <- plot_data %>%
  group_by(CellType, Pathway) %>%
  summarise(
    mean_AD = mean(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    mean_Control = mean(EnrichmentScore[Group == "ctl"], na.rm = TRUE),  # 注意：您的Control组标签是"ctl"
    n_AD = sum(Group == "AD"),
    n_Control = sum(Group == "ctl"),
    sd_AD = sd(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    sd_Control = sd(EnrichmentScore[Group == "ctl"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Delta = mean_AD - mean_Control,
    se = sqrt((sd_AD^2 / n_AD) + (sd_Control^2 / n_Control)),
    CI_lower = Delta - 1.96 * se,
    CI_upper = Delta + 1.96 * se
  )

# 步骤7：执行统计检验（可选，用于标注显著性）---------------------------
# 注意：这里使用Wilcoxon检验，但您也可以选择其他方法
stats_results <- plot_data %>%
  group_by(CellType, Pathway) %>%
  do(tidy(wilcox.test(EnrichmentScore ~ Group, data = .))) %>%
  ungroup() %>%
  mutate(
    Significance = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 合并效应量和显著性
effect_data <- effect_data %>%
  left_join(stats_results %>% select(CellType, Pathway, Significance), 
            by = c("CellType", "Pathway"))

# 步骤8：创建决定性的点图 ----------------------------------------------
delta_plot <- ggplot(effect_data, 
                     aes(x = CellType, y = Delta, color = Pathway)) +
  geom_point(size = 4, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(
    aes(label = Significance, y = CI_upper + 0.05 * max(Delta, na.rm = TRUE)),
    position = position_dodge(width = 0.7),
    size = 5,
    fontface = "bold",
    vjust = 0
  ) +
  scale_color_manual(
    values = c("Early Response" = "#1f77b4", "Late Response" = "#ff7f0e")
  ) +
  labs(
    title = "Effect Size of Estrogen Response Pathway Activity",
    subtitle = "AD vs Control (Δ = AD mean - Control mean)",
    x = "Cell Type",
    y = "Effect Size (Δ)",
    color = "Pathway",
    caption = "*** p < 0.001; ** p < 0.01; * p < 0.05; ns: not significant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

# 步骤9：保存和展示结果 ------------------------------------------------
ggsave("Estrogen_Response_Effect_Size_Plot_Final.png", delta_plot,
       width = 12, height = 8, dpi = 300, bg = "white")

# 保存数据
write.csv(effect_data, "Estrogen_Response_Effect_Sizes_Final.csv", row.names = FALSE)
write.csv(plot_data, "Estrogen_Response_Plot_Data_Final.csv", row.names = FALSE)

# 显示图表
print(delta_plot)

