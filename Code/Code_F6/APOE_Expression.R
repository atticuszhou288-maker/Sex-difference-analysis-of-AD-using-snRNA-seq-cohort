# 0. 准备：提取目标细胞并创建数据框
library(dplyr)
library(tidyr)

# 0.1 提取星形胶质细胞 (Astrocytes) 数据
ast_cells <- colnames(mzl)[mzl$cell_type_ident == "ast"]
# 获取APOE及其他关键基因在ast中的表达量（使用标准化数据）
ast_expr_data <- FetchData(mzl, 
                           vars = c("APOE", "MEGF10", "LRP1", "CLU"), 
                           cells = ast_cells)
# 合并元数据
ast_meta_data <- mzl@meta.data[ast_cells, c("orig.ident", "diagnosis", "sex")]
ast_df <- cbind(ast_meta_data, ast_expr_data)
ast_df$cell_type <- "Astrocyte"

# 0.2 提取小胶质细胞 (Microglia) 数据
mic_cells <- colnames(mzl)[mzl$cell_type_ident == "mic"]
# 获取APOE表达及功能评分（假设你已计算并存入元数据，如未存入请参考备注）
mic_expr_data <- FetchData(mzl, 
                           vars = c("APOE"), 
                           cells = mic_cells)
mic_meta_data <- mzl@meta.data[mic_cells, c("orig.ident", "diagnosis", "sex")]
# 请确保你的Metabolic_Score和Clearance_Score已作为元数据列存在，否则需要重新计算
# 假设它们已在mzl$Metabolic_Score和mzl$Clearance_Score中
mic_score_data <- mzl@meta.data[mic_cells, c("Metabolic_Score", "Clearance_Score")]
mic_df <- cbind(mic_meta_data, mic_expr_data, mic_score_data)
mic_df$cell_type <- "Microglia"

# 0.3 创建合并数据框用于跨细胞分析
# 按样本（orig.ident）汇总星形胶质细胞的平均APOE
ast_sample_avg <- ast_df %>%
  group_by(orig.ident, diagnosis, sex) %>%
  summarise(ast_APOE_mean = mean(APOE, na.rm = TRUE), .groups = 'drop')
# 按样本汇总小胶质细胞的平均APOE及功能评分
mic_sample_avg <- mic_df %>%
  group_by(orig.ident, diagnosis, sex) %>%
  summarise(
    mic_APOE_mean = mean(APOE, na.rm = TRUE),
    mic_Metabolic_mean = mean(Metabolic_Score, na.rm = TRUE),
    mic_Clearance_mean = mean(Clearance_Score, na.rm = TRUE),
    .groups = 'drop'
  )
# 合并为样本层面数据
sample_level_df <- inner_join(ast_sample_avg, mic_sample_avg, by = c("orig.ident", "diagnosis", "sex"))

# 4. APOE在所有细胞类型中的表达分布（小提琴图）。
# 方法A：使用Seurat的VlnPlot（针对特定基因较快）
library(Seurat)
# 提取所有细胞的APOE表达和元数据，用于ggplot2自定义绘图（更灵活）
expr_apoe <- FetchData(mzl, vars = "APOE")
meta_for_plot <- mzl@meta.data[, c("cell_type_ident", "diagnosis", "sex")]
plot_df <- cbind(meta_for_plot, expr_apoe)

# 为了绘图清晰，可以按细胞类型、诊断、性别分组后展示
# 使用ggplot2绘制分面小提琴图
p_global_vln <- ggplot(plot_df, aes(x = cell_type_ident, y = APOE, fill = diagnosis)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  facet_grid(sex ~ ., scales = "free_y") + # 按性别分面
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(x = "Cell Type", y = "APOE Expression", 
       title = "APOE Expression Across Major Cell Types",
       fill = "Diagnosis") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_global_vln)

ggsave(filename = "APOE_Global_Expression_Violin.png", plot = p_global_vln,
       width = 10, height = 7, dpi = 300)

library(Seurat)
library(ggplot2)
library(ggpubr) # 用于添加统计显著性标记
library(dplyr)

# 1. 准备数据（与你的代码一致）
expr_apoe <- FetchData(mzl, vars = "APOE")
meta_for_plot <- mzl@meta.data[, c("cell_type_ident", "diagnosis", "sex")]
plot_df <- cbind(meta_for_plot, expr_apoe)

# 2. 创建绘图对象并添加显著性标记
p_global_vln <- ggplot(plot_df, aes(x = cell_type_ident, y = APOE, fill = diagnosis)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7) +
  # 核心：添加组间比较的显著性标记
  stat_compare_means(
    aes(group = diagnosis),           # 指定要比较的分组
    method = "wilcox.test",           # 使用非参数Wilcoxon秩和检验
    label = "p.signif",               # 显示符号标记 (ns, *, **, ***, ****)
    # label = "p.format",             # 或者，显示具体p值 (例如 p=0.03)
    hide.ns = FALSE,                  # TRUE将隐藏不显著(ns)的标记，这里我们选择显示
    show.legend = FALSE,              # 不显示图例
    size = 4,                         # 标记字体大小
    vjust = 0.5,                      # 垂直微调标记位置
    # 调整标记的y轴位置，这里设定在每组数据最大值的上方
    # bracket.nudge.y = 0.05           # 可以微调标记的垂直位置
  ) +
  facet_grid(sex ~ ., scales = "free_y") + # 按性别分面
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(x = "Cell Type", y = "APOE Expression", 
       title = "APOE Expression Across Major Cell Types",
       fill = "Diagnosis",
       subtitle = "Statistical test: Wilcoxon test; * p<0.05, ** p<0.01, *** p<0.001, **** p<0.0001") + # 添加副标题说明
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. 显示图形
print(p_global_vln)

# 4. 保存图形
ggsave(filename = "APOE_Global_Expression_Violin_with_Stats.png", 
       plot = p_global_vln,
       width = 10, 
       height = 7, 
       dpi = 300)