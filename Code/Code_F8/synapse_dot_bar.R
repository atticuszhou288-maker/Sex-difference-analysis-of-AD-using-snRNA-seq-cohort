# 1. 点图改进：只标注不显著的基因
expression_dot_plot_complete <- ggplot(
  gene_expression_summary, 
  aes(x = Log2FC, y = reorder(Gene, Log2FC))
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(aes(color = Category, size = -log10(P_Value + 1e-300)), 
             alpha = 0.8) +
  # 为不显著的基因添加ns标签，显著基因只显示数值
  geom_text_repel(
    aes(label = ifelse(P_Value >= 0.05, "ns", sprintf("%.2f", Log2FC))), 
    size = 3.5, 
    box.padding = 0.3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c(
    "Postsynaptic" = "#1B9E77",
    "Presynaptic" = "#D95F02",
    "Synaptic Adhesion" = "#7570B3"
  )) +
  scale_size_continuous(range = c(3, 8), name = "-log10(p-value)") +
  labs(
    title = "Synaptic Genes Expression Changes in AD vs Control",
    x = "Log2 Fold Change (AD vs Control)",
    y = "Gene",
    subtitle = "Postsynaptic (green), Presynaptic (orange), Adhesion (purple)",
    color = "Gene Category"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )

print(expression_dot_plot_complete)

# 保存
ggsave("synaptic_genes_complete_dotplot_v2.png", 
       expression_dot_plot_complete, 
       width = 10, height = 8, dpi = 300)

# 2. 柱状图改进：调整y轴范围，让所有柱子都可见
functional_summary <- gene_expression_summary %>%
  group_by(Category) %>%
  summarise(
    n_genes = n(),
    n_up = sum(Direction == "Up"),
    n_down = sum(Direction == "Down"),
    mean_log2fc = mean(Log2FC, na.rm = TRUE),
    sd_log2fc = sd(Log2FC, na.rm = TRUE)
  )

print(functional_summary)

# 计算y轴范围
y_min <- min(functional_summary$mean_log2fc - functional_summary$sd_log2fc, na.rm = TRUE)
y_max <- max(functional_summary$mean_log2fc + functional_summary$sd_log2fc, na.rm = TRUE)
y_range <- max(abs(y_min), abs(y_max))

functional_bar <- ggplot(functional_summary, aes(x = Category, y = mean_log2fc)) +
  geom_bar(stat = "identity", aes(fill = Category), alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_log2fc - sd_log2fc, 
                    ymax = mean_log2fc + sd_log2fc), 
                width = 0.2) +
  geom_text(aes(label = sprintf("%.2f", mean_log2fc)), 
            vjust = -0.5, size = 4) +
  geom_text(aes(y = y_min - 0.1, label = paste0("n=", n_genes)), 
            size = 3.5, color = "gray40") +
  scale_fill_manual(values = c(
    "Postsynaptic" = "#1B9E77",
    "Presynaptic" = "#D95F02",
    "Synaptic Adhesion" = "#7570B3"
  )) +
  labs(
    title = "Average Expression Changes by Synaptic Category",
    x = "Synaptic Category",
    y = "Average Log2 Fold Change (AD vs Control)",
    fill = "Category"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none"
  ) +
  # 设置对称的y轴范围
  ylim(y_min - 0.15, y_max + 0.15)

print(functional_bar)

ggsave("synaptic_functional_summary_bar_v2.png", 
       functional_bar, 
       width = 8, height = 6, dpi = 300)

# 3. 热图改进：颜色梯度调整为2,0,-2，0为白色
# 重新提取表达数据
exact_gene_names <- rownames(expr_matrix)[toupper(rownames(expr_matrix)) %in% 
                                            toupper(gene_expression_summary$Gene)]

hc_expr_data <- expr_matrix[exact_gene_names, ]

# 随机选择部分细胞（为了可视化清晰）
set.seed(123)
n_cells <- min(2000, ncol(hc_expr_data))
sampled_cells <- sample(colnames(hc_expr_data), n_cells)

# 创建注释
annotation_df <- data.frame(
  Diagnosis = ex_neu@meta.data[sampled_cells, "diagnosis"],
  Sex = ex_neu@meta.data[sampled_cells, "sex"],
  row.names = sampled_cells
)

# 设置颜色
ann_colors <- list(
  Diagnosis = c(AD = "#E41A1C", ctl = "#377EB8"),
  Sex = c(female = "#F781BF", male = "#4DAF4A")
)

# 绘制热图，颜色梯度为2,0,-2，0为白色
library(pheatmap)
heatmap_v2 <- pheatmap(
  hc_expr_data[, sampled_cells],
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  scale = "row",  # 按行标准化
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2, 2, length.out = 101),  # 设置颜色断点
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  fontsize_row = 10,
  main = "Synaptic Genes Expression - Scaled by Row",
  border_color = NA
)

# 保存热图
png("synaptic_genes_heatmap_v2.png", width = 1000, height = 800, res = 150)
print(heatmap_v2)
dev.off()