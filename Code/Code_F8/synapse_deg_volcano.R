# 修正火山图 - 处理p值为0的情况
volcano_data <- synaptic_degs

# 检查并修正p值为0的情况
if(any(volcano_data$p_val == 0)) {
  cat("发现p值为0的基因，进行修正...\n")
  
  # 找到最小非零p值
  min_nonzero_p <- min(volcano_data$p_val[volcano_data$p_val > 0], na.rm = TRUE)
  replacement_p <- min_nonzero_p / 10  # 使用最小非零p值的1/10
  
  # 修正p值为0的情况
  volcano_data$p_val[volcano_data$p_val == 0] <- replacement_p
  
  # 重新计算调整后的p值（如果需要）
  volcano_data$p_val_adj <- p.adjust(volcano_data$p_val, method = "BH")
}

# 添加显著性标签
volcano_data$Significance <- ifelse(
  volcano_data$p_val_adj < 0.05 & abs(volcano_data$avg_log2FC) > 0.25,
  "Significant",
  "Non-significant"
)

volcano_data$Direction <- ifelse(
  volcano_data$avg_log2FC > 0, "Upregulated", "Downregulated"
)

# 限制y轴范围
max_y_val <- 320
current_max <- max(-log10(volcano_data$p_val_adj), na.rm = TRUE)
if(current_max > max_y_val) {
  cat("调整y轴范围: 当前最大值 =", round(current_max, 2), 
      "，限制为", max_y_val, "\n")
  volcano_data$log10_padj <- pmin(-log10(volcano_data$p_val_adj), max_y_val)
} else {
  volcano_data$log10_padj <- -log10(volcano_data$p_val_adj)
}

# 创建火山图
volcano_plot <- ggplot(volcano_data, 
                       aes(x = avg_log2FC, y = log10_padj)) +
  # 添加网格线
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed", color = "gray60", alpha = 0.7) +
  geom_vline(xintercept = c(-0.25, 0, 0.25), 
             linetype = "dashed", color = "gray60", alpha = 0.7) +
  
  # 添加点
  geom_point(aes(color = Direction, size = abs(avg_log2FC)), 
             alpha = 0.8, shape = 19) +
  
  # 添加基因标签
  geom_text_repel(
    aes(label = gene),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.color = "gray50",
    segment.alpha = 0.5
  ) +
  
  # 设置颜色和大小
  scale_color_manual(
    values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")
  ) +
  scale_size_continuous(
    range = c(3, 6),
    name = "|Log2FC|"
  ) +
  
  # 坐标轴和标签
  labs(
    title = "Volcano Plot: Synaptic-Related DEGs in AD vs Control",
    x = "Log2 Fold Change (AD vs Control)",
    y = "-Log10(Adjusted P-value)",
    color = "Expression Change"
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  ) +
  
  # 坐标轴范围
  xlim(min(volcano_data$avg_log2FC) - 0.1, max(volcano_data$avg_log2FC) + 0.1) +
  ylim(0, max(volcano_data$log10_padj) * 1.1)

print(volcano_plot)

# 保存火山图
ggsave("synaptic_volcano_plot.png", volcano_plot, 
       width = 12, height = 9, dpi = 300)