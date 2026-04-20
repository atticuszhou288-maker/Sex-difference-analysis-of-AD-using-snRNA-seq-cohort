library(Seurat)
library(tidyverse)
female_cells <- mzl[, mzl$sex == "female"]

# 验证细胞数量
print(paste("女性细胞数量:", ncol(female_cells)))
print(paste("占总数比例:", round(ncol(female_cells)/ncol(mzl)*100, 1), "%"))

# 2. 提取标准化表达数据（RNA@data）
female_expr <- GetAssayData(female_cells, slot = "data")

# 3. 确认诊断分组
print("诊断分组分布:")
table(female_cells@meta.data$diagnosis)

# 4. 准备分析基因列表
target_genes <- c("ESR1", "BDNF", "SIRT1", "PGR", "NRCAM", "TFF1", "GREB1", 
                  "CCND1", "MYC", "BCL2", "VEGFA", "ESR2")

# 筛选实际存在的基因
valid_genes <- rownames(female_expr)[rownames(female_expr) %in% target_genes]
print(paste("有效分析基因:", paste(valid_genes, collapse = ", ")))

# 加载必要库
library(Seurat)
library(tidyverse)
library(corrplot)
library(igraph)
library(ggraph)
library(patchwork)
library(ggrepel)

# 1. 有效基因列表（根据实际存在调整）
valid_genes <- c("ESR1", "BDNF", "SIRT1", "PGR", "NRCAM", 
                 "GREB1", "VEGFA", "MYC", "CCND1", "ESR2", "BCL2")

# 2. 提取标准化表达数据
female_expr <- GetAssayData(female_cells, slot = "data")

# 3. 创建结果存储列表
cor_results <- list()
p_results <- list()

# 4. 计算斯皮尔曼相关系数
for (diag_group in c("AD", "ctl")) {
  # 获取当前组的细胞ID
  diag_cells <- colnames(female_cells)[female_cells$diagnosis == diag_group]
  
  # 提取基因表达矩阵
  expr_subset <- t(as.matrix(female_expr[valid_genes, diag_cells]))
  
  # 初始化相关性和p值矩阵
  cor_matrix <- matrix(NA, nrow = length(valid_genes), ncol = length(valid_genes), 
                       dimnames = list(valid_genes, valid_genes))
  p_matrix <- matrix(NA, nrow = length(valid_genes), ncol = length(valid_genes), 
                     dimnames = list(valid_genes, valid_genes))
  
  # 计算每对基因的相关性和p值
  for (i in 1:length(valid_genes)) {
    for (j in 1:length(valid_genes)) {
      if (i != j) {  # 跳过自相关
        gene1 <- valid_genes[i]
        gene2 <- valid_genes[j]
        
        # 确保两组数据都有表达
        if (sum(expr_subset[, gene1] > 0) > 10 && sum(expr_subset[, gene2] > 0) > 10) {
          test_result <- cor.test(expr_subset[, gene1], expr_subset[, gene2], 
                                  method = "spearman", exact = FALSE)
          cor_matrix[i, j] <- test_result$estimate
          p_matrix[i, j] <- test_result$p.value
        }
      }
    }
  }
  
  # 存储结果
  cor_results[[diag_group]] <- cor_matrix
  p_results[[diag_group]] <- p_matrix
}

# 5. 可视化：分组相关性热图
# 创建热图函数
create_corr_heatmap <- function(cor_matrix, p_matrix, title) {
  # 标记不显著的相关性
  cor_to_plot <- cor_matrix
  cor_to_plot[p_matrix > 0.05] <- 0
  
  corrplot(
    cor_to_plot,
    method = "color",
    type = "upper",
    tl.col = "black",
    tl.cex = 0.8,
    title = title,
    mar = c(0, 0, 2, 0),
    col = colorRampPalette(c("blue", "white", "red"))(200),
    diag = FALSE,
    cl.pos = "r",
    cl.ratio = 0.1,
    number.cex = 0.7  # 减小数字大小
  )
}

# 创建并排热图
png("ESR1_Target_Correlations_Heatmaps.png", width = 16, height = 10, units = "in", res = 300)
par(mfrow = c(1, 2))
create_corr_heatmap(cor_results$AD, p_results$AD, "Female AD: ESR1-Target Correlations")
create_corr_heatmap(cor_results$ctl, p_results$ctl, "Female Control: ESR1-Target Correlations")
par(mfrow = c(1, 1))
dev.off()

# 6. 可视化：相关性网络图
# 创建网络图函数
create_corr_network <- function(cor_matrix, p_matrix, title, min_cor = 0.2) {
  # 创建边列表
  edges <- data.frame()
  genes <- rownames(cor_matrix)
  
  for (i in 1:length(genes)) {
    for (j in 1:length(genes)) {
      if (i != j) {  # 跳过自相关
        cor_val <- cor_matrix[i, j]
        p_val <- p_matrix[i, j]
        
        if (!is.na(p_val) && p_val < 0.05 && abs(cor_val) > min_cor) {
          edges <- rbind(edges, data.frame(
            from = genes[i],
            to = genes[j],
            correlation = cor_val,
            p.value = p_val
          ))
        }
      }
    }
  }
  
  # 如果边数为0，创建空图
  if (nrow(edges) == 0) {
    plot.new()
    title(main = title, sub = "No significant correlations found")
    return(NULL)
  }
  
  # 创建节点列表
  nodes <- data.frame(
    id = genes,
    label = genes,
    type = ifelse(genes == "ESR1", "core", "target"),
    size = ifelse(genes == "ESR1", 15, 8)
  )
  
  # 创建网络图
  net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
  
  # 设置边的属性
  E(net)$width <- abs(E(net)$correlation) * 5
  E(net)$color <- ifelse(E(net)$correlation > 0, "red", "blue")
  
  # 创建ggraph图
  ggraph(net, layout = "fr") +  # 使用力导向布局
    geom_edge_link(
      aes(width = width, color = color),
      alpha = 0.7,
      show.legend = TRUE
    ) +
    geom_node_point(
      aes(size = size, fill = type),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    geom_node_text(
      aes(label = label),
      repel = TRUE,
      size = 4,
      fontface = "bold"
    ) +
    scale_size_identity() +
    scale_fill_manual(
      values = c("core" = "#FF6B6B", "target" = "#4ECDC4"),
      guide = "none"
    ) +
    scale_edge_color_identity() +
    scale_edge_width_continuous(
      name = "|Correlation|",
      range = c(0.5, 3)
    ) +
    theme_void() +
    labs(
      title = title,
      subtitle = "Significant correlations (p < 0.05) shown"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "bottom"
    )
}

# 创建网络图
net_AD <- create_corr_network(
  cor_results$AD, 
  p_results$AD, 
  "Female AD: ESR1-Target Network"
)

net_CTL <- create_corr_network(
  cor_results$ctl, 
  p_
  results$ctl, 
  "Female Control: ESR1-Target Network"
)

# 并排显示网络图
if (!is.null(net_AD) && !is.null(net_CTL)) {
  combined_nets <- net_AD + net_CTL +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = "Estrogen Signaling Compensation Networks",
      subtitle = "Comparison between AD and Control Groups",
      theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
    )
  
  # 保存网络图
  ggsave("ESR1_Target_Networks.png", combined_nets, 
         width = 16, height = 8, dpi = 300)
} else {
  # 如果网络图有问题，尝试保存单图
  if (!is.null(net_AD)) ggsave("ESR1_Network_AD.png", net_AD, width = 8, height = 8, dpi = 300)
  if (!is.null(net_CTL)) ggsave("ESR1_Network_CTL.png", net_CTL, width = 8, height = 8, dpi = 300)
}

# 7. 结果分析表格
# 创建结果数据框
results_df <- data.frame()

for (gene in valid_genes) {
  if (gene != "ESR1") {
    # 获取ESR1与该基因的相关性
    ad_cor <- cor_results$AD["ESR1", gene]
    ad_p <- p_results$AD["ESR1", gene]
    
    ctl_cor <- cor_results$ctl["ESR1", gene]
    ctl_p <- p_results$ctl["ESR1", gene]
    
    results_df <- rbind(results_df, data.frame(
      Target_Gene = gene,
      AD_Correlation = round(ad_cor, 3),
      AD_Pvalue = round(ad_p, 4),
      Control_Correlation = round(ctl_cor, 3),
      Control_Pvalue = round(ctl_p, 4),
      Correlation_Difference = round(ad_cor - ctl_cor, 3)
    ))
  }
}

# 添加显著性标记
results_df <- results_df %>%
  mutate(
    AD_Significance = case_when(
      AD_Pvalue < 0.001 ~ "***",
      AD_Pvalue < 0.01 ~ "**",
      AD_Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    Control_Significance = case_when(
      Control_Pvalue < 0.001 ~ "***",
      Control_Pvalue < 0.01 ~ "**",
      Control_Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    AD_Correlation_Label = paste0(AD_Correlation, AD_Significance),
    Control_Correlation_Label = paste0(Control_Correlation, Control_Significance)
  )

# 保存结果
write.csv(results_df, "ESR1_Target_Correlation_Results.csv", row.names = FALSE)

# 8. 关键可视化：相关性变化点图
cor_change_plot <- ggplot(results_df, 
                          aes(x = AD_Correlation, y = Control_Correlation,
                              color = Correlation_Difference,
                              size = abs(Correlation_Difference))) +
  geom_point(alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_text_repel(aes(label = Target_Gene), size = 4.5, box.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray40") +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", 
                        midpoint = 0, name = "Correlation\nDifference") +
  scale_size_continuous(range = c(4, 12), name = "|Correlation\nDifference|") +
  labs(
    title = "ESR1-Target Gene Correlation Changes",
    subtitle = "AD vs Control Comparison in Female Cells",
    x = "AD Group Correlation",
    y = "Control Group Correlation",
    caption = "Dotted line indicates no correlation (r=0)\nDiagonal line indicates equal correlation in both groups"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 14),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1  # 保持正方形比例
  )

# 保存点图
ggsave("ESR1_Correlation_Changes.png", cor_change_plot, 
       width = 10, height = 8, dpi = 300)

# 9. 统计结果报告
# 计算组间差异的统计显著性
cor_diff_test <- cor.test(results_df$AD_Correlation, results_df$Control_Correlation, 
                          method = "spearman", exact = FALSE)

# 创建报告文本
report_text <- paste(
  "## Estrogen Signaling Compensation Analysis Report\n\n",
  "### Data Summary\n",
  "- Female AD cells: ", sum(female_cells$diagnosis == "AD"), "\n",
  "- Female Control cells: ", sum(female_cells$diagnosis == "ctl"), "\n",
  "- Genes analyzed: ", paste(valid_genes, collapse = ", "), "\n\n",
  
  "### Key Findings\n",
  "1. **Overall correlation pattern**:\n",
  "   - Spearman correlation between AD and Control group correlations: ρ = ", 
  round(cor_diff_test$estimate, 3), " (p = ", round(cor_diff_test$p.value, 4), ")\n\n",
  
  "2. **Genes with strongest correlation changes**:\n"
)

# 添加变化最大的基因
top_changed <- results_df[order(-abs(results_df$Correlation_Difference)), ]
for (i in 1:min(3, nrow(top_changed))) {
  gene <- top_changed$Target_Gene[i]
  diff <- top_changed$Correlation_Difference[i]
  report_text <- paste0(
    report_text,
    "   - ", gene, ": Δr = ", round(diff, 3), 
    " (AD: ", top_changed$AD_Correlation[i], 
    ", Control: ", top_changed$Control_Correlation[i], ")\n"
  )
}

# 添加显著相关基因信息
significant_AD <- results_df[results_df$AD_Pvalue < 0.05, ]
if (nrow(significant_AD) > 0) {
  report_text <- paste0(
    report_text,
    "\n3. **Genes significantly correlated with ESR1 in AD**:\n",
    paste("   - ", significant_AD$Target_Gene, ": r = ", significant_AD$AD_Correlation, 
          " (p = ", significant_AD$AD_Pvalue, ")", collapse = "\n"),
    "\n"
  )
} else {
  report_text <- paste0(report_text, "\n3. No significant correlations in AD group\n")
}

significant_CTL <- results_df[results_df$Control_Pvalue < 0.05, ]
if (nrow(significant_CTL) > 0) {
  report_text <- paste0(
    report_text,
    "\n4. **Genes significantly correlated with ESR1 in Control**:\n",
    paste("   - ", significant_CTL$Target_Gene, ": r = ", significant_CTL$Control_Correlation, 
          " (p = ", significant_CTL$Control_Pvalue, ")", collapse = "\n"),
    "\n"
  )
} else {
  report_text <- paste0(report_text, "\n4. No significant correlations in Control group\n")
}

# 添加解释和建议
report_text <- paste0(
  report_text,
  "\n### Biological Interpretation\n",
  "1. Positive correlation difference (Δr > 0) indicates stronger estrogen signaling connectivity in AD\n",
  "2. Negative correlation difference (Δr < 0) suggests disrupted estrogen signaling in AD\n\n",
  
  "### Recommendations\n",
  "1. Validate top changed genes (", paste(top_changed$Target_Gene[1:3], collapse = ", "), 
  ") in independent dataset\n",
  "2. Perform pathway analysis on correlated gene sets\n",
  "3. Investigate cell-type specific effects of estrogen signaling\n"
)

# 保存报告
writeLines(report_text, "ESR1_Correlation_Analysis_Report.md")

# 10. 完成所有分析
print("===================================================")
print("雌激素信号代偿分析已完成！生成结果包括：")
print("1. 热图对比：ESR1_Target_Correlations_Heatmaps.png")
print("2. 网络图：ESR1_Target_Networks.png")
print("3. 相关性变化点图：ESR1_Correlation_Changes.png")
print("4. 结果表格：ESR1_Target_Correlation_Results.csv")
print("5. 分析报告：ESR1_Correlation_Analysis_Report.md")
print("===================================================")

# 定义基因列表
genes <- c("ESR1", "ESR2", "PGR", "GREB1", "BDNF", "SIRT1", "NRCAM", "VEGFA", "MYC", "CCND1", "BCL2")

# 定义基因功能关系
edges <- data.frame(
  from = c("ESR1", "ESR1", "ESR1", "ESR1", "ESR1", "Estrogen Response", "Estrogen Response"),
  to = c("ESR2", "PGR", "GREB1", "SIRT1", "BDNF", "Neuroprotection", "Cell Proliferation"),
  relationship = c("Coreceptor", "Primary Target", "Primary Target", "Secondary Target", "Secondary Target", "Functional Pathway", "Functional Pathway")
)

# 添加二级关系
edges <- rbind(edges, data.frame(
  from = c("Estrogen Response", "Estrogen Response", "Neuroprotection", "Neuroprotection", "Cell Proliferation"),
  to = c("NRCAM", "VEGFA", "BCL2", "CCND1", "MYC"),
  relationship = c("Downstream Effector", "Downstream Effector", "Downstream Effector", "Downstream Effector", "Downstream Effector")
))

# 获取所有唯一节点ID
node_ids <- unique(c(edges$from, edges$to))

# 创建节点数据框 - 修复错误
nodes <- data.frame(
  id = node_ids,
  type = ifelse(node_ids %in% genes, "Gene", "Function"),
  group = ifelse(node_ids %in% genes, "Gene", "Function")
)

# 创建网络
net <- graph_from_data_frame(d = edges, vertices = nodes)

# 设置布局
set.seed(123)
layout <- create_layout(net, layout = "fr")

# 可视化
gene_network <- ggraph(layout) +
  geom_edge_link(
    aes(color = relationship),
    arrow = arrow(length = unit(2, 'mm')),
    end_cap = circle(5, 'mm'),
    start_cap = circle(5, 'mm'),
    alpha = 0.8,
    width = 1
  ) +
  geom_node_point(
    aes(size = ifelse(type == "Gene", 15, 10), fill = type),
    shape = 21,
    color = "black",
    stroke = 1.5
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,
    size = 5,
    fontface = "bold"
  ) +
  scale_edge_color_manual(
    values = c(
      "Coreceptor" = "#8c96c6",
      "Primary Target" = "#8c6bb1",
      "Secondary Target" = "#88419d",
      "Functional Pathway" = "#6e016b",
      "Downstream Effector" = "#4d004b"
    ),
    name = "Relationship"
  ) +
  scale_fill_manual(
    values = c("Gene" = "#fdcc8a", "Function" = "#b3cde3"),
    name = "Node Type"
  ) +
  scale_size_identity() +
  theme_void() +
  labs(
    title = "Estrogen Signaling Pathway Gene Relationships",
    subtitle = "Rationale for gene selection in correlation analysis",
    caption = "Gene relationships based on established biological pathways and estrogen signaling mechanisms"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 16),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

# 显示基因关系网络图
print(gene_network)

# 保存基因关系网络图
ggsave("Estrogen_Gene_Relationship_Network.png", gene_network, 
       width = 14, height = 10, dpi = 300, bg = "white")



