# 提取兴奋性神经元
ex_neu <- subset(mzl, subset = cell_type_ident == "ex.neu")

# 设置Idents为diagnosis
Idents(ex_neu) <- "diagnosis"

# 寻找AD vs ctl的差异表达基因
markers <- FindMarkers(
  ex_neu,
  ident.1 = "AD",  # AD组
  ident.2 = "ctl", # control组
  logfc.threshold = 0.25,  # 稍微宽松一点
  min.pct = 0.1,  # 在至少10%的细胞中表达
  test.use = "wilcox"  # Wilcoxon秩和检验
)

# 查看top DEGs
top_genes <- markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  arrange(desc(abs(avg_log2FC)))

print(head(top_genes, 20))

# 突触后密度蛋白复合物基因
postsynaptic_density <- c(
  "DLG4", "DLG1", "DLG2", "DLG3",  # PSD-95家族
  "GRIN1", "GRIN2A", "GRIN2B",     # NMDA受体
  "GRIA1", "GRIA2", "GRIA3",       # AMPA受体
  "GABRA1", "GABRA2", "GABRB1",    # GABA受体
  "SYNGAP1", "SHANK1", "SHANK2", "SHANK3",
  "HOMER1", "HOMER2", "HOMER3",
  "CASK", "MAGI2", "LRRTM1", "LRRTM2"
)

# 突触前基因
presynaptic <- c(
  "SNAP25", "SYT1", "SYT2", "STX1A", "STX1B",
  "VAMP1", "VAMP2", "SYN1", "SYN2", "SYN3",
  "RAB3A", "RIMS1", "RIMS2", "UNC13A", "UNC13B",
  "SYPH", "SYP", "NRXN1", "NRXN2", "NRXN3"
)

# 突触囊泡相关基因
synaptic_vesicle <- c(
  "SV2A", "SV2B", "SV2C",
  "SYT1", "SYT2", "SYT3", "SYT4", "SYT5",
  "SYN1", "SYN2", "SYN3",
  "VAMP1", "VAMP2", "VAMP3",
  "SNAP25", "SNAP23", "SNAP29"
)

# 合并所有突触相关基因
synaptic_genes <- unique(c(postsynaptic_density, presynaptic, synaptic_vesicle))

# 转换为小写以匹配Seurat基因名
synaptic_genes <- toupper(synaptic_genes)

# 确保所有必要的包都已加载
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)  # 包含rownames_to_column函数

# 检查包是否加载成功
sessionInfo()

# 使用tibble包的rownames_to_column
synaptic_degs <- top_genes %>%
  tibble::rownames_to_column("gene") %>%
  filter(gene %in% synaptic_genes) %>%
  arrange(desc(avg_log2FC))

print(paste("找到", nrow(synaptic_degs), "个突触相关DEGs"))
print(synaptic_degs)

# 1. 查看详细的DEGs统计
cat("=== 突触相关DEGs详细信息 ===\n")
print(synaptic_degs)

# 2. 可视化这些基因的表达变化
library(ggplot2)
library(ggrepel)

# 创建火山图
volcano_plot <- ggplot(synaptic_degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(aes(color = ifelse(avg_log2FC > 0, "Up", "Down")), size = 3) +
  geom_text_repel(aes(label = gene), size = 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "gray") +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue")) +
  theme_classic() +
  labs(
    title = " (AD vs Control)",
    x = "Log2 Fold Change",
    y = "-Log10(adjusted p-value)",
    color = "Expression Change"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
max_y <- min(400, max(-log10(volcano_data$p_val_adj), na.rm = TRUE))
print(volcano_plot)

# 3. 创建热图展示这8个基因
library(pheatmap)

# 提取这8个基因的表达数据
expr_matrix <- GetAssayData(ex_neu, slot = "data")[synaptic_degs$gene, ]

# 创建注释信息
annotation_df <- data.frame(
  Diagnosis = ex_neu$diagnosis,
  Sex = ex_neu$sex,
  row.names = colnames(ex_neu)
)

# 为了加快计算，可以随机选择部分细胞展示
set.seed(123)
sample_cells <- sample(colnames(ex_neu), size = min(5000, ncol(ex_neu)))

# 绘制热图
pheatmap(expr_matrix[, sample_cells],
         annotation_col = annotation_df[sample_cells, ],
         show_colnames = FALSE,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "突触相关DEGs表达热图",
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         fontsize_row = 10,
         fontsize_col = 8)

# 4. 功能分组分析
# 按基因功能分类
functional_groups <- data.frame(
  gene = synaptic_degs$gene,
  function_group = c(
    "NMDA受体",      # GRIN1
    "PSD支架蛋白",    # DLG2
    "突触黏附分子",   # NRXN2
    "突触支架蛋白",   # MAGI2
    "突触囊泡蛋白",   # SYP
    "SNARE复合体",    # SNAP25
    "突触小泡调节",   # RAB3A
    "突触囊泡蛋白"    # VAMP2
  ),
  direction = ifelse(synaptic_degs$avg_log2FC > 0, "上调", "下调")
)

print(functional_groups)

# 5. 绘制功能分组条形图
functional_summary <- ggplot(functional_groups, aes(x = function_group, fill = direction)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("上调" = "red", "下调" = "blue")) +
  theme_classic() +
  labs(
    title = "突触相关DEGs功能分类",
    x = "功能分类",
    y = "基因数量",
    fill = "表达变化"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14))

print(functional_summary)

# 6. 继续之前的小提琴图分析（添加新的关键基因）
key_genes <- c("SNAP25", "SYT1", "DLG4", "GRIN2B", "SYN1", "GRIN1", "DLG2", "VAMP2")

# 检查哪些基因在数据中存在
available_genes <- key_genes[key_genes %in% rownames(ex_neu)]
cat("可用的基因:", available_genes, "\n")

# 绘制小提琴图
vln_plots <- VlnPlot(
  ex_neu,
  features = available_genes,
  split.by = "diagnosis",
  group.by = "sex",
  pt.size = 0,
  ncol = 4,
  combine = FALSE
)

# 自定义图形
library(patchwork)

for (i in seq_along(vln_plots)) {
  vln_plots[[i]] <- vln_plots[[i]] +
    theme(legend.position = "right",
          plot.title = element_text(size = 12)) +
    scale_fill_manual(values = c("female" = "#FF6B6B", "male" = "#4ECDC4")) +
    xlab("") +
    ylab("Expression")
}

# 组合图形
combined_plot <- wrap_plots(vln_plots, ncol = 4) +
  plot_annotation(
    title = "突触关键基因在男女AD神经元中的表达",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(combined_plot)

# 7. 统计分析
cat("\n=== 统计分析 ===\n")

# 对每个基因进行t检验
for (gene in available_genes) {
  cat("\n基因:", gene, "\n")
  
  expr_data <- FetchData(ex_neu, vars = c(gene, "diagnosis", "sex"))
  colnames(expr_data)[1] <- "expression"
  
  # 按性别分别分析
  for (s in c("female", "male")) {
    subset_data <- expr_data %>% filter(sex == s)
    t_test <- t.test(expression ~ diagnosis, data = subset_data)
    
    mean_ad <- mean(subset_data$expression[subset_data$diagnosis == "AD"], na.rm = TRUE)
    mean_ctl <- mean(subset_data$expression[subset_data$diagnosis == "ctl"], na.rm = TRUE)
    
    cat(s, ": AD平均表达 =", round(mean_ad, 3), 
        ", Control平均表达 =", round(mean_ctl, 3),
        ", Fold Change =", round(mean_ad - mean_ctl, 3),
        ", p-value =", format(t_test$p.value, scientific = TRUE, digits = 3), "\n")
  }
}

# 8. 保存结果
# 创建分析报告目录
dir.create("synaptic_dysfunction_analysis", showWarnings = FALSE)

# 保存数据
write.csv(synaptic_degs, 
          file = "synaptic_dysfunction_analysis/synaptic_related_DEGs.csv",
          row.names = FALSE)

write.csv(functional_groups,
          file = "synaptic_dysfunction_analysis/functional_groups.csv",
          row.names = FALSE)

# 保存图形
ggsave("synaptic_dysfunction_analysis/volcano_plot.png", volcano_plot, width = 10, height = 8)
ggsave("synaptic_dysfunction_analysis/functional_summary.png", functional_summary, width = 10, height = 6)
ggsave("synaptic_dysfunction_analysis/violin_plots.png", combined_plot, width = 16, height = 10)

# 9. 生成分析摘要
cat("\n=== 突触失能假说分析摘要 ===\n")
cat("1. 发现了", nrow(synaptic_degs), "个突触相关差异表达基因\n")
cat("2. 上调基因 (", sum(synaptic_degs$avg_log2FC > 0), "个):", 
    paste(synaptic_degs$gene[synaptic_degs$avg_log2FC > 0], collapse = ", "), "\n")
cat("3. 下调基因 (", sum(synaptic_degs$avg_log2FC < 0), "个):", 
    paste(synaptic_degs$gene[synaptic_degs$avg_log2FC < 0], collapse = ", "), "\n")
cat("4. 变化最显著的基因:", synaptic_degs$gene[which.max(abs(synaptic_degs$avg_log2FC))], 
    "(log2FC =", round(max(abs(synaptic_degs$avg_log2FC)), 3), ")\n")

# 10. 突触功能通路分析
cat("\n=== 突触功能通路受影响情况 ===\n")
cat("• NMDA受体信号: GRIN1 (上调)\n")
cat("• 突触后支架: DLG2 (上调), MAGI2 (上调)\n")
cat("• SNARE复合体: SNAP25 (下调), VAMP2 (下调)\n")
cat("• 突触囊泡循环: RAB3A (下调)\n")
cat("• 突触黏附: NRXN2 (上调)\n")
cat("• 突触囊泡蛋白: SYP (下调)\n")

