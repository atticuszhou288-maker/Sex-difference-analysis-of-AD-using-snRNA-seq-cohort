# 创建性别特异性的热图，显示每个基因在男性和女性中的Log2FC

# 1. 重新计算性别特异性的Log2FC
# 提取每个基因在男性和女性中的AD vs Control的Log2FC
sex_specific_fc <- data.frame()

for(gene in gene_expression_summary$Gene) {
  # 找到确切的基因名
  exact_name <- rownames(expr_matrix)[toupper(rownames(expr_matrix)) == gene]
  
  if(length(exact_name) > 0) {
    # 提取表达数据
    expr_values <- expr_matrix[exact_name[1], ]
    
    # 计算女性中的log2FC
    female_ad <- mean(expr_values[ex_neu$sex == "female" & ex_neu$diagnosis == "AD"])
    female_ctl <- mean(expr_values[ex_neu$sex == "female" & ex_neu$diagnosis == "ctl"])
    female_fc <- log2((female_ad + 0.001) / (female_ctl + 0.001))
    
    # 计算男性中的log2FC
    male_ad <- mean(expr_values[ex_neu$sex == "male" & ex_neu$diagnosis == "AD"])
    male_ctl <- mean(expr_values[ex_neu$sex == "male" & ex_neu$diagnosis == "ctl"])
    male_fc <- log2((male_ad + 0.001) / (male_ctl + 0.001))
    
    sex_specific_fc <- rbind(sex_specific_fc, data.frame(
      Gene = gene,
      Female_FC = female_fc,
      Male_FC = male_fc,
      Category = gene_expression_summary$Category[gene_expression_summary$Gene == gene]
    ))
  }
}

print(sex_specific_fc)

# 2. 转换为矩阵格式
heatmap_data <- sex_specific_fc %>%
  select(Gene, Female_FC, Male_FC) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# 重命名列，使其更清晰
colnames(heatmap_data) <- c("Female", "Male")

# 3. 创建行注释（基因类别）
row_annotation <- sex_specific_fc %>%
  select(Gene, Category) %>%
  distinct() %>%
  column_to_rownames("Gene")

# 4. 设置颜色
# 行注释颜色
row_colors <- list(
  Category = c(
    "Postsynaptic" = "#1B9E77",
    "Presynaptic" = "#D95F02",
    "Synaptic Adhesion" = "#7570B3"
  )
)

# 热图颜色：从蓝色(-2)到白色(0)到红色(2)
# 创建101个颜色的梯度
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(101)
breaks <- seq(-2, 2, length.out = 102)  # 102个断点，对应101个颜色

# 5. 绘制热图
library(pheatmap)

sex_specific_heatmap <- pheatmap(
  heatmap_data,
  annotation_row = row_annotation,
  annotation_colors = row_colors,
  color = heatmap_colors,
  breaks = breaks,
  cluster_rows = TRUE,   # 聚类行
  cluster_cols = FALSE,  # 不聚类列，保持Female和Male顺序
  scale = "none",        # 不进行标准化，因为我们已经是log2FC
  main = "Sex-Specific Expression Changes of Synaptic Genes\n(AD vs Control, Log2 Fold Change)",
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize = 10,
  angle_col = 0,
  border_color = NA,
  cellwidth = 40,
  cellheight = 20,
  display_numbers = TRUE,  # 在单元格中显示数值
  number_format = "%.2f",  # 保留两位小数
  number_color = "black",  # 数值颜色
  fontsize_number = 8      # 数值字体大小
)

# 保存热图
png("sex_specific_synaptic_changes_v2.png", width = 800, height = 900, res = 300)
print(sex_specific_heatmap)
dev.off()

# 6. 创建另一个版本：按类别排序
# 首先按类别排序，然后在类别内按平均FC排序
ordered_genes <- sex_specific_fc %>%
  mutate(Avg_FC = (Female_FC + Male_FC) / 2) %>%
  arrange(Category, desc(Avg_FC)) %>%
  pull(Gene)

# 按排序重新组织数据
heatmap_data_ordered <- heatmap_data[ordered_genes, ]

# 绘制按类别排序的热图
sex_specific_heatmap_ordered <- pheatmap(
  heatmap_data_ordered,
  annotation_row = row_annotation[ordered_genes, , drop = FALSE],
  annotation_colors = row_colors,
  color = heatmap_colors,
  breaks = breaks,
  cluster_rows = FALSE,   # 不聚类行，使用我们定义的顺序
  cluster_cols = FALSE,   # 不聚类列
  scale = "none",
  main = "Sex-Specific Expression Changes by Category",
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize = 10,
  angle_col = 0,
  border_color = NA,
  cellwidth = 40,
  cellheight = 20,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 8,
  gaps_row = cumsum(table(row_annotation[ordered_genes, "Category"]))  # 在类别之间添加间隔
)

# 保存排序版本的热图
png("sex_specific_synaptic_changes_by_category.png", width = 800, height = 900, res = 300)
print(sex_specific_heatmap_ordered)
dev.off()

# 7. 创建差异热图：展示性别差异
# 计算男性FC - 女性FC，显示性别差异
gender_diff <- sex_specific_fc %>%
  mutate(Gender_Difference = Male_FC - Female_FC) %>%
  select(Gene, Gender_Difference, Category) %>%
  column_to_rownames("Gene")

# 设置颜色梯度：从蓝色(女性更高)到白色(无差异)到红色(男性更高)
diff_colors <- colorRampPalette(c("blue", "white", "red"))(101)
diff_breaks <- seq(-1, 1, length.out = 102)  # 限制范围在-1到1之间

# 绘制性别差异热图
gender_diff_heatmap <- pheatmap(
  as.matrix(gender_diff[, "Gender_Difference", drop = FALSE]),
  annotation_row = gender_diff[, "C", drop = FALSE],
  annotation_colors = row_colors,
  color = diff_colors,
  breaks = diff_breaks,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  main = "Gender Differences in Synaptic Gene Expression Changes\n(Male FC - Female FC)",
  fontsize_row = 10,
  fontsize_col = 12,
  fontsize = 10,
  angle_col = 0,
  border_color = NA,
  cellwidth = 40,
  cellheight = 20,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  fontsize_number = 8
)

# 保存性别差异热图
png("gender_differences_synaptic_changes.png", width = 600, height = 900, res = 150)
print(gender_diff_heatmap)
dev.off()

# 8. 创建总结图：综合展示性别特异性结果
library(ggplot2)
library(tidyr)

# 将数据转换为长格式
sex_fc_long <- sex_specific_fc %>%
  pivot_longer(cols = c(Female_FC, Male_FC), 
               names_to = "Sex", 
               values_to = "Log2FC") %>%
  mutate(Sex = gsub("_FC", "", Sex))

# 创建点图展示性别差异
gender_comparison_plot <- ggplot(sex_fc_long, aes(x = Sex, y = Log2FC, group = Gene)) +
  geom_line(aes(color = Category), alpha = 0.5) +
  geom_point(aes(color = Category), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ Category, ncol = 1, scales = "free_y") +
  scale_color_manual(values = c(
    "Postsynaptic" = "#1B9E77",
    "Presynaptic" = "#D95F02",
    "Synaptic Adhesion" = "#7570B3"
  )) +
  labs(
    title = "Gender Comparison of Synaptic Gene Expression Changes",
    x = "Sex",
    y = "Log2 Fold Change (AD vs Control)",
    color = "Category"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(gender_comparison_plot)
ggsave("gender_comparison_plot.png", gender_comparison_plot, 
       width = 8, height = 10, dpi = 300)

# 9. 统计性别特异性结果
cat("\n=== 性别特异性分析结果 ===\n\n")

# 统计每个基因的性别差异
significant_gender_diff <- gender_diff %>%
  mutate(Significant_Diff = abs(Gender_Difference) > 0.1)  # 设置阈值

cat("基因性别差异统计:\n")
cat(sprintf("• 总基因数: %d\n", nrow(significant_gender_diff)))
cat(sprintf("• 性别差异 > 0.1的基因数: %d\n", sum(significant_gender_diff$Significant_Diff)))

# 查看性别差异最大的基因
top_gender_diff <- gender_diff %>%
  arrange(desc(abs(Gender_Difference))) %>%
  head(5)

cat("\n性别差异最大的基因:\n")
for(i in 1:nrow(top_gender_diff)) {
  cat(sprintf("  %d. %s: 男性FC - 女性FC = %.3f\n", 
              i, 
              rownames(top_gender_diff)[i],
              top_gender_diff$Gender_Difference[i]))
}

# 按类别统计性别差异
cat("\n按类别统计性别差异:\n")
for(cat in unique(sex_specific_fc$Category)) {
  cat_genes <- sex_specific_fc %>% filter(Category == cat)
  mean_female_fc <- mean(cat_genes$Female_FC)
  mean_male_fc <- mean(cat_genes$Male_FC)
  mean_diff <- mean_male_fc - mean_female_fc
  
  cat(sprintf("  %s:\n", cat))
  cat(sprintf("    • 平均女性FC: %.3f\n", mean_female_fc))
  cat(sprintf("    • 平均男性FC: %.3f\n", mean_male_fc))
  cat(sprintf("    • 平均性别差异: %.3f\n", mean_diff))
  cat(sprintf("    • 基因数: %d\n\n", nrow(cat_genes)))
}

# 10. 保存所有性别特异性分析结果
gender_analysis_dir <- "gender_specific_analysis"
dir.create(gender_analysis_dir, showWarnings = FALSE)

# 保存数据
write.csv(sex_specific_fc, 
          file.path(gender_analysis_dir, "sex_specific_fc_data.csv"),
          row.names = FALSE)

write.csv(gender_diff, 
          file.path(gender_analysis_dir, "gender_differences_data.csv"),
          row.names = TRUE)

# 保存图形
png(file.path(gender_analysis_dir, "sex_specific_heatmap.png"), 
    width = 800, height = 900, res = 150)
print(sex_specific_heatmap)
dev.off()

png(file.path(gender_analysis_dir, "sex_specific_heatmap_by_category.png"), 
    width = 800, height = 900, res = 150)
print(sex_specific_heatmap_ordered)
dev.off()

png(file.path(gender_analysis_dir, "gender_differences_heatmap.png"), 
    width = 600, height = 900, res = 150)
print(gender_diff_heatmap)
dev.off()

ggsave(file.path(gender_analysis_dir, "gender_comparison_plot.png"), 
       gender_comparison_plot, width = 8, height = 10, dpi = 300)

# 保存分析报告
sink(file.path(gender_analysis_dir, "gender_analysis_report.txt"))
cat("=== 性别特异性分析报告 ===\n")
cat("分析日期:", date(), "\n\n")

cat("数据概述:\n")
cat(sprintf("• 分析的基因数: %d\n", nrow(sex_specific_fc)))
cat(sprintf("• 基因类别: %s\n", paste(unique(sex_specific_fc$Category), collapse = ", ")))
cat("\n")

cat("性别差异总结:\n")
print(gender_diff)

cat("\n\n性别特异性变化模式:\n")
# 按基因显示
for(i in 1:nrow(sex_specific_fc)) {
  gene_info <- sex_specific_fc[i, ]
  cat(sprintf("%s (%s):\n", gene_info$Gene, gene_info$Category))
  cat(sprintf("  女性: log2FC = %.3f\n", gene_info$Female_FC))
  cat(sprintf("  男性: log2FC = %.3f\n", gene_info$Male_FC))
  cat(sprintf("  性别差异(男-女): %.3f\n\n", gene_info$Male_FC - gene_info$Female_FC))
}

cat("\n关键发现:\n")
cat("1. 性别差异模式:\n")
cat("   • 大多数基因在男性和女性中显示相似的变化方向\n")
cat("   • 但变化幅度存在性别差异\n")
cat("2. 类别特异性模式:\n")
cat("   • 突触前基因在男性和女性中都普遍下调\n")
cat("   • 突触后基因在男性和女性中都普遍上调\n")
cat("3. 最大的性别差异见于:\n")
for(i in 1:nrow(top_gender_diff)) {
  cat(sprintf("   • %s (差异 = %.3f)\n", 
              rownames(top_gender_diff)[i],
              top_gender_diff$Gender_Difference[i]))
}

cat("\n生物学意义:\n")
cat("这些性别特异性差异可能反映:\n")
cat("• 男性和女性AD患者中突触病理的不同机制\n")
cat("• 性激素对突触基因表达的调节作用\n")
cat("• 解释AD患病率和进展的性别差异\n")

cat("\n生成的图形文件:\n")
cat("1. sex_specific_heatmap.png - 性别特异性热图\n")
cat("2. sex_specific_heatmap_by_category.png - 按类别排序的热图\n")
cat("3. gender_differences_heatmap.png - 性别差异热图\n")
cat("4. gender_comparison_plot.png - 性别比较点图\n")

sink()

cat("\n=== 性别特异性分析完成 ===\n")
cat("所有结果已保存到 '", gender_analysis_dir, "' 目录\n")
cat("主要输出:\n")
cat("1. sex_specific_fc_data.csv - 性别特异性log2FC数据\n")
cat("2. gender_differences_data.csv - 性别差异数据\n")
cat("3. 四种不同的热图和比较图\n")
cat("4. 详细的分析报告\n")