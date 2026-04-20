library(pheatmap)
library(RColorBrewer)

# 准备数据
deg_genes <- synaptic_degs$gene

# 确保使用正确的基因名
available_genes <- deg_genes[deg_genes %in% rownames(ex_neu)]
cat("可用的基因数量:", length(available_genes), "\n")
print(available_genes)

# 获取表达矩阵（使用log-normalized数据）
expr_data <- GetAssayData(ex_neu, slot = "data")[available_genes, ]

# 创建元数据用于排序
metadata <- ex_neu@meta.data
metadata$cell_id <- rownames(metadata)

# 按性别和诊断排序
metadata_sorted <- metadata %>%
  arrange(sex, diagnosis)  # 先按性别，再按诊断排序

# 按排序后的顺序重排列
expr_sorted <- expr_data[, metadata_sorted$cell_id]

# 创建注释信息
annotation_col <- data.frame(
  Sex = metadata_sorted$sex,
  Diagnosis = metadata_sorted$diagnosis,
  row.names = metadata_sorted$cell_id
)

# 设置颜色
sex_colors <- c("female" = "#F781BF", "male" = "#4DAF4A")
diagnosis_colors <- c("AD" = "#E41A1C", "ctl" = "#377EB8")

ann_colors <- list(
  Sex = sex_colors,
  Diagnosis = diagnosis_colors
)

# 绘制分组热图（不聚类）
heatmap_grouped <- pheatmap(
  expr_sorted,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  cluster_cols = FALSE,  # 不进行列聚类
  cluster_rows = TRUE,   # 只对行聚类
  fontsize_row = 10,
  fontsize_col = 8,
  main = "Synaptic DEGs Expression - Grouped by Sex and Diagnosis",
  gaps_col = cumsum(table(metadata_sorted$sex)),  # 在性别之间添加分隔线
  border_color = NA,
  silent = FALSE
)

# 保存热图
png("synaptic_heatmap_grouped.png", width = 1400, height = 900, res = 150)
print(heatmap_grouped)
dev.off()

# 创建更清晰的分组热图 - 只显示部分细胞
set.seed(123)
# 每个组别抽取相同数量的细胞
sample_per_group <- 500  # 每个性别-诊断组合抽取500个细胞

# 获取每个组的细胞
group_cells <- list()
for(s in c("female", "male")) {
  for(d in c("AD", "ctl")) {
    group_idx <- which(metadata$sex == s & metadata$diagnosis == d)
    if(length(group_idx) > sample_per_group) {
      selected <- sample(group_idx, sample_per_group)
    } else {
      selected <- group_idx
    }
    group_cells[[paste(s, d, sep = "_")]] <- metadata$cell_id[selected]
  }
}

# 合并所有选择的细胞
selected_cells <- unlist(group_cells)

# 按组别排序
selected_metadata <- metadata[selected_cells, ]
selected_metadata <- selected_metadata %>%
  arrange(sex, diagnosis)

# 获取表达数据
expr_selected <- expr_data[available_genes, selected_metadata$cell_id]

# 创建注释
annotation_selected <- data.frame(
  Sex = selected_metadata$sex,
  Diagnosis = selected_metadata$diagnosis,
  row.names = selected_metadata$cell_id
)

# 绘制更清晰的热图
heatmap_clear <- pheatmap(
  expr_selected,
  annotation_col = annotation_selected,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  cluster_cols = FALSE,  # 不聚类
  cluster_rows = TRUE,
  fontsize_row = 10,
  fontsize = 8,
  main = "Synaptic DEGs Expression (Sampled Cells)",
  gaps_col = cumsum(table(selected_metadata$sex)),  # 在性别之间添加分隔线
  border_color = NA,
  cellwidth = NA,
  cellheight = 15
)

# 保存清晰热图
png("synaptic_heatmap_clear.png", width = 1200, height = 800, res = 150)
print(heatmap_clear)
dev.off()