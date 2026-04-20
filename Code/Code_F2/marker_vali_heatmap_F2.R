library(Seurat)
library(pheatmap)
library(RColorBrewer)

# 1. 定义细胞类型顺序和验证基因
cell_type_order <- c("ex.neu", "in.neu", "ast", "mic", "oli", "opc")
marker_genes <- list(
  ex.neu = c("SLC17A7", "CAMK2A"),
  in.neu = c("GAD1", "GAD2"),
  ast = c("GFAP", "AQP4"),
  mic = c("CSF1R", "CD68"),
  oli = c("MBP", "PLP1"),
  opc = c("PDGFRA", "SOX10")
)

# 2. 按定义顺序重组基因名
gene_order <- unlist(marker_genes[cell_type_order])
gene_labels <- rep(names(marker_genes), sapply(marker_genes, length))

# 3. 手动计算平均表达
expr_matrix <- GetAssayData(mzl, assay = "RNA", slot = "data")
cell_types <- mzl$cell_type_ident

# 初始化平均表达矩阵
avg_expr <- matrix(0, nrow = length(gene_order), ncol = length(cell_type_order))
rownames(avg_expr) <- gene_order
colnames(avg_expr) <- cell_type_order

# 填充矩阵
for (gene in gene_order) {
  for (ct in cell_type_order) {
    cells_in_ct <- which(cell_types == ct)
    avg_expr[gene, ct] <- mean(expr_matrix[gene, cells_in_ct])
  }
}

# 4. Z-score标准化
z_score <- t(scale(t(avg_expr)))
z_score[is.na(z_score)] <- 0

# 5. 创建简洁配色方案
## 细胞类型配色（Set2调色板）
cell_type_colors <- brewer.pal(6, "Set2")
names(cell_type_colors) <- cell_type_order

## 热图配色（经典蓝红渐变）
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)  # 蓝-白-红

# 6. 创建注释数据框
annotation_row <- data.frame(
  CellType = factor(gene_labels, levels = cell_type_order),
  row.names = gene_order
)

# 7. 简洁热图绘制
pheatmap(
  mat = z_score,
  color = heatmap_colors,  # 蓝红配色
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_row = cumsum(sapply(marker_genes, length))[-length(marker_genes)],  # 行分组线
  annotation_row = annotation_row,  # 行注释
  annotation_colors = list(CellType = cell_type_colors),  # 注释颜色
  main = "Cell Type Marker Validation",
  fontsize = 10,
  angle_col = 45,  # 列名倾斜
  border_color = "grey90",  # 浅灰色边框
  cellwidth = 20,
  cellheight = 12,
  legend = TRUE
)

# 8. 保存高清图
ggsave("clean_validation_heatmap.png", 
       width = 9, 
       height = 7, 
       dpi = 300)
