# 1. 检查Seurat对象的基本信息
cat("=== Seurat对象基本信息 ===\n")
cat("总基因数:", nrow(mzl), "\n")
cat("总细胞数:", ncol(mzl), "\n")
cat("已使用的内存:", format(object.size(mzl), units = "MB"), "\n")

# 2. 检查基因过滤情况
cat("\n=== 基因过滤情况 ===\n")
cat("原始基因数:", nrow(mzl@assays$RNA@counts), "\n")
cat("当前基因数:", nrow(mzl), "\n")
cat("过滤掉的基因比例:", 
    round((1 - nrow(mzl)/nrow(mzl@assays$RNA@counts)) * 100, 1), "%\n")

# 3. 检查基因表达情况
gene_stats <- data.frame(
  gene = rownames(mzl),
  mean_expression = rowMeans(mzl@assays$RNA@data),
  expressed_cells = rowSums(mzl@assays$RNA@data > 0),
  total_cells = ncol(mzl)
) %>% arrange(desc(mean_expression))

# 显示表达量最高的10个基因
cat("\n表达量最高的10个基因:\n")
print(head(gene_stats, 10))

# 显示表达量最低的10个基因
cat("\n表达量最低的10个基因:\n")
print(tail(gene_stats, 10))

# 4. 检查细胞类型中的基因表达情况
cell_types <- unique(mzl$cell_type_ident)
cell_gene_stats <- list()

for (ct in cell_types) {
  # 获取该细胞类型的细胞
  ct_cells <- colnames(mzl)[mzl$cell_type_ident == ct]
  
  # 计算该细胞类型中的基因表达
  ct_data <- mzl@assays$RNA@data[, ct_cells]
  
  cell_gene_stats[[ct]] <- data.frame(
    cell_type = ct,
    gene = rownames(ct_data),
    mean_expression = rowMeans(ct_data),
    expressed_cells = rowSums(ct_data > 0),
    total_cells = length(ct_cells)
  )
}

# 合并所有细胞类型的数据
all_cell_gene_stats <- do.call(rbind, cell_gene_stats)

# 5. 可视化基因表达分布
ggplot(gene_stats, aes(x = mean_expression)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  scale_x_log10() +
  labs(title = "基因平均表达量分布",
       x = "平均表达量(log10尺度)",
       y = "基因数量") +
  theme_minimal()

# 6. 检查差异表达分析的输入数据
cat("\n=== 差异表达分析输入检查 ===\n")
for (ct in cell_types) {
  ct_cells <- colnames(mzl)[mzl$cell_type_ident == ct]
  ad_cells <- ct_cells[mzl$diagnosis[ct_cells] == "AD"]
  ctl_cells <- ct_cells[mzl$diagnosis[ct_cells] == "ctl"]
  
  cat("细胞类型:", ct, "\n")
  cat("  总细胞数:", length(ct_cells), "\n")
  cat("  AD组细胞数:", length(ad_cells), "\n")
  cat("  对照组细胞数:", length(ctl_cells), "\n")
  
  if (length(ad_cells) < 3 || length(ctl_cells) < 3) {
    cat("  警告: 某一组的细胞数不足3个!\n")
  }
}

# 7. 检查差异表达分析函数
cat("\n=== 差异表达分析函数检查 ===\n")
test_result <- run_de_analysis(mzl, cell_type = "ex.neu", sex_group = "male")
cat("测试分析返回的基因数:", nrow(test_result), "\n")

# 重新创建对象
mzl_new <- CreateSeuratObject(
  counts = original_counts,
  min.cells = 3,
  min.features = 200,
  meta.data = original_metadata  # 直接添加元数据
)

# 预处理
mzl_new <- NormalizeData(mzl_new)
mzl_new <- FindVariableFeatures(mzl_new, nfeatures = 3000)
mzl_new <- ScaleData(mzl_new)
mzl_new <- RunPCA(mzl_new)

# 验证元数据存在
cat("元数据列:", colnames(mzl_new@meta.data), "\n")

# 运行差异分析
test_result <- run_de_analysis(mzl_new, cell_type = "ex.neu", sex_group = "male")

