# ===== 1. 数据统一修复 =====
cat("\n===== 执行数据统一修复 =====")

# 获取对象细胞标识符（最可靠来源）
cell_ids <- colnames(mzl)

# 修复元数据行名
rownames(mzl@meta.data) <- cell_ids
cat("\n✅ 元数据行名已统一")

# 修复RNA检测数据的列名
colnames(mzl@assays$RNA@data) <- cell_ids
cat("\n✅ RNA检测数据列名已统一")

# 修复RNA计数数据（如果存在）
if (!is.null(mzl@assays$RNA@counts)) {
  colnames(mzl@assays$RNA@counts) <- cell_ids
  cat("\n✅ RNA计数数据列名已统一")
}

# 修复降维数据（如果存在）
if (length(mzl@reductions) > 0) {
  for (red_name in names(mzl@reductions)) {
    colnames(mzl@reductions[[red_name]]@cell.embeddings) <- NULL
    rownames(mzl@reductions[[red_name]]@cell.embeddings) <- cell_ids
  }
  cat("\n✅ 降维数据行名已统一")
}

# ===== 2. 简化但高效的手动通信分析 =====
compute_communication_robust <- function(grp) {
  cat("\n===== 计算分组:", grp, "=====")
  
  # 1. 获取分组细胞索引
  grp_cells <- which(mzl$analysis_group == grp)
  
  if (length(grp_cells) == 0) {
    cat("\n❌ 分组中没有细胞:", grp)
    return(NULL)
  }
  
  # 2. 提取表达数据（直接使用矩阵格式）
  expr_data <- as.matrix(mzl@assays$RNA@data[, grp_cells, drop = FALSE])
  
  # 3. 获取细胞类型标识
  cell_types <- mzl@meta.data$cell_type_ident[grp_cells]
  
  # 4. 计算每个细胞类型的平均表达
  type_avg <- matrix(0, nrow = nrow(expr_data), ncol = length(unique(cell_types)),
                     dimnames = list(rownames(expr_data), unique(cell_types)))
  
  for (ct in colnames(type_avg)) {
    ct_cells <- expr_data[, cell_types == ct, drop = FALSE]
    type_avg[, ct] <- rowMeans(ct_cells, na.rm = TRUE)
  }
  
  # 5. 加载配体-受体数据库
  lr_db <- CellChatDB.human$interaction
  lr_pairs <- unique(lr_db[, c("ligand", "receptor")])
  
  # 6. 创建通信矩阵
  comm_matrix <- matrix(0, nrow = ncol(type_avg), ncol = ncol(type_avg),
                        dimnames = list(colnames(type_avg), colnames(type_avg)))
  
  # 7. 计算通信强度
  for (i in 1:nrow(lr_pairs)) {
    ligand <- lr_pairs[i, "ligand"]
    receptor <- lr_pairs[i, "receptor"]
    
    if (ligand %in% rownames(type_avg) && receptor %in% rownames(type_avg)) {
      # 通信强度 = 配体表达 × 受体表达
      comm_strength <- outer(
        type_avg[ligand, ], 
        type_avg[receptor, ],
        FUN = "*"
      )
      comm_matrix <- comm_matrix + comm_strength
    }
  }
  
  # 8. 保存结果
  result <- list(
    group = grp,
    communication_matrix = comm_matrix,
    cell_types = colnames(type_avg),
    timestamp = Sys.time()
  )
  
  saveRDS(result, paste0("CommResult_", grp, ".rds"))
  
  # 9. 生成热图
  pdf(paste0("CommHeatmap_", grp, ".pdf"), width = 10, height = 8)
  heatmap(comm_matrix, 
          main = paste(grp, "Cell-Cell Communication"),
          col = colorRampPalette(c("blue", "white", "red"))(100),
          margins = c(10, 10),
          cexRow = 0.8, cexCol = 0.8)
  dev.off()
  
  cat("\n✅ 分组计算完成:", grp)
  return(result)
}

# ===== 3. 执行分析 =====
# 分析所有分组
groups <- c("male_AD", "male_ctl", "female_AD", "female_ctl")
results <- list()

for (grp in groups) {
  results[[grp]] <- compute_communication_robust(grp)
}

# ===== 4. 合并结果 =====
# 创建合并矩阵
merged_matrix <- matrix(0, 
                        nrow = nrow(results[[1]]$communication_matrix),
                        ncol = length(groups),
                        dimnames = list(rownames(results[[1]]$communication_matrix), groups))

for (grp in groups) {
  if (!is.null(results[[grp]])) {
    # 取矩阵的对角线和作为该组的总体通信强度
    diag_sum <- diag(results[[grp]]$communication_matrix)
    merged_matrix[, grp] <- diag_sum
  }
}

# 保存合并结果
saveRDS(merged_matrix, "Combined_Communication_Results.rds")

# 生成比较热图
pdf("Combined_Communication.pdf", width = 12, height = 9)
heatmap(merged_matrix,
        main = "Cell-Cell Communication Comparison",
        xlab = "Groups",
        ylab = "Cell Types",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        margins = c(10, 8),
        cexRow = 0.7, cexCol = 0.8)
dev.off()

# ===== 5. 最终报告 =====
cat("\n\n===== 分析完成! =====")
cat("\n生成结果:")
cat("\n- 各分组通信矩阵: CommResult_*.rds")
cat("\n- 各分组热图: CommHeatmap_*.pdf")
cat("\n- 合并结果: Combined_Communication_Results.rds")
cat("\n- 比较热图: Combined_Communication.pdf")
cat("\n\n请在当前目录检查生成的文件!")

