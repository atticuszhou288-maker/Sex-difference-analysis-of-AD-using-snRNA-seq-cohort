cat("测试运行'Overall'细胞类型...\n")
test_de <- run_de_analysis(mzl_new, cell_type = NULL, sex_group = "male")

if (nrow(test_de) > 0) {
  cat("成功运行差异分析! 发现", nrow(test_de), "个基因\n")
  print(head(test_de))
} else {
  cat("警告: 未发现差异基因\n")
}

if (nrow(test_de) > 0) {
  cell_types <- c("Overall", "ex.neu", "in.neu", "ast", "oli", "mic", "opc")
  all_de_genes <- data.frame()
  
  for (ct in cell_types) {
    cat("\n处理细胞类型:", ct, "\n")
    
    # 男性样本
    cat(" - 男性样本...")
    de_male <- run_de_analysis(mzl_new, 
                               cell_type = if(ct != "Overall") ct else NULL,
                               sex_group = "male")
    cat("完成 (", nrow(de_male), "个基因)\n")
    
    # 女性样本
    cat(" - 女性样本...")
    de_female <- run_de_analysis(mzl_new, 
                                 cell_type = if(ct != "Overall") ct else NULL,
                                 sex_group = "female")
    cat("完成 (", nrow(de_female), "个基因)\n")
    
    # 添加元数据
    if (nrow(de_male) > 0) {
      de_male$CellType <- ct
      de_male$Sex <- "male"
      all_de_genes <- bind_rows(all_de_genes, de_male)
    }
    
    if (nrow(de_female) > 0) {
      de_female$CellType <- ct
      de_female$Sex <- "female"
      all_de_genes <- bind_rows(all_de_genes, de_female)
    }
  }

  saveRDS(all_de_genes, "all_de_genes_results.rds")
  cat("\n差异分析完成! 保存结果到 all_de_genes_results.rds\n")
  
  sig_genes <- all_de_genes %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
    arrange(CellType, Sex, desc(avg_log2FC))
  
  cat("发现显著差异基因数量:", nrow(sig_genes), "\n")
}

# 常用阈值标准 (可根据需求调整)
p_val_adj_threshold <- 0.05     # 调整后p值
log2FC_threshold <- 0.58        # 约1.5倍变化 (|log2FC|>0.58)
library(dplyr)

significant_degs <- all_de_gene_results %>%
  filter(
    p_val_adj < p_val_adj_threshold,  # 显著性过滤
    abs(avg_log2FC) > log2FC_threshold # 表达量变化过滤
  )
significant_degs <- significant_degs %>%
  mutate(Direction = ifelse(avg_log2FC > 0, "Up", "Down"))
# 查看概况
cat("显著DEG数量:", nrow(significant_degs), "\n")
print(table(significant_degs$CellType, significant_degs$Sex))

# 检查示例行
head(significant_degs[, c("gene", "CellType", "Sex", "avg_log2FC", "p_val_adj")])
# 保存为CSV (推荐)
write.csv(significant_degs, "significant_degs.csv", row.names = FALSE)

# 保存为RDS (保留数据类型)
saveRDS(significant_degs, "significant_degs.rds")
