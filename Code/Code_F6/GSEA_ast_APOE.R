# Part 1.1: GSEA on genes correlated with APOE in Female AD Astrocytes (修正版)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(enrichplot)

# 1. 准备女性AD星形胶质细胞的表达矩阵
female_ad_ast_cells <- colnames(mzl)[mzl$cell_type_ident == "ast" & 
                                       mzl$sex == "female" & 
                                       mzl$diagnosis == "AD"]
cat("Number of female AD astrocyte cells:", length(female_ad_ast_cells), "\n")

expr_matrix <- GetAssayData(mzl, slot = "data")[, female_ad_ast_cells]

# 2. 计算所有基因与APOE的Spearman相关系数
apoe_expr <- as.numeric(expr_matrix["APOE", ])

# 只计算在足够多细胞中表达的基因（减少计算量和噪音）
min_cells <- ceiling(0.1 * length(female_ad_ast_cells))
genes_to_keep <- rownames(expr_matrix)[rowSums(expr_matrix > 0) >= min_cells]
expr_matrix <- expr_matrix[genes_to_keep, ]

cat("Analyzing", nrow(expr_matrix), "genes expressed in ≥10% of cells\n")

# 3. 构建排序基因列表 - 使用更稳健的方法
cor_df <- data.frame(
  gene = character(0),
  rho = numeric(0),
  pval = numeric(0),
  stringsAsFactors = FALSE
)

# 逐基因计算相关性，避免apply中的警告
pb <- txtProgressBar(min = 0, max = nrow(expr_matrix), style = 3)
for(i in 1:nrow(expr_matrix)) {
  gene_expr <- as.numeric(expr_matrix[i, ])
  
  # 跳过表达量恒定的基因
  if(sd(gene_expr) == 0) next
  
  # 计算Spearman相关性
  cor_test <- tryCatch({
    cor.test(gene_expr, apoe_expr, method = "spearman", exact = FALSE)  # exact = FALSE避免大量数据时的计算问题
  }, error = function(e) {
    return(NULL)
  })
  
  if(!is.null(cor_test)) {
    cor_df <- rbind(cor_df, data.frame(
      gene = rownames(expr_matrix)[i],
      rho = cor_test$estimate,
      pval = cor_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
  
  setTxtProgressBar(pb, i)
}
close(pb)

cat("Successfully calculated correlations for", nrow(cor_df), "genes\n")

# 检查数据框结构
cat("Column names in cor_df:", colnames(cor_df), "\n")
cat("First few rows of cor_df:\n")
print(head(cor_df))

# 按相关系数从大到小排序（递减顺序！）
cor_df <- cor_df %>% arrange(desc(rho))

# 创建排序基因列表：必须是递减排序的命名向量
ranked_gene_list <- setNames(cor_df$rho, cor_df$gene)

# 检查排序是否正确
cat("\nGene list is sorted (decreasing):", is.unsorted(rev(ranked_gene_list)), "\n")
cat("First 5 genes (highest positive correlation with APOE):\n")
print(head(ranked_gene_list, 5))
cat("\nLast 5 genes (highest negative correlation with APOE):\n")
print(tail(ranked_gene_list, 5))

# 检查是否有NA或Inf值
cat("\nNumber of NA values in ranked_gene_list:", sum(is.na(ranked_gene_list)), "\n")
cat("Number of Inf/-Inf values:", sum(is.infinite(ranked_gene_list)), "\n")

# 移除任何NA或Inf值
ranked_gene_list <- ranked_gene_list[!is.na(ranked_gene_list) & !is.infinite(ranked_gene_list)]
cat("After removing NA/Inf, remaining genes:", length(ranked_gene_list), "\n")

# 4. 获取基因集 - 使用更相关的基因集
# 可以尝试多种基因集，这里先使用Hallmark
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

# 选择与星形胶质细胞功能和AD相关的通路
selected_pathways <- c(
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_FATTY_ACID_METABOLISM", 
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)
hallmark_subset <- hallmark_sets %>% filter(gs_name %in% selected_pathways)

# 5. 运行GSEA（使用更宽松的参数，先查看趋势）
cat("\nRunning GSEA...\n")
set.seed(123)

# 先检查基因列表是否有效
if(length(ranked_gene_list) < 100) {
  cat("Warning: Gene list is too short for meaningful GSEA\n")
} else {
  # 运行GSEA
  gsea_result <- tryCatch({
    GSEA(
      geneList = ranked_gene_list,
      exponent = 1,
      minGSSize = 10,      # 最小基因集大小
      maxGSSize = 500,     # 最大基因集大小
      pvalueCutoff = 0.2,  # 放宽p值阈值查看趋势
      pAdjustMethod = "BH",
      TERM2GENE = hallmark_subset,
      seed = TRUE,
      verbose = TRUE       # 显示进度
    )
  }, error = function(e) {
    cat("GSEA error:", e$message, "\n")
    return(NULL)
  })
  
  # 6. 可视化结果
  if (!is.null(gsea_result) && nrow(gsea_result@result) > 0) {
    # 显示结果摘要
    cat("\nGSEA Results Summary:\n")
    print(head(gsea_result@result, 10))
    
    # 绘制GSEA点图
    p1 <- dotplot(gsea_result, 
                  showCategory = 15, 
                  color = "p.adjust",
                  title = "GSEA: Pathways correlated with APOE expression\n(Female AD Astrocytes)",
                  font.size = 10)
    print(p1)
    
    # 保存图片
    ggsave("./DEG_Results/GSEA_Dotplot.png", p1, width = 10, height = 7, dpi = 300)
    
    # 如果有显著富集的通路，绘制单个通路的富集图
    sig_pathways <- gsea_result@result %>% 
      filter(p.adjust < 0.2) %>% 
      pull(ID)
    
    if (length(sig_pathways) > 0) {
      cat("\nSignificantly enriched pathways (p.adjust < 0.2):\n")
      print(sig_pathways)
      
      # 绘制第一个显著通路的富集图
      p2 <- gseaplot2(gsea_result, 
                      geneSetID = sig_pathways[1],
                      title = paste("Enrichment Plot:", sig_pathways[1]),
                      color = "firebrick",
                      base_size = 11)
      print(p2)
      ggsave("./DEG_Results/GSEA_Enrichment_Plot.png", p2, width = 8, height = 6, dpi = 300)
    }
    
    # 保存GSEA结果表
    write.csv(gsea_result@result, 
              file = "./DEG_Results/GSEA_APOE_Female_AD_Astrocytes.csv", 
              row.names = FALSE)
    cat("\nGSEA results saved to: ./DEG_Results/GSEA_APOE_Female_AD_Astrocytes.csv\n")
    
  } else {
    cat("GSEA did not find any significantly enriched pathways.\n")
    
    # 尝试使用更广泛的基因集
    cat("Trying with all Hallmark gene sets...\n")
    
    gsea_result_all <- tryCatch({
      GSEA(
        geneList = ranked_gene_list,
        exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        pvalueCutoff = 0.3,  # 进一步放宽阈值
        TERM2GENE = hallmark_sets,
        seed = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      cat("GSEA with all Hallmark error:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(gsea_result_all) && nrow(gsea_result_all@result) > 0) {
      cat("Found", nrow(gsea_result_all@result), "potentially enriched pathways\n")
      print(head(gsea_result_all@result[order(gsea_result_all@result$p.adjust), ], 10))
    }
  }
}

# 7. 保存中间结果
saveRDS(cor_df, file = "./DEG_Results/APOE_Correlation_Results.rds")
saveRDS(ranked_gene_list, file = "./DEG_Results/Ranked_Gene_List.rds")
cat("\nIntermediate results saved.\n")

# 8. 绘制相关性分布图
p_cor_dist <- ggplot(cor_df, aes(x = rho)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Spearman Correlation with APOE",
       subtitle = paste("Female AD Astrocytes (n =", nrow(cor_df), "genes)"),
       x = "Spearman Correlation Coefficient (ρ)",
       y = "Number of Genes") +
  theme_minimal()

print(p_cor_dist)
ggsave("./DEG_Results/APOE_Correlation_Distribution.png", p_cor_dist, width = 8, height = 6, dpi = 300)

# 如果有显著富集的通路，绘制多个通路的富集图
sig_pathways <- gsea_result@result %>% 
  filter(p.adjust < 0.2) %>% 
  arrange(p.adjust) %>%  # 按p值排序
  pull(ID)

if (length(sig_pathways) > 0) {
  cat("\n显著富集通路 (p.adjust < 0.2):\n")
  print(sig_pathways)
  
  # 选择要绘制的通路数量（最多6个，避免过于拥挤）
  n_to_plot <- min(6, length(sig_pathways))
  pathways_to_plot <- sig_pathways[1:n_to_plot]
  
  # 绘制多个通路的富集图
  cat(paste0("\n绘制前", n_to_plot, "个最显著的通路...\n"))
  
  # 方法2：为每个通路单独绘制并保存
  cat("为每个显著通路单独绘制富集图...\n")
  for (i in 1:length(sig_pathways)) {
    pathway <- sig_pathways[i]
    
    # 为每个通路创建单独的图
    single_plot <- gseaplot2(gsea_result, 
                             geneSetID = pathway,
                             title = paste("Enrichment Plot:", pathway),
                             color = "firebrick",
                             base_size = 11)
    
    # 保存单个图
    file_name <- paste0("./DEG_Results/GSEA_", 
                        gsub("HALLMARK_", "", pathway),  # 移除HALLMARK_前缀
                        "_Enrichment.png")
    ggsave(file_name, single_plot, width = 8, height = 6, dpi = 300)
    
    # 打印前3个通路
    if (i <= 3) {
      print(single_plot)
    }
  }
  
  cat(paste0("\n已保存", length(sig_pathways), "个显著通路的富集图到DEG_Results文件夹\n"))
  
  # 方法3：创建点图显示所有显著通路（替代可视化）
  if (length(sig_pathways) > 1) {
    sig_data <- gsea_result@result %>% 
      filter(p.adjust < 0.2) %>% 
      arrange(NES)  # 按NES排序
    
    # 创建点图
    dot_plot <- ggplot(sig_data, aes(x = NES, y = reorder(Description, NES))) +
      geom_point(aes(size = -log10(p.adjust), color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", 
                            midpoint = 0) +
      labs(title = "Significantly Enriched Pathways",
           x = "Normalized Enrichment Score (NES)",
           y = "Pathway",
           size = "-log10(adj. p-value)",
           color = "NES") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 10),
            plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(dot_plot)
    ggsave("./DEG_Results/GSEA_Significant_Pathways_DotPlot.png", 
           dot_plot, width = 10, height = 6 + length(sig_pathways)*0.3, dpi = 300)
  }
  
} else {
  cat("\n没有显著富集的通路 (p.adjust < 0.2)\n")
}