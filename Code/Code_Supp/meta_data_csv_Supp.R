# 提取并导出临床元数据
clinical_metadata <- mzl@meta.data[, c("orig.ident", "diagnosis", "sex", "age", "apoE", "braak", "cerad")]

# 去除重复行（每个样本保留唯一记录）
clinical_metadata_unique <- clinical_metadata[!duplicated(clinical_metadata$orig.ident), ]

# 导出为CSV文件
write.csv(clinical_metadata_unique, 
          file = "clinical_metadata.csv", 
          row.names = FALSE)

# 显示导出成功信息
message("临床元数据已成功导出为 clinical_metadata.csv")
print(paste("包含", nrow(clinical_metadata_unique), "个样本的临床数据"))

