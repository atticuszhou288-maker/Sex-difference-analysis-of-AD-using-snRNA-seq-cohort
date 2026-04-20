library(Seurat)
library(dplyr)
library(tidyr)

# 计算每个样本中每种细胞类型的细胞数量
cell_type_counts <- mzl@meta.data %>%
  filter(diagnosis %in% c("AD", "ctl") & 
           sex %in% c("female", "male") &
           cell_type_ident %in% c("ast", "ex.neu", "in.neu", "mic", "oli", "opc")) %>%
  count(orig.ident, diagnosis, sex, cell_type_ident, name = "cell_count") %>%
  arrange(orig.ident, cell_type_ident)

# 转换为宽格式以便查看
cell_type_counts_wide <- cell_type_counts %>%
  pivot_wider(
    names_from = cell_type_ident,
    values_from = cell_count,
    values_fill = 0
  ) %>%
  # 计算每个样本的总细胞数
  mutate(total_cells = rowSums(across(c(ast, ex.neu, in.neu, mic, oli, opc)))) %>%
  select(orig.ident, diagnosis, sex, total_cells, everything())

# 添加百分比列
cell_type_counts_wide <- cell_type_counts_wide %>%
  mutate(
    ast_pct = round(ast / total_cells * 100, 1),
    ex.neu_pct = round(ex.neu / total_cells * 100, 1),
    in.neu_pct = round(in.neu / total_cells * 100, 1),
    mic_pct = round(mic / total_cells * 100, 1),
    oli_pct = round(oli / total_cells * 100, 1),
    opc_pct = round(opc / total_cells * 100, 1)
  )

# 输出结果（显示前10个样本）
cat("=== SAMPLE-LEVEL CELL TYPE COUNTS (FIRST 10 SAMPLES) ===\n")
print(cell_type_counts_wide, n=75)

# 保存完整结果到CSV
write.csv(cell_type_counts_wide, "cell_type_counts_per_sample.csv", row.names = FALSE)

# 计算各细胞类型在不同分组中的平均比例
avg_proportions <- cell_type_counts_wide %>%
  group_by(diagnosis, sex) %>%
  summarise(
    n_samples = n(),
    avg_ast = mean(ast_pct),
    avg_ex.neu = mean(ex.neu_pct),
    avg_in.neu = mean(in.neu_pct),
    avg_mic = mean(mic_pct),
    avg_oli = mean(oli_pct),
    avg_opc = mean(opc_pct),
    .groups = "drop"
  )

cat("\n=== AVERAGE CELL PROPORTIONS BY DIAGNOSIS AND SEX ===\n")
print(avg_proportions)
