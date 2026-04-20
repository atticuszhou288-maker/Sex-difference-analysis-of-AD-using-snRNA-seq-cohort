# 加载必要库
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# 1. 定义关键基因列表 - 按功能分组严格排序
key_genes_ordered <- c(
  "TREM2", "TYROBP", "C1QA", "APOE", "CLU", 
  "IL1B", "ESR1", "SNAP25", "SYT1", "DLG4"
)

# 2. 定义表达量阈值 - 低于此值视为不表达
expression_threshold <- 0.001  # 可根据数据调整

# 3. 数据准备与log2(AD/CTL)计算（添加表达量过滤）
log2fc_data <- stats_df %>%
  select(sex, cell_type_ident, diagnosis, ends_with("_mean")) %>%
  pivot_longer(
    cols = ends_with("_mean"),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(gene = str_remove(gene, "_mean")) %>%
  filter(gene %in% key_genes_ordered) %>%
  # 按组计算AD/CTL比值，添加表达量检查
  group_by(sex, cell_type_ident, gene) %>%
  summarize(
    AD_expression = expression[diagnosis == "AD"],
    CTL_expression = expression[diagnosis == "ctl"],
    # 检查是否表达量足够
    is_expressed = (AD_expression > expression_threshold) | (CTL_expression > expression_threshold),
    # 计算log2FC，对不表达的设为NA
    log2FC = if (is_expressed) {
      # 处理分母为0的情况
      if (CTL_expression == 0) {
        if (AD_expression == 0) 0 else 10  # 极大上调
      } else {
        log2(AD_expression/CTL_expression)
      }
    } else {
      NA_real_  # 不表达设为NA
    }
  ) %>%
  ungroup() %>%
  # 处理特殊值
  mutate(log2FC = case_when(
    is.infinite(log2FC) & log2FC > 0 ~ 3,   # 极大上调设为3
    is.infinite(log2FC) & log2FC < 0 ~ -3,  # 极大下调设为-3
    TRUE ~ log2FC
  )) %>%
  # 创建组合标签
  mutate(group = paste(sex, cell_type_ident, sep = " - ")) %>%
  # 应用严格的顺序控制
  mutate(
    sex = factor(sex, levels = c("female", "male")),
    cell_type_ident = factor(cell_type_ident, levels = c("ast", "mic", "oli", "opc", "ex.neu", "in.neu")),
    gene = factor(gene, levels = key_genes_ordered)
  ) %>%
  # 按预定顺序排列
  arrange(sex, cell_type_ident, gene)

# 4. 转换为热图矩阵（保持严格顺序）
heatmap_matrix <- log2fc_data %>%
  select(gene, group, log2FC) %>%
  pivot_wider(names_from = group, values_from = log2FC) %>%
  mutate(gene = factor(gene, levels = key_genes_ordered)) %>%
  arrange(gene) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# 5. 创建列顺序
col_order <- log2fc_data %>%
  distinct(sex, cell_type_ident, group) %>%
  arrange(sex, cell_type_ident) %>%
  pull(group)
heatmap_matrix <- heatmap_matrix[, col_order]

# 6. 创建热图 - 空白格子表示不表达
p <- pheatmap(
  heatmap_matrix,
  # ====== 颜色和断点 ======
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-1.5, 1.5, length.out = 101),
  na_col = "white",  # NA值显示为白色空白格子
  
  # ====== 严格禁用聚类 ======
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  # ====== 布局优化 ======
  cellwidth = 26,
  cellheight = 24,
  fontsize_row = 13,
  fontsize_col = 12,
  
  # ====== 标注优化 ======
  display_numbers = matrix(ifelse(!is.na(heatmap_matrix) & abs(heatmap_matrix) > 0.3, 
                                  sprintf("%.1f", heatmap_matrix), ""),
                           nrow = nrow(heatmap_matrix)),
  number_color = "black",
  fontsize_number = 11,
  
  # ====== 图例优化 ======
  legend = TRUE,
  legend_breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
  legend_labels = c("-1.5", "-1", "-0.5", "0", "0.5", "1", "1.5"),
  
  # ====== 分组注释 ======
  annotation_col = log2fc_data %>%
    select(group, sex, cell_type_ident) %>%
    distinct() %>%
    column_to_rownames("group"),
  annotation_colors = list(
    sex = c(female = "#F8766D", male = "#00BFC4"),
    cell_type_ident = c(ast = "#1B9E77", mic = "#D95F02", oli = "#7570B3", 
                        opc = "#E7298A", ex.neu = "#66A61E", in.neu = "#E6AB02")
  ),
  
  # ====== 标题和边框 ======
  main = "Log2(AD/Control) Expression Change of Key Genes\n(Blank cells indicate low/no expression)",
  border_color = "grey70",
  gaps_col = c(6)
)

# 7. 保存最终热图
ggsave("Final_AD_Expression_Heatmap_With_NA.png", p, 
       width = 20, height = 14, dpi = 300)

# 8. 添加图例说明 - 修复版本
# 创建图例说明文本
legend_text <- textGrob(
  "Blank cells indicate genes with expression < 0.001 in both AD and Control groups",
  gp = gpar(fontsize = 12, col = "grey40"),
  x = 0.5, y = 0.5, 
  just = "centre"
)

# 创建空白绘图区域作为图例容器
legend_plot <- grid.rect(gp = gpar(col = NA, fill = NA))  # 透明背景

# 将图例添加到热图下方
# 创建整体布局
final_layout <- grid.arrange(
  p$gtable,             # 热图部分
  legend_plot,          # 空白区域（用于放置图例）
  legend_text,          # 图例文本
  nrow = 3,             # 垂直排列
  heights = c(0.90, 0.02, 0.08)  # 各组件高度比例
)

# 保存带图例说明的图像
ggsave("Final_AD_Expression_Heatmap_With_Legend.png", 
       plot = final_layout, 
       width = 20, height = 16, dpi = 300)  # 增加高度以适应图例