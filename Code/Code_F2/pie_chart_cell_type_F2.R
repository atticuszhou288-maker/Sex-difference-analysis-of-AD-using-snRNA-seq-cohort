library(Seurat)
library(ggplot2)
library(dplyr)

# 1. 确保元数据列存在
if (!"cell_type_ident" %in% colnames(mzl@meta.data)) stop("缺少'cell_type_ident'列")
if (!all(c("diagnosis", "sex") %in% colnames(mzl@meta.data))) stop("缺少'diagnosis'或'sex'列")

# 2. 创建分组标签（AD/ctl + female/male）
mzl$group <- interaction(mzl$diagnosis, mzl$sex, sep = "_")

# 3. 计算比例数据（仅保留所需分组）
prop_data <- mzl@meta.data %>% 
  filter(group %in% c("AD_female", "AD_male", "ctl_female", "ctl_male")) %>%
  count(group, cell_type_ident) %>%
  group_by(group) %>% 
  mutate(percent = n / sum(n) * 100) %>% 
  ungroup()

# 4. 标准神经细胞配色
cell_colors <- c(
  "ex.neu" = "#1F77B4",  # 蓝色
  "in.neu" = "#FF7F0E",  # 橙色
  "ast"    = "#2CA02C",  # 绿色
  "mic"    = "#D62728",  # 红色
  "oli"    = "#9467BD",  # 紫色
  "opc"    = "#8C564B"   # 棕色
)

# 5. 创建饼图函数（简化可靠版）
create_pie <- function(group_name) {
  group_data <- filter(prop_data, group == group_name)
  
  ggplot(group_data, aes(x = "", y = percent, fill = cell_type_ident)) +
    geom_col(color = "white", linewidth = 0.3, width = 1) +
    coord_polar("y") +
    geom_label(
      aes(label = ifelse(percent > 1, sprintf("%.1f%%", percent), "")),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 3.5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cell_colors) +
    labs(title = group_name) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "none"
    )
}

# 6. 生成四个饼图（确保顺序正确）
pie_list <- list(
  "AD_female" = create_pie("AD_female"),
  "AD_male"   = create_pie("AD_male"),
  "ctl_female" = create_pie("ctl_female"),
  "ctl_male"  = create_pie("ctl_male")
)

# 7. 创建共享图例（独立图形）
legend_data <- data.frame(
  cell_type = names(cell_colors),
  color = cell_colors
)

legend_plot <- ggplot(legend_data, aes(x = cell_type, y = 1, fill = cell_type)) +
  geom_tile() +
  scale_fill_manual(values = cell_colors, name = "Cell Type") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.8, "lines")
  )

# 8. 保存所有图形（三个文件）
## 8.1 保存四个独立饼图
for (name in names(pie_list)) {
  ggsave(
    paste0("pie_", name, ".png"), 
    pie_list[[name]],
    width = 5, 
    height = 5,
    dpi = 300
  )
}

## 8.2 保存图例
ggsave("cell_legend.png", legend_plot, width = 6, height = 1, dpi = 300)

## 8.3 保存组合PDF（2×2 + 图例）
pdf("cell_proportions.pdf", width = 10, height = 10)
layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE), 
       heights = c(1,1,0.2))
for (i in 1:4) {
  print(pie_list[[i]])
}
print(legend_plot)
dev.off()


