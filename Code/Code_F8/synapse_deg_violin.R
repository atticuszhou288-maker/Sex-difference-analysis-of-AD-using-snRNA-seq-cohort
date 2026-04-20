# 方法1：使用Seurat的VlnPlot，但修正颜色映射
# 创建分组变量
ex_neu$group <- paste(ex_neu$sex, ex_neu$diagnosis, sep = ".")

# 为每个组分配颜色
group_colors <- c(
  "female.AD" = "#F781BF",
  "female.ctl" = "#C6DBEF",
  "male.AD" = "#4DAF4A",
  "male.ctl" = "#C7E9C0"
)

# 使用Seurat VlnPlot，手动设置颜色
violin_plot_seurat <- VlnPlot(
  ex_neu,
  features = available_genes,
  group.by = "group",  # 按组分组
  cols = group_colors, # 手动指定颜色
  pt.size = 0,
  ncol = 4,
  combine = TRUE
) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    title = "Synaptic DEGs Expression by Group",
    x = "Group (Sex.Diagnosis)",
    y = "Expression Level"
  )

print(violin_plot_seurat)

# 方法2：使用ggplot2完全控制
library(ggplot2)
library(dplyr)
library(tidyr)

# 准备数据
violin_data <- FetchData(ex_neu, vars = c(available_genes, "sex", "diagnosis"))

# 转换为长格式
violin_long <- violin_data %>%
  pivot_longer(
    cols = all_of(available_genes),
    names_to = "Gene",
    values_to = "Expression"
  )

# 创建分组变量
violin_long$Group <- paste(violin_long$sex, violin_long$diagnosis, sep = ".")

# 设置因子顺序，确保图形顺序一致
violin_long$Group <- factor(violin_long$Group, 
                            levels = c("female.ctl", "female.AD", 
                                       "male.ctl", "male.AD"))
violin_long$Gene <- factor(violin_long$Gene, levels = available_genes)

# 创建完整的小提琴图
violin_ggplot <- ggplot(violin_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c(
      "female.ctl" = "#C6DBEF",  # 浅蓝色
      "female.AD" = "#F781BF",   # 粉色
      "male.ctl" = "#C7E9C0",    # 浅绿色
      "male.AD" = "#4DAF4A"      # 绿色
    ),
    labels = c(
      "female.ctl" = "Female Control",
      "female.AD" = "Female AD",
      "male.ctl" = "Male Control", 
      "male.AD" = "Male AD"
    )
  ) +
  labs(
    title = "Synaptic DEG Expression by Sex and Diagnosis",
    x = "Group (Sex - Diagnosis)",
    y = "Expression Level",
    fill = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    strip.background = element_rect(fill = "lightgray", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = guide_legend(ncol = 1))

print(violin_ggplot)

# 方法3：分开绘制每个性别，然后组合
# 女性数据
female_data <- violin_long %>% filter(sex == "female")
male_data <- violin_long %>% filter(sex == "male")

# 女性小提琴图
violin_female <- ggplot(female_data, aes(x = diagnosis, y = Expression, fill = diagnosis)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c("AD" = "#F781BF", "ctl" = "#C6DBEF")
  ) +
  labs(
    title = "Female: Synaptic DEG Expression in AD vs Control",
    x = "Diagnosis",
    y = "Expression Level",
    fill = "Diagnosis"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    strip.background = element_rect(fill = "lightgray")
  )

# 男性小提琴图
violin_male <- ggplot(male_data, aes(x = diagnosis, y = Expression, fill = diagnosis)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.5, outlier.shape = NA) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
  scale_fill_manual(
    values = c("AD" = "#4DAF4A", "ctl" = "#C7E9C0")
  ) +
  labs(
    title = "Male: Synaptic DEG Expression in AD vs Control",
    x = "Diagnosis",
    y = "Expression Level",
    fill = "Diagnosis"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    strip.background = element_rect(fill = "lightgray")
  )

# 组合图形
library(patchwork)
combined_sex_violin <- (violin_female / violin_male) +
  plot_annotation(
    title = "Synaptic DEG Expression by Sex",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

print(combined_sex_violin)

# 保存所有小提琴图
ggsave("synaptic_violin_all_groups.png", violin_ggplot, 
       width = 16, height = 12, dpi = 300)
ggsave("synaptic_violin_separate_sex.png", combined_sex_violin, 
       width = 16, height = 16, dpi = 300)
ggsave("synaptic_violin_seurat.png", violin_plot_seurat, 
       width = 16, height = 12, dpi = 300)