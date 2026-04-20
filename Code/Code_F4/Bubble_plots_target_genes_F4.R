# 加载必要库
library(tidyverse)
library(ggplot2)
library(scales) # 用于对数刻度

# 1. 计算每个基因在CTL组的平均表达量（作为本底规模）
base_expression <- stats_df %>%
  filter(diagnosis == "ctl") %>% # 只考虑对照组
  select(ends_with("_mean")) %>% # 选择表达量列
  summarise(across(ends_with("_mean"), mean, na.rm = TRUE)) %>% # 计算每列的平均值
  pivot_longer(
    everything(),
    names_to = "gene",
    values_to = "mean_expression"
  ) %>%
  mutate(gene = str_remove(gene, "_mean")) %>%
  # 过滤关键基因
  filter(gene %in% key_genes_ordered)

# 2. 合并log2FC数据和本底表达数据
bubble_data <- log2fc_data %>%
  # 注意：log2fc_data已经包含每个组合的log2FC
  left_join(base_expression, by = "gene") %>%
  # 对于本图，我们可以将每个基因的本底表达视为该基因在所有细胞类型和性别中的平均（CTL组）
  # 如果需要，也可以按细胞类型和性别计算，但这里为了简化，使用全局平均
  mutate(log10_mean_expression = log10(mean_expression + 1e-5)) # 避免log(0)，加一个小常数

ggplot(bubble_data, aes(x = group, y = gene)) +
  geom_point(
    aes(size = mean_expression, fill = log2FC), 
    shape = 21, color = "black", stroke = 0.3
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-1.5, 1.5),
    oob = scales::squish,   # 修改点1：极值裁剪
    na.value = "white",     # 修改点2：低表达设为白色
    name = "Log2 FC"
  ) +

  # 大小梯度（表达量）
  scale_size_continuous(range = c(1, 10), 
                        trans = "log10",  # 对数转换，因为表达量差异可能很大
                        breaks = trans_breaks("log10", function(x) 10^x),
                        labels = trans_format("log10", math_format(10^.x)),
                        name = "Mean Expression (Control)") +
  # 主题和布局
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position = "right",
    legend.box = "vertical"
  ) +
  labs(
    title = "Gene Expression Changes in AD vs Control",
    subtitle = "Bubble size: average expression in Control; Color: Log2(AD/Control) FC"
  ) +
  # 添加分隔线（按性别）
  geom_vline(xintercept = 6.5, linetype = "dashed", color = "grey40") +
  # 添加性别标签
  annotate("text", x = 3.5, y = length(key_genes_ordered) + 0.5, 
           label = "Female", size = 5, fontface = "bold", color = "#F8766D") +
  annotate("text", x = 9.5, y = length(key_genes_ordered) + 0.5, 
           label = "Male", size = 5, fontface = "bold", color = "#00BFC4") +
  # 调整坐标轴范围，留出标签空间
  coord_cartesian(ylim = c(0.5, length(key_genes_ordered) + 0.8), clip = "off")

# 4. 保存气泡图
ggsave("AD_Expression_Bubble_Plot.png", bubble_data, 
       width = 16, height = 10, dpi = 300)
