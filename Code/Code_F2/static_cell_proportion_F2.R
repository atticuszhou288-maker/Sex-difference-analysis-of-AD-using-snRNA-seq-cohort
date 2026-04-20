library(Seurat)
library(tidyverse)

# 1. 计算整体比例数据
cell_prop <- mzl@meta.data %>%
  filter(diagnosis %in% c("AD", "ctl") & sex %in% c("female", "male")) %>%
  count(diagnosis, sex, cell_type_ident) %>%
  group_by(diagnosis, sex) %>%
  mutate(total = sum(n), prop = n/total) %>%
  ungroup()

# 2. 比例差异检验函数
proportion_diff_test <- function(data, cell_type, sex_group) {
  # 提取特定细胞类型和性别的数据
  group_data <- data %>%
    filter(cell_type_ident == cell_type & sex == sex_group)
  
  # 获取AD组数据
  ad_group <- group_data %>% filter(diagnosis == "AD")
  n_ad <- ad_group$n
  n_total_ad <- ad_group$total
  
  # 获取Control组数据
  ctl_group <- group_data %>% filter(diagnosis == "ctl")
  n_ctl <- ctl_group$n
  n_total_ctl <- ctl_group$total
  
  # 计算比例差异
  p_ad <- n_ad / n_total_ad
  p_ctl <- n_ctl / n_total_ctl
  delta <- p_ad - p_ctl
  
  # 计算合并比例
  p_pooled <- (n_ad + n_ctl) / (n_total_ad + n_total_ctl)
  
  # 计算z统计量
  se <- sqrt(p_pooled * (1 - p_pooled) * (1/n_total_ad + 1/n_total_ctl))
  z <- delta / se
  
  # 计算p值 (双侧检验)
  p_value <- 2 * pnorm(-abs(z))
  
  # 返回结果
  data.frame(
    cell_type = cell_type,
    sex = sex_group,
    delta = delta,
    z = z,
    p_value = p_value
  )
}

# 3. 对所有细胞类型和性别组合进行检验
cell_types <- unique(cell_prop$cell_type_ident)
sex_groups <- c("female", "male")

results <- expand.grid(cell_type = cell_types, sex = sex_groups) %>%
  rowwise() %>%
  do(proportion_diff_test(cell_prop, .$cell_type, .$sex)) %>%
  ungroup()

# 4. 多重检验校正
results$adj_p_value <- p.adjust(results$p_value, method = "BH")
results$significance <- ifelse(results$adj_p_value < 0.05, "Significant", "Not significant")

# 5. 可视化
ggplot(results, aes(x = cell_type, y = delta, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = ifelse(significance == "Significant", "*", "")),
            position = position_dodge(width = 0.7), vjust = -0.5, size = 6) +
  scale_fill_manual(values = c("#E64B35", "#4DBBD5")) +
  labs(x = "Cell Type",
       y = "Proportion Difference (AD - Control)",
       title = "Proportion Differences by Diagnosis and Sex",
       fill = "Sex") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("cell_proportion_differences.png", width = 10, height = 8, dpi = 300)

