# 创建样本级别的汇总数据（每个orig.ident代表一个个体）
library(dplyr)
sample_meta <- mzl@meta.data %>%
  distinct(orig.ident, .keep_all = TRUE) %>% # 每个样本只保留一行
  select(orig.ident, diagnosis, sex, age, braak, cerad) %>%
  mutate(
    braak = factor(braak, levels = sort(unique(braak)), ordered = TRUE),
    cerad = factor(cerad, levels = sort(unique(cerad)), ordered = TRUE)
  )

# 检查数据结构
head(sample_meta)

library(ggplot2)
library(ggpubr)

# 创建样本级别的箱线图（修正颜色映射）
p_age <- ggplot(sample_meta, aes(x = diagnosis, y = age, fill = sex)) +
  geom_boxplot(
    outlier.size =1, 
    color = "black", 
    lwd = 0.5,
    alpha = 0.8
  ) +
  # 修正点：使用小写性别标签匹配数据
  scale_fill_manual(
    values = c("female" = "#E64B35", "male" = "#4DBBD5"), # 改为小写
    name = "Sex"
  ) +
  geom_jitter(
    width = 0.2, 
    height = 0, 
    aes(color = sex),
    size = 1,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("female" = "#B03A26", "male" = "#2874A6") # 改为小写
  ) +
  stat_compare_means(
    aes(group = sex), 
    label = "p.format",
    method = "t.test",
    label.y = max(sample_meta$age) + 3
  ) +
  geom_text(
    data = sample_meta %>% group_by(diagnosis, sex) %>% tally(),
    aes(x = diagnosis, y = min(sample_meta$age) - 5, 
        label = paste0("n=", n)),
    position = position_dodge(width = 0.8),
    size = 2
  ) +
  labs(
    x = "diagnosis", 
    y = "Age", 
    title = "Age Distribution"
  ) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "right",
    aspect.ratio = 1
  )


print(p_age)



# 创建Braak评分图（精细坐标轴）
p_braak <- ggplot(sample_meta %>% mutate(diagnosis = ifelse(diagnosis == "Control", "ctl", diagnosis)), 
                  aes(x = age, y = as.numeric(braak))) +
  geom_point(
    aes(color = diagnosis),
    size = 2.5,
    alpha = 0.7
  ) +
  geom_smooth(
    aes(color = diagnosis, fill = diagnosis),
    method = "lm",
    se = TRUE,
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c("AD" = "#D55E00", "ctl" = "#0072B2"),
    name = "Diagnosis"
  ) +
  scale_fill_manual(
    values = c("AD" = "#D55E00", "ctl" = "#0072B2"),
    name = "Diagnosis"
  ) +
  labs(
    x = "Age", 
    y = "Braak Stage",
    title = "Braak Staging by Age"
  ) +
  theme_classic(base_size = 8) +  # 基础文字大小改为8pt
  theme(
    legend.position = "right",
    # 精细坐标轴设置
    axis.line = element_line(color = "black", size = 0.2),  # 极细坐标轴线 (0.2pt)
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),  # 细框线 (0.3pt)
    # 文字优化
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  # 标题10pt
    axis.title = element_text(size = 8),  # 轴标题8pt
    axis.text = element_text(size = 7, color = "black"),  # 轴刻度文字7pt
    legend.title = element_text(size = 8),  # 图例标题8pt
    legend.text = element_text(size = 7),  # 图例文字7pt
    legend.key.size = unit(0.6, "lines"),  # 更紧凑的图例
    aspect.ratio = 1
  )

print(p_braak)

# 创建CERAD评分图（精细坐标轴）
p_cerad <- ggplot(sample_meta %>% mutate(diagnosis = ifelse(diagnosis == "Control", "ctl", diagnosis)), 
                  aes(x = age, y = as.numeric(cerad))) +
  geom_point(
    aes(color = diagnosis),
    size = 2.5,
    alpha = 0.7
  ) +
  geom_smooth(
    aes(color = diagnosis, fill = diagnosis),
    method = "lm",
    se = TRUE,
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c("AD" = "#D55E00", "ctl" = "#0072B2"),
    name = "Diagnosis"
  ) +
  scale_fill_manual(
    values = c("AD" = "#D55E00", "ctl" = "#0072B2"),
    name = "Diagnosis"
  ) +
  labs(
    x = "Age", 
    y = "CERAD Score",
    title = "CERAD Scoring by Age"
  ) +
  theme_classic(base_size = 8) +  # 基础文字大小改为8pt
  theme(
    legend.position = "right",
    # 精细坐标轴设置
    axis.line = element_line(color = "black", size = 0.2),  # 极细坐标轴线 (0.2pt)
    panel.border = element_rect(color = "black", fill = NA, size = 0.3),  # 细框线 (0.3pt)
    # 文字优化
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),  # 标题10pt
    axis.title = element_text(size = 8),  # 轴标题8pt
    axis.text = element_text(size = 7, color = "black"),  # 轴刻度文字7pt
    legend.title = element_text(size = 8),  # 图例标题8pt
    legend.text = element_text(size = 7),  # 图例文字7pt
    legend.key.size = unit(0.6, "lines"),  # 更紧凑的图例
    aspect.ratio = 1
  )

print(p_cerad)
