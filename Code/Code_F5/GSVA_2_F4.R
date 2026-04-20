# 步骤1：分离女性样本数据 ------------------------------------------
female_data <- plot_data[plot_data$Gender == "Female", ]

# 步骤2：女性样本效应量重新计算（只关注AD vs Control）----------------
female_effect <- female_data %>%
  group_by(CellType, Pathway) %>%
  summarise(
    mean_AD = mean(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    mean_Control = mean(EnrichmentScore[Group == "ctl"], na.rm = TRUE),
    n_AD = sum(Group == "AD"),
    n_Control = sum(Group == "ctl"),
    sd_AD = sd(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    sd_Control = sd(EnrichmentScore[Group == "ctl"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Delta = mean_AD - mean_Control,
    se = sqrt((sd_AD^2 / n_AD) + (sd_Control^2 / n_Control)),
    CI_lower = Delta - 1.96 * se,
    CI_upper = Delta + 1.96 * se
  )

# 步骤3：女性样本统计检验 ------------------------------------------
female_stats <- female_data %>%
  group_by(CellType, Pathway) %>%
  do({
    if(n_distinct(.$Group) < 2) {
      data.frame(p.value = NA_real_)
    } else {
      tidy(wilcox.test(EnrichmentScore ~ Group, data = .))
    }
  }) %>%
  ungroup() %>%
  mutate(
    Significance = case_when(
      is.na(p.value) ~ "NA",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# 步骤4：修正女性专属图表（消除错误）-------------------------------
female_plot <- ggplot(female_effect, 
                      aes(x = CellType, y = Delta, color = Pathway)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.7),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  geom_text(
    data = female_stats,
    aes(label = Significance, y = female_effect$CI_upper + 0.05),
    position = position_dodge(width = 0.7),
    size = 5,
    fontface = "bold",
    vjust = 0
  ) +
  scale_color_manual(
    values = c("Early Response" = "#e63946", "Late Response" = "#2a9d8f")
  ) +
  labs(
    title = "Estrogen Pathway Activity in Female AD vs Control",
    subtitle = "Effect Size (Δ) = AD Mean - Control Mean",
    x = "Cell Type",
    y = "Effect Size (Δ)",
    color = "Pathway",
    caption = "*** p < 0.001; ** p < 0.01; * p < 0.05; ns: not significant"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

# 步骤5：创建男性对照图表 -----------------------------------------
male_data <- plot_data[plot_data$Gender == "Male", ]

male_effect <- male_data %>%
  group_by(CellType, Pathway) %>%
  summarise(
    mean_AD = mean(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    mean_Control = mean(EnrichmentScore[Group == "ctl"], na.rm = TRUE),
    n_AD = sum(Group == "AD"),
    n_Control = sum(Group == "ctl"),
    sd_AD = sd(EnrichmentScore[Group == "AD"], na.rm = TRUE),
    sd_Control = sd(EnrichmentScore[Group == "ctl"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Delta = mean_AD - mean_Control,
    se = sqrt((sd_AD^2 / n_AD) + (sd_Control^2 / n_Control)),
    CI_lower = Delta - 1.96 * se,
    CI_upper = Delta + 1.96 * se
  )

male_plot <- ggplot(male_effect, 
                    aes(x = CellType, y = Delta, color = Pathway)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.7),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  labs(
    title = "Male Control: Estrogen Pathway Activity (AD vs Control)",
    subtitle = "Effect Size (Δ) = AD Mean - Control Mean",
    x = "Cell Type",
    y = "Effect Size (Δ)",
    color = "Pathway"
  ) +
  scale_color_manual(
    values = c("Early Response" = "#1f77b4", "Late Response" = "#ff7f0e")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "grey50"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey60"),
    legend.position = "bottom"
  )

# 步骤6：组合图表并保存 ------------------------------------------
library(patchwork)

# 垂直组合：女性在上，男性在下
combined_plot <- female_plot / male_plot +
  plot_layout(heights = c(2, 1)) +  # 女性图表占2/3空间
  plot_annotation(title = "Estrogen Pathway Activity: Female Focus with Male Control",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16)))

# 保存组合图表
ggsave("Estrogen_Response_Female_Focus_with_Male_Control.png", 
       combined_plot, 
       width = 12, height = 14, dpi = 300, bg = "white")

# 步骤7：直接对比性别差异（效应量差异）---------------------------
# 计算男女效应量差异
gender_delta_diff <- female_effect %>%
  rename(Female_Delta = Delta) %>%
  left_join(male_effect %>% select(CellType, Pathway, Male_Delta = Delta), 
            by = c("CellType", "Pathway")) %>%
  mutate(
    Delta_Diff = Female_Delta - Male_Delta,
    Fold_Change = ifelse(abs(Male_Delta) > 0.001, 
                         Female_Delta / Male_Delta, 
                         NA)
  )

# 创建性别差异热图
gender_diff_heatmap <- ggplot(gender_delta_diff, 
                              aes(x = CellType, y = Pathway, fill = Delta_Diff)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", Delta_Diff)), 
            color = "black", size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#377eb8", 
    mid = "white",
    high = "#e41a1c",
    midpoint = 0,
    name = "Δ(Female - Male)"
  ) +
  labs(
    title = "Sex Difference in Estrogen Pathway Response",
    subtitle = "Δ Difference = (Female AD-ctl) - (Male AD-ctl)",
    x = "Cell Type",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    panel.grid = element_blank()
  )

# 保存热图
ggsave("Sex_Difference_Heatmap.png", gender_diff_heatmap,
       width = 10, height = 6, dpi = 300, bg = "white")

# ====================== 数据保存方案 ======================

# 步骤1：创建专用目录（如果不存在）
if(!dir.exists("GSVA_Results")) {
  dir.create("GSVA_Results")
  cat("创建目录: GSVA_Results\n")
}

# 步骤2：保存原始GSVA结果 --------------------------------
# 保存完整GSVA分数矩阵
write.csv(as.data.frame(gsva_results), 
          "GSVA_Results/Full_GSVA_Scores.csv", 
          row.names = TRUE)

# 保存RDS格式（保留原始结构）
saveRDS(gsva_results, "GSVA_Results/Full_GSVA_Scores.rds")

# 步骤3：保存带元数据的结果 ------------------------------
# 完整元数据合并结果
write.csv(gsva_df, 
          "GSVA_Results/GSVA_Results_with_Metadata.csv", 
          row.names = FALSE)

# 步骤4：保存性别分层数据 --------------------------------
# 女性样本GSVA结果
female_gsva <- gsva_results[, plot_data$CellID[plot_data$Gender == "Female"]]
write.csv(as.data.frame(female_gsva), 
          "GSVA_Results/Female_GSVA_Scores.csv", 
          row.names = TRUE)

# 男性样本GSVA结果
male_gsva <- gsva_results[, plot_data$CellID[plot_data$Gender == "Male"]]
write.csv(as.data.frame(male_gsva), 
          "GSVA_Results/Male_GSVA_Scores.csv", 
          row.names = TRUE)

# 步骤5：保存效应量数据 --------------------------------
# 女性效应量
write.csv(female_effect, 
          "GSVA_Results/Female_Effect_Sizes.csv", 
          row.names = FALSE)

# 男性效应量
write.csv(male_effect, 
          "GSVA_Results/Male_Effect_Sizes.csv", 
          row.names = FALSE)

# 性别差异热图数据
gender_diff_data <- female_effect %>%
  select(CellType, Pathway, Female_Delta = Delta) %>%
  left_join(male_effect %>% select(CellType, Pathway, Male_Delta = Delta),
            by = c("CellType", "Pathway")) %>%
  mutate(Delta_Diff = Female_Delta - Male_Delta)

write.csv(gender_diff_data, 
          "GSVA_Results/Sex_Difference_Data.csv", 
          row.names = FALSE)

# 步骤6：保存统计结果 --------------------------------
# 女性统计结果
write.csv(female_stats, 
          "GSVA_Results/Female_Statistical_Results.csv", 
          row.names = FALSE)

# 男性统计结果
write.csv(male_stats, 
          "GSVA_Results/Male_Statistical_Results.csv", 
          row.names = FALSE)

# 步骤7：保存绘图数据 --------------------------------
# 热图数据
write.csv(gender_diff_data, 
          "GSVA_Results/Sex_Difference_Heatmap_Data.csv", 
          row.names = FALSE)

# 女性点图数据
write.csv(female_effect, 
          "GSVA_Results/Female_Effect_Plot_Data.csv", 
          row.names = FALSE)

# 男性点图数据
write.csv(male_effect, 
          "GSVA_Results/Male_Effect_Plot_Data.csv", 
          row.names = FALSE)

# 步骤8：保存分析日志 --------------------------------
analysis_log <- data.frame(
  AnalysisStep = c("GSVA Calculation", 
                   "Gender Stratification",
                   "Effect Size Calculation",
                   "Statistical Testing"),
  Timestamp = c(Sys.time(), Sys.time(), Sys.time(), Sys.time()),
  Parameters = c(
    paste("Method: ssgsea | Gene Sets: 2 | Cells:", ncol(gsva_results)),
    paste("Female Cells:", sum(plot_data$Gender == "Female"),
          "| Male Cells:", sum(plot_data$Gender == "Male")),
    paste("Cell Types:", length(unique(plot_data$CellType))),
    paste("Test: Wilcoxon | Significance Levels: p<0.05, p<0.01, p<0.001")
  )
)

write.csv(analysis_log, 
          "GSVA_Results/Analysis_Log.csv", 
          row.names = FALSE)

# 步骤9：创建ZIP压缩包（可选）------------------------
# 将所有结果打包为单个文件
zip::zipr("GSVA_Results_Full.zip", 
          list.files("GSVA_Results", full.names = TRUE))

# 步骤10：输出保存摘要 --------------------------------
cat("\n===== GSVA分析结果保存完成 =====\n")
cat("保存目录: GSVA_Results/\n")
cat("包含文件:\n")
cat("1. Full_GSVA_Scores.csv - 完整GSVA分数矩阵\n")
cat("2. GSVA_Results_with_Metadata.csv - 带元数据的GSVA结果\n")
cat("3. Female_GSVA_Scores.csv - 女性样本GSVA分数\n")
cat("4. Male_GSVA_Scores.csv - 男性样本GSVA分数\n")
cat("5. Female_Effect_Sizes.csv - 女性效应量数据\n")
cat("6. Male_Effect_Sizes.csv - 男性效应量数据\n")
cat("7. Sex_Difference_Data.csv - 性别差异数据\n")
cat("8. Female_Statistical_Results.csv - 女性统计结果\n")
cat("9. Male_Statistical_Results.csv - 男性统计结果\n")
cat("10. Sex_Difference_Heatmap_Data.csv - 热图数据\n")
cat("11. Analysis_Log.csv - 分析日志\n")
cat("12. Full_GSVA_Scores.rds - R格式原始数据\n")
cat("\n压缩包: GSVA_Results_Full.zip (包含所有结果)\n")

# ====================== 精确匹配的显著性标记 ======================

# 步骤1：修正组合图（添加男性显著性标记）------------------------
# 创建组合图（含男性和女性的显著性标记）
combined_plot_corrected <- female_plot_corrected / male_plot_corrected +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(title = "Estrogen Pathway Activity: Female Focus with Male Control",
                  theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18)))

# 保存修正后的组合图
ggsave("Estrogen_Response_Female_Focus_with_Male_Control_Corrected.png", 
       combined_plot_corrected, 
       width = 12, height = 14, dpi = 300, bg = "white")

# 步骤2：精确女性图表（标记颜色匹配通路）------------------------
female_plot_precise <- ggplot(female_effect, 
                              aes(x = CellType, y = Delta, color = Pathway)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.7),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  
  # 为每个点添加匹配颜色的显著性标记
  geom_text(
    data = female_effect %>% 
      left_join(female_stats, by = c("CellType", "Pathway")) %>%
      mutate(LabelY = CI_upper + 0.05 * max(CI_upper - CI_lower, na.rm = TRUE)),
    aes(y = LabelY, label = Significance, color = Pathway),  # 关键：颜色与点匹配
    position = position_dodge(width = 0.7),
    size = 5,
    fontface = "bold",
    vjust = 0
  ) +
  
  scale_color_manual(
    values = c("Early Response" = "#e63946", "Late Response" = "#2a9d8f")
  ) +
  labs(
    title = "Estrogen Pathway Activity in Female AD vs Control",
    subtitle = "Effect Size (Δ) = AD Mean - Control Mean",
    x = "Cell Type",
    y = "Effect Size (Δ)",
    color = "Pathway",
    caption = "Significance: *** p < 0.001; ** p < 0.01; * p < 0.05; ns: not significant"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, color = "#e63946"),
    plot.subtitle = element_text(color = "grey40", size = 16, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  ) +
  # 扩展Y轴空间确保标记可见
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# 步骤3：精确男性图表（标记颜色匹配通路）------------------------
male_plot_precise <- ggplot(male_effect, 
                            aes(x = CellType, y = Delta, color = Pathway)) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  geom_errorbar(
    aes(ymin = CI_lower, ymax = CI_upper),
    width = 0.2,
    position = position_dodge(width = 0.7),
    linewidth = 1
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
  
  # 为每个点添加匹配颜色的显著性标记
  geom_text(
    data = male_effect %>% 
      left_join(male_stats, by = c("CellType", "Pathway")) %>%
      mutate(LabelY = CI_upper + 0.05 * max(CI_upper - CI_lower, na.rm = TRUE)),
    aes(y = LabelY, label = Significance, color = Pathway),  # 关键：颜色与点匹配
    position = position_dodge(width = 0.7),
    size = 5,
    fontface = "bold",
    vjust = 0
  ) +
  
  scale_color_manual(
    values = c("Early Response" = "#1f77b4", "Late Response" = "#ff7f0e")
  ) +
  labs(
    title = "Estrogen Pathway Activity in Male AD vs Control",
    subtitle = "Effect Size (Δ) = AD Mean - Control Mean",
    x = "Cell Type",
    y = "Effect Size (Δ)",
    color = "Pathway",
    caption = "Significance: *** p < 0.001; ** p < 0.01; * p < 0.05; ns: not significant"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5, color = "#1f77b4"),
    plot.subtitle = element_text(color = "grey40", size = 16, hjust = 0.5),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom"
  ) +
  # 扩展Y轴空间确保标记可见
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# 步骤4：保存精确匹配的图表 -----------------------------------
# 女性图表（精确匹配）
ggsave("Estrogen_Response_Female_Focus_Precise.png", 
       female_plot_precise,
       width = 14, height = 8, dpi = 300, bg = "white")

# 男性图表（精确匹配）
ggsave("Estrogen_Response_Male_Focus_Precise.png", 
       male_plot_precise,
       width = 14, height = 8, dpi = 300, bg = "white")