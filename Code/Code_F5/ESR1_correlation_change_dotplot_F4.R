# 创建点图展示相关性的变化
dot_diff_plot <- ggplot(results_df, 
                        aes(x = reorder(Target_Gene, Correlation_Difference), 
                            y = Correlation_Difference)) +
  geom_segment(aes(x = reorder(Target_Gene, Correlation_Difference), 
                   xend = reorder(Target_Gene, Correlation_Difference),
                   y = 0, yend = Correlation_Difference),
               color = "gray50", size = 1) +
  geom_point(aes(fill = ifelse(Correlation_Difference > 0, "Increase", "Decrease"),
                 size = abs(Correlation_Difference) * 300),
             shape = 21, color = "black", stroke = 1) +
  geom_text(aes(label = sprintf("%+.3f", Correlation_Difference)), 
            vjust = -0.8, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(
    values = c("Increase" = "#e41a1c", "Decrease" = "#377eb8"),
    name = "Change Direction"
  ) +
  scale_size_identity() +
  labs(
    title = "ESR1-Target Gene Correlation Differences: AD vs Control",
    subtitle = "Size represents magnitude of change | Color indicates direction",
    x = "Target Gene",
    y = "ΔCorrelation (AD - Control)",
    caption = "Positive values: Stronger correlation in AD group\nNegative values: Weaker correlation in AD group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 14),
    legend.position = "bottom"
  ) +
  coord_flip()

# 保存点图
ggsave("ESR1_Correlation_Differences_DotPlot.png", dot_diff_plot, 
       width = 10, height = 7, dpi = 300, bg = "white")

