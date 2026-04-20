# ==================================
# Cell-Type Specific DEG Analysis (English Version)
# ==================================

# Load required packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(knitr)

# Prepare data: identify cell-type specific genes
celltype_specific_df <- sig_deg %>%
  filter(CellType != "Overall") %>%  # Remove Overall category
  group_by(gene) %>%
  mutate(
    is_specific = n_distinct(CellType) == 1  # TRUE if gene appears in only one cell type
  ) %>%
  ungroup()

# Calculate cell-type specific DEG statistics
specificity_stats <- celltype_specific_df %>%
  group_by(CellType) %>%
  summarise(
    Total_DEGs = n(),
    Non_Specific = sum(!is_specific),
    Specific_DEGs = sum(is_specific),
    .groups = "drop"
  ) %>%
  # Add overall summary
  bind_rows(
    summarise(.,
              across(where(is.numeric), sum),
              CellType = "Overall"
    ) %>%
      mutate(
        Specific_Ratio = Specific_DEGs / Total_DEGs * 100
      )
  )

# Create formatted statistics table
formatted_stats <- specificity_stats %>%
  mutate(
    Specific_Ratio = sprintf("%.1f%%", Specific_Ratio)
  ) %>%
  rename(
    `Cell Type` = CellType,
    `Total DEGs` = Total_DEGs,
    `Non-Specific DEGs` = Non_Specific,
    `Specific DEGs` = Specific_DEGs,
    `Specificity Ratio` = Specific_Ratio
  )

# Print statistics table
cat("===== Cell-Type Specific DEG Statistics =====\n")
kable(formatted_stats, align = "c", caption = "Summary of Cell-Type Specific DEGs")

# ==================================
# Visualization: Cell-Type Specific DEG Distribution
# ==================================

# Prepare data for plotting (exclude Overall for visualization)
plot_data <- specificity_stats %>%
  filter(CellType != "Overall") %>%
  pivot_longer(
    cols = c(Specific_DEGs, Non_Specific),
    names_to = "Gene_Type",
    values_to = "Count"
  ) %>%
  mutate(
    Gene_Type = factor(Gene_Type,
                       levels = c("Non_Specific", "Specific_DEGs"),
                       labels = c("Non-Specific DEGs", "Cell-Type Specific DEGs"))
  )

# Create distribution plot
specificity_plot <- ggplot(plot_data, aes(x = CellType, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.5, color = "black") +
  scale_fill_manual(values = c("Cell-Type Specific DEGs" = "#3498db", 
                               "Non-Specific DEGs" = "#e74c3c")) +
  labs(
    title = "Distribution of Cell-Type Specific DEGs",
    x = "Cell Type",
    y = "Number of Genes",
    fill = "Gene Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Display the plot
print(specificity_plot)

# ==================================
# Visualization: Specificity Ratio by Cell Type
# ==================================

ratio_plot_data <- specificity_stats %>%
  filter(CellType != "Overall") %>%
  mutate(
    Label = sprintf("%.1f%%", Specific_Ratio)
  )

specificity_ratio_plot <- ggplot(ratio_plot_data, 
                                 aes(x = reorder(CellType, Specific_Ratio), 
                                     y = Specific_Ratio)) +
  geom_bar(stat = "identity", fill = "#2ecc71", alpha = 0.8) +
  geom_text(aes(label = Label), vjust = -0.5, size = 4, color = "black") +
  labs(
    title = "Proportion of Cell-Type Specific DEGs",
    x = "Cell Type",
    y = "Specificity Ratio (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, max(ratio_plot_data$Specific_Ratio) * 1.2))

# Display the ratio plot
print(specificity_ratio_plot)

# ==================================
# Save Results
# ==================================

# Save statistics
write.csv(specificity_stats, "celltype_specific_deg_statistics.csv", row.names = FALSE)

# Save plots
ggsave("specific_deg_distribution.png", specificity_plot, 
       width = 10, height = 7, dpi = 300)
ggsave("specificity_ratio.png", specificity_ratio_plot, 
       width = 8, height = 6, dpi = 300)
