library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

deg_stats <- all_de_genes %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  mutate(
    Direction = ifelse(avg_log2FC > 0, "Upregulated", "Downregulated")
  ) %>%
  group_by(CellType, Sex, Direction) %>%
  summarise(Count = n(), .groups = "drop") %>%

  mutate(SignedCount = ifelse(Direction == "Downregulated", -Count, Count)) %>%

  mutate(Group = paste(Sex, Direction, sep = "_"))

heatmap_data <- deg_stats %>%
  select(CellType, Group, SignedCount) %>%

  complete(CellType, Group, fill = list(SignedCount = 0)) %>%

  pivot_wider(
    names_from = Group,
    values_from = SignedCount
  ) %>%

print(heatmap_data)

heatmap_data_renamed <- heatmap_data %>%
  rename(
    "Female_Up" = "female_Upregulated",
    "Female_Down" = "female_Downregulated",
    "Male_Up" = "male_Upregulated",
    "Male_Down" = "male_Downregulated"
  ) %>%

  select(CellType, Female_Up, Male_Up, Female_Down, Male_Down)

print(heatmap_data_renamed)

rownames(heatmap_data_renamed) <- heatmap_data_renamed$CellType
heatmap_matrix <- as.matrix(heatmap_data_renamed[, -1])

data_range <- range(heatmap_matrix)
if (diff(data_range) == 0) data_range <- c(-1, 1)  
breaks <- seq(data_range[1], data_range[2], length.out = 101)

p <- pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("#3182bd", "#f7fbff", "#de2d26"))(100),
  breaks = breaks,
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  fontsize_row = 12,
  fontsize_col = 12,
  cellwidth = 40,
  cellheight = 30,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 45,
  border_color = "grey80",
  main = "Differential Expressed Genes",
  legend = TRUE
)

png("DEG_final_corrected_heatmap.png", width = 12, height = 10, units = "in", res = 300)
p
dev.off()

stats_table <- heatmap_data_renamed %>%
  rename("Cell Type" = CellType)

write.csv(stats_table, "DEG_final_corrected_stats.csv", row.names = FALSE)

print(stats_table)


library(pheatmap)
library(RColorBrewer)

heatmap_matrix <- matrix(
  c(
    # ast
    137, 143, -474, -123,
    # ex.neu
    303, 68, -393, -138,
    # in.neu
    80, 63, -187, -145,
    # mic
    151, 150, -556, -157,
    # oli
    446, 154, -585, -96,
    # opc
    42, 51, -117, -48,
    # Overall
    239, 67, -371, -71
  ),
  ncol = 4,
  byrow = TRUE
)

rownames(heatmap_matrix) <- c("ast", "ex.neu", "in.neu", "mic", "oli", "opc", "Overall")
colnames(heatmap_matrix) <- c("Female_Up", "Male_Up", "Female_Down", "Male_Down")

print(heatmap_matrix)

data_range <- range(heatmap_matrix)
if (diff(data_range) == 0) data_range <- c(-600, 600)

p <- pheatmap(
  heatmap_matrix,
  color = colorRampPalette(c("#3182bd", "#f7fbff", "#de2d26"))(100),
  breaks = seq(data_range[1], data_range[2], length.out = 101),
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "black",
  fontsize_row = 12,
  fontsize_col = 12,
  cellwidth = 40,
  cellheight = 30,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  angle_col = 45,
  border_color = "grey80",
  main = "Differential Expressed Genes",
  legend = TRUE,
  show_rownames = TRUE  
)

png("DEG_manual_heatmap.png", width = 12, height = 10, units = "in", res = 300)
p
dev.off()

stats_table <- as.data.frame(heatmap_matrix)
stats_table$CellType <- rownames(stats_table)
stats_table <- stats_table[, c("CellType", "Female_Up", "Male_Up", "Female_Down", "Male_Down")]

write.csv(stats_table, "DEG_manual_stats.csv", row.names = FALSE)

print(stats_table)

