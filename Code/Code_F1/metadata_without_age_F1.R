library(ggplot2)
library(dplyr)
library(patchwork)
library(RColorBrewer)

sample_data <- mzl@meta.data %>%
  distinct(orig.ident, .keep_all = TRUE) %>%
  select(orig.ident, diagnosis, sex, apoE, braak, cerad)

sex_dist <- ggplot(sample_data, aes(x = diagnosis, fill = sex)) +
  geom_bar(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("male" = "#4D8BD6", "female" = "#FF6B6B")) +
  labs(title = "Sex Distribution by Diagnosis", 
       x = "Diagnosis Group", 
       y = "Number of Samples", 
       fill = "Sex") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.major.x = element_blank())

braak_dist <- ggplot(sample_data, aes(x = factor(braak), fill = diagnosis)) +
  geom_bar(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(title = "Braak Stage Distribution", 
       x = "Braak Stage", 
       y = "Number of Samples", 
       fill = "Diagnosis") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

cerad_dist <- ggplot(sample_data, aes(x = factor(cerad), fill = diagnosis)) +
  geom_bar(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(title = "CERAD Score Distribution", 
       x = "CERAD Score", 
       y = "Number of Samples", 
       fill = "Diagnosis") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10))

apoE_levels <- sort(unique(sample_data$apoE))
complete_data <- expand.grid(apoE = apoE_levels, diagnosis = unique(sample_data$diagnosis), sex = unique(sample_data$sex))
count_data <- sample_data %>%
  group_by(apoE, diagnosis, sex) %>%
  summarise(count = n(), .groups = "drop")
apoE_data <- left_join(complete_data, count_data, by = c("apoE", "diagnosis", "sex")) %>%
  mutate(count = ifelse(is.na(count), 0, count))

apoE_dist <- ggplot(apoE_data, aes(x = factor(apoE), y = count, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = ifelse(count > 0, count, "")), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("male" = "#4D8BD6", "female" = "#FF6B6B")) +
  facet_wrap(~ diagnosis, nrow = 1) +
  labs(title = "apoE Genotype Distribution by Diagnosis and Sex", 
       x = "apoE Genotype", 
       y = "Number of Samples", 
       fill = "Sex") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(size = 10, face = "bold"))

braak_box <- ggplot(sample_data, aes(x = diagnosis, y = braak, fill = diagnosis)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(title = "Braak Stage by Diagnosis", 
       x = "Diagnosis Group", 
       y = "Braak Stage") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none")

cerad_box <- ggplot(sample_data, aes(x = diagnosis, y = cerad, fill = diagnosis)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("AD" = "#E41A1C", "ctl" = "#377EB8")) +
  labs(title = "CERAD Score by Diagnosis", 
       x = "Diagnosis Group", 
       y = "CERAD Score") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "none")

layout <- "
AABB
CCDD
EEFF
GGHH
"

combined_plot <- sex_dist + braak_dist + cerad_dist + apoE_dist + 
  braak_box + cerad_box + plot_spacer() + plot_spacer() +
  plot_layout(design = layout, heights = c(1, 1, 1, 0.2)) +
  plot_annotation(
    title = "Sample Clinical Characteristics",
    subtitle = "Based on available metadata fields",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 14, hjust = 0.5))
  )

ggsave("sample_clinical_features.png", combined_plot, width = 16, height = 14, dpi = 300)

combined_plot
