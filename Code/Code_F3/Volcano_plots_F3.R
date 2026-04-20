# 1. 加载必要的包
library(ggplot2)
library(dplyr)

# 2. 加载数据
all_de_genes <- readRDS("all_de_genes_results.rds")

# 3. 创建输出目录
if (!dir.exists("Individual_Volcanoes")) dir.create("Individual_Volcanoes")

# 4. 获取所有细胞类型
cell_types <- unique(all_de_genes$CellType)
print(paste("发现细胞类型:", paste(cell_types, collapse = ", ")))

# 5. 定义精调火山图函数
create_volcano_plot <- function(de_data, title) {
  # 处理没有数据的情况
  if (nrow(de_data) == 0) {
    return(
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "无显著基因", size = 4) + 
        labs(title = title) +
        theme_void()
    )
  }
  
  # 准备数据
  plot_data <- de_data %>%
    mutate(
      Significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.25 ~ "上调",
        p_val_adj < 0.05 & avg_log2FC < -0.25 ~ "下调",
        TRUE ~ "不显著"
      ),
      neg_log10_p = -log10(p_val_adj)
    )
  
  # 处理p值为0的情况
  plot_data$neg_log10_p[is.infinite(plot_data$neg_log10_p)] <- 350
  
  # 创建火山图
  p <- ggplot(plot_data, aes(x = avg_log2FC, y = neg_log10_p, color = Significance)) +
    geom_point(
      alpha = 0.6, 
      size = 0.6,  # 点大小
      shape = 16    # 实心圆点
    ) +
    geom_hline(
      yintercept = -log10(0.05), 
      linetype = "dashed", 
      color = "black", 
      linewidth = 0.3,
      alpha = 0.6
    ) +
    geom_vline(
      xintercept = c(-0.25, 0.25), 
      linetype = "dashed", 
      color = "black", 
      linewidth = 0.3,
      alpha = 0.6
    ) +
    scale_color_manual(
      values = c("上调" = "#e41a1c", "下调" = "#377eb8", "不显著" = "gray70")
    ) +
    coord_cartesian(ylim = c(0, 300)) +
    labs(
      title = title,
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10]("Adjusted p-value"))
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 1,  # 1:1长宽比
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 5, 5, "mm")
    )
  
  return(p)
}

# 6. 循环生成并保存所有火山图
for (ct in cell_types) {
  # 处理男性样本
  male_results <- all_de_genes %>% 
    filter(CellType == ct, Sex == "male")
  
  p_male <- create_volcano_plot(
    male_results, 
    paste(ct, "(Male)")
  )
  
  male_file <- file.path("Individual_Volcanoes", paste0("Volcano_", gsub(" ", "_", ct), "_Male.png"))
  
  ggsave(
    male_file,
    plot = p_male,
    width = 4,    # 4英寸宽度
    height = 4,   # 4英寸高度
    dpi = 300,    # 300 DPI
    bg = "white"  # 白色背景
  )
  
  cat("已保存男性火山图:", male_file, "\n")
  
  # 处理女性样本
  female_results <- all_de_genes %>% 
    filter(CellType == ct, Sex == "female")
  
  p_female <- create_volcano_plot(
    female_results, 
    paste(ct, "(Female)")
  )
  
  female_file <- file.path("Individual_Volcanoes", paste0("Volcano_", gsub(" ", "_", ct), "_Female.png"))
  
  ggsave(
    female_file,
    plot = p_female,
    width = 4,    # 4英寸宽度
    height = 4,   # 4英寸高度
    dpi = 300,    # 300 DPI
    bg = "white"  # 白色背景
  )
  
  cat("已保存女性火山图:", female_file, "\n")
}

# 7. 完成消息
cat("\n已完成所有火山图的生成和保存！\n")
cat("请检查目录: ", normalizePath("Individual_Volcanoes"), "\n")

