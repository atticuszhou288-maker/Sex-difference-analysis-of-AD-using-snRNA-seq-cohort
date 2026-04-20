names(mzl@reductions)
DimPlot(mzl, reduction = "umap", group.by = "sex")

png("gender_umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(mzl, reduction = "umap", group.by = "sex")
dev.off()

png("split_sex_umap.png", width = 12, height = 6, units = "in", res = 300)
DimPlot(mzl, reduction = "umap", group.by = "sex", split.by = "sex", pt.size = 0.3)
dev.off()

png("celltype_by_sex_umap.png", width = 14, height = 10, units = "in", res = 300)
DimPlot(mzl, reduction = "umap", group.by = "cell_type_ident", split.by = "sex", pt.size = 0.3)
dev.off()

mzl$group <- paste(mzl$diagnosis, mzl$sex, sep = "_")

mzl$group <- factor(mzl$group, levels = c("AD_F", "AD_M", "CON_F", "CON_M"))

table(mzl$group)

gender_colors <- c("male" = "#0066CC",
                   "female" = "#FF6B6B")

png("celltype_by_sex_umap.png", width = 14, height = 10, units = "in", res = 300)
DimPlot(mzl, 
        reduction = "umap", 
        group.by = "cell_type_ident", 
        split.by = "sex", 
        pt.size = 0.3,
        cols = celltype_colors) +  
  ggtitle("Cell Type Distribution by Sex")  
dev.off()

DimPlot(mzl, 
        reduction = "umap", 
        group.by = "cell_type_ident", 
        split.by = "sex", 
        pt.size = 0.3,
        cols = celltype_colors) + 
  ggtitle("Cell Type Distribution by Sex")

png("gender_overlaid_umap.png", width = 12, height = 8, units = "in", res = 300)
DimPlot(mzl, 
        reduction = "umap", 
        group.by = "sex", 
        pt.size = 0.3,
        cols = gender_colors) +  
  ggtitle("Overall Cell Distribution by Sex") 
dev.off()

DimPlot(mzl, 
        reduction = "umap", 
        group.by = "sex", 
        pt.size = 0.3,
        cols = gender_colors) + 
  ggtitle("Overall Cell Distribution by Sex")

p1 <- DimPlot(mzl, 
              reduction = "umap", 
              group.by = "cell_type_ident", 
              split.by = "sex", 
              pt.size = 0.3,
              cols = celltype_colors) + 
  ggtitle("Cell Type Distribution by Sex")

p2 <- DimPlot(mzl, 
              reduction = "umap", 
              group.by = "sex", 
              pt.size = 0.3,
              cols = gender_colors) + 
  ggtitle("Overall Cell Distribution by Sex")

combined_plot <- p1 / p2  
combined_plot

png("combined_gender_umaps.png", width = 14, height = 16, units = "in", res = 300)
combined_plot
dev.off()

ad_data <- subset(mzl, diagnosis == "AD")

control_data <- subset(mzl, diagnosis == "Control")

DimPlot(ad_data, reduction = "umap", group.by = "sex", cols = c("deepskyblue", "red")) +
  ggtitle("AD Group")

DimPlot(control_data, reduction = "umap", group.by = "sex", cols = c("deepskyblue", "red")) +
  ggtitle("Control Group")

png("disease_sex_umap.png", width = 14, height = 10, units = "in", res = 300)
DimPlot(mzl, 
        reduction = "umap", 
        group.by = "cell_type_ident", 
        split.by = "disease_sex", 
        ncol = 2,  
        pt.size = 0.3,
        cols = celltype_colors) + 
  ggtitle("Cell Type Distribution by Disease and Sex")
dev.off()

png("AD_sex_overlaid_umap.png", width = 12, height = 8, units = "in", res = 300)

DimPlot(ad_data, 
        reduction = "umap", 
        group.by = "sex",  
        pt.size = 0.3,     
        cols = c("male" = "#0066CC", "female" = "#FF6B6B"), 
        order = c("female", "male")) +  
  ggtitle("AD Patients: Male vs Female Cell Distribution") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

dev.off()

png("AD_vs_CTL_overlaid_vector.png", width = 12, height = 9, units = "in", res = 300)

DimPlot(mzl, 
        reduction = "umap", 
        group.by = "diagnosis",
        pt.size = 0.4,           
        cols = c("AD" = "#FF5252", "ctl" = "#4FC3F7"),  
        order = c("AD"),        
        raster = FALSE) +        
  ggtitle("AD vs Control Distribution (Vector)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "right")

dev.off()

