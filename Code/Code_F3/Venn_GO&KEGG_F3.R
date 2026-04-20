SATURATED_COLORS <- list(
  male = "#1E88E5",   
  female = "#D81B60"
)

create_optimized_venn <- function(ct) {

  ct_deg <- sig_deg %>% filter(CellType == ct)
  
  female_genes <- ct_deg %>% filter(Sex == "female") %>% pull(gene) %>% unique()
  male_genes <- ct_deg %>% filter(Sex == "male") %>% pull(gene) %>% unique()

  venn_list <- list(
    Female = female_genes,
    Male = male_genes
  )
  
  venn_plot <- ggVennDiagram(venn_list, label_alpha = 0) +
    scale_fill_gradient(low = "#F0F0F0", high = SATURATED_COLORS$female) +
    scale_color_manual(values = c(SATURATED_COLORS$female, SATURATED_COLORS$male)) +
    labs(title = paste(ct, "Sex-Specific DEGs"),
         caption = paste("Total DEGs:", length(unique(ct_deg$gene)))) +
    theme_classic() +
    theme(
      axis.line = element_blank(),   
      axis.text = element_blank(),     
      axis.ticks = element_blank(),     
      axis.title = element_blank(),       
      panel.grid = element_blank()        
    )
  

  venn_filename <- paste0("venn_", gsub(" ", "_", ct), ".png")
  ggsave(venn_filename, venn_plot, width = 8, height = 8)

  
  return(venn_filename)
}


BAR_COLORS <- list(
  male_low = "#C6E2FF",   
  male_high = "#1E62D9",   
  female_low = "#FFC0CB",  
  female_high = "#C00000"  
)


create_enrichment_barplot <- function(enrich_result, plot_type, sex, cell_type) {

  if (is.null(enrich_result) || 
      !methods::is(enrich_result, "enrichResult") ||
      is.null(enrich_result@result) || 
      nrow(enrich_result@result) == 0) {
    cat(cell_type, sex, plot_type)
    return(NULL)
  }
  
  if (sex == "male") {
    color_low <- BAR_COLORS$male_low
    color_high <- BAR_COLORS$male_high
    bar_color <- BAR_COLORS$male_high
  } else {
    color_low <- BAR_COLORS$female_low
    color_high <- BAR_COLORS$female_high
    bar_color <- BAR_COLORS$female_high
  }

  display_n <- min(15, nrow(enrich_result@result))
  top_terms <- head(enrich_result[order(enrich_result$p.adjust), ], display_n)

  p <- ggplot(top_terms, aes(x = reorder(Description, -log10(p.adjust)), 
                             y = -log10(p.adjust), 
                             fill = Count)) +
    geom_col(width = 0.7) +  
    scale_fill_gradient(low = color_low, high = color_high) + 
    labs(title = paste(cell_type, "-", sex, "-", plot_type),
         x = "Pathway",
         y = "-log10(Adjusted p-value)",
         fill = "Gene Count") + 
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.text.x = element_text(angle = 0, hjust = 1), 
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    coord_flip() 
  

  filename <- paste0("bar_", gsub(" ", "_", cell_type), "_", sex, "_", plot_type, ".png")
  ggsave(filename, p, width = 12, height = max(8, display_n * 0.6))

  
  return(filename)
}


for (cell_type in cell_types) {

  enrich_data <- celltype_enrichments[[cell_type]]
  

  if (!is.null(enrich_data$female) && !is.null(enrich_data$female$GO)) {
    create_enrichment_barplot(
      enrich_data$female$GO, 
      "GO", 
      "female", 
      cell_type
    )
  }
  

  if (!is.null(enrich_data$female) && !is.null(enrich_data$female$KEGG)) {
    create_enrichment_barplot(
      enrich_data$female$KEGG, 
      "KEGG", 
      "female", 
      cell_type
    )
  }
  

  if (!is.null(enrich_data$male) && !is.null(enrich_data$male$GO)) {
    create_enrichment_barplot(
      enrich_data$male$GO, 
      "GO", 
      "male", 
      cell_type
    )
  }
  

  if (!is.null(enrich_data$male) && !is.null(enrich_data$male$KEGG)) {
    create_enrichment_barplot(
      enrich_data$male$KEGG, 
      "KEGG", 
      "male", 
      cell_type
    )
  }
}


