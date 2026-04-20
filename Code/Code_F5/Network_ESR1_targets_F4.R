# 定义基因列表
genes <- c("ESR1", "ESR2", "PGR", "GREB1", "BDNF", "SIRT1", "NRCAM", "VEGFA", "MYC", "CCND1", "BCL2")

# 定义基因功能关系
edges <- data.frame(
  from = c("ESR1", "ESR1", "ESR1", "ESR1", "ESR1", "Estrogen Response", "Estrogen Response"),
  to = c("ESR2", "PGR", "GREB1", "SIRT1", "BDNF", "Neuroprotection", "Cell Proliferation"),
  relationship = c("Coreceptor", "Primary Target", "Primary Target", "Secondary Target", "Secondary Target", "Functional Pathway", "Functional Pathway")
)

# 添加二级关系
edges <- rbind(edges, data.frame(
  from = c("Estrogen Response", "Estrogen Response", "Neuroprotection", "Neuroprotection", "Cell Proliferation"),
  to = c("NRCAM", "VEGFA", "BCL2", "CCND1", "MYC"),
  relationship = c("Downstream Effector", "Downstream Effector", "Downstream Effector", "Downstream Effector", "Downstream Effector")
))

# 获取所有唯一节点ID
node_ids <- unique(c(edges$from, edges$to))

# 创建节点数据框
nodes <- data.frame(
  name = node_ids,
  type = ifelse(node_ids %in% genes, "Gene", "Function")
)

# 创建网络 - 确保包含关系属性
net <- graph_from_data_frame(d = edges, vertices = nodes)

# 使用自动布局避免错误
layout <- create_layout(net, layout = "fr")

# 添加节点分类信息
layout$category <- ifelse(layout$name %in% genes, "Gene", 
                          ifelse(layout$name %in% c("Estrogen Response", "Neuroprotection", "Cell Proliferation"), 
                                 "Functional Module", "Pathway Component"))

# 创建边标签数据
edge_labels <- data.frame(
  from = edges$from,
  to = edges$to,
  relationship = edges$relationship,
  x = (layout$x[match(edges$from, layout$name)] + layout$x[match(edges$to, layout$name)]) / 2,
  y = (layout$y[match(edges$from, layout$name)] + layout$y[match(edges$to, layout$name)]) / 2
)

# 可视化 - 简化版确保不报错
gene_network <- ggraph(layout) +
  # 边设置
  geom_edge_link(
    aes(color = relationship),
    arrow = arrow(length = unit(3, 'mm'), type = "closed"),
    end_cap = circle(7, 'mm'),
    start_cap = circle(7, 'mm'),
    alpha = 0.8,
    width = 1.2
  ) +
  # 边标签（关系类型）
  geom_label(
    data = edge_labels,
    aes(x = x, y = y, label = relationship, fill = relationship),
    size = 3,
    color = "white",
    alpha = 0.9,
    label.size = 0,
    label.padding = unit(0.15, "lines"),
    show.legend = FALSE
  ) +
  # 节点设置
  geom_node_point(
    aes(fill = category),
    size = 15,
    shape = 21,
    color = "black",
    stroke = 1.5
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE,  # 自动调整标签位置避免重叠
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    max.overlaps = 100
  ) +
  # 颜色标度
  scale_edge_color_manual(
    values = c(
      "Coreceptor" = "#4daf4a",
      "Primary Target" = "#984ea3",
      "Secondary Target" = "#ff7f00",
      "Functional Pathway" = "#377eb8",
      "Downstream Effector" = "#a65628"
    ),
    name = "Relationship"
  ) +
  scale_fill_manual(
    values = c(
      "Functional Module" = "#8da0cb",
      "Pathway Component" = "#66c2a5",
      "Gene" = "#fdcc8a"
    ),
    name = "Node Category"
  ) +
  # 主题和图例
  theme_void() +
  labs(
    title = "Estrogen Signaling Pathway Gene Relationships",
    subtitle = "Rationale for gene selection in correlation analysis",
    caption = "Relationship types labeled on edges | Automatic layout"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 16),
    plot.caption = element_text(hjust = 0.5, color = "grey50", size = 12),
    legend.position = "right",
    legend.box = "vertical",
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold", size = 12)
  )

# 显示基因关系网络图
print(gene_network)

# 保存基因关系网络图
ggsave("Estrogen_Gene_Relationship_Network.png", gene_network, 
       width = 13, height = 10, dpi = 300, bg = "white")
