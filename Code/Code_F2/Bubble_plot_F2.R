ggplot(delta_prop, aes(x = sex, y = cell_type_ident)) +
  geom_point(aes(size = abs(delta), color = delta)) +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red", 
    midpoint = 0,
    name = "Δ Value"
  ) +
  scale_size_continuous(range = c(3, 10), name = "|Δ| Magnitude") +
  geom_text(aes(label = delta_label), size = 3.5, vjust = -1.5) +
  labs(
    x = "Sex",
    y = "Cell Type",
    title = "Cell Type Proportion Changes in AD vs Control"
  ) +
  theme_minimal()
