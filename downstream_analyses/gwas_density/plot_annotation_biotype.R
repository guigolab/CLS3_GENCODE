library(ggplot2)
library(dplyr)
library(svglite)

exonic_data <- data.frame(
  subgroup = c("bigTranscriptome", "CMfinderCRSs", "fantomCat", "fantomEnhancers", "GENCODE20+", "miTranscriptome", "non targeted", "NONCODE", "phyloCSF", "refSeq", "UCE", "VISTAenhancers"),
  value = c(10.50, 11.30, 10.71, 12.94, 10.06, 9.80, 8.36, 10.03, 13.73, 8.42, 12.65, 5.22),
  group = c("lncRNA", "Putative", "lncRNA", "Enhancers", "lncRNA", "lncRNA", "lncRNA", "lncRNA", "Putative", "lncRNA", "Putative", "Enhancers"),
  set = "Exonic")

gene_body_data <- data.frame(
  subgroup = c("bigTranscriptome", "CMfinderCRSs", "fantomCat", "fantomEnhancers", "GENCODE20+", "miTranscriptome", "non targeted", "NONCODE", "phyloCSF", "refSeq", "UCE", "VISTAenhancers"),
  value = c(5.08, 5.15, 4.84, 6.02, 4.88, 4.65, 4.43, 4.78, 5.90, 4.40, 4.58, 4.87),
  group = c("lncRNA", "Putative", "lncRNA", "Enhancers", "lncRNA", "lncRNA", "lncRNA", "lncRNA", "Putative", "lncRNA", "Putative", "Enhancers"),
  set = "Gene Body")

combined_data <- bind_rows(exonic_data, gene_body_data)

# Prepare the data for plotting
plot_df_combined <- combined_data %>%
  mutate(
    sum_length = value,
    value_label = round(value, 1), 
    text_y = sum_length + 0.5
  )

# Reorder subgroups by group and within each group by value
plot_df_combined <- plot_df_combined %>%
  arrange(group, -value) %>%  
  mutate(subgroup = factor(subgroup, levels = unique(subgroup)))  


color_palette_sets <- c(
  "Exonic" = "#c17ba0",
  "Gene Body" = "#C83737"
)

plt_combined <- ggplot(plot_df_combined) +
  geom_hline(yintercept = seq(0, ceiling(max(plot_df_combined$sum_length)), by = 1), color = "lightgrey") +
  geom_col(aes(x = subgroup, y = sum_length, fill = set,group = set),
    position = position_dodge(width = 0.9), width = 0.8, show.legend = TRUE, alpha = 0.9) +
  geom_text(aes(x = subgroup, y = text_y, label = value_label, group = set),
    position = position_dodge(width = 0.9), color = "black", size = 6, hjust = -0.15) + coord_polar(start = 0) +  theme_minimal() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 15, hjust = 1),
    legend.title = element_blank(), legend.text = element_text(size = 12), plot.margin = margin(0, 0, 0, 0)) + scale_fill_manual(values = color_palette_sets) +
  labs(title = "")

svglite("exonic_vs_gene_body.svg", width = 12, height = 10)
print(plt_combined)
dev.off()

# Plot display has then be polished using InkScape
