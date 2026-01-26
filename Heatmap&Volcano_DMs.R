# Load necessary packages
library(ggplot2)
library(pheatmap)
library(dplyr)

# Read the filtered metabolite data table
data <- read.csv("metabolite_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Preview the data
head(data)

# Filter metabolites with VIP > 1 and p < 0.05
filtered_data <- data %>%
  filter(VIP_pred_(TQN) > 1, Pvalue < 0.05)

# Preview filtered metabolites
head(filtered_data)

# 1. Plot Heatmap
# Extract the required metabolites and data
heatmap_data <- filtered_data %>%
  select(Metabolite, contains("TG")) %>%
  column_to_rownames("Metabolite")

# Ensure the expression matrix is numeric
heatmap_data <- as.data.frame(sapply(heatmap_data, as.numeric))

# Define group labels (here it is assumed the data contains a 'Group' column, modify as needed)
group_labels <- data.frame(
  Group = c("TS", "TG", "NS", "TG", "MS", "MG"),
  row.names = colnames(heatmap_data)
)

# Plot the heatmap
pheatmap(heatmap_data,
         scale = "row",                  # Scale by rows (standardize)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "average",
         annotation_col = group_labels,  # Sample group information
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of VIP > 1 and p < 0.05 Metabolites")

# Save the heatmap
ggsave("heatmap_VIP1_pvalue0.05.png")

# 2. Plot Volcano Plot
# Calculate log2FC (here assuming there is a log2FC column, adjust if needed)
filtered_data$Log2FC <- log2(filtered_data$FC.Solid.G)

# Plot the volcano plot
volcano_plot <- ggplot(filtered_data, aes(x = Log2FC, y = -log10(Pvalue))) +
  geom_point(aes(color = Regulate), alpha = 0.7) + 
  scale_color_manual(values = c("red", "blue")) +   # Red for upregulated, blue for downregulated
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10(p-value)") +
  theme(plot.title = element_text(hjust = 0.5))

# Save the volcano plot
ggsave("volcano_plot_VIP1_pvalue0.05.png", plot = volcano_plot)

# Export the filtered metabolites
write.csv(filtered_data, "filtered_metabolites.csv", row.names = FALSE)
