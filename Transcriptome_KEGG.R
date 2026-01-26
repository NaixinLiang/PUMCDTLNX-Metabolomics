# Load necessary packages
library(clusterProfiler)
library(KEGGREST)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)  # Human annotation for KEGG IDs (for human species)

# Read the transcriptome data (assuming it has a 'KEGG ID' column)
transcriptome_data <- read.csv("transcriptome_data.csv", header = TRUE, stringsAsFactors = FALSE)

# Preview the data
head(transcriptome_data)

# Assuming the KEGG ID column is named "KEGG ID"
kegg_ids_transcriptome <- transcriptome_data$`KEGG ID`

# Remove NA values
valid_kegg_ids_transcriptome <- na.omit(kegg_ids_transcriptome)

# Perform KEGG enrichment analysis using clusterProfiler
kegg_enrichment_transcriptome <- enrichKEGG(gene = valid_kegg_ids_transcriptome, 
                                            organism = 'hsa',   # Human species
                                            pvalueCutoff = 0.05,
                                            qvalueCutoff = 0.1)

# View the enrichment result
summary(kegg_enrichment_transcriptome)

# Plot the KEGG enrichment bubble plot
kegg_dotplot_transcriptome <- dotplot(kegg_enrichment_transcriptome)
ggsave("KEGG_Enrichment_Bubble_Plot_Transcriptome.png", kegg_dotplot_transcriptome)

# Plot the KEGG enrichment bar plot
kegg_barplot_transcriptome <- barplot(kegg_enrichment_transcriptome, showCategory = 10)
ggsave("KEGG_Enrichment_Bar_Plot_Transcriptome.png", kegg_barplot_transcriptome)

# Save the KEGG enrichment result as a CSV file
write.csv(as.data.frame(kegg_enrichment_transcriptome), "KEGG_Enrichment_Transcriptome_Results.csv", row.names = FALSE)
