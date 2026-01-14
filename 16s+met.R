library(Hmisc)
library(reshape2)
library(dplyr)
library(pheatmap)

# 1. Reading Data
genus <- read.csv("16S_Genus_Abundance.csv", row.names = 1, check.names = FALSE)
metab <- read.csv("Metabolite_Processed_Data.csv", row.names = 1, check.names = FALSE)

# 2. Data Preprocessing
# Filter out genus names containing "norank" or "unclassified", remove "g__" prefix
genus <- genus[!grepl("norank|unclassified", rownames(genus), ignore.case = TRUE), ]
rownames(genus) <- sub("^g__", "", rownames(genus))

# 3. Transpose Data (samples as rows, variables as columns)
genus_t <- t(as.matrix(genus))
rownames(genus_t) <- colnames(genus)
colnames(genus_t) <- rownames(genus)

metab_t <- t(as.matrix(metab))
rownames(metab_t) <- colnames(metab)
colnames(metab_t) <- rownames(metab)

# 4. Sample Matching
common_samples <- intersect(rownames(genus_t), rownames(metab_t))
genus_t <- genus_t[common_samples, ]
metab_t <- metab_t[common_samples, ]
head(rownames(genus_t), 10)
head(rownames(metab_t), 10)

# 5. Spearman Correlation Calculation
res <- rcorr(as.matrix(metab_t), as.matrix(genus_t), type = "spearman")

n_metab <- ncol(metab_t)
cor_matrix <- res$r[1:n_metab, (n_metab + 1):ncol(res$r)]
p_matrix <- res$P[1:n_metab, (n_metab + 1):ncol(res$P)]

# 6. Reshape to Long Format and Merge P-values
cor_df <- melt(cor_matrix, varnames = c("Metabolite", "Genus"), value.name = "SpearmanR")
p_df <- melt(p_matrix, varnames = c("Metabolite", "Genus"), value.name = "Pvalue")
cor_df$Pvalue <- p_df$Pvalue

# 7. Filter Significant Correlations
sig_df <- cor_df %>%
  filter(abs(SpearmanR) > 0.65, Pvalue < 0.05) %>%
  arrange(desc(abs(SpearmanR)))

write.csv(sig_df, "Significant_Spearman_Correlation_byPvalue.csv", row.names = FALSE)

# 8. Get Significant Metabolites and Genera
sig_metabs <- unique(sig_df$Metabolite)
sig_genera <- unique(sig_df$Genus)

# 9. Subset Correlation Matrix
heatmap_mat <- cor_matrix[sig_metabs, sig_genera, drop = FALSE]

# 10. Handle NA/NaN/Inf Values
heatmap_mat[is.na(heatmap_mat)] <- 0
heatmap_mat[is.nan(heatmap_mat)] <- 0
heatmap_mat[is.infinite(heatmap_mat)] <- 0

# 11. Remove Rows and Columns with No Variance
keep_rows <- apply(heatmap_mat, 1, function(x) sd(x) != 0)
keep_cols <- apply(heatmap_mat, 2, function(x) sd(x) != 0)
heatmap_mat <- heatmap_mat[keep_rows, keep_cols, drop = FALSE]

# 12. Get Valid Row and Column Names
valid_metabs <- rownames(heatmap_mat)
valid_genera <- colnames(heatmap_mat)

# 13. Filter sig_df to Keep Valid Combinations
sig_df_filt <- sig_df %>%
  filter(Metabolite %in% valid_metabs, Genus %in% valid_genera)

# 14. Build Annotation Matrix with Stars (fix index out of bounds error)
annotation_mat <- matrix("", 
                         nrow = length(valid_metabs), 
                         ncol = length(valid_genera),
                         dimnames = list(valid_metabs, valid_genera))

# Use row and column position indices instead of name indices
for (i in seq_len(nrow(sig_df_filt))) {
  m <- as.character(sig_df_filt$Metabolite[i])
  g <- as.character(sig_df_filt$Genus[i])
  
  row_idx <- which(rownames(annotation_mat) == m)
  col_idx <- which(colnames(annotation_mat) == g)
  
  if (length(row_idx) > 0 && length(col_idx) > 0) {
    annotation_mat[row_idx, col_idx] <- "*"
  }
}

# 15. Clustering Switches
row_cluster <- length(valid_metabs) >= 2
col_cluster <- length(valid_genera) >= 2

# 16. Generate Heatmap and Save PDF
pdf("Spearman_Heatmap_withStars.pdf", width = 12, height = 10)
pheatmap(heatmap_mat,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_rows = row_cluster,
         cluster_cols = col_cluster,
         display_numbers = annotation_mat,
         number_color = "black",
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize_number = 20,
         angle_col = 45,
         main = "Spearman Correlations with Significant Stars")
dev.off()