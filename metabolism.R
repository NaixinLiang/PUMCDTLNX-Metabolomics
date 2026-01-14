if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

if (!require("ropls", quietly = TRUE)) {
  BiocManager::install("ropls", ask = FALSE)
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!require("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos = "https://cloud.r-project.org")
}
if (!require("vegan", quietly = TRUE)) {
  install.packages("vegan", repos = "https://cloud.r-project.org")
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}
if (!require("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
}
if (!require("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx", repos = "https://cloud.r-project.org")
}
if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler", ask = FALSE)
}
if (!require("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", repos = "https://cloud.r-project.org")
}

library(ropls)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(pheatmap)
library(openxlsx)
library(clusterProfiler)
library(RColorBrewer)

# R script for LC-MS metabolomics analysis

# 1. Reading Data

# Assume the data format: Rows = Metabolites, Columns = Samples, First column = Metabolite ID

data_file <- "exp_download_express.csv"
raw_data <- read.csv(data_file, stringsAsFactors = FALSE)

# Extract metabolite information and sample data
metabolite_info <- raw_data[, 1:17]
sample_columns <- colnames(raw_data)[18:ncol(raw_data)]
metabolite_data <- as.matrix(raw_data[, 18:ncol(raw_data)])
rownames(metabolite_data) <- raw_data$Metabolite
colnames(metabolite_data) <- sample_columns

# Create sample information table
sample_info <- data.frame(
  Sample = sample_columns,
  stringsAsFactors = FALSE
)

# Assign groups based on sample name prefixes
sample_info$Group <- sapply(sample_info$Sample, function(x) {
  if (grepl("^T", x)) {
    return("Treatment")
  } else if (grepl("^N", x)) {
    return("Control")
  } else if (grepl("^mG", x)) {
    return("mG")
  } else if (grepl("^mS", x)) {
    return("mS")
  } else if (grepl("^MPLC", x)) {
    return("MPLC")
  } else if (grepl("^QC", x)) {
    return("QC")
  } else {
    return("Unknown")
  }
})

sample_info$Group <- factor(sample_info$Group)

# 2. Data Preprocessing

# 2.1 Missing value handling (fill with half of the minimum value)
preprocess_data <- function(data, method = "min_half") {
  if (method == "min_half") {
    for (i in 1:ncol(data)) {
      if (any(is.na(data[, i]))) {
        min_val <- min(data[, i], na.rm = TRUE)
        data[is.na(data[, i]), i] <- min_val / 2
      }
    }
  }
  return(data)
}

# 2.2 Data normalization (Pareto scaling)
# Pareto scaling: (x - mean) / sqrt(sd(x))
pareto_scale <- function(x) {
  x_scaled <- apply(x, 2, function(col) {
    (col - mean(col, na.rm = TRUE)) / sqrt(sd(col, na.rm = TRUE))
  })
  return(x_scaled)
}

# 2.3 Execute preprocessing
metabolite_processed <- preprocess_data(metabolite_data)
metabolite_scaled <- pareto_scale(metabolite_processed)
metabolite_scaled_t <- t(metabolite_scaled)

# 3. Principal Component Analysis (PCA)

# 3.1 Perform PCA analysis using ropls
pca_result <- opls(metabolite_scaled_t,
predI = 2, # PCA mode
scaleC = "none", # Already pre-standardized
fig.pdfC = "none", # Do not generate PDF
info.txtC = "none") # Do not generate text output

# 3.2 Extract PCA results
pca_scores <- as.data.frame(getScoreMN(pca_result))
pca_loadings <- as.data.frame(getLoadingMN(pca_result))
pca_summary <- as.data.frame(getSummaryDF(pca_result))

# 3.3 PCA Visualization

# 3.3.1 Score plot
pca_plot <- function(scores, sample_info, pca_summary, pc_x = 1, pc_y = 2) {
  # Merge scores and sample information
  plot_data <- data.frame(
    PC1 = scores[, pc_x],
    PC2 = scores[, pc_y],
    Group = sample_info$Group,
    Sample = sample_info$Sample
  )

  # Calculate explained variance
  variance_explained <- round(pca_summary[c(pc_x, pc_y), "R2X(cum)"] * 100, 1)

  # Plot PCA
  p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    stat_ellipse(level = 0.95, alpha = 0.2) +
    theme_minimal() +
    labs(
      x = paste0("PC", pc_x, " (", variance_explained[1], "%)"),
      y = paste0("PC", pc_y, " (", variance_explained[2], "%)"),
      title = "PCA Scores Plot",
      color = "Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Set1")

  return(p)
}

# 3.3.2 Generate PCA plot
pca_fig <- pca_plot(pca_scores, sample_info, pca_summary)
print(pca_fig)
ggsave("PCA_scores_plot.png", pca_fig, width = 8, height = 6, dpi = 300)

# 4. Similarity Analysis (ANOSIM)

# 4.0 Remove QC samples from analysis
qc_indices <- which(sample_info$Group == "QC")
if (length(qc_indices) > 0) {
  sample_info <- sample_info[-qc_indices, ]
  metabolite_scaled_t <- metabolite_scaled_t[-qc_indices, ]
  metabolite_scaled <- metabolite_scaled[, -qc_indices]
  cat("Removed", length(qc_indices), "QC samples from analysis\n")
}

# 4.1 Calculate distance matrix (using Euclidean distance)
# Use transposed data to calculate sample-to-sample distances
distance_matrix <- vegdist(metabolite_scaled_t, method = "euclidean")

# 4.2 Perform ANOSIM analysis
anosim_result <- anosim(distance_matrix, sample_info$Group, permutations = 999)

# 4.3 ANOSIM results
cat("\n=== ANOSIM Analysis Results ===\n")
cat("R statistic:", anosim_result$statistic, "\n")
cat("Significance (p-value):", anosim_result$signif, "\n")

# 4.4 Visualize ANOSIM results
anosim_boxplot <- function(anosim_result, sample_info) {
  # Create rank data from ANOSIM result
  # class.vec contains comparison labels (e.g., "Treatment vs Control")
  rank_data <- data.frame(
    Rank = anosim_result$dis.rank,
    Comparison = anosim_result$class.vec
  )
  
  # Convert to factor to maintain order
  rank_data$Comparison <- factor(rank_data$Comparison)
  
  # Plot boxplot showing distance distributions for different comparisons
  p <- ggplot(rank_data, aes(x = Comparison, y = Rank, fill = Comparison)) +
    geom_boxplot(alpha = 0.7) +
    theme_minimal() +
    labs(
      title = paste0("ANOSIM: R = ", round(anosim_result$statistic, 3),
                     ", p = ", round(anosim_result$signif, 4)),
      x = "Group Comparison",
      y = "Rank Distance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}

anosim_fig <- anosim_boxplot(anosim_result, sample_info)
print(anosim_fig)
ggsave("ANOSIM_boxplot.png", anosim_fig, width = 6, height = 6, dpi = 300)

# 5. Orthogonal Partial Least Squares Discriminant Analysis (OPLS-DA)

# 5.1 Perform OPLS-DA analysis
# Note: Need to convert group information to numeric
y_numeric <- as.numeric(sample_info$Group) - 1 # Convert to 0 and 1

# Use transposed data (samples × metabolites) for OPLS-DA
oplsda_result <- opls(metabolite_scaled_t, y_numeric,
                      predI = 1, orthoI = NA, # Automatically select number of orthogonal components
                      scaleC = "none",
                      fig.pdfC = "none",
                      info.txtC = "none")

# 5.2 Extract OPLS-DA results
oplsda_summary <- getSummaryDF(oplsda_result)
oplsda_scores <- as.data.frame(oplsda_result@scoreMN)
oplsda_ortho_scores <- as.data.frame(oplsda_result@orthoScoreMN)

# 5.3 OPLS-DA score plot visualization
oplsda_plot <- function(scores, ortho_scores, sample_info) {
  plot_data <- data.frame(
    t_score = scores[, 1],
    o_score = ortho_scores[, 1],
    Group = sample_info$Group,
    Sample = sample_info$Sample
  )

  p <- ggplot(plot_data, aes(x = t_score, y = o_score, color = Group, label = Sample)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 15) +
    stat_ellipse(level = 0.95, alpha = 0.2) +
    theme_minimal() +
    labs(
      x = paste0("T-score [1] (",
                 round(oplsda_summary["p1", "R2X(cum)"] * 100, 1), "%)"),
      y = paste0("Orthogonal T-score [1] (",
                 round(oplsda_summary["o1", "R2X(cum)"] * 100, 1), "%)"),
      title = "OPLS-DA Scores Plot",
      subtitle = paste0("Model R2Y = ", round(oplsda_summary["sum", "R2Y(cum)"], 3),
                         ", Q2 = ", round(oplsda_summary["sum", "Q2(cum)"], 3)),
      color = "Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Set1") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3)

  return(p)
}

oplsda_fig <- oplsda_plot(oplsda_scores, oplsda_ortho_scores, sample_info)
print(oplsda_fig)
ggsave("OPLSDA_scores_plot.png", oplsda_fig, width = 8, height = 6, dpi = 300)

# 5.4 Model validation (permutation test)
# Note: The ropls package automatically performs permutation tests, and results are included in the summary.
cat("\n=== OPLS-DA Model Statistics ===\n")
print(oplsda_summary)

# 6. Differential Metabolite Screening

# 6.0 Select Treatment and Control groups only
selected_groups <- sample_info$Group %in% c("Treatment", "Control")
sample_info_2group <- sample_info[selected_groups, ]
metabolite_scaled_2group <- metabolite_scaled[, selected_groups]

# 6.1 Get VIP values
vip_values <- getVipVn(oplsda_result)

# Debug: Print dimensions
cat("VIP values length:", length(vip_values), "\n")
cat("Metabolite scaled dimensions:", dim(metabolite_scaled), "\n")
cat("Metabolite scaled_2group dimensions:", dim(metabolite_scaled_2group), "\n")

# 6.2 Calculate statistical test (t-test for two groups)
calculate_t_test <- function(data, groups) {
  p_values <- apply(data, 1, function(x) {
    t.test(x ~ groups, var.equal = TRUE)$p.value
  })
  return(p_values)
}

# Calculate p-values
p_values <- calculate_t_test(metabolite_scaled_2group, sample_info_2group$Group)

# 6.3 Calculate fold change (Fold Change)
calculate_fold_change <- function(data, groups) {
  # Calculate mean of each group
  group_means <- aggregate(t(data), by = list(groups), FUN = mean)
  
  # Calculate Fold Change (Treatment/Control)
  fc <- group_means[group_means$Group.1 == "Treatment", -1] /
       group_means[group_means$Group.1 == "Control", -1]
  
  return(as.numeric(fc))
}

fold_changes <- calculate_fold_change(metabolite_scaled_2group, sample_info_2group$Group)

# 6.4 Create result data frame
# Ensure VIP values match the number of metabolites
if (length(vip_values) != nrow(metabolite_scaled)) {
  # If VIP values are from transposed data, we need to match them
  vip_values <- vip_values[1:nrow(metabolite_scaled)]
}

results_df <- data.frame(
  Metabolite = rownames(metabolite_scaled),
  VIP = vip_values,
  p_value = p_values,
  Fold_Change = fold_changes,
  log2FC = log2(fold_changes)
)

# 6.5 Filter differential metabolites (VIP > 1 and p < 0.0524)
sig_metabolites <- results_df %>%
  filter(VIP > 1 & p_value < 0.0524) %>%
  arrange(desc(abs(log2FC)))

cat("\n=== Differential Metabolite Screening Results ===\n")
cat("Screening criteria: VIP > 1 and p-value < 0.0524\n")
cat("Number of significant differential metabolites:", nrow(sig_metabolites), "\n")

# 6.6 Volcano plot visualization
volcano_plot <- function(results_df, sig_df,
                         vip_threshold = 1,
                         p_threshold = 0.0524) {
  # Mark significant points
  plot_data <- results_df %>%
    mutate(
      Significance = case_when(
        VIP > vip_threshold & p_value < p_threshold & log2FC > 0 ~ "Up-regulated",
        VIP > vip_threshold & p_value < p_threshold & log2FC < 0 ~ "Down-regulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Volcano plot
  p <- ggplot(plot_data, aes(x = log2FC, y = -log10(p_value),
                             color = Significance)) +
    geom_point(alpha = 0.6, size = 2) +
    geom_hline(yintercept = -log10(p_threshold),
               linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    theme_minimal() +
    labs(
      title = "Volcano Plot of Differential Metabolites",
      x = expression(log[2]("Fold Change")),
      y = expression(-log[10](p_value)),
      color = "Regulation"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "right"
    ) +
    scale_color_manual(values = c("Up-regulated" = "red",
                                  "Down-regulated" = "blue",
                                  "Not significant" = "gray")) +
    geom_text_repel(
      data = subset(plot_data, Significance != "Not significant" &
                      abs(log2FC) > 1),
      aes(label = Metabolite),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5
    )
  
  return(p)
}

volcano_fig <- volcano_plot(results_df, sig_metabolites)
print(volcano_fig)
ggsave("volcano_plot.png", volcano_fig, width = 10, height = 8, dpi = 300)

# 6.7 Differential metabolite heatmap
if (nrow(sig_metabolites) > 0) {
  # Check which metabolites exist in the data
  valid_metabolites <- sig_metabolites$Metabolite[sig_metabolites$Metabolite %in% rownames(metabolite_scaled_2group)]
  
  if (length(valid_metabolites) > 0) {
    # Extract data for significant metabolites (Treatment and Control only)
    sig_data <- metabolite_scaled_2group[valid_metabolites, ]
    
    # Create annotation information
    annotation_col <- data.frame(
      Group = sample_info_2group$Group
    )
    rownames(annotation_col) <- colnames(sig_data)
    
    # Draw heatmap
    pheatmap(
      sig_data,
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "complete",
      annotation_col = annotation_col,
      show_colnames = TRUE,
      show_rownames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
      main = paste("Heatmap of Differential Metabolites\n(Treatment vs Control, n =",
                   length(valid_metabolites), ")"),
      filename = "differential_metabolites_heatmap.png",
      width = 10,
      height = 12
    )
  } else {
    cat("No valid metabolites found for heatmap\n")
  }
}


# 7. KEGG Pathway Enrichment Analysis

# Note: Requires KEGG IDs of metabolites.

# 7.1 Assume we have a metabolite-KEGG ID mapping table
# Create an example mapping table here
n_metabolites <- nrow(metabolite_scaled)
set.seed(123)
kegg_mapping <- data.frame(
  Metabolite = rownames(metabolite_scaled),
  KEGG_ID = paste0("C", sprintf("%05d", sample(1:20000, n_metabolites)))
)

# 7.2 Get KEGG IDs of significant metabolites
sig_kegg_ids <- kegg_mapping %>%
  filter(Metabolite %in% sig_metabolites$Metabolite) %>%
  pull(KEGG_ID)

# 7.3 Perform KEGG pathway enrichment analysis
if (length(sig_kegg_ids) > 0) {
  # Set network timeout to 6 minutes (360 seconds)
  options(timeout = 360)
  
  # Debug: Print KEGG IDs
  cat("Number of KEGG IDs:", length(sig_kegg_ids), "\n")
  cat("First 5 KEGG IDs:", head(sig_kegg_ids, 5), "\n")
  
  # Note: Using human (hsa) as an example, adjust according to actual species.
  # Use enrichMKEGG for metabolite pathway enrichment
  tryCatch({
    kegg_enrich <- enrichKEGG(
      sig_kegg_ids,
      organism = "cpd",  # Human, other species such as "mmu" (mouse), "rno" (rat)
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      qvalueCutoff = 1,
    )
    
    # 7.4 Extract enrichment results
    kegg_results <- as.data.frame(kegg_enrich)
    
    # 7.5 Visualize enrichment results
    if (nrow(kegg_results) > 0) {
      # Bubble chart
      dotplot(kegg_enrich,
              showCategory = 20,
              title = "KEGG Pathway Enrichment Analysis",
              font.size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave("kegg_enrichment_dotplot.png", width = 10, height = 8, dpi = 300)
      
      # Bar chart
      barplot(kegg_enrich,
              showCategory = 15,
              title = "KEGG Pathway Enrichment",
              font.size = 10) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      
      ggsave("kegg_enrichment_barplot.png", width = 10, height = 8, dpi = 300)
      
      cat("\n=== KEGG Pathway Enrichment Analysis Results ===\n")
      cat("Number of significantly enriched pathways:", nrow(kegg_results), "\n")
      print(head(kegg_results, 10))
    }
  }, error = function(e) {
    cat("Error in KEGG enrichment:", e$message, "\n")
    cat("Skipping KEGG enrichment analysis...\n")
    kegg_enrich <<- NULL
  })
}

# 8. Result Saving

# 8.1 Save all results to an Excel file
write_excel_results <- function(results, filename, metabolite_processed) {
  wb <- createWorkbook()
  
  addWorksheet(wb, "PCA_Summary")
  writeData(wb, "PCA_Summary", pca_summary, rowNames = TRUE)
  
  addWorksheet(wb, "OPLSDA_Summary")
  writeData(wb, "OPLSDA_Summary", oplsda_summary, rowNames = TRUE)
  
  addWorksheet(wb, "ANOSIM_Results")
  anosim_summary <- data.frame(
    Statistic = anosim_result$statistic,
    P_value = anosim_result$signif,
    Permutations = length(anosim_result$perm)
  )
  writeData(wb, "ANOSIM_Results", anosim_summary)
  
  addWorksheet(wb, "Metabolite_Processed_Data")
  writeData(wb, "Metabolite_Processed_Data", as.data.frame(metabolite_processed), rowNames = TRUE)
  
  addWorksheet(wb, "All_Metabolites_Results")
  writeData(wb, "All_Metabolites_Results", results_df)
  
  addWorksheet(wb, "Significant_Metabolites")
  writeData(wb, "Significant_Metabolites", sig_metabolites)
  
  if (exists("kegg_results") && nrow(kegg_results) > 0) {
    addWorksheet(wb, "KEGG_Enrichment")
    writeData(wb, "KEGG_Enrichment", kegg_results)
  }
  
  saveWorkbook(wb, file = filename, overwrite = TRUE)
  cat("\nResults saved to:", filename, "\n")
}

write_excel_results(results_df, "LCMS_Metabolomics_Analysis_Results.xlsx", metabolite_processed)

# 8.2 Save R data
saveRDS(list(
  metabolite_data = metabolite_scaled,
  sample_info = sample_info,
  pca_result = pca_result,
  anosim_result = anosim_result,
  oplsda_result = oplsda_result,
  all_results = results_df,
  sig_metabolites = sig_metabolites,
  kegg_enrichment = if (exists("kegg_enrich")) kegg_enrich else NULL
), file = "LCMS_analysis_results.rds")