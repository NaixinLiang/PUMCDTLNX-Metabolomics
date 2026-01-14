if (!require("vegan", quietly = TRUE)) {
  install.packages("vegan", repos = "https://cloud.r-project.org")
}
if (!require("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2", repos = "https://cloud.r-project.org")
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}
if (!require("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel", repos = "https://cloud.r-project.org")
}
if (!require("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer", repos = "https://cloud.r-project.org")
}
if (!require("boot", quietly = TRUE)) {
  install.packages("boot", repos = "https://cloud.r-project.org")
}

library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(boot)

cat("=== 16S Data Analysis: PCoA and ANOSIM ===\n\n")

data_file <- "ASV_Taxon_Depth_asv.full.xls"
raw_data <- read.delim(data_file, sep = "\t", stringsAsFactors = FALSE)

cat("Data loaded successfully. Dimensions:", nrow(raw_data), "rows x", ncol(raw_data), "columns\n")

sample_columns <- colnames(raw_data)[10:(ncol(raw_data)-3)]

cat("Total samples found:", length(sample_columns), "\n")
cat("Samples:", paste(sample_columns, collapse = ", "), "\n")

sample_info <- data.frame(
  Sample = sample_columns,
  stringsAsFactors = FALSE
)

sample_info$Group <- sapply(sample_info$Sample, function(x) {
  if (grepl("^GGO", x)) {
    return("GGO")
  } else if (grepl("^mG", x)) {
    return("mG")
  } else if (grepl("^mS", x)) {
    return("mS")
  } else if (grepl("^S", x)) {
    return("S")
  } else {
    return("Unknown")
  }
})

sample_info$Group <- factor(sample_info$Group)

cat("\nSample groups:\n")
print(table(sample_info$Group))

genus_data <- raw_data %>%
  filter(!is.na(genus) & genus != "unclassified_c__Alphaproteobacteria" & 
         genus != "unclassified_p__Pseudomonadota" & 
         genus != "unclassified_k__norank_d__Bacteria" &
         genus != "unclassified_p__Pseudomonadota" &
         genus != "unclassified_c__Alphaproteobacteria" &
         genus != "unclassified_c__Alphaproteobacteria" &
         !grepl("^unclassified", genus))

genus_abundance <- genus_data %>%
  group_by(genus) %>%
  summarise(
    across(all_of(sample_columns), sum, .names = "{.col}"),
    .groups = "drop"
  )

abundance_matrix <- as.matrix(genus_abundance[, sample_columns])
rownames(abundance_matrix) <- genus_abundance$genus

cat("\nGenus abundance matrix created:", nrow(abundance_matrix), "genera x", ncol(abundance_matrix), "samples\n")

abundance_relative <- sweep(abundance_matrix, 2, colSums(abundance_matrix), "/") * 100

distance_matrix <- vegdist(t(abundance_relative), method = "bray")

cat("\nDistance matrix calculated using Bray-Curtis method\n")

pcoa_result <- cmdscale(distance_matrix, k = 2, eig = TRUE)

pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("Axis1", "Axis2")
pcoa_scores$Sample <- sample_info$Sample
pcoa_scores$Group <- sample_info$Group

cat("\n=== PCoA Results ===\n")
cat("Eigenvalues:\n")
print(round(pcoa_result$eig[1:10], 4))
cat("\nVariance explained by first two axes:", 
    round(sum(abs(pcoa_result$eig[1:2])) / sum(abs(pcoa_result$eig)) * 100, 2), "%\n")

cat("\n=== ANOSIM Results ===\n")
group_counts <- table(sample_info$Group)
cat("Sample counts per group:\n")
print(group_counts)

if (all(group_counts == 1)) {
  cat("\nWARNING: Each group has only one sample. ANOSIM requires replicates within groups.\n")
  cat("Skipping ANOSIM analysis.\n")
  anosim_result <- NULL
  anosim_plot <- NULL
} else {
  anosim_result <- anosim(distance_matrix, sample_info$Group, permutations = 999)
  
  cat("R statistic:", round(anosim_result$statistic, 4), "\n")
  cat("Significance (p-value):", round(anosim_result$signif, 4), "\n")
  cat("Number of permutations:", length(anosim_result$perm), "\n")
  
  if (anosim_result$signif < 0.05) {
    cat("Result: Significant differences between groups (p < 0.05)\n")
  } else {
    cat("Result: No significant differences between groups (p >= 0.05)\n")
  }
}

variance_explained <- round(abs(pcoa_result$eig[1:2]) / sum(abs(pcoa_result$eig)) * 100, 2)

pcoa_plot <- ggplot(pcoa_scores, aes(x = Axis1, y = Axis2, color = Group, label = Sample)) +
  geom_point(size = 6, alpha = 0.8) +
  geom_text_repel(size = 4, max.overlaps = 20, box.padding = 0.5) +
  stat_ellipse(level = 0.95, alpha = 0.2, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  labs(
    title = "PCoA Analysis (Bray-Curtis Distance)",
    subtitle = paste0("Variance explained: Axis1 = ", variance_explained[1], "%, Axis2 = ", variance_explained[2], "%"),
    x = paste0("Axis 1 (", variance_explained[1], "%)"),
    y = paste0("Axis 2 (", variance_explained[2], "%)"),
    color = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1") +
  geom_hline(yintercept = 0, linetype = "solid", alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", alpha = 0.3)

print(pcoa_plot)
ggsave("16S_PCoA_bray_curtis.png", pcoa_plot, width = 10, height = 8, dpi = 300)

if (!is.null(anosim_result)) {
  anosim_data <- data.frame(
    Rank = anosim_result$dis.rank,
    Comparison = anosim_result$class.vec
  )
  
  anosim_data$Comparison <- factor(anosim_data$Comparison)
  
  anosim_plot <- ggplot(anosim_data, aes(x = Comparison, y = Rank, fill = Comparison)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("ANOSIM Analysis\nR = ", round(anosim_result$statistic, 3), 
                     ", p = ", round(anosim_result$signif, 4)),
      x = "Group Comparison",
      y = "Rank Distance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  print(anosim_plot)
  ggsave("16S_ANOSIM_bray_curtis.png", anosim_plot, width = 8, height = 6, dpi = 300)
}

cat("\n=== Wilcoxon Rank Sum Test: S vs GGO ===\n")

S_samples <- sample_info$Sample[sample_info$Group == "S"]
GGO_samples <- sample_info$Sample[sample_info$Group == "GGO"]

cat("S group samples:", paste(S_samples, collapse = ", "), "\n")
cat("GGO group samples:", paste(GGO_samples, collapse = ", "), "\n")

S_abundance <- abundance_relative[, S_samples]
GGO_abundance <- abundance_relative[, GGO_samples]

wilcoxon_results <- data.frame(
  Genus = rownames(abundance_relative),
  S_mean = rowMeans(S_abundance),
  GGO_mean = rowMeans(GGO_abundance),
  S_median = apply(S_abundance, 1, median),
  GGO_median = apply(GGO_abundance, 1, median),
  stringsAsFactors = FALSE
)

wilcoxon_results$W_statistic <- NA
wilcoxon_results$p_value <- NA
wilcoxon_results$CI_lower <- NA
wilcoxon_results$CI_upper <- NA

for (i in 1:nrow(wilcoxon_results)) {
  S_data <- as.numeric(S_abundance[i, ])
  GGO_data <- as.numeric(GGO_abundance[i, ])
  
  wilcox_test <- tryCatch({
    wilcox.test(S_data, GGO_data, alternative = "two.sided", exact = FALSE)
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(wilcox_test)) {
    wilcoxon_results$W_statistic[i] <- wilcox_test$statistic
    wilcoxon_results$p_value[i] <- wilcox_test$p.value
    
    combined_data <- c(S_data, GGO_data)
    n_S <- length(S_data)
    
    if (length(unique(combined_data)) > 1 && n_S > 1 && length(GGO_data) > 1) {
      boot_fn <- function(data, idx) {
        S_boot <- data[idx[1:n_S]]
        GGO_boot <- data[idx[(n_S+1):length(idx)]]
        return(mean(S_boot) - mean(GGO_boot))
      }
      
      boot_result <- tryCatch({
        boot(data = combined_data, statistic = boot_fn, R = 9999)
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(boot_result)) {
        ci_result <- tryCatch({
          boot.ci(boot_result, type = "perc", conf = 0.95)
        }, error = function(e) {
          return(NULL)
        })
        
        if (!is.null(ci_result)) {
          wilcoxon_results$CI_lower[i] <- ci_result$percent[4]
          wilcoxon_results$CI_upper[i] <- ci_result$percent[5]
        }
      }
    }
  }
}

wilcoxon_results$adj_p_value <- p.adjust(wilcoxon_results$p_value, method = "BH")
wilcoxon_results$significant <- wilcoxon_results$adj_p_value < 0.05

unclassified_mask <- !grepl("unclassified", wilcoxon_results$Genus, ignore.case = TRUE)
wilcoxon_results_filtered <- wilcoxon_results[unclassified_mask, ]

cat("\nWilcoxon test completed for", nrow(wilcoxon_results_filtered), "genera (excluding unclassified)\n")
cat("Significant genera:", sum(wilcoxon_results_filtered$significant), "\n")

write.csv(wilcoxon_results_filtered, "16S_Wilcoxon_S_vs_GGO.csv", row.names = FALSE)

cat("\n=== Results Summary ===\n")

results_summary <- data.frame(
  Metric = c("ANOSIM R statistic", "ANOSIM p-value", "Number of permutations",
             "Variance explained by Axis1 (%)", "Variance explained by Axis2 (%)",
             "Total variance explained by first 2 axes (%)",
             "Number of genera tested", "Number of significant genera (S vs GGO)"),
  Value = c(if (!is.null(anosim_result)) round(anosim_result$statistic, 4) else NA,
            if (!is.null(anosim_result)) round(anosim_result$signif, 4) else NA,
            if (!is.null(anosim_result)) length(anosim_result$perm) else NA,
            variance_explained[1],
            variance_explained[2],
            sum(variance_explained),
            nrow(wilcoxon_results_filtered),
            sum(wilcoxon_results_filtered$significant))
)

write.csv(results_summary, "16S_Analysis_Summary.csv", row.names = FALSE)

write.csv(pcoa_scores, "16S_PCoA_Scores.csv", row.names = FALSE)

genus_summary <- as.data.frame(abundance_relative)
genus_summary$Genus <- rownames(abundance_relative)
genus_summary <- genus_summary[, c("Genus", sample_columns)]
genus_summary[, sample_columns] <- lapply(genus_summary[, sample_columns], round, 2)
write.csv(genus_summary, "16S_Genus_Abundance.csv", row.names = FALSE)

cat("\n=== Analysis Complete ===\n")
cat("Output files:\n")
cat("  - 16S_PCoA_bray_curtis.png (PCoA plot)\n")
cat("  - 16S_ANOSIM_bray_curtis.png (ANOSIM plot)\n")
cat("  - 16S_Analysis_Summary.csv (Analysis summary)\n")
cat("  - 16S_PCoA_Scores.csv (PCoA coordinates)\n")
cat("  - 16S_Genus_Abundance.csv (Genus abundance data)\n")
cat("  - 16S_Wilcoxon_S_vs_GGO.csv (Wilcoxon test results: S vs GGO)\n")
