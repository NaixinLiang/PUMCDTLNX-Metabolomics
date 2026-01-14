library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

countData <- read.table("gene.count.matrix.xls", header = TRUE, sep = "\t")
rownames(countData) <- countData[, 1]
countData <- countData[, -1]
countData <- round(countData)

sampleNames <- colnames(countData)

colData <- data.frame(
  sample = sampleNames,
  group = factor(c(rep("GGO", 17), rep("Solid", 16)),
                 levels = c("GGO", "Solid")),
  row.names = sampleNames
)

all(colnames(countData) == rownames(colData))

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ group)

dds$group <- relevel(dds$group, ref = "GGO")

dds <- DESeq(dds)

res <- results(dds,
               contrast = c("group", "Solid", "GGO"),
               alpha = 0.05)

write.csv(as.data.frame(res),
          file = "DESeq2_GGO_vs_Solid_results.csv")

sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 2)]

vsd <- vst(dds, blind = TRUE)

expr <- assay(vsd)[sig_genes, ]

annotation_col <- data.frame(
  Group = colData(dds)$group
)
rownames(annotation_col) <- colnames(expr)

ann_colors <- list(
  Group = c(GGO = "#1F77B4", Solid = "#FF7F0E")
)

pheatmap(
  expr,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  breaks = seq(-3, 3, length.out = 101),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  treeheight_row = 100,
  treeheight_col = 50,
  main = "Heatmap of Differentially Expressed Genes\n(Solid vs GGO, padj<0.05 & |log2FC|>1)",
  filename = "DEGs_heatmap_vst_scaled.pdf",
  width = 10,
  height = 12
)
