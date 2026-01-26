expr_raw <- read.csv(file_path,
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

library(dplyr)
library(tibble)

expr_matrix <- expr_raw %>%
  rownames_to_column("symbol") %>%
  group_by(symbol) %>%
  summarise(across(where(is.numeric), max)) %>%
  column_to_rownames("symbol")

samples <- colnames(expr_matrix)
group <- ifelse(grepl("01A", samples), "Tumor",
                ifelse(grepl("11A", samples), "Normal", NA))
keep <- !is.na(group)
expr_matrix <- expr_matrix[, keep]
group <- group[keep]
group <- factor(group, levels = c("Normal", "Tumor"))

library(edgeR)

dge <- DGEList(counts = expr_matrix, group = group)
dge <- calcNormFactors(dge)
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)

res <- topTags(lrt, n = Inf)$table
res$FC <- 2^res$logFC

res_filtered <- res %>%
  filter(abs(logFC) >= 1, FDR < 0.05)

write.csv(res_filtered, "TCGA_LUAD_DEG_result_edgeR.csv", row.names = TRUE)
