## =========================
## DEG heatmap (P < 0.05 & log2FC > 1)
## =========================

## -------- User settings (EDIT HERE) --------
deg_file   <- "deg_result.csv"         # differential result table (CSV/TSV)
expr_file  <- "expr_matrix.csv"        # expression matrix (CSV/TSV): rows=genes, cols=samples
sample_file <- "sample_info.csv"       # optional: sample annotation (CSV/TSV) with columns: Sample, Group
outdir     <- "DEG_heatmap_P0.05_log2FC1"
p_cutoff   <- 0.05
lfc_cutoff <- 1
top_n      <- NA   # set e.g. 50 for top 50 genes by log2FC; NA means keep all
log2_expr  <- TRUE # if expression matrix is not log2-scaled, keep TRUE to log2(x+1)
## ------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## -------- Helpers --------
read_any <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv")) {
    df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    stop("Unsupported file extension: ", ext, " (use .csv/.tsv/.txt)")
  }
  df
}

pick_col <- function(df, candidates_regex, fallback_first = FALSE) {
  cn <- colnames(df)
  hit <- cn[grepl(candidates_regex, cn, ignore.case = TRUE)]
  if (length(hit) >= 1) return(hit[1])
  if (fallback_first) return(cn[1])
  return(NA_character_)
}

## -------- Load packages --------
suppressPackageStartupMessages({
  library(dplyr)
  library(pheatmap)
})

## =========================
## 1) Read DEG table & filter genes
## =========================
deg <- read_any(deg_file)

gene_col <- pick_col(deg, "^(gene[_ ]?id|geneid|id)$", fallback_first = TRUE)
p_col    <- pick_col(deg, "^(pvalue|p[_ ]?value|pval|p\\.value)$")
lfc_col  <- pick_col(deg, "log2fc|log2foldchange|log2_foldchange|log2\\s*fc")

if (is.na(p_col))   stop("Cannot find P-value column. Expected something like: Pvalue / pvalue / p.value")
if (is.na(lfc_col)) stop("Cannot find log2FC column. Expected something like: Log2FC / log2FoldChange")

deg2 <- deg %>%
  mutate(
    .gene = as.character(.data[[gene_col]]),
    .p    = suppressWarnings(as.numeric(.data[[p_col]])),
    .lfc  = suppressWarnings(as.numeric(.data[[lfc_col]]))
  ) %>%
  filter(!is.na(.gene), .gene != "", !is.na(.p), !is.na(.lfc)) %>%
  filter(.p < p_cutoff, .lfc > lfc_cutoff) %>%
  arrange(desc(.lfc))

if (nrow(deg2) == 0) {
  stop("No genes passed filters: P < ", p_cutoff, " and log2FC > ", lfc_cutoff)
}

if (!is.na(top_n)) {
  deg2 <- deg2 %>% slice_head(n = top_n)
}

write.csv(deg2, file.path(outdir, "DEG_filtered_P_lt_0.05_log2FC_gt_1.csv"), row.names = FALSE)

sel_genes <- deg2$.gene

## =========================
## 2) Read expression matrix & subset
## =========================
expr <- read_any(expr_file)

expr_gene_col <- pick_col(expr, "^(gene[_ ]?id|geneid|id)$", fallback_first = TRUE)

expr_mat <- expr
rownames(expr_mat) <- make.unique(as.character(expr_mat[[expr_gene_col]]))
expr_mat[[expr_gene_col]] <- NULL

## Make numeric matrix
expr_mat <- as.matrix(expr_mat)
mode(expr_mat) <- "numeric"

## Keep only genes in DEG list and in expression matrix
keep_genes <- intersect(sel_genes, rownames(expr_mat))
if (length(keep_genes) == 0) stop("No overlap between DEG genes and expression matrix rownames.")

expr_sub <- expr_mat[keep_genes, , drop = FALSE]

## Reorder rows by DEG log2FC (descending)
expr_sub <- expr_sub[match(deg2$.gene, rownames(expr_sub), nomatch = 0), , drop = FALSE]
expr_sub <- expr_sub[rownames(expr_sub) != "", , drop = FALSE]

## Optional log2 transform
if (isTRUE(log2_expr)) {
  expr_sub <- log2(expr_sub + 1)
}

## =========================
## 3) Optional sample annotation (Group)
## =========================
anno_col <- NULL
if (!is.null(sample_file) && file.exists(sample_file)) {
  si <- read_any(sample_file)
  if (!all(c("Sample", "Group") %in% colnames(si))) {
    stop("sample_info must contain columns: Sample, Group")
  }
  si <- si %>% mutate(Sample = as.character(Sample), Group = as.character(Group))
  si <- si[match(colnames(expr_sub), si$Sample), , drop = FALSE]
  if (any(is.na(si$Sample))) {
    stop("Some samples in expression matrix are missing in sample_info.csv")
  }
  rownames(si) <- si$Sample
  anno_col <- data.frame(Group = si$Group, row.names = si$Sample, stringsAsFactors = FALSE)
}

## =========================
## 4) Plot & export heatmap
## =========================
## Figure size heuristic
n_genes <- nrow(expr_sub)
n_samp  <- ncol(expr_sub)
png_w <- max(1600, 250 + n_samp * 120)
png_h <- max(1200, 250 + n_genes * 18)

## Draw heatmap (row Z-score)
png(file.path(outdir, "Heatmap_DEG_P_lt_0.05_log2FC_gt_1.png"),
    width = png_w, height = png_h, res = 200)
pheatmap(
  expr_sub,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  annotation_col = anno_col,
  show_rownames = (n_genes <= 80),
  fontsize_row = 7,
  border_color = NA
)
dev.off()

pdf(file.path(outdir, "Heatmap_DEG_P_lt_0.05_log2FC_gt_1.pdf"),
    width = max(8, 2 + n_samp * 0.35),
    height = max(6, 2 + n_genes * 0.12))
pheatmap(
  expr_sub,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "average",
  annotation_col = anno_col,
  show_rownames = (n_genes <= 80),
  fontsize_row = 7,
  border_color = NA
)
dev.off()

message("Done. Outputs saved in: ", normalizePath(outdir, winslash = "/"))
message("Filtered gene list: ", file.path(outdir, "DEG_filtered_P_lt_0.05_log2FC_gt_1.csv"))
message("Heatmap PNG/PDF exported.")
