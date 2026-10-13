# Locate the project root (parent of scripts/) regardless of how the script is run:
#   - Rscript: the script path is passed via --file=
#   - RStudio (Source/Run): ask rstudioapi for the active document's path
#   - otherwise: fall back to the working directory
args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
if (!(length(script_file) == 1 && nzchar(script_file)) &&
    requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  script_file <- rstudioapi::getActiveDocumentContext()$path
}
project_root <- if (length(script_file) == 1 && nzchar(script_file)) {
  normalizePath(file.path(dirname(script_file), ".."))
} else {
  getwd()
}
data_dir <- file.path(project_root, "data")
out_dir  <- file.path(project_root, "results")
if (!dir.exists(data_dir)) {
  stop("Cannot find the 'data' folder (looked under: ", data_dir,
       "). Set the working directory to the project root, or run via Rscript.")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load required packages
library(pheatmap)   # heatmap plotting
library(tidyverse)  # data manipulation
library(Cairo)      # high-quality PDF output

# ==================== Read data ====================
data <- read.csv(file.path(data_dir, "01_human_phospho_5_delete_con_at_least1.csv"),
                 fileEncoding = "UTF-8-BOM")

# Set rownames to protein identifiers (first column)
rownames(data) <- data[, 1]

# Keep expression data (columns 11 onward)
expression_data <- data[, 11:ncol(data)]

# Convert "NaN" to NA
expression_data[expression_data == "NaN"] <- NA

# Convert to numeric
expression_data <- as.data.frame(
  lapply(expression_data, function(x) as.numeric(as.character(x)))
)

# ==================== Missing value imputation ====================
# Fill NA with tiny random values near zero
fill_with_small_random <- function(data, min_val = 1e-10, max_val = 1e-8) {
  data_filled <- apply(data, 2, function(col) {
    col[is.na(col)] <- runif(sum(is.na(col)), min = min_val, max = max_val)
    return(col)
  })
  return(as.data.frame(data_filled))
}

filled_data <- fill_with_small_random(expression_data)

# ==================== Z-score standardization ====================
scaled_data <- t(scale(t(filled_data)))

# ==================== Define function: order clusters by specified columns ====================
get_clustered_row_order <- function(scaled_data, k = NULL, h = NULL, focus_cols = NULL) {
  # Hierarchical clustering
  hc <- hclust(dist(scaled_data), method = "ward.D2")

  # Cut the tree
  if (!is.null(k)) {
    cluster_assign <- cutree(hc, k = k)
  } else if (!is.null(h)) {
    cluster_assign <- cutree(hc, h = h)
  } else {
    stop("Please specify at least one of k (number of clusters) or h (tree height)")
  }

  # Default to all columns; if focus_cols is given, use only those
  if (is.null(focus_cols)) {
    focus_cols <- 1:ncol(scaled_data)
  }

  # Mean of each cluster (based on focus_cols)
  cluster_means <- tapply(
    rowMeans(scaled_data[, focus_cols, drop = FALSE]),
    cluster_assign, mean
  )

  # Order clusters
  ordered_clusters <- order(cluster_means)

  # Row order
  row_order <- order(match(cluster_assign, ordered_clusters))

  return(list(
    hc = hc,
    cluster_assign = cluster_assign,
    cluster_means = cluster_means,
    row_order = row_order
  ))
}

# ==================== Cluster and order ====================
# Example: cut into 5 clusters, order by the first 3 sample columns
res <- get_clustered_row_order(scaled_data, k = 5, focus_cols = 1:3)

# ==================== Plot heatmap ====================
heatmap_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

p2 <- pheatmap(
  scaled_data[res$row_order, ],
  color = heatmap_colors,
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = FALSE,   # row order already fixed
  cluster_cols = FALSE,   # whether to cluster columns can be set to FALSE
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  fontsize_col = 8,
  border_color = NA,
  main = "Phosphorylated Protein Quantity Heatmap",
  labels_col = c("Oocytes-1", "Oocytes-2", "Oocytes-3",
                 "1C-2", "1C-3", "4-6C-1", "4-6C-2", "4-6C-3",
                 "8C-1", "8C-2", "8C-3", "Morula-1", "Morula-2", "Morula-3",
                 "Blastocyst-1", "Blastocyst-2", "Blastocyst-3")
)

# ==================== Save PDF ====================
#install.packages("Cairo")
library(Cairo)
CairoPDF(file.path(out_dir, "01_phospho_5_delete_con&1_atleast1_heatmap.1.pdf"),
         width = 12, height = 10, pointsize = 12)
print(p2)
dev.off()

# # ==================== Export corresponding matrix (optional) ====================
# # Standardized matrix (in heatmap row order)
# heatmap_matrix <- scaled_data[res$row_order, ]
# write.csv(as.data.frame(heatmap_matrix),
#           file.path(out_dir, "01_phospho_5_delete_con&1_atleast1_heatmap_matrix.csv"),
#           row.names = TRUE)
#
# # Raw imputed matrix (in heatmap row order)
# heatmap_matrix_raw <- filled_data[res$row_order, ]
# write.csv(as.data.frame(heatmap_matrix_raw),
#           file.path(out_dir, "01_phospho_5_delete_con&1_atleast1_heatmap_raw.csv"),
#           row.names = TRUE)

# ==================== Export matrices with annotation info ====================

# Extract the first 10 annotation columns
annotation_info <- data[, 1:10]

# Standardized matrix (in heatmap order)
heatmap_matrix <- scaled_data[res$row_order, ]
heatmap_matrix_df <- as.data.frame(heatmap_matrix)

# Raw imputed matrix (in heatmap order)
heatmap_matrix_raw <- filled_data[res$row_order, ]
heatmap_matrix_raw_df <- as.data.frame(heatmap_matrix_raw)

# Matching annotation info (same row order)
annotation_info_ordered <- annotation_info[res$row_order, ]

# Combine annotation + data matrix
heatmap_matrix_out <- cbind(annotation_info_ordered, heatmap_matrix_df)
heatmap_matrix_raw_out <- cbind(annotation_info_ordered, heatmap_matrix_raw_df)

# Write CSV
write.csv(heatmap_matrix_out,
          file.path(out_dir, "20250901-phospho-5-delete_con&1_atleast1_heatmap_matrix.csv"),
          row.names = FALSE)

write.csv(heatmap_matrix_raw_out,
          file.path(out_dir, "20250901-phospho-5-delete_con&1_atleast1_heatmap_raw.csv"),
          row.names = FALSE)
