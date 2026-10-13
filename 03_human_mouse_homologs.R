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
library(dplyr)
library(stringdist)

# 1. Read data
human <- read.csv(file.path(data_dir, "01_human_phospho_5_delete_con_at_least1.csv"),
                  stringsAsFactors = FALSE)
mouse <- read.csv(file.path(data_dir, "03_mouse_phospho.csv"),
                  stringsAsFactors = FALSE)

# 2. Trim human peptides (15 aa -> 13 aa)
human$Peptide_trimmed <- substr(human$PTM.FlankingRegion, 2,
                                nchar(human$PTM.FlankingRegion) - 1)

# 3. Generate human-mouse peptide combinations
pairs <- expand.grid(
  Mouse = mouse$Peptide,
  Human = human$Peptide_trimmed,
  stringsAsFactors = FALSE
)

# 4. Compute pairwise distances
pairs$diff <- stringdist(pairs$Mouse, pairs$Human, method = "hamming")

# 5. Keep homolog pairs with distance <= 4
homologs <- pairs %>% filter(diff <= 4)

# Inspect the data structure to be used downstream
table(nchar(human$Peptide_trimmed))
table(nchar(mouse$Peptide))
colnames(mouse)

# 6. Merge human table + mouse SYMBOL, and reorder columns
homologs_info <- homologs %>%
  left_join(human %>% select(PTM.FlankingRegion, Peptide_trimmed, everything()),
            by = c("Human" = "Peptide_trimmed")) %>%
  left_join(mouse %>% select(Peptide, SYMBOL),
            by = c("Mouse" = "Peptide")) %>%
  select(SYMBOL, everything()) %>%   # put SYMBOL first
  arrange(diff)                      # sort ascending by distance

# 7. Preview first 10 rows
head(homologs_info, 10)

# 8. Write result
write.csv(homologs_info, file.path(out_dir, "03_blast.csv"), row.names = FALSE)
