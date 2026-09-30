# =============================================================================
# Supplementary Figure 7d preprocessing 2 - miRNA differential expression (edgeR)
# =============================================================================
#
# Description: Runs edgeR quasi-likelihood F-test differential expression on
#   cleaned miRNA CPM counts for the 4 ADE groups, using the same 6 pairwise
#   contrast design as the proteomics DEA. edgeR is used rather than DESeq2
#   or limma: a field benchmark of miRNA-seq DE methods found DESeq2 and
#   limma-based approaches tend to systematically underestimate log2FC for
#   miRNA-seq specifically, while edgeR performs comparably to purpose-built
#   miRNA-seq tools.
#
# Prerequisites: supp_figure_7_output_data/supp_figure_7d_wide_counts.csv,
#   supp_figure_7_output_data/supp_figure_7d_group_id_mapping.csv (from preprocessing step 1)
#
# Output: supp_figure_7d_edgeR_results_*.csv, supp_figure_7d_DEA_significant_counts_summary.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "supp_figure_7") {
  if (file.exists("supp_figure_7")) setwd("supp_figure_7")
  else if (file.exists(file.path("..", "supp_figure_7"))) setwd(file.path("..", "supp_figure_7"))
}

OUTPUT_DIR <- "supp_figure_7_output_data"
counts_file <- file.path(OUTPUT_DIR, "supp_figure_7d_wide_counts.csv")
mapping_file <- file.path(OUTPUT_DIR, "supp_figure_7d_group_id_mapping.csv")
if (!file.exists(counts_file)) {
  stop("Run supp_figure_7d_preprocessing_1_data_cleaning.R first. Expected: ", counts_file)
}

library(tidyverse)
library(limma)
library(edgeR)

mature_normalized_miRNA_R4A <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
group_id_df <- read.csv(mapping_file, stringsAsFactors = FALSE)

# Set miRNA as rownames, build a miRNA x sample count matrix
counts_matrix <- mature_normalized_miRNA_R4A %>%
  column_to_rownames("miRNA") %>%
  as.matrix()

# Create group factor that matches exactly with the counts matrix
group_vector <- setNames(group_id_df$Group, group_id_df$ID)
group_factor <- factor(group_vector[colnames(counts_matrix)])

missing_ids <- colnames(counts_matrix)[is.na(group_factor)]
if (length(missing_ids) > 0) {
  stop("Sample columns with no group mapping: ", paste(missing_ids, collapse = ", "))
}

print("Dimensions of counts matrix:")
print(dim(counts_matrix))
print("Group factor:")
print(table(group_factor))

# Create DGEList object
dge <- DGEList(counts = counts_matrix, group = group_factor)

# Filter low expressed miRNAs
keep <- filterByExpr(dge)
cat("miRNAs retained after expression filtering:", sum(keep), "of", length(keep), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Set up design matrix
design <- model.matrix(~0 + group_factor)
colnames(design) <- levels(group_factor)

# Estimate dispersion and fit the model
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Define contrasts: 6 pairwise comparisons among the 4 ADE groups
contrasts_mirna <- makeContrasts(
  COV_TBI_vs_Control = COV_TBI_Exo - Cont_Exo,
  TBI_Only_vs_Control = TBI_Only_Exo - Cont_Exo,
  COV_Only_vs_Control = COV_Only_Exo - Cont_Exo,
  TBI_Only_vs_COV_TBI = TBI_Only_Exo - COV_TBI_Exo,
  COV_Only_vs_COV_TBI = COV_Only_Exo - COV_TBI_Exo,
  TBI_Only_vs_COV_Only = TBI_Only_Exo - COV_Only_Exo,
  levels = colnames(design)
)

# Perform differential expression tests
results_list <- lapply(colnames(contrasts_mirna), function(contrast_name) {
  qlf <- glmQLFTest(fit, contrast = contrasts_mirna[, contrast_name])
  topTags(qlf, n = Inf)$table %>% rownames_to_column("miRNA")
})
names(results_list) <- colnames(contrasts_mirna)

# Function to get significant miRNAs
get_significant_mirnas <- function(results, contrast_name, p_threshold = 0.05, fc_threshold = 1) {
  significant_mirnas <- results %>%
    filter(FDR < p_threshold & abs(logFC) > fc_threshold)

  print(paste("Contrast:", contrast_name))
  print(paste("Number of significant miRNAs:", nrow(significant_mirnas)))
  print(paste("Up-regulated:", sum(significant_mirnas$logFC > 0)))
  print(paste("Down-regulated:", sum(significant_mirnas$logFC < 0)))
  print("-------------------")

  return(significant_mirnas)
}

significant_mirnas_list <- lapply(names(results_list), function(contrast) {
  get_significant_mirnas(results_list[[contrast]], contrast)
})
names(significant_mirnas_list) <- names(results_list)

sapply(significant_mirnas_list, nrow)

# Save results to output folder
dir.create(OUTPUT_DIR, showWarnings = FALSE)
write.csv(as.data.frame(design), file.path(OUTPUT_DIR, "supp_figure_7d_design_matrix.csv"), row.names = FALSE)
for (nm in names(results_list)) {
  safe <- gsub("[^A-Za-z0-9_]", "_", nm)
  write.csv(results_list[[nm]], file.path(OUTPUT_DIR, sprintf("supp_figure_7d_edgeR_results_%s.csv", safe)), row.names = FALSE)
  write.csv(significant_mirnas_list[[nm]], file.path(OUTPUT_DIR, sprintf("supp_figure_7d_significant_%s.csv", safe)), row.names = FALSE)
}
sig_counts <- data.frame(contrast = names(significant_mirnas_list), n_significant = sapply(significant_mirnas_list, nrow))
write.csv(sig_counts, file.path(OUTPUT_DIR, "supp_figure_7d_DEA_significant_counts_summary.csv"), row.names = FALSE)

cat("DEA completed. All outputs saved as CSV to ", OUTPUT_DIR, ".\n")
