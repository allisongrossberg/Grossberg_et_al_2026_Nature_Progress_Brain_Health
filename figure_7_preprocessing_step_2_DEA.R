# =============================================================================
# Figure 7 – Preprocessing step 2: Differential Expression Analysis
# =============================================================================
#
# Description: Runs limma differential expression for all design contrasts on
#   cleaned proteomics data. Writes limma results and significant protein
#   counts to figure_7_output_data.
#
# REPRODUCIBILITY NOTE: This script queries BioMart (Ensembl) online for
#   UniProt-to-Entrez ID mapping. Ensembl updates these mappings periodically,
#   so a small number of IDs may map differently across runs. The limma DEA
#   itself is fully deterministic. Original analysis used Ensembl accessed
#   Aug 2024; minor mapping differences (~6 proteins) are expected.
#
# Prerequisites: figure_7_output_data/pro_HA_cleaned_long_data.csv (from preprocessing step 1)
#
# Output: pro_HA_limma_results_*.csv, pro_HA_significant_protein_counts.csv, etc.
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

OUTPUT_DIR <- "figure_7_output_data"
cleaned_file <- file.path(OUTPUT_DIR, "pro_HA_cleaned_long_data.csv")
if (!file.exists(cleaned_file)) {
  stop("Run figure_7_preprocessing_1_data_cleaning.R first. Expected: ", cleaned_file)
}
protein_sample_HA_long_filtered_sc_log <- read.csv(cleaned_file, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Load packages and build design matrix for limma
# -----------------------------------------------------------------------------
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(stringr)
library(igraph)
library(biomaRt)
library(tidyverse)

#Convert to wide format
protein_HA_wide_filtered_sc_log <- protein_sample_HA_long_filtered_sc_log %>%
  dplyr::select(Entry_Name, Condition, log2_Intensity) %>%
  pivot_wider(names_from = Condition, values_from = log2_Intensity)

# Convert tibble to data frame
protein_HA_wide_filtered_sc_log_df <- as.data.frame(protein_HA_wide_filtered_sc_log)

# Set row names to be the Entry_Name
rownames(protein_HA_wide_filtered_sc_log_df) <- protein_HA_wide_filtered_sc_log_df$Entry_Name

# Remove the Entry_Name column
protein_HA_wide_filtered_sc_log_df <- protein_HA_wide_filtered_sc_log_df[, -1]

print(colnames(protein_HA_wide_filtered_sc_log_df))

# Extract condition names without replicate numbers
clean_conditions <- unique(gsub("_\\d+$", "", trimws(colnames(protein_HA_wide_filtered_sc_log_df))))
print(clean_conditions)

# Create a factor for the conditions
condition_factor <- factor(gsub("_\\d+$", "", trimws(colnames(protein_HA_wide_filtered_sc_log_df))))

# Create the design matrix
design <- model.matrix(~0 + condition_factor)
colnames(design) <- clean_conditions

print(colnames(design))

dim(design)

# Define contrasts
contrasts <- makeContrasts(
  COV_Only_vs_Control = COV_Only_Exo - Cont_Exo,
  COV_TBI_vs_Control = COV_TBI_Exo - Cont_Exo,
  TBI_Only_vs_Control = TBI_Only_Exo - Cont_Exo,
  IL1b_vs_Control = IL1b - Cont_Exo,
  TNFa_vs_Control = TNFa - Cont_Exo,
  LPS_vs_Control = LPS - Cont_Exo,
  ODN_D_SLO3_vs_Control = ODN_D_SLO3 - Cont_Exo,
  Poly_A_U_vs_Control = Poly_A_U - Cont_Exo,
  COV_Only_vs_COV_TBI = COV_Only_Exo - COV_TBI_Exo,
  TBI_Only_vs_COV_TBI = TBI_Only_Exo - COV_TBI_Exo,
  IL1b_vs_COV_TBI = IL1b - COV_TBI_Exo,
  TNFa_vs_COV_TBI = TNFa - COV_TBI_Exo,
  LPS_vs_COV_TBI = LPS - COV_TBI_Exo,
  ODN_D_SLO3_vs_COV_TBI = ODN_D_SLO3 - COV_TBI_Exo,
  Poly_A_U_vs_COV_TBI = Poly_A_U - COV_TBI_Exo,
  COV_Only_vs_NegCTRL = COV_Only_Exo - Neg_CTRL,
  COV_TBI_vs_NegCTRL = COV_TBI_Exo - Neg_CTRL,
  TBI_Only_vs_NegCTRL = TBI_Only_Exo - Neg_CTRL,
  Cont_Exo_vs_NegCTRL = Cont_Exo - Neg_CTRL,
  IL1b_vs_NegCTRL = IL1b - Neg_CTRL,
  TNFa_vs_NegCTRL = TNFa - Neg_CTRL,
  LPS_vs_NegCTRL = LPS - Neg_CTRL,
  ODN_D_SLO3_vs_NegCTRL = ODN_D_SLO3 - Neg_CTRL,
  Poly_A_U_vs_NegCTRL = Poly_A_U - Neg_CTRL,
  TBI_Only_vs_COV_Only = TBI_Only_Exo - COV_Only_Exo,
  levels = design
)

# Fit the linear model
fit <- lmFit(protein_HA_wide_filtered_sc_log_df, design)

# Apply contrasts
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

# Get results
results_list <- lapply(colnames(contrasts), function(contrast) {
  topTable(fit2, coef = contrast, n = Inf) %>%
    rownames_to_column("Entry_Name") %>%
    arrange(adj.P.Val)
})
names(results_list) <- colnames(contrasts)

head(results_list[[1]])

# Function to get significant proteins and print the result
get_significant_proteins <- function(results, contrast_name) {
  significant_proteins <- results %>% 
    filter(adj.P.Val < 0.05 & abs(logFC) > 1)
  
  print(paste("Contrast:", contrast_name))
  print(paste("Number of significant proteins:", nrow(significant_proteins)))
  print(paste("Up-regulated:", sum(significant_proteins$logFC > 0)))
  print(paste("Down-regulated:", sum(significant_proteins$logFC < 0)))
  print("-------------------")
  
  return(significant_proteins)
}

# Apply the function to all contrasts
significant_proteins_list <- lapply(names(results_list), function(contrast) {
  get_significant_proteins(results_list[[contrast]], contrast)
})

names(significant_proteins_list) <- names(results_list)

sapply(significant_proteins_list, nrow)

######
# Entrez ID Mapping: BioMart (matching original old_code/COAST_Pro_HA_2_DEA.R)

uniprot_ids <- protein_sample_HA_long_filtered_sc_log$Protein_ID
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
id_conversion <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id"),
                      filters = "uniprot_gn_id",
                      values = uniprot_ids,
                      mart = ensembl)
# No collapse: match old_code exactly (one row per UniProt–Entrez pair from BioMart).
protein_sample_HA_long_filtered_sc_log_ep <- protein_sample_HA_long_filtered_sc_log %>%
  left_join(id_conversion, by = c("Protein_ID" = "uniprot_gn_id")) %>%
  mutate(has_entrez_id = !is.na(entrezgene_id))

entry_name_entrez_ids <- protein_sample_HA_long_filtered_sc_log_ep %>%
  mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  dplyr::select(Entry_Name, entrezgene_id) %>%
  filter(!is.na(entrezgene_id)) %>%
  distinct()

######
# Map Entrez IDs to Limma results

# Function to map Entrez IDs to Limma results
map_entrez_to_limma_results <- function(limma_results, mapping_df) {
  # Join the limma results with the mapping dataframe based on Entry_Name
  results_with_entrez <- limma_results %>%
    left_join(mapping_df, by = "Entry_Name") %>%
    # Select relevant columns
    dplyr::select(Entry_Name, logFC, AveExpr, t, P.Value, adj.P.Val, B, entrezgene_id) %>%
    # drop duplicates
    distinct(Entry_Name, .keep_all = TRUE)
  
  return(results_with_entrez)
}

# Apply the mapping to all contrasts
results_with_entrez_list <- lapply(names(results_list), function(contrast) {
  map_entrez_to_limma_results(results_list[[contrast]], protein_sample_HA_long_filtered_sc_log_ep)
})

names(results_with_entrez_list) <- names(results_list)

# Function to get significant Entrez IDs
get_significant_entrez <- function(results_with_entrez, contrast_name) {
  significant_genes <- results_with_entrez %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
    pull(entrezgene_id) %>%
    na.omit() %>%
    unique()
  
  if (length(significant_genes) == 0) {
    message(paste("No significant Entrez IDs found for contrast:", contrast_name))
  }
  
  return(significant_genes)
}

# Get significant Entrez IDs for all contrasts
entrez_lists <- lapply(names(results_with_entrez_list), function(contrast) {
  get_significant_entrez(results_with_entrez_list[[contrast]], contrast)
})
names(entrez_lists) <- names(results_with_entrez_list)

# Function to get significant Entrez IDs dataframes
get_significant_entrez_df <- function(results_with_entrez, contrast_name) {
  significant_genes_df <- results_with_entrez %>%
    filter(adj.P.Val < 0.05 & abs(logFC) > 1)
  return(significant_genes_df)
}

# Get significant Entrez IDs dataframes for all contrasts
significant_entrez_dfs <- lapply(names(results_with_entrez_list), function(contrast) {
  get_significant_entrez_df(results_with_entrez_list[[contrast]], contrast)
})
names(significant_entrez_dfs) <- names(results_with_entrez_list)

# Save intermediate DEA results and all dataframes to output folder
dir.create(OUTPUT_DIR, showWarnings = FALSE)
# Wide format matrix used for limma
protein_HA_wide_for_limma <- protein_HA_wide_filtered_sc_log_df
protein_HA_wide_for_limma$Entry_Name <- rownames(protein_HA_wide_for_limma)
protein_HA_wide_for_limma <- protein_HA_wide_for_limma[, c("Entry_Name", setdiff(colnames(protein_HA_wide_for_limma), "Entry_Name"))]
write.csv(protein_HA_wide_for_limma, file.path(OUTPUT_DIR, "pro_HA_wide_for_limma.csv"), row.names = FALSE)
# Design matrix
write.csv(as.data.frame(design), file.path(OUTPUT_DIR, "pro_HA_design_matrix.csv"), row.names = FALSE)
# BioMart Entrez mapping
write.csv(id_conversion, file.path(OUTPUT_DIR, "pro_HA_biomart_entrez_mapping.csv"), row.names = FALSE)
# entry_name_entrez_ids: Entry_Name, entrezgene_id (already minimal)
write.csv(entry_name_entrez_ids, file.path(OUTPUT_DIR, "pro_HA_entry_name_entrez_ids.csv"), row.names = FALSE)
# limma results: figure_7b needs Entry_Name, logFC, adj.P.Val; keep standard topTable columns
limma_cols <- c("Entry_Name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
for (nm in names(results_list)) {
  safe <- gsub("[^A-Za-z0-9_]", "_", nm)
  r <- results_list[[nm]]
  keep_limma <- limma_cols[limma_cols %in% names(r)]
  write.csv(r %>% dplyr::select(dplyr::all_of(keep_limma)), file.path(OUTPUT_DIR, sprintf("pro_HA_limma_results_%s.csv", safe)), row.names = FALSE)
  # results_with_entrez: step 3/4/5 need entrezgene_id, logFC; figure scripts need Entry_Name, logFC, adj.P.Val
  rwe <- results_with_entrez_list[[nm]]
  keep_rwe <- c("Entry_Name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "entrezgene_id")
  keep_rwe <- keep_rwe[keep_rwe %in% names(rwe)]
  write.csv(rwe %>% dplyr::select(dplyr::all_of(keep_rwe)), file.path(OUTPUT_DIR, sprintf("pro_HA_results_with_entrez_%s.csv", safe)), row.names = FALSE)
}
sig_counts <- data.frame(contrast = names(significant_proteins_list), n_significant = sapply(significant_proteins_list, nrow))
write.csv(sig_counts, file.path(OUTPUT_DIR, "pro_HA_significant_protein_counts.csv"), row.names = FALSE)

# significant_entrez: figure_7d needs entrezgene_id; keep Entry_Name, logFC, adj.P.Val for reproducibility
sig_entrez_cols <- c("Entry_Name", "logFC", "adj.P.Val", "entrezgene_id")
for (nm in names(significant_entrez_dfs)) {
  safe <- gsub("[^A-Za-z0-9_]", "_", nm)
  s <- significant_entrez_dfs[[nm]]
  keep_sig <- sig_entrez_cols[sig_entrez_cols %in% names(s)]
  write.csv(s %>% dplyr::select(dplyr::all_of(keep_sig)), file.path(OUTPUT_DIR, sprintf("pro_HA_significant_entrez_%s.csv", safe)), row.names = FALSE)
}

cat("DEA completed. All outputs saved as CSV to ", OUTPUT_DIR, ".\n")
