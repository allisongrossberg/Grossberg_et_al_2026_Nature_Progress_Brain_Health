# =============================================================================
# Figure 7 – Preprocessing step 3: Pathway Enrichment Analysis
# =============================================================================
#
# Description: Over-representation analysis (GO BP, KEGG) for significant
#   proteins per contrast. Writes enrichment tables to figure_7_output_data.
#
# REPRODUCIBILITY NOTE: enrichKEGG() queries the KEGG REST API at runtime.
#   KEGG updates pathway-gene memberships between releases, so the set of
#   enriched pathways may differ slightly across runs. The original analysis
#   (Aug 2024) used KEGG Release 111.0; the current KEGG release is 117.0
#   (Jan 2026). Biological conclusions are unchanged, but exact pathway
#   counts in downstream Venn diagrams (figure_7e) may shift.
#
# Prerequisites: figure_7_output_data (run preprocessing step 2 first)
#
# Output: pro_HA_GO_BP_enrichment_*.csv, pro_HA_KEGG_enrichment_*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

OUTPUT_DIR <- "figure_7_output_data"
library(dplyr)

# Load results_with_entrez from DEA CSVs (reproducible)
entrez_files <- list.files(OUTPUT_DIR, pattern = "^pro_HA_results_with_entrez_.*\\.csv$", full.names = TRUE)
if (length(entrez_files) == 0) {
  stop("Run figure_7_preprocessing_step_2_DEA.R first. Expected pro_HA_results_with_entrez_*.csv in ", OUTPUT_DIR)
}
results_with_entrez_list <- setNames(
  lapply(entrez_files, function(f) {
    d <- read.csv(f, stringsAsFactors = FALSE)
    d$entrezgene_id <- as.character(d$entrezgene_id)
    d
  }),
  gsub("^pro_HA_results_with_entrez_|\\.csv$", "", basename(entrez_files)))

# Same get_significant_entrez logic as old DEA / COAST_Pro_HA_2_DEA.R
if (!exists("entrez_lists")) {
  get_significant_entrez <- function(results_with_entrez, contrast_name) {
    significant_genes <- results_with_entrez %>%
      filter(adj.P.Val < 0.05 & abs(logFC) > 1) %>%
      pull(entrezgene_id) %>% na.omit() %>% unique()
    if (length(significant_genes) == 0) {
      message(paste("No significant Entrez IDs found for contrast:", contrast_name))
    }
    return(significant_genes)
  }
  entrez_lists <- lapply(names(results_with_entrez_list), function(contrast) {
    get_significant_entrez(results_with_entrez_list[[contrast]], contrast)
  })
  names(entrez_lists) <- names(results_with_entrez_list)
}

######
#Load required libraries
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)
library(DOSE)
library(ggplot2)
library(factoextra)
library(STRINGdb)
library(igraph)
library(UpSetR)

####Pathway Enrichment Analysis - GO (BF) and KEGG 

# Function to perform GO enrichment analysis for a condition
perform_enrichment_KEGG_GO_BP <- function(entrez_ids, contrast) {
  if (length(entrez_ids) == 0) {
    return(NULL)
  }
  go_enrichment <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  
  kegg_enrichment <- enrichKEGG(
    gene = entrez_ids,
    organism = 'hsa',
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )
  return(list(GO = go_enrichment, KEGG = kegg_enrichment))
}

# Perform GO enrichment analysis for each condition
enrichment_by_contrast_1 <- lapply(names(entrez_lists), function(contrast) {
  perform_enrichment_KEGG_GO_BP(entrez_lists[[contrast]], contrast)
})
names(enrichment_by_contrast_1) <- names(entrez_lists)

#drop NULL values
enrichment_by_contrast_1 <- Filter(Negate(is.null), enrichment_by_contrast_1)

# Save pathway enrichment results (only columns needed downstream: figure_7e uses Description)
dir.create(OUTPUT_DIR, showWarnings = FALSE)
enrich_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
for (nm in names(enrichment_by_contrast_1)) {
  enc <- enrichment_by_contrast_1[[nm]]
  safe_name <- gsub("[^A-Za-z0-9_]", "_", nm)
  if (!is.null(enc$GO) && nrow(enc$GO@result) > 0) {
    go_res <- enc$GO@result
    keep_go <- enrich_cols[enrich_cols %in% names(go_res)]
    write.csv(go_res %>% dplyr::select(dplyr::all_of(keep_go)), file.path(OUTPUT_DIR, sprintf("pro_HA_GO_BP_enrichment_%s.csv", safe_name)), row.names = FALSE)
  }
  if (!is.null(enc$KEGG) && nrow(enc$KEGG@result) > 0) {
    kegg_res <- enc$KEGG@result
    keep_kegg <- enrich_cols[enrich_cols %in% names(kegg_res)]
    write.csv(kegg_res %>% dplyr::select(dplyr::all_of(keep_kegg)), file.path(OUTPUT_DIR, sprintf("pro_HA_KEGG_enrichment_%s.csv", safe_name)), row.names = FALSE)
  }
}

cat("Pathway enrichment completed. All results saved as CSV to ", OUTPUT_DIR, ".\n")
