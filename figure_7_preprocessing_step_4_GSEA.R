# =============================================================================
# Figure 7 – Preprocessing step 4: Gene Set Enrichment Analysis (GSEA)
# =============================================================================
#
# Description: GSEA on ranked protein list per contrast. Writes GSEA KEGG
#   (and optionally GO) results to figure_7_output_data.
#
# NOTE: gseKEGG() and enrichKEGG() query the KEGG REST API (rest.kegg.jp)
#   at runtime, so results depend on the KEGG release available at the time
#   of execution. Pathway-gene memberships change between releases, which can
#   shift borderline p-values and alter which pathways are significant.
#   Original analysis (Aug 2024): KEGG Release 111.0 (2024-07-01)
#   Last modified (Sep 2025):     KEGG Release 115.1 (2025-08-01)
#   Most recent run (Mar 2026):   KEGG Release 117.0 (2026-01-01)
#
# Prerequisites: figure_7_output_data (run preprocessing step 2 first)
#
# Output: pro_HA_GSEA_KEGG_*.csv (and optional GO)
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

OUTPUT_DIR <- "figure_7_output_data"

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

######
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(ggplot2)
library(dplyr)
library(ReactomePA)

set.seed(123)
perform_GSEA_Proteomics <- function(entrez_ids, contrast) {
  entrez_ids <- sort(entrez_ids, decreasing = TRUE)
  
  go_gsea <- gseGO(geneList = entrez_ids,
                   ont = "MF",
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   minGSSize = 15,
                   maxGSSize = 2000,
                   pvalueCutoff = 0.05,
                   verbose = FALSE)
  
  kegg_gsea <- gseKEGG(geneList = entrez_ids,
                       organism = 'hsa',
                       minGSSize = 15,
                       maxGSSize = 2000,
                       pvalueCutoff = 0.05,
                       verbose = FALSE)
  
  return(list(GO = go_gsea, KEGG = kegg_gsea))
}

contrasts <- lapply(names(results_with_entrez_list), function(contrast_name) {
  df <- results_with_entrez_list[[contrast_name]]
  df <- df[!is.na(df$entrezgene_id), ]
  df <- df[!duplicated(df$entrezgene_id), ]
  
  entrez_col <- "entrezgene_id"
  logfc_col <- "logFC"
  pvalue_col <- "adj.P.Val"
  
  ranking_metric <- with(df, -log10(get(pvalue_col)) * sign(get(logfc_col)))
  
  entrez_ids <- df[[entrez_col]]
  
  valid_indices <- !is.na(entrez_ids) & !is.na(ranking_metric)
  entrez_ids <- entrez_ids[valid_indices]
  ranking_metric <- ranking_metric[valid_indices]
  
  unique_indices <- !duplicated(entrez_ids)
  entrez_ids <- entrez_ids[unique_indices]
  ranking_metric <- ranking_metric[unique_indices]
  
  list(entrez_ids = entrez_ids, metric = ranking_metric)
})
names(contrasts) <- names(results_with_entrez_list)

# Perform GSEA for each contrast
gsea_results_proteomics <- list()
for (contrast_name in names(contrasts)) {
  contrast_data <- contrasts[[contrast_name]]
  gene_list <- setNames(contrast_data$metric, contrast_data$entrez_ids)
  gsea_results_proteomics[[contrast_name]] <- perform_GSEA_Proteomics(gene_list, contrast_name)
}

for (contrast_name in names(gsea_results_proteomics)) {
  cat("\nContrast:", contrast_name, "\n")
  cat("GO results:", nrow(gsea_results_proteomics[[contrast_name]]$GO@result), "\n")
  cat("KEGG results:", nrow(gsea_results_proteomics[[contrast_name]]$KEGG@result), "\n")
}

# Save GSEA results (only columns needed downstream: figure_7c uses Description; step 5 uses KEGG result + genesets)
dir.create(OUTPUT_DIR, showWarnings = FALSE)
gsea_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge")
for (nm in names(gsea_results_proteomics)) {
  gsea <- gsea_results_proteomics[[nm]]
  safe_name <- gsub("[^A-Za-z0-9_]", "_", nm)
  if (!is.null(gsea$GO) && nrow(gsea$GO@result) > 0) {
    go_res <- gsea$GO@result
    keep_go <- gsea_cols[gsea_cols %in% names(go_res)]
    write.csv(go_res %>% dplyr::select(dplyr::all_of(keep_go)), file.path(OUTPUT_DIR, sprintf("pro_HA_GSEA_GO_%s.csv", safe_name)), row.names = FALSE)
  }
  if (!is.null(gsea$KEGG) && nrow(gsea$KEGG@result) > 0) {
    kegg_res <- gsea$KEGG@result
    keep_kegg <- gsea_cols[gsea_cols %in% names(kegg_res)]
    write.csv(kegg_res %>% dplyr::select(dplyr::all_of(keep_kegg)), file.path(OUTPUT_DIR, sprintf("pro_HA_GSEA_KEGG_%s.csv", safe_name)), row.names = FALSE)
    # KEGG gene sets: step 5 needs kegg_id, entrez_id only
    gs <- gsea$KEGG@geneSets
    if (length(gs) > 0) {
      geneset_df <- do.call(rbind, lapply(names(gs), function(kid) data.frame(kegg_id = kid, entrez_id = gs[[kid]], stringsAsFactors = FALSE)))
      write.csv(geneset_df, file.path(OUTPUT_DIR, sprintf("pro_HA_GSEA_KEGG_genesets_%s.csv", safe_name)), row.names = FALSE)
    }
  }
}

cat("GSEA completed. All outputs saved as CSV to ", OUTPUT_DIR, ".\n")
