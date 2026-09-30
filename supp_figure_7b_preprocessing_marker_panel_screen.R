# =============================================================================
# Supplementary Figure 7b preprocessing - Systematic literature marker panel screen
# =============================================================================
#
# Description: Screens the full measured proteome against every marker
#   listed in Table 1 ("Potential markers of reactive astrocytes") of
#   Escartin et al. 2021, Nat Neurosci (PMID 33589835).
#   supp_figure_7b_marker_panel_reference.csv is a direct transcription of
#   that table (25 markers). For every marker quantified in this
#   proteomics dataset, this reports its COV_TBI_vs_Control (and
#   COV_Only/TBI_Only) logFC and significance.
#
# Prerequisites: supp_figure_7b_input_data/supp_figure_7b_marker_panel_reference.csv;
#   supp_figure_7b_input_data/pro_HA_results_with_entrez_{COV_TBI,COV_Only,TBI_Only}_vs_Control.csv
#
# Output: supp_figure_7b_output_data/supp_figure_7b_marker_panel_screen_results.csv
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")

input_dir  <- "supp_figure_7/supp_figure_7_input_data"
output_dir <- "supp_figure_7/supp_figure_7_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

panel <- read.csv(file.path(input_dir, "supp_figure_7b_marker_panel_reference.csv"), stringsAsFactors = FALSE)

# Map gene symbols to Entrez IDs (same ID space used by the figure_7 DEA
# pipeline). Joining on Entrez ID rather than the raw proteomics file's
# Entry_Name avoids a symbol collision: Clusterin (Entrez 1191) is listed
# under Entry_Name "CLUS" there, while Entry_Name "CLU" is actually CLUH
# (Entrez 23277), a different gene.
sym_map <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = panel$gene_symbol, keytype = "SYMBOL", columns = "ENTREZID")
) %>%
  dplyr::rename(gene_symbol = SYMBOL, entrezgene_id = ENTREZID) %>%
  distinct(gene_symbol, .keep_all = TRUE)

panel_mapped <- panel %>% left_join(sym_map, by = "gene_symbol")

read_entrez <- function(path, suffix) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$entrezgene_id <- as.character(d$entrezgene_id)
  d <- d %>% dplyr::select(entrezgene_id, Entry_Name, logFC, adj.P.Val) %>% distinct(entrezgene_id, .keep_all = TRUE)
  names(d)[names(d) != "entrezgene_id"] <- paste0(names(d)[names(d) != "entrezgene_id"], "_", suffix)
  d
}

cov_tbi  <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_COV_TBI_vs_Control.csv"),  "COV_TBI")
cov_only <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_COV_Only_vs_Control.csv"), "COV_Only")
tbi_only <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_TBI_Only_vs_Control.csv"), "TBI_Only")

panel_mapped$entrezgene_id <- as.character(panel_mapped$entrezgene_id)

screened <- panel_mapped %>%
  left_join(cov_tbi,  by = "entrezgene_id") %>%
  left_join(cov_only, by = "entrezgene_id") %>%
  left_join(tbi_only, by = "entrezgene_id") %>%
  mutate(
    detected_in_proteomics = !is.na(Entry_Name_COV_TBI),
    significant_COV_TBI = detected_in_proteomics & adj.P.Val_COV_TBI < 0.05,
    significant_COV_Only = detected_in_proteomics & !is.na(adj.P.Val_COV_Only) & adj.P.Val_COV_Only < 0.05,
    significant_TBI_Only = detected_in_proteomics & !is.na(adj.P.Val_TBI_Only) & adj.P.Val_TBI_Only < 0.05
  ) %>%
  arrange(desc(detected_in_proteomics), table1_section, gene_symbol)

write.csv(screened, file.path(output_dir, "supp_figure_7b_marker_panel_screen_results.csv"), row.names = FALSE)

n_detected <- sum(screened$detected_in_proteomics)
n_total <- nrow(screened)
cat(sprintf("Marker panel coverage: %d / %d Table 1 markers (%.0f%%) quantified in this proteomics dataset.\n",
            n_detected, n_total, 100 * n_detected / n_total))

hits <- screened %>% filter(significant_COV_TBI)
cat(sprintf("\n%d panel markers are significant in COV_TBI_vs_Control:\n", nrow(hits)))
if (nrow(hits) > 0) {
  print(hits %>% dplyr::select(gene_symbol, table1_section, logFC_COV_TBI, adj.P.Val_COV_TBI))
}
