# =============================================================================
# Supplementary Figure 7e preprocessing 1 - miRNA-target overlap with COV_TBI DE proteins
# =============================================================================
#
# Description: Pulls each of the 7 significant COV_TBI miRNAs'
#   experimentally validated target genes (multiMiR, aggregating
#   miRTarBase/TarBase/miRecords; predicted-only targets are excluded as
#   high false-positive), then intersects those targets against (1) the
#   full COV_TBI_vs_Control significant DE protein list (501 proteins) and
#   (2) the Escartin et al. 2021 marker panel genes significant in COV_TBI.
#   For every overlap, checks whether the miRNA/target directions are
#   inverse (consistent with direct repression) or same-direction (kept
#   and flagged, not dropped, since it could reflect indirect regulation).
#
#   Caveat: this is validated-target overlap with a DE gene list, not
#   proof these miRNAs regulate these proteins in this dataset.
#
# Prerequisites: supp_figure_7_input_data/mirna_significant_COV_TBI_vs_Control.csv;
#   supp_figure_7_input_data/pro_HA_results_with_entrez_COV_TBI_vs_Control.csv;
#   supp_figure_7_output_data/supp_figure_7b_marker_panel_screen_results.csv
#
# Output: supp_figure_7_output_data/supp_figure_7e_mirna_target_overlap_all.csv
#         supp_figure_7_output_data/supp_figure_7e_mirna_target_overlap_summary.csv
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")

input_dir  <- "supp_figure_7/supp_figure_7_input_data"
output_dir <- "supp_figure_7/supp_figure_7_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(dplyr)
  library(multiMiR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# -----------------------------------------------------------------------------
# 1. The 7 significant miRNAs (COV_TBI_vs_Control)
# -----------------------------------------------------------------------------
mirnas <- read.csv(file.path(input_dir, "mirna_significant_COV_TBI_vs_Control.csv"), stringsAsFactors = FALSE)
cat(sprintf("%d significant miRNAs (COV_TBI_vs_Control)\n", nrow(mirnas)))

# -----------------------------------------------------------------------------
# 2. Target gene set A: full COV_TBI_vs_Control significant DE protein list
#    (adj.P.Val < 0.05 & |logFC| > 2 -- same threshold as figure_7b / the
#    top25 up/down heatmap)
# -----------------------------------------------------------------------------
dea <- read.csv(file.path(input_dir, "pro_HA_results_with_entrez_COV_TBI_vs_Control.csv"), stringsAsFactors = FALSE)
dea$entrezgene_id <- as.character(dea$entrezgene_id)
dea_sig <- dea %>% filter(!is.na(adj.P.Val), adj.P.Val < 0.05, abs(logFC) > 2, !is.na(entrezgene_id))

sym_map <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = unique(dea_sig$entrezgene_id), keytype = "ENTREZID", columns = "SYMBOL")
)
names(sym_map) <- c("entrezgene_id", "symbol")
sym_map <- sym_map %>% distinct(entrezgene_id, .keep_all = TRUE)

dea_sig <- dea_sig %>% left_join(sym_map, by = "entrezgene_id") %>%
  mutate(symbol = ifelse(is.na(symbol), Entry_Name, symbol)) %>%
  distinct(symbol, .keep_all = TRUE)
cat(sprintf("Target set A: %d significant COV_TBI DE proteins (%d up, %d down)\n",
            nrow(dea_sig), sum(dea_sig$logFC > 0), sum(dea_sig$logFC < 0)))

# -----------------------------------------------------------------------------
# 3. Target gene set B: Escartin panel genes significant in COV_TBI
# -----------------------------------------------------------------------------
panel_path <- file.path(output_dir, "supp_figure_7b_marker_panel_screen_results.csv")
if (!file.exists(panel_path)) stop("Run supp_figure_7b_preprocessing_marker_panel_screen.R first. Expected: ", panel_path)
panel <- read.csv(panel_path, stringsAsFactors = FALSE)
panel_sig <- panel %>% filter(significant_COV_TBI == TRUE)
cat(sprintf("Target set B: %d Escartin panel genes significant in COV_TBI\n", nrow(panel_sig)))

# -----------------------------------------------------------------------------
# 4. Pull validated targets per miRNA (multiMiR), intersect with both sets
# -----------------------------------------------------------------------------
all_hits <- list()
for (i in seq_len(nrow(mirnas))) {
  m <- mirnas$miRNA[i]
  mirna_logFC <- mirnas$logFC[i]
  res <- tryCatch(
    multiMiR::get_multimir(mirna = m, table = "validated", summary = TRUE),
    error = function(e) { message("multiMiR failed for ", m, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(res) || nrow(res@data) == 0) {
    cat(sprintf("%-18s: no validated targets returned\n", m)); next
  }
  targets <- res@data %>% distinct(target_symbol, database, .keep_all = FALSE) %>%
    group_by(target_symbol) %>% summarise(databases = paste(sort(unique(database)), collapse = ";"), .groups = "drop")

  hit_A <- targets %>% inner_join(dea_sig %>% dplyr::select(symbol, logFC_target = logFC, adj.P.Val_target = adj.P.Val),
                                   by = c("target_symbol" = "symbol")) %>% mutate(target_source = "COV_TBI DE protein")
  hit_B <- targets %>% inner_join(panel_sig %>% dplyr::select(gene_symbol, logFC_target = logFC_COV_TBI, adj.P.Val_target = adj.P.Val_COV_TBI, category),
                                   by = c("target_symbol" = "gene_symbol")) %>% mutate(target_source = "Escartin panel")

  hit <- bind_rows(hit_A, hit_B) %>%
    mutate(miRNA = m, miRNA_logFC = mirna_logFC,
           directionality_consistent = sign(miRNA_logFC) != sign(logFC_target))

  n_A <- nrow(hit_A); n_B <- nrow(hit_B)
  cat(sprintf("%-18s: %d validated targets checked -> %d hit DE-protein set, %d hit Escartin panel\n",
              m, length(unique(targets$target_symbol)), n_A, n_B))
  all_hits[[m]] <- hit
}

combined <- bind_rows(all_hits)

if (nrow(combined) == 0) {
  cat("\nNo overlaps found between any miRNA's validated targets and either gene set.\n")
} else {
  combined <- combined %>%
    dplyr::select(miRNA, miRNA_logFC, target_symbol, target_source, category,
                   logFC_target, adj.P.Val_target, directionality_consistent, databases) %>%
    arrange(miRNA, desc(directionality_consistent), target_source)
  write.csv(combined, file.path(output_dir, "supp_figure_7e_mirna_target_overlap_all.csv"), row.names = FALSE)

  summary_tab <- combined %>% group_by(miRNA, target_source) %>%
    summarise(n_targets = n(), n_direction_consistent = sum(directionality_consistent), .groups = "drop")
  write.csv(summary_tab, file.path(output_dir, "supp_figure_7e_mirna_target_overlap_summary.csv"), row.names = FALSE)

  cat(sprintf("\n%d total overlap rows (%d with directionality consistent with repression)\n",
              nrow(combined), sum(combined$directionality_consistent)))
  cat("\nFull hit table:\n")
  print(combined, n = 100)
}
