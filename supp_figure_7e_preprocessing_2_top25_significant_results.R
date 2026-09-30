# =============================================================================
# Supplementary Figure 7e preprocessing 2 - Top up/down-regulated proteins (COV_TBI)
# =============================================================================
#
# Description: Feeder script for supp_figure_7e.R, which shows this gene
#   list's COV_TBI logFC as lollipop bars alongside the miRNA-target
#   network. Selects the top 13 up-regulated + top 12 down-regulated
#   proteins in COV_TBI_vs_Control (score = -log10(adj.P.Val) * logFC,
#   filtered to adj.P.Val < 0.05 & |logFC| > 2), the same convention as
#   figure_7b.R's volcano plot. Only COV_TBI has enough significant
#   proteins for a top25 list (COV_Only has none, TBI_Only has 2); their
#   values are still saved to the output CSV but not plotted.
#
# Prerequisites: supp_figure_7_input_data/pro_HA_results_with_entrez_{COV_TBI,
#   COV_Only,TBI_Only}_vs_Control.csv
#
# Output: supp_figure_7_output_data/supp_figure_7e_top25_significant_results.csv
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

read_entrez <- function(path, suffix) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  d$entrezgene_id <- as.character(d$entrezgene_id)
  d <- d %>% dplyr::select(entrezgene_id, Entry_Name, logFC, adj.P.Val) %>% distinct(entrezgene_id, .keep_all = TRUE)
  names(d)[names(d) != "entrezgene_id"] <- paste0(names(d)[names(d) != "entrezgene_id"], "_", suffix)
  d
}

cov_tbi_raw <- read.csv(file.path(input_dir, "pro_HA_results_with_entrez_COV_TBI_vs_Control.csv"), stringsAsFactors = FALSE)
cov_tbi_raw$entrezgene_id <- as.character(cov_tbi_raw$entrezgene_id)

cov_tbi  <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_COV_TBI_vs_Control.csv"),  "COV_TBI")
cov_only <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_COV_Only_vs_Control.csv"), "COV_Only")
tbi_only <- read_entrez(file.path(input_dir, "pro_HA_results_with_entrez_TBI_Only_vs_Control.csv"), "TBI_Only")

# Selection matches figure_7b.R: top N up + top N down by score, gated on
# adj.P.Val < 0.05 & |logFC| > 2 -- not a single pooled significance ranking.
top_upregulated <- cov_tbi_raw %>%
  filter(adj.P.Val < 0.05 & logFC > 2) %>%
  mutate(score = -log10(adj.P.Val) * logFC) %>%
  slice_max(score, n = 13, with_ties = FALSE)

top_downregulated <- cov_tbi_raw %>%
  filter(adj.P.Val < 0.05 & logFC < -2) %>%
  mutate(score = -log10(adj.P.Val) * abs(logFC)) %>%
  slice_max(score, n = 12, with_ties = FALSE)

top25_ids <- bind_rows(top_upregulated, top_downregulated)$entrezgene_id

top25 <- cov_tbi %>% filter(entrezgene_id %in% top25_ids)

sym <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = top25$entrezgene_id, keytype = "ENTREZID", columns = "SYMBOL")
)
names(sym) <- c("entrezgene_id", "gene_symbol")
sym <- sym %>% distinct(entrezgene_id, .keep_all = TRUE)

screened <- top25 %>%
  left_join(sym, by = "entrezgene_id") %>%
  mutate(gene_symbol = ifelse(is.na(gene_symbol), Entry_Name_COV_TBI, gene_symbol)) %>%
  left_join(cov_only, by = "entrezgene_id") %>%
  left_join(tbi_only, by = "entrezgene_id") %>%
  mutate(
    significant_COV_TBI  = adj.P.Val_COV_TBI < 0.05,
    significant_COV_Only = !is.na(adj.P.Val_COV_Only) & adj.P.Val_COV_Only < 0.05,
    significant_TBI_Only = !is.na(adj.P.Val_TBI_Only) & adj.P.Val_TBI_Only < 0.05
  ) %>%
  arrange(adj.P.Val_COV_TBI)

write.csv(screened, file.path(output_dir, "supp_figure_7e_top25_significant_results.csv"), row.names = FALSE)
message("Saved: ", file.path(output_dir, "supp_figure_7e_top25_significant_results.csv"))
