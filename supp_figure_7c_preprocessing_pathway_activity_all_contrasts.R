# =============================================================================
# Supplementary Figure 7c preprocessing – Pathway activity (KEGG, Hallmark,
#   Reactome, PROGENy), all 8 contrasts
# =============================================================================
#
# Note: queries the live KEGG and Reactome APIs, so results will drift from
#   the published figure as those databases update. supp_figure_7c.R
#   therefore reads from a frozen snapshot rather than this script's live
#   output; this script is kept for methods transparency.
#
# Description: Feeder script for supp_figure_7c.R. Runs Hallmark GSEA,
#   Reactome GSEA, KEGG GSEA, and PROGENy activity across all 8 contrasts
#   vs Control. Unlike figure_7's KEGG/Reactome pipeline (which only keeps
#   pathways passing p < 0.05), this uses pvalueCutoff = 1 so every
#   pathway's score is kept for every contrast, needed for a heatmap that
#   shows each selected pathway's score regardless of significance.
#   Hallmark and PROGENy already return every pathway's score by default.
#
#   leading_edge column: the gene symbols GSEA identifies as driving each
#   pathway's enrichment score (semicolon-joined), feeding supp_figure_7c.R's
#   gene x pathway grid.
#
# Prerequisites: supp_figure_7_input_data/pro_HA_results_with_entrez_{contrast}_
#   vs_Control.csv for all 8 contrasts; supp_figure_7_input_data/
#   supp_figure_7c_hallmark_gene_sets.gmt
#
# Output: supp_figure_7_output_data/supp_figure_7c_pathway_activity_all_contrasts.csv
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")

input_dir  <- "supp_figure_7/supp_figure_7_input_data"
output_dir <- "supp_figure_7/supp_figure_7_output_data"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(fgsea)
  library(ReactomePA)
  library(clusterProfiler)
  library(decoupleR)
})

set.seed(123)

contrasts <- c("COV_Only", "TBI_Only", "COV_TBI", "IL1b", "LPS", "Poly_A_U", "ODN_D_SLO3", "TNFa")

# -----------------------------------------------------------------------------
# 1. Load all 8 contrasts, build ranking metric = -log10(adj.P.Val)*sign(logFC)
#    (same convention as figure_7_preprocessing_step_4_GSEA.R), keyed by
#    both Entrez ID and symbol.
# -----------------------------------------------------------------------------
raw <- lapply(contrasts, function(c) {
  d <- read.csv(file.path(input_dir, sprintf("pro_HA_results_with_entrez_%s_vs_Control.csv", c)), stringsAsFactors = FALSE)
  d$entrezgene_id <- as.character(d$entrezgene_id)
  d <- d %>% filter(!is.na(entrezgene_id)) %>% distinct(entrezgene_id, .keep_all = TRUE)
  d %>% mutate(ranking_metric = -log10(adj.P.Val) * sign(logFC))
})
names(raw) <- contrasts

sym_map <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = unique(raw[["COV_TBI"]]$entrezgene_id), keytype = "ENTREZID", columns = "SYMBOL")
)
names(sym_map) <- c("entrezgene_id", "symbol")
sym_map <- sym_map %>% distinct(entrezgene_id, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# 2. Hallmark GSEA (fgsea) -- returns every gene set's NES regardless of
#    significance already, no pvalueCutoff filtering to worry about.
# -----------------------------------------------------------------------------
hallmark_sets <- fgsea::gmtPathways(file.path(input_dir, "supp_figure_7c_hallmark_gene_sets.gmt"))

hallmark_all <- bind_rows(lapply(contrasts, function(c) {
  d <- raw[[c]] %>% left_join(sym_map, by = "entrezgene_id") %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE)
  ranked <- sort(setNames(d$ranking_metric, d$symbol), decreasing = TRUE)
  res <- fgsea::fgsea(pathways = hallmark_sets, stats = ranked, minSize = 15, maxSize = 2000) %>% as.data.frame()
  cat(sprintf("Hallmark %-11s: %d gene sets scored\n", c, nrow(res)))
  data.frame(method = "Hallmark", contrast = c, pathway = res$pathway, score = res$NES,
             pvalue = res$pval, padj = res$padj, size = res$size,
             leading_edge = sapply(res$leadingEdge, paste, collapse = ";"), stringsAsFactors = FALSE)
}))

# core_enrichment (clusterProfiler/ReactomePA) is a "/"-delimited string of the
# ENTREZ IDs driving each pathway's enrichment score -- translate to symbols
entrez_string_to_symbols <- function(core_enrichment) {
  sapply(strsplit(core_enrichment, "/"), function(ids) {
    syms <- sym_map$symbol[match(ids, sym_map$entrezgene_id)]
    paste(syms[!is.na(syms)], collapse = ";")
  })
}

# -----------------------------------------------------------------------------
# 3. Reactome GSEA (gsePathway) -- pvalueCutoff = 1 so every pathway is kept
# -----------------------------------------------------------------------------
reactome_all <- bind_rows(lapply(contrasts, function(c) {
  d <- raw[[c]]
  ranked <- sort(setNames(d$ranking_metric, d$entrezgene_id), decreasing = TRUE)
  gsea <- tryCatch(
    ReactomePA::gsePathway(geneList = ranked, organism = "human",
                            minGSSize = 15, maxGSSize = 2000,
                            pvalueCutoff = 1, verbose = FALSE),
    error = function(e) { message("Reactome failed for ", c, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(gsea)) return(NULL)
  res <- as.data.frame(gsea)
  cat(sprintf("Reactome %-11s: %d pathways scored\n", c, nrow(res)))
  if (nrow(res) == 0) return(NULL)
  data.frame(method = "Reactome", contrast = c, pathway = res$Description, score = res$NES,
             pvalue = res$pvalue, padj = res$p.adjust, size = res$setSize,
             leading_edge = entrez_string_to_symbols(res$core_enrichment), stringsAsFactors = FALSE)
}))

# -----------------------------------------------------------------------------
# 4. KEGG GSEA (gseKEGG) -- pvalueCutoff = 1, otherwise identical to
#    figure_7_preprocessing_step_4_GSEA.R (minGSSize=15, maxGSSize=2000, seed 123)
# -----------------------------------------------------------------------------
kegg_all <- bind_rows(lapply(contrasts, function(c) {
  d <- raw[[c]]
  ranked <- sort(setNames(d$ranking_metric, d$entrezgene_id), decreasing = TRUE)
  gsea <- tryCatch(
    clusterProfiler::gseKEGG(geneList = ranked, organism = "hsa",
                              minGSSize = 15, maxGSSize = 2000,
                              pvalueCutoff = 1, verbose = FALSE),
    error = function(e) { message("KEGG failed for ", c, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(gsea)) return(NULL)
  res <- as.data.frame(gsea)
  cat(sprintf("KEGG     %-11s: %d pathways scored\n", c, nrow(res)))
  if (nrow(res) == 0) return(NULL)
  data.frame(method = "KEGG", contrast = c, pathway = res$Description, score = res$NES,
             pvalue = res$pvalue, padj = res$p.adjust, size = res$setSize,
             leading_edge = entrez_string_to_symbols(res$core_enrichment), stringsAsFactors = FALSE)
}))

# -----------------------------------------------------------------------------
# 5. PROGENy activity (decoupleR::run_mlm) -- single call across all 8
#    contrasts at once (t-statistic matrix, symbol-keyed, one column per
#    contrast).
# -----------------------------------------------------------------------------
sym_t <- lapply(contrasts, function(c) {
  d <- raw[[c]] %>% left_join(sym_map, by = "entrezgene_id") %>% filter(!is.na(symbol)) %>% distinct(symbol, .keep_all = TRUE)
  setNames(d$t, d$symbol)
})
names(sym_t) <- contrasts
common_symbols <- Reduce(intersect, lapply(sym_t, names))
progeny_mat <- sapply(sym_t, function(x) x[common_symbols])
rownames(progeny_mat) <- common_symbols
cat(sprintf("PROGENy input matrix: %d proteins x %d contrasts\n", nrow(progeny_mat), ncol(progeny_mat)))

net <- get_progeny(organism = "human", top = 500)
progeny_res <- run_mlm(mat = progeny_mat, net = net, .source = "source", .target = "target", .mor = "weight", minsize = 5)
progeny_res$padj <- ave(progeny_res$p_value, progeny_res$condition, FUN = function(p) p.adjust(p, method = "BH"))

progeny_all <- data.frame(
  method = "PROGENy", contrast = progeny_res$condition, pathway = progeny_res$source,
  score = progeny_res$score, pvalue = progeny_res$p_value, padj = progeny_res$padj,
  size = NA_real_, stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# 6. Combine and save
# -----------------------------------------------------------------------------
combined <- bind_rows(hallmark_all, reactome_all, kegg_all, progeny_all)
write.csv(combined, file.path(output_dir, "supp_figure_7c_pathway_activity_all_contrasts.csv"), row.names = FALSE)

cat(sprintf("\nSaved %d rows (method x pathway x contrast) to supp_figure_7c_pathway_activity_all_contrasts.csv\n", nrow(combined)))
print(table(combined$method, combined$contrast))
