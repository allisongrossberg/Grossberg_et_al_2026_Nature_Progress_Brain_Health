# =============================================================================
# Supplementary Figure 7c – Pathway enrichment heatmap + shared leading-edge
#   driver genes (Hallmark/KEGG/Reactome, COVID-19+mTBI+ ADEs vs. control)
# =============================================================================
#
# Note: supp_figure_7c_preprocessing_pathway_activity_all_contrasts.R queries
#   the live KEGG and Reactome APIs, so re-running it will not exactly
#   reproduce this panel (pathway-gene memberships change between database
#   releases). This script instead loads the frozen pathway-activity table
#   from the original 2024 analysis (supp_figure_7_output_data/
#   original_2024_objects/).
#
# Description: Selects pathways reaching adjusted p < 0.05 across Hallmark,
#   KEGG, and Reactome GSEA (Reactome capped at the top 15 by adjusted
#   p-value), then filters to 19 pathways relevant to infection/injury/
#   inflammation/CNS/neurodegeneration. A second panel to the right shows
#   the leading-edge genes driving each pathway's enrichment, colored by
#   each gene's own Log2FC (COV_TBI vs Control); genes sharing an identical
#   set of pathways are collapsed into one column ("<gene> +N"). Both axes
#   are hierarchically clustered by shared pathway/gene membership so
#   related pathways and genes sit next to each other.
#
# Prerequisites: supp_figure_7c_output_data/supp_figure_7c_pathway_activity_all_contrasts.csv
#   (from supp_figure_7c_pathway_activity_all_contrasts.R)
#
# Output: supp_figure_7c_pathway_driver_genes_heatmap.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
if (basename(getwd()) == "supp_figure_7") setwd("..")

input_dir  <- "supp_figure_7/supp_figure_7_input_data"
output_dir <- "supp_figure_7/supp_figure_7_output_data"
tif_dir    <- "supp_figure_7"

suppressMessages({
  library(dplyr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
  library(stringr)
})

d <- read.csv(file.path(output_dir, "original_2024_objects/supp_figure_7c_pathway_activity_all_contrasts_original.csv"), stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Same pathway selection as supp_figure_7c_combined_pathways_heatmap.R
# -----------------------------------------------------------------------------
cns_keywords <- "neuro|synap|axon|glia|astrocyt|myelin|interferon|complement|interleukin|chemokine|toll.like|nf.kappab|nf-kb|nrf2|oxidative|antioxidant"

cov_tbi_sig <- d %>% filter(contrast == "COV_TBI", !is.na(padj), padj < 0.05)

hallmark_rows <- cov_tbi_sig %>% filter(method == "Hallmark") %>% pull(pathway)
kegg_rows     <- cov_tbi_sig %>% filter(method == "KEGG") %>% pull(pathway)

reactome_sig <- cov_tbi_sig %>% filter(method == "Reactome") %>%
  mutate(cns_relevant = grepl(cns_keywords, pathway, ignore.case = TRUE)) %>%
  arrange(desc(cns_relevant), padj)
reactome_rows <- reactome_sig %>% head(15) %>% pull(pathway)

selected <- bind_rows(
  cov_tbi_sig %>% filter(method == "Hallmark", pathway %in% hallmark_rows),
  cov_tbi_sig %>% filter(method == "KEGG",     pathway %in% kegg_rows),
  cov_tbi_sig %>% filter(method == "Reactome", pathway %in% reactome_rows)
) %>% mutate(method = factor(method, levels = c("Hallmark", "KEGG", "Reactome")))

# Manually curated scope filter (same as supp_figure_7c_combined_pathways_
# heatmap.R -- see that script's comment for the full rationale and
# borderline calls): drop pathways not related to infection, injury,
# inflammation (especially CNS), or neurodegenerative disease. Oxidative
# Phosphorylation, Chemical carcinogenesis - reactive oxygen species, and
# Apoptosis are kept despite not being disease names themselves -- they're
# mechanistically central to the neurodegenerative cluster (same driver
# genes, see the ATP5F1A/VDAC/SDH block below); Diabetic cardiomyopathy is
# dropped since it shares that gene set by GSEA construction but isn't
# itself CNS/infection/injury/neurodegeneration-relevant.
excluded_pathways <- c(
  "Myc Targets V1", "Protein Secretion", "Xenobiotic Metabolism",
  "Ribosome", "Glutathione metabolism", "Carbon metabolism", "Spliceosome",
  "Aminoacyl-tRNA biosynthesis", "Biosynthesis of amino acids", "Diabetic cardiomyopathy",
  "Eukaryotic Translation Termination", "Response of EIF2AK4 (GCN2) to amino acid deficiency",
  "Eukaryotic Translation Elongation", "L13a-mediated translational silencing of Ceruloplasmin expression",
  "GTP hydrolysis and joining of the 60S ribosomal subunit",
  "Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)",
  "Salmonella infection", "Shigellosis", "Allograft Rejection"
)
selected <- selected %>% filter(!pathway %in% excluded_pathways)

all_pathways <- selected$pathway
n_pathways <- length(all_pathways)

cap <- 2

fmt_val <- function(x) {
  r <- round(x, 1)
  ifelse(r == 0, "0.0", sprintf("%.1f", r))
}
sig_stars <- function(p) {
  dplyr::case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}
selected <- selected %>% mutate(value_label = fmt_val(score), star_label = sig_stars(padj))

# -----------------------------------------------------------------------------
# Leading-edge driver genes, collapsed by identical pathway signature (see
# header). No minimum-pathway-count threshold: an earlier version required
# a gene to be in >=8 of the selected pathways' leading edges, tuned for the
# original 38-pathway set -- but that number doesn't generalize (it's a
# fraction-of-total heuristic in disguise), and after the infection/injury/
# inflammation/neurodegeneration scope filter cut the pathway count to 22,
# the same fixed >=8 nearly emptied the panel. Instead, every unique
# signature is ranked by how many pathways it spans and the top N shown --
# a rank-based cutoff that doesn't need re-tuning if the pathway count
# changes again.
# -----------------------------------------------------------------------------
gene_pathway_pairs <- selected %>%
  transmute(pathway = as.character(pathway), method, leading_edge) %>%
  tidyr::separate_rows(leading_edge, sep = ";") %>%
  mutate(gene = trimws(leading_edge)) %>%
  filter(gene != "") %>%
  distinct(pathway, method, gene)

COLLAPSE_THRESHOLD <- 4
gene_signature <- gene_pathway_pairs %>%
  group_by(gene) %>%
  summarise(sig = paste(sort(pathway), collapse = "|"), n_pathways = n(), .groups = "drop")

# Every unique signature shown, no top-N cap -- matches the standard GSEA
# Leading Edge Analysis / clusterProfiler cnetplot convention of displaying
# the full leading-edge set for the selected pathways rather than
# restricting to genes that are widely shared; only exact duplicates (genes
# with an identical pathway signature) are collapsed, which is a display
# simplification, not a relevance filter -- every gene from every selected
# pathway's leading edge is represented by some column.
sig_groups <- gene_signature %>%
  group_by(sig) %>%
  summarise(genes = list(sort(gene)), n_pathways = dplyr::first(n_pathways), n_genes = dplyr::n(), .groups = "drop") %>%
  arrange(desc(n_pathways)) %>%
  rowwise() %>%
  mutate(representative = genes[[1]]) %>%
  ungroup()

driver_columns <- bind_rows(lapply(seq_len(nrow(sig_groups)), function(i) {
  g <- sig_groups$genes[[i]]
  if (sig_groups$n_genes[i] > COLLAPSE_THRESHOLD) {
    data.frame(gene = sig_groups$representative[i],
               display_label = paste0(sig_groups$representative[i], " +", sig_groups$n_genes[i] - 1),
               stringsAsFactors = FALSE)
  } else {
    data.frame(gene = g, display_label = g, stringsAsFactors = FALSE)
  }
}))
driver_genes <- driver_columns$gene
all_columns <- driver_columns

collapsed <- sig_groups %>% filter(n_genes > COLLAPSE_THRESHOLD)
cat(sprintf("%d driver columns shown (every unique pathway-signature, no top-N cap; %d collapsed groups of >%d genes sharing an identical signature; %d selected pathways)\n",
            nrow(driver_columns), nrow(collapsed), COLLAPSE_THRESHOLD, n_pathways))
collapsed_caption_lines <- character(0)
if (nrow(collapsed) > 0) {
  for (i in seq_len(nrow(collapsed))) {
    label <- paste0(collapsed$representative[i], " +", collapsed$n_genes[i] - 1)
    members <- paste(collapsed$genes[[i]], collapse = ", ")
    cat(sprintf("  %s: %s\n", label, members))
    # Wrap long membership lists to a fixed character width -- at 25+ gene
    # symbols, one unwrapped line runs off the canvas regardless of plot width
    wrapped <- strwrap(sprintf("%s = %s", label, members), width = 85, exdent = 4)
    collapsed_caption_lines <- c(collapsed_caption_lines, paste(wrapped, collapse = "\n"))
  }
}

# -----------------------------------------------------------------------------
# Each displayed column's own Log2FC (COV_TBI vs Control) and significance --
# cell fill is the gene's OWN differential-expression value (direction +
# magnitude), not just a flat "is a member" indicator. Collapsed driver
# columns (e.g. "ATP5F1A +20") use the representative gene's own value, same
# as their display label. Symbol -> Entrez mapping (not a direct symbol join)
# because supp_figure_7c_marker_panel_screen.R found this dataset's raw Entry_Name
# field doesn't always match the gene symbol 1:1 (e.g. Clusterin is "CLUS").
# -----------------------------------------------------------------------------
sym_map <- suppressWarnings(
  AnnotationDbi::select(org.Hs.eg.db, keys = unique(all_columns$gene), keytype = "SYMBOL", columns = "ENTREZID")
) %>%
  dplyr::rename(gene = SYMBOL, entrezgene_id = ENTREZID) %>%
  dplyr::distinct(gene, .keep_all = TRUE)
cov_tbi_de <- read.csv(file.path(input_dir, "pro_HA_results_with_entrez_COV_TBI_vs_Control.csv"), stringsAsFactors = FALSE) %>%
  dplyr::mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  dplyr::distinct(entrezgene_id, .keep_all = TRUE)
gene_de <- sym_map %>%
  dplyr::left_join(cov_tbi_de, by = "entrezgene_id") %>%
  dplyr::select(gene, logFC, adj.P.Val)

not_sig <- gene_de %>% dplyr::filter(is.na(adj.P.Val) | adj.P.Val >= 0.05)
if (nrow(not_sig) > 0) {
  cat("NOTE: these displayed genes are NOT individually significant (adj.P.Val>=0.05) despite being a GSEA leading-edge member:\n")
  print(not_sig)
  sig_caption_line <- sprintf(
    "%d of %d displayed genes are NOT individually significant (adj.P.Val>=0.05) in COV_TBI_vs_Control: %s.",
    nrow(not_sig), nrow(gene_de), paste(not_sig$gene, collapse = ", ")
  )
} else {
  cat(sprintf("All %d displayed genes are individually significant (adj.P.Val<0.05) in COV_TBI_vs_Control.\n", nrow(gene_de)))
  sig_caption_line <- sprintf("All %d displayed genes are individually significant (adj.P.Val<0.05) in COV_TBI_vs_Control.", nrow(gene_de))
}
gene_logfc <- setNames(gene_de$logFC, gene_de$gene)

# -----------------------------------------------------------------------------
# Presence pairs for every displayed column (driver + cross-reference)
# -----------------------------------------------------------------------------
present_pairs <- gene_pathway_pairs %>% filter(gene %in% all_columns$gene) %>% distinct(pathway, gene)

# -----------------------------------------------------------------------------
# Hierarchical clustering, both axes, on the binary pathway x gene presence
# matrix -- see header note
# -----------------------------------------------------------------------------
mat <- matrix(0L, nrow = n_pathways, ncol = nrow(all_columns),
              dimnames = list(all_pathways, all_columns$gene))
for (i in seq_len(nrow(present_pairs))) mat[present_pairs$pathway[i], present_pairs$gene[i]] <- 1L

# Columns (genes): every displayed gene is present in >=1 pathway by
# construction, so no all-zero columns to worry about
gene_hc <- hclust(dist(t(mat), method = "binary"), method = "average")
gene_cluster_order <- colnames(mat)[gene_hc$order]

# Rows (pathways), clustered WITHIN each method facet (keeps the Hallmark/
# KEGG/Reactome grouping). All-zero rows (no displayed gene present) have an
# undefined 0/0 binary distance to each other, so they're pulled out and
# appended at the end of their method block ordered by score instead
cluster_rows_for_method <- function(m) {
  pathways_m <- selected$pathway[selected$method == m]
  sub <- mat[as.character(pathways_m), , drop = FALSE]
  has_any <- rowSums(sub) > 0
  clustered <- character(0)
  if (sum(has_any) > 1) {
    hc <- hclust(dist(sub[has_any, , drop = FALSE], method = "binary"), method = "average")
    clustered <- rownames(sub)[has_any][hc$order]
  } else if (sum(has_any) == 1) {
    clustered <- rownames(sub)[has_any]
  }
  blank <- pathways_m[!(pathways_m %in% clustered)]
  blank_order <- selected %>% filter(pathway %in% blank) %>% arrange(desc(score)) %>% pull(pathway)
  c(clustered, as.character(blank_order))
}
pathway_final_order <- unlist(lapply(levels(selected$method), cluster_rows_for_method))

pathway_method_lookup <- setNames(as.character(selected$method), as.character(selected$pathway))

n_genes <- length(gene_cluster_order)
gene_labels <- all_columns$display_label[match(gene_cluster_order, all_columns$gene)]

# Fixed (not data-max) cap with oob=squish: the driver-block genes are all
# modest (0.2-1.7 Log2FC), but a few cross-reference genes are large outliers
# (USP18 = 8.3) -- capping at the data max would wash the driver genes out to
# near-white, since they'd all sit near the middle of a -8..8 scale.
# Saturating the outliers instead preserves color resolution for the block
# structure that's the actual point of this panel.
gene_cap <- 2.5

method_colors <- c(
  "Hallmark" = "#8A5A1E",
  "KEGG"     = "#1E6B6B",
  "Reactome" = "#6B3FA0"
)

# -----------------------------------------------------------------------------
# Build display matrices for ComplexHeatmap (both already in final row/column
# order -- pathway_final_order top-down, gene_cluster_order left-right --
# from the clustering above, so cluster_rows/cluster_columns = FALSE below)
# -----------------------------------------------------------------------------
nes_idx <- match(pathway_final_order, as.character(selected$pathway))
nes_mat <- matrix(selected$score[nes_idx], ncol = 1, dimnames = list(pathway_final_order, "NES"))
nes_value_label <- selected$value_label[nes_idx]
nes_star_label  <- selected$star_label[nes_idx]

fill_mat <- matrix(NA_real_, nrow = length(pathway_final_order), ncol = n_genes,
                    dimnames = list(pathway_final_order, gene_cluster_order))
for (i in seq_len(nrow(present_pairs))) {
  p <- present_pairs$pathway[i]; g <- present_pairs$gene[i]
  if (g %in% gene_cluster_order) fill_mat[p, g] <- gene_logfc[[g]]
}

pathway_method <- factor(unname(pathway_method_lookup[pathway_final_order]),
                          levels = c("Hallmark", "KEGG", "Reactome"))
# Full pathway names as row labels, wrapped at 50 chars -- with this
# pathway set, only the single 82-character JAK-STAT/IL-12 name wraps (to
# 2 lines); every other name stays on one line.
pathway_row_labels <- str_wrap(pathway_final_order, width = 50)

nes_col_fun  <- circlize::colorRamp2(c(-cap, -cap / 2, 0, cap / 2, cap),
                                      c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D"))
gene_col_fun <- circlize::colorRamp2(c(-gene_cap, -gene_cap / 2, 0, gene_cap / 2, gene_cap),
                                      c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D"))

row_ann <- rowAnnotation(
  Database = pathway_method, col = list(Database = method_colors),
  show_annotation_name = FALSE, simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(Database = list(
    title = "Database", title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 10), grid_height = unit(3.5, "mm"),
    grid_width = unit(3.5, "mm"), direction = "horizontal", nrow = 1
  ))
)

nes_cell_fun <- function(j, i, x, y, width, height, fill) {
  grid.text(nes_value_label[i], x - unit(0.7, "mm"), y, hjust = 1, gp = gpar(fontsize = 7.5, fontface = "bold"))
  grid.text(nes_star_label[i], x + unit(0.7, "mm"), y, hjust = 0, gp = gpar(fontsize = 7.5, fontface = "bold"))
}
# Full pathway names as a dedicated anno_text() column (not Heatmap's own
# row_names_side="left") -- when a rowAnnotation sits to a heatmap's left in
# an ht_list, row_names_side="left" on that heatmap silently renders no text
# at all (reproduced in isolation; not specific to this script/data), so the
# row-name mechanism can't be used here. anno_text() sidesteps it entirely,
# and its multi-line "\n" text (from str_wrap above) renders fine as long as
# the shared row height (set via `height` below) is tall enough for 2 lines.
name_ann <- rowAnnotation(
  Pathway = anno_text(pathway_row_labels, gp = gpar(fontsize = 8.5, fontface = "bold"),
                       just = "right", location = 1),
  show_annotation_name = FALSE
)

# Explicit shared row height (applies to every component in ht_list) sized
# for up to 2 lines of wrapped pathway-name text per row -- without this,
# ComplexHeatmap divides whatever body height the device happens to leave
# evenly across rows, which is fine for 1-line codes but not tall enough for
# wrapped 2-line names.
row_height <- unit(0.4 * length(pathway_final_order), "in")

ht_nes <- Heatmap(
  nes_mat, name = "NES", col = nes_col_fun, rect_gp = gpar(col = "white", lwd = 1.2),
  width = unit(1.5, "cm"), height = row_height, cluster_rows = FALSE, cluster_columns = FALSE,
  show_column_names = FALSE, show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = pathway_method, row_title = NULL,
  cell_fun = nes_cell_fun,
  heatmap_legend_param = list(title = "NES", title_gp = gpar(fontsize = 11, fontface = "bold"),
                               labels_gp = gpar(fontsize = 10), at = seq(-cap, cap, by = 1),
                               legend_width = unit(2.4, "cm"), direction = "horizontal")
)

ht_genes <- Heatmap(
  fill_mat, name = "Gene Log2FC", col = gene_col_fun, na_col = "grey93",
  rect_gp = gpar(col = "white", lwd = 0.3),
  width = unit(0.34 * n_genes, "cm"), height = row_height,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_row_dend = FALSE, show_column_dend = FALSE,
  row_split = pathway_method, row_title = NULL,
  column_labels = gene_labels, column_names_side = "top",
  column_names_rot = 90, column_names_gp = gpar(fontsize = 6.5),
  heatmap_legend_param = list(title = "Gene Log2FC", title_gp = gpar(fontsize = 11, fontface = "bold"),
                               labels_gp = gpar(fontsize = 10), at = seq(-gene_cap, gene_cap, by = 1.25),
                               legend_width = unit(2.4, "cm"), direction = "horizontal")
)

ht_list <- row_ann + name_ann + ht_nes + ht_genes

# No caption on the figure itself; methods/caveat text and the collapsed-
# group full membership belong in the manuscript figure legend instead
# (console output above already has both).

dir.create(tif_dir, showWarnings = FALSE, recursive = TRUE)
pdf_path <- file.path(tif_dir, "supp_figure_7c_pathway_driver_genes_heatmap.pdf")
cairo_pdf(pdf_path, width = 20, height = 10)
draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     merge_legend = TRUE,
     ht_gap = unit(c(1, 1, 3), "mm"), padding = unit(c(2, 4, 2, 2), "mm"),
     column_title = NULL)
dev.off()
message("Saved: ", pdf_path)

tif_path <- sub("\\.pdf$", ".tif", pdf_path)
tiff(tif_path, width = 20, height = 10, units = "in", res = 600, compression = "lzw")
draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     merge_legend = TRUE,
     ht_gap = unit(c(1, 1, 3), "mm"), padding = unit(c(2, 4, 2, 2), "mm"),
     column_title = NULL)
dev.off()
message("Saved: ", tif_path)
