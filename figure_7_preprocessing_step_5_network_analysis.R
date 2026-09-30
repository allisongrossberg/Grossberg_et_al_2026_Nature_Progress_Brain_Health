# =============================================================================
# Figure 7 – Preprocessing step 5: Network Analysis
# =============================================================================
#
# Description: Builds protein interaction networks per contrast from
#   significant proteins and interaction DB. Writes edge lists for figures 7g–7j.
#
# NOTE: Network plots depend on GSEA KEGG results (step 4) which query
#   the KEGG REST API at runtime. Pathway lists may differ between KEGG
#   releases. STRING DB version is pinned at v11.5 (interaction scores
#   unchanged), so differences come solely from KEGG pathway definitions.
#   Original analysis (Aug 2024): KEGG Release 111.0 (2024-07-01)
#   Last modified (Sep 2025):     KEGG Release 115.1 (2025-08-01)
#   Most recent run (Mar 2026):   KEGG Release 117.0 (2026-01-01)
#
# Prerequisites: figure_7_output_data (run preprocessing steps 2, 4 first)
#
# Output: pro_HA_network_edges_*.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

OUTPUT_DIR <- "figure_7_output_data"

# Load all inputs from CSVs (reproducible)
entry_name_entrez_csv <- file.path(OUTPUT_DIR, "pro_HA_entry_name_entrez_ids.csv")
if (!file.exists(entry_name_entrez_csv)) {
  stop("Run figure_7_preprocessing_step_2_DEA.R first. Expected: ", entry_name_entrez_csv)
}
entry_name_entrez_ids <- read.csv(entry_name_entrez_csv, stringsAsFactors = FALSE)
entry_name_entrez_ids$entrezgene_id <- as.character(entry_name_entrez_ids$entrezgene_id)

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

# Build gsea_results_proteomics from GSEA KEGG result + geneset CSVs (no RDS)
gsea_results_proteomics <- list()
kegg_result_files <- list.files(OUTPUT_DIR, pattern = "^pro_HA_GSEA_KEGG_[^g].*\\.csv$", full.names = TRUE)
kegg_result_files <- kegg_result_files[!grepl("genesets_", kegg_result_files)]
for (f in kegg_result_files) {
  nm <- gsub("^pro_HA_GSEA_KEGG_|\\.csv$", "", basename(f))
  result_df <- read.csv(f, stringsAsFactors = FALSE)
  gs_file <- file.path(OUTPUT_DIR, sprintf("pro_HA_GSEA_KEGG_genesets_%s.csv", nm))
  if (file.exists(gs_file)) {
    gs_df <- read.csv(gs_file, stringsAsFactors = FALSE)
    gs_df$entrez_id <- as.character(gs_df$entrez_id)
    geneSets <- split(gs_df$entrez_id, gs_df$kegg_id)
  } else {
    geneSets <- setNames(as.list(rep(list(character(0)), nrow(result_df))), result_df$ID)
  }
  gsea_results_proteomics[[nm]] <- list(KEGG = list(result = result_df, geneSets = geneSets))
}
cat("Loaded GSEA KEGG results and genesets from CSV for ", length(gsea_results_proteomics), " contrasts.\n")

#####
library(STRINGdb)
library(tidyverse)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(reshape2)
library(viridis)
library(scales)
library(ggvenn)

generate_proteomics_network_plot <- function(contrast, min_degree_percentile = 0.9, edge_score_percentile = 0.65, pathway_symbol_size = 5, regulation_symbol_size = 5) {
  results_df <- gsea_results_proteomics[[contrast]]$KEGG$result
  kegg_ids <- results_df$ID
  geneSets <- gsea_results_proteomics[[contrast]]$KEGG$geneSets
  kegg_column <- list()
  entrez_column <- list()
  for (kegg_id in kegg_ids) {
    for (entrez_id in geneSets[[kegg_id]]) {
      kegg_column <- append(kegg_column, kegg_id)
      entrez_column <- append(entrez_column, entrez_id)
    }
  }
  kegg_entrez_df <- data.frame(unlist(kegg_column), unlist(entrez_column))
  names(kegg_entrez_df) <- c("kegg_id", "entrezgene_id")
  kegg_results_df <- kegg_entrez_df %>%
    left_join(results_df, by = c("kegg_id" = "ID"))
  
  kegg_results_df <- kegg_results_df %>%
    left_join(entry_name_entrez_ids, by = "entrezgene_id") %>%
    filter(!is.na(Entry_Name))
  
  kegg_results_df <- kegg_results_df[!kegg_results_df$Description %in% c("Salmonella infection", "Shigellosis", "Spliceosome", "Diabetic cardiomyopathy", "Aminoacyl-tRNA biosynthesis", "Ribosome"), ]
  
  kegg_results_df_final <- kegg_results_df %>%
    left_join(results_with_entrez_list[[contrast]] %>%
                mutate(entrezgene_id = as.character(entrezgene_id)) %>%
                dplyr::select(entrezgene_id, logFC), by = "entrezgene_id")
  
  kegg_results_df_final$direction <- ifelse(kegg_results_df_final$logFC > 0, "up", "down")
  
  string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 250, protocol = "http")
  string_mapped <- string_db$map(kegg_results_df_final, c("entrezgene_id"))
  string_ids <- string_mapped$STRING_id
  string_interactions <- string_db$get_interactions(string_ids)
  string_interactions_filtered <- string_interactions %>%
    filter(combined_score > quantile(combined_score, 0.10))
  
  string_interactions_filtered <- string_interactions_filtered %>%
    left_join(string_mapped %>% dplyr::select(STRING_id, Entry_Name), by = c("from" = "STRING_id")) %>%
    dplyr::rename(from_gene = Entry_Name) %>%
    left_join(string_mapped %>% dplyr::select(STRING_id, Entry_Name), by = c("to" = "STRING_id")) %>%
    dplyr::rename(to_gene = Entry_Name) %>%
    left_join(string_mapped %>% dplyr::select(STRING_id, entrezgene_id), by = c("from" = "STRING_id")) %>%
    dplyr::rename(from_entrezgene_id = entrezgene_id) %>%
    left_join(string_mapped %>% dplyr::select(STRING_id, entrezgene_id), by = c("to" = "STRING_id")) %>%
    dplyr::rename(to_entrezgene_id = entrezgene_id) %>%
    distinct()
  
  g <- graph_from_data_frame(d = string_interactions_filtered, directed = FALSE)
  V(g)$degree <- degree(g)
  
  high_degree_nodes <- V(g)[degree(g) > quantile(degree(g), min_degree_percentile)]
  g_filtered <- induced_subgraph(g, high_degree_nodes)
  
  if (igraph::vcount(g_filtered) == 0) {
    return(list(network = ggplot2::ggplot() + ggplot2::theme_void(), g = g_filtered))
  }
  
  E(g_filtered)$weight <- E(g_filtered)$combined_score
  g_filtered <- delete_edges(g_filtered, E(g_filtered)[weight < quantile(E(g_filtered)$weight, edge_score_percentile)])
  
  if (igraph::vcount(g_filtered) == 0) {
    return(list(network = ggplot2::ggplot() + ggplot2::theme_void(), g = g_filtered))
  }
  
  communities <- cluster_louvain(g_filtered)
  V(g_filtered)$community <- communities$membership
  
  V(g_filtered)$gene_name <- V(g_filtered)$name
  V(g_filtered)$gene_name[V(g_filtered)$name %in% string_interactions_filtered$from] <-
    string_interactions_filtered$from_gene[match(V(g_filtered)$name[V(g_filtered)$name %in% string_interactions_filtered$from], string_interactions_filtered$from)]
  V(g_filtered)$gene_name[V(g_filtered)$name %in% string_interactions_filtered$to] <-
    string_interactions_filtered$to_gene[match(V(g_filtered)$name[V(g_filtered)$name %in% string_interactions_filtered$to], string_interactions_filtered$to)]
  
  V(g_filtered)$pathway <- string_mapped$Description[match(V(g_filtered)$name, string_mapped$STRING_id)]
  
  V(g_filtered)$direction <- string_mapped$direction[match(V(g_filtered)$name, string_mapped$STRING_id)]
  
  V(g_filtered)$node_shape <- ifelse(V(g_filtered)$direction == "up", "Up-regulated", "Down-regulated")
  
  custom_colors <- c("#1E3221", "#395724", "#739D52", "#A2B9FC",
                     "#7583B7", "#DF4416", "#E86D3D", "#FFD35A",
                     "#FCDC94", "#FEFBD8", "#EECEB9", "#BB9AB1",
                     "#D66EB7", "#8C3061", "#522258", "#201E43", "#8B4513")
  set.seed(123)
  network <- ggraph(g_filtered, layout = "stress") +
    geom_edge_link(aes(edge_alpha = combined_score, edge_width = combined_score),
                   edge_colour = "#B4B4B3") +
    geom_node_point(aes(size = degree, color = pathway, shape = node_shape)) +
    geom_node_text(aes(label = gene_name), repel = FALSE, size = 4.5) +
    scale_edge_width_continuous(range = c(0.1, 0.7), name = "Interaction Strength") +
    scale_size_continuous("Degree", range = c(2, 8)) +
    scale_color_manual(values = custom_colors, name = "Pathway") +
    theme_graph(base_family = "Helvetica") +
    theme(legend.position = "right",
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 22)) +
    guides(size = guide_legend(order = 1, title = "Node Degree"),
           edge_width = guide_legend(order = 2, title = "Interaction Strength"),
           edge_alpha = guide_none(),
           color = guide_legend(order = 3, title = "Pathway",
                                override.aes = list(size = pathway_symbol_size)),
           shape = guide_legend(order = 4, title = "Regulation",
                                override.aes = list(size = regulation_symbol_size)))
  
  return(list(network = network, g = g_filtered))
}

# Generate and save network edges for contrasts used by figure_7f and figure_7g–6j
CONTRASTS_FOR_NETWORKS <- c("COV_TBI_vs_Control", "COV_TBI_vs_NegCTRL", "Poly_A_U_vs_NegCTRL", "IL1b_vs_NegCTRL")
dir.create(OUTPUT_DIR, showWarnings = FALSE)
for (contrast in CONTRASTS_FOR_NETWORKS) {
  if (contrast %in% names(gsea_results_proteomics) &&
      !is.null(gsea_results_proteomics[[contrast]]$KEGG) &&
      nrow(gsea_results_proteomics[[contrast]]$KEGG$result) > 0) {
    pl <- generate_proteomics_network_plot(contrast)
    g <- pl$g
    if (igraph::vcount(g) > 0) {
      ed <- as.data.frame(igraph::as_edgelist(g, names = TRUE))
      names(ed) <- c("V1", "V2")
      # Save only edge list columns (V1, V2) for figure_7f/7g–7j
      write.csv(ed %>% dplyr::select(V1, V2), file.path(OUTPUT_DIR, sprintf("pro_HA_network_edges_%s.csv", gsub("[^A-Za-z0-9_]", "_", contrast))), row.names = FALSE)
    }
  }
}
