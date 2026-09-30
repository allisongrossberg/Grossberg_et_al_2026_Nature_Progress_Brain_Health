# =============================================================================
# Figure 7f – Jaccard similarity heatmap (network comparison)
# =============================================================================
#
# Description: Heatmap of Jaccard similarity between protein interaction
#   networks across contrasts. Requires preprocessing steps 1–5.
#
# REPRODUCIBILITY NOTE: Jaccard scores depend on network structure, which
#   depends on GSEA KEGG results (step 4) querying the KEGG REST API.
#   Pathway-gene memberships change between KEGG releases, so similarity
#   scores may differ slightly across runs. The original published figures
#   used KEGG Release 111.0 (Aug 2024). See step 4/5 headers for details.
#
# Prerequisites: figure_7_output_data (run preprocessing steps 1–5 first)
#
# Output: figure_7f_jaccard_similarity_heatmap.pdf
# =============================================================================
cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI

suppressMessages({
  library(showtext)
  sysfonts::font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf", bold = "/System/Library/Fonts/Supplemental/Arial Bold.ttf", italic = "/System/Library/Fonts/Supplemental/Arial Italic.ttf", bolditalic = "/System/Library/Fonts/Supplemental/Arial Bold Italic.ttf")
  showtext_auto()
  showtext_opts(dpi = 600)
})
if (basename(getwd()) != "figure_7") {
  if (file.exists("figure_7")) setwd("figure_7")
  else if (file.exists(file.path("..", "figure_7"))) setwd(file.path("..", "figure_7"))
}

library(ggplot2)
library(reshape2)
library(viridis)
library(scales)
library(tidyverse)
library(igraph)

OUTPUT_DIR <- "figure_7_output_data"
# Original heatmap uses 4 comparisons (including TNF-alpha, whose network has
# zero significant proteins/edges, giving 0% similarity by construction)
ORIGINAL_EDGES_DIR <- "figure_7_input_data/original_2024_objects/network_edges_original"
CONTRASTS_FOR_7F <- c("COV_TBI_vs_Control", "COV_TBI_vs_NegCTRL", "Poly_A_U_vs_NegCTRL", "IL1b_vs_NegCTRL", "TNFa_vs_NegCTRL")

# Load graph from the true original 2024 network edge and node data (extracted
# from the original analysis project's cached network graphs -- see
# figure_7_input_data/original_2024_objects/README for provenance). The
# reproducible pipeline's own preprocessing step 5 output
# (figure_7_output_data/pro_HA_network_edges_*.csv) has fewer edges than these
# originals for 3 of 4 contrasts (STRING/KEGG version drift since 2024), so
# this figure specifically uses the frozen original edge lists to match what
# was actually published. The node list is loaded separately (not just
# derived from the edge list) because some networks include isolated nodes
# -- significant proteins with no STRING interaction above the confidence
# threshold, kept as vertices with no edges -- which a plain edge-list graph
# would silently drop, undercounting node similarity (e.g. IL1b_vs_NegCTRL
# has 6 such isolated nodes).
load_network_from_file <- function(contrast) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", contrast)
  f <- file.path(ORIGINAL_EDGES_DIR, paste0("pro_HA_network_edges_", safe_name, "_original.csv"))
  nf <- file.path(ORIGINAL_EDGES_DIR, paste0("pro_HA_network_nodes_", safe_name, "_original.csv"))
  if (!file.exists(f) || !file.exists(nf)) {
    stop("Missing ", f, " or ", nf, ".")
  }
  ed <- read.csv(f, stringsAsFactors = FALSE)
  nodes <- read.csv(nf, stringsAsFactors = FALSE)
  if (nrow(nodes) == 0) return(igraph::make_empty_graph())
  if (nrow(ed) == 0 || ncol(ed) < 2) {
    return(igraph::graph_from_data_frame(data.frame(from = character(0), to = character(0)), directed = FALSE, vertices = nodes))
  }
  igraph::graph_from_data_frame(ed[, 1:2], directed = FALSE, vertices = nodes)
}

graphs <- setNames(lapply(CONTRASTS_FOR_7F, load_network_from_file), CONTRASTS_FOR_7F)

# Original COAST_Pro_HA_6 Jaccard formula (no edge canonicalization)
calculate_overlap <- function(g1, g2) {
  nodes1 <- V(g1)$name
  nodes2 <- V(g2)$name
  edges1 <- as.data.frame(igraph::as_edgelist(g1, names = TRUE))
  edges2 <- as.data.frame(igraph::as_edgelist(g2, names = TRUE))
  names(edges1) <- c("V1", "V2")
  names(edges2) <- c("V1", "V2")
  edges1[] <- lapply(edges1, as.character)
  edges2[] <- lapply(edges2, as.character)
  node_jaccard <- length(intersect(nodes1, nodes2)) / length(union(nodes1, nodes2))
  edge_union_n <- nrow(dplyr::bind_rows(edges1, edges2) %>% dplyr::distinct())
  edge_jaccard <- if (edge_union_n > 0) nrow(dplyr::inner_join(edges1, edges2, by = c("V1", "V2"))) / edge_union_n else 0
  return(list(node_similarity = node_jaccard, edge_similarity = edge_jaccard))
}

cov_tbi_poly_au <- calculate_overlap(graphs[["COV_TBI_vs_NegCTRL"]], graphs[["Poly_A_U_vs_NegCTRL"]])
cov_tbi_il1b <- calculate_overlap(graphs[["COV_TBI_vs_NegCTRL"]], graphs[["IL1b_vs_NegCTRL"]])
cov_tbi_conts <- calculate_overlap(graphs[["COV_TBI_vs_Control"]], graphs[["COV_TBI_vs_NegCTRL"]])
cov_tbi_tnfa <- calculate_overlap(graphs[["COV_TBI_vs_NegCTRL"]], graphs[["TNFa_vs_NegCTRL"]])

# Original COAST_Pro_HA_6: 4 rows, including TNF-alpha (its network has zero
# significant proteins, so similarity to it is 0% by construction)
similarity_scores <- data.frame(
  comparison = c(
    "COVID-19 (+) mTBI (+) ADEs vs PBS Control vs Poly A:U vs PBS Control ",
    "COVID-19 (+) mTBI (+) ADEs vs PBS Control vs IL-1β vs PBS Control",
    "COVID-19 (+) mTBI (+) ADEs vs PBS Control vs COVID-19 (+) mTBI (+) ADEs vs COVID-19 (-) mTBI (-) ADEs",
    "COVID-19 (+) mTBI (+) ADEs vs PBS Control vs TNF-α vs PBS Control"
  ),
  node_similarity = c(cov_tbi_poly_au$node_similarity, cov_tbi_il1b$node_similarity, cov_tbi_conts$node_similarity, cov_tbi_tnfa$node_similarity),
  edge_similarity = c(cov_tbi_poly_au$edge_similarity, cov_tbi_il1b$edge_similarity, cov_tbi_conts$edge_similarity, cov_tbi_tnfa$edge_similarity)
)

# Original: pivot node_similarity, edge_similarity; order by mean descending
plot_data <- similarity_scores %>%
  tidyr::pivot_longer(cols = c(node_similarity, edge_similarity),
                      names_to = "similarity_type",
                      values_to = "score") %>%
  group_by(comparison) %>%
  mutate(mean_similarity = mean(score)) %>%
  ungroup() %>%
  mutate(comparison = fct_reorder(comparison, mean_similarity, .desc = TRUE))

# Define the function to add line breaks
add_line_breaks <- function(x, max_length = 30) {
  sapply(x, function(y) {
    words <- str_split(y, "\\s+")[[1]]
    lines <- c()
    current_line <- words[1]
    
    for (i in 2:length(words)) {
      if (nchar(current_line) + nchar(words[i]) + 1 <= max_length) {
        current_line <- paste(current_line, words[i])
      } else {
        lines <- c(lines, current_line)
        current_line <- words[i]
      }
    }
    
    lines <- c(lines, current_line)
    paste(lines, collapse = "\n")
  })
}

# Create the plot
jaccard_similarity_heatmap <- ggplot(plot_data, aes(x = comparison, y = similarity_type, fill = score)) +
  geom_tile() +
  geom_text(aes(label = percent(score, accuracy = 0.1)), color = "white", size = 8.5) +
  scale_fill_gradient(low = "#F0E7DB", high = "#000000", labels = percent)+
  theme_minimal() +
  coord_flip()+
  scale_y_discrete(labels = c("node_similarity" = "Node Similarity", 
                              "edge_similarity" = "Edge Similarity")) +
  scale_x_discrete(labels = add_line_breaks) +
  labs(x = "Comparison", y = "Similarity Type", fill = "Jaccard\nSimilarity") +
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        legend.key.size = unit(1, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.key.width = unit(0.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 22, face = "bold"),
  ) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 12))  
jaccard_similarity_heatmap

dir.create(OUTPUT_DIR, showWarnings = FALSE)
write.csv(similarity_scores, file.path(OUTPUT_DIR, "figure_7f_jaccard_similarity_data.csv"), row.names = FALSE)

path_output <- "figure_7f_jaccard_similarity_heatmap.pdf"
ggsave(path_output, plot = jaccard_similarity_heatmap, width = 12, height = 20, dpi = 600)
