# Original 2024 analysis objects

Frozen objects extracted from the original analysis workspace 
that produced the published figure_7 KEGG/STRING-dependent
panels: the jaccard similarity heatmap (7f), both venn diagrams (7d, 7e), and
the four protein-interaction network plots (7g-7j).

These figures depend on live KEGG/STRING API queries (see the reproducibility
notes in figure_7_preprocessing_step_4_GSEA.R and
figure_7_preprocessing_step_5_network_analysis.R); pathway membership and
interaction scores drift with database updates, so re-querying today does not
reproduce the exact published figures. This folder freezes the actual data
and rendered objects used, so the pipeline can regenerate the published
figures exactly rather than a same-but-different recomputation.

## Contents

- `network_edges_original/pro_HA_network_edges_*_original.csv` -- true
  original protein-interaction edge lists per contrast (5 contrasts,
  including the empty TNFa_vs_NegCTRL network). Used by figure_7f.R.
- `venn_plot_dea_pro_HA.rds`, `venn_plot_pathway_pro_HA.rds` -- cached ggplot
  objects for the two venn diagrams. Used directly by figure_7d.R/7e.R
  (re-exported as vector PDF) rather than recomputed, since venn circle
  packing isn't guaranteed to reproduce identically across package versions.
- `cov_tbi_control_plot.rds`, `cov_tbi_neg_control_plot.rds`,
  `poly_a_u_neg_control_plot.rds`, `il1b_neg_control_plot.rds` -- cached
  network plot objects (each a list with `$network` ggplot + `$g` igraph).
  Used directly by figure_7g.R-7j.R for the same reason (force-directed
  layout isn't guaranteed reproducible).
- `venn_data.rds`/`figure_7d_venn_gene_sets_original.csv`,
  `venn_data_2.rds`/`figure_7e_venn_pathway_sets_original.csv` -- the raw
  gene/pathway ID sets behind the venn diagrams (verified: set sizes and
  pairwise overlaps match the published diagrams).
- `circos_matrix.rds`/`figure_7c_circos_matrix_original.csv` -- the
  condition x pathway enrichment-score matrix behind the circos plot
  (figure_7c currently uses a direct TIFF conversion; not yet wired into a
  script).
- `kegg_list.rds`, `kegg_df.rds`, `kegg_category.rds`,
  `network_list.rds` -- supporting KEGG pathway metadata and the full
  network object list, kept for reference.
