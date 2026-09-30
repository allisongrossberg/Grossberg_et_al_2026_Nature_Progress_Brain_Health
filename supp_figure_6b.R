# =============================================================================
# Supplementary Figure 6b – Reactive-astrocyte morphology heatmap, HA TLR
#   Inhibitor experiment (process thickness, length, solidity, circularity)
# =============================================================================
# Reads per-cell shape/skeleton/thickness (from the process-preserving
# branch mask) and hypertrophy area (from the area mask) in
# morphology_v2_dapiseeded_results.csv, plus per-image whole-field GFAP/
# Vimentin MFI normalized by DAPI count (mfi_per_image_both_channels.csv,
# the same Normed_MFI used in the published Fig 3-5).
#
# The 4 panels describe the reactive process transition (thin/diffuse ->
# thick/consolidated) seen in the images:
#   1. Process thickness = median per-cell process_median_thickness_um (caliber)
#   2. Process length    = median per-cell skeleton_length_um (extent)
#   3. Solidity          = median per-cell solidity (compactness/density)
#   4. Circularity       = median per-cell circularity (roundness)
# All four have large effects vs. PBS (Cohen's d ~1.0-1.3). Total cell area
# is not a panel: diffuse resting cells spread thin over a large footprint,
# so area doesn't track hypertrophy (d ~0.5 only); area, whole-cell
# thickness, branch count, and MFI are written to the stats CSV as
# supplementary instead.
#
# Output: supp_figure_6b_hypertrophy_heatmap_<TAG>.pdf  (TAG default v2corrected)
#         supp_figure_6b_output_data/supp_figure_6b_<TAG>_HA_TLR_means.csv
#         supp_figure_6b_output_data/supp_figure_6b_<TAG>_HA_TLR_vs_control_holm.csv
#         supp_figure_6b_output_data/supp_figure_6b_<TAG>_HA_TLR_consolidation_redundancy.csv
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))
if (basename(getwd()) == "supp_figure_6") setwd("..")
output_dir <- "supp_figure_6/supp_figure_6_output_data"
input_dir  <- "supp_figure_6/supp_figure_6_input_data"
tif_dir    <- "supp_figure_6"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(rstatix); library(seriation); library(patchwork)
})

# Env overrides let the SAME script make both the ungated (primary) and the
# soma-gated (supplement) versions. Defaults = ungated primary.
cells_path <- Sys.getenv("COMMENT_D_CELLS_CSV", unset = file.path(input_dir, "morphology_v2_dapiseeded_results.csv"))
TAG        <- Sys.getenv("COMMENT_D_TAG", unset = "v2corrected")
mfi_path   <- file.path(input_dir, "mfi_per_image_both_channels.csv")

cells <- read_csv(cells_path, show_col_types = FALSE)
mfi   <- read_csv(mfi_path,   show_col_types = FALSE)

# QC: drop physically implausible objects (failed segmentations). A real
# astrocyte's local thickness (largest inscribed circle = soma diameter) is
# at most ~15-20 um; a handful of images have merged-blob segmentation
# failures with per-cell median thickness up to ~160-200 um, which inflate
# the mean thickness of a few conditions without being statistically
# significant themselves. Capping at 20 um (just above the GFAP p99 of 18
# um) removes <1% (GFAP) / <2% (Vimentin) of cells. NA-thickness cells
# (real thin cells with no resolvable skeleton) are kept and contribute to
# area only.
THICKNESS_CAP_UM <- 20
n0 <- nrow(cells)
cells <- cells %>% filter(is.na(whole_cell_median_thickness_um) |
                          whole_cell_median_thickness_um <= THICKNESS_CAP_UM)
cat(sprintf("QC thickness cap %g um: dropped %d of %d cells (%.2f%%)\n",
            THICKNESS_CAP_UM, n0 - nrow(cells), n0, 100 * (n0 - nrow(cells)) / n0))
cat(sprintf("Corrected HA_TLR: %d images, %d cells\n", n_distinct(cells$image), nrow(cells)))

# per-image aggregation (median across cells), then join whole-field MFI
per_image <- cells %>%
  group_by(image, channel) %>%
  summarise(
    proc_thickness = median(process_median_thickness_um, na.rm = TRUE),  # caliber
    skeleton    = median(skeleton_length_um, na.rm = TRUE),   # process length
    solidity    = median(solidity, na.rm = TRUE),             # compactness
    circularity = median(circularity, na.rm = TRUE),          # roundness
    # kept for the stats CSV only (not heatmap panels): area is confounded by
    # diffuse spread in controls; whole-cell thickness / branch count / MFI
    area        = median(area_um2, na.rm = TRUE),
    thickness   = median(whole_cell_median_thickness_um, na.rm = TRUE),
    branches    = median(n_branches, na.rm = TRUE),
    n_cells     = n(), .groups = "drop"
  ) %>%
  left_join(mfi, by = c("image", "channel")) %>%
  rename(mfi = mfi_normed)

# ---- condition parsing / exclusions ----
assign_cond_ha_tlr <- function(df) {
  df <- df %>%
    mutate(Image_File_Clean_2 = str_replace(image, "TH1020_Exo_Cont", "TH1020_Cont")) %>%
    mutate(Image_File_Clean_2 = str_remove(Image_File_Clean_2, "^COAST_HA_ADE_"))
  df %>% mutate(Condition = sub("(_[^_]+){2}$", "", Image_File_Clean_2),
                Condition = gsub("_[0-9]+$", "", Condition))
}
REFERENCE_CONDITION <- "Neg_Cont"
EXCLUDED_PREFIXES <- c("TL2_C29", "CLI_095", "TH1020")
excluded_regex <- paste0("^(", paste(EXCLUDED_PREFIXES, collapse = "|"), ")")
per_image <- per_image %>% assign_cond_ha_tlr() %>%
  filter(Condition != "PEG", !str_detect(Condition, "_Veh$"), !str_detect(Condition, excluded_regex))

ALL_METRICS     <- c("proc_thickness", "skeleton", "solidity", "circularity",
                     "area", "thickness", "mfi", "branches")
# See header for panel selection rationale.
HEATMAP_METRICS <- c("proc_thickness", "skeleton", "solidity", "circularity")

# ---- aggregate technical replicates to the experimental unit ----
# The 10 images in a well are technical replicates (fields); the WELL is the
# independent unit (n = 6 wells per condition). Statistics are computed on
# well-level means, NOT on images, to avoid pseudoreplication. Well id =
# condition + the first trailing file-name number (..._<well>_<field>_AG);
# the second number is the field within the well.
per_image <- per_image %>%
  mutate(block = str_match(image, "_(\\d+)_\\d+_AG$")[, 2],
         well  = ifelse(is.na(block), image, paste0(Condition, "_", block)))
per_well <- per_image %>%
  group_by(Condition, well, channel) %>%
  summarise(across(all_of(ALL_METRICS), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# consolidation-axis redundancy (well level)
redundancy <- per_well %>% group_by(channel) %>%
  summarise(sol_circ  = cor(solidity, circularity, use = "complete.obs"),
            circ_skel = cor(circularity, skeleton, use = "complete.obs"),
            sol_skel  = cor(solidity, skeleton, use = "complete.obs"),
            n_wells = n(), .groups = "drop")
write.csv(redundancy, file.path(output_dir, paste0("supp_figure_6b_", TAG, "_HA_TLR_consolidation_redundancy.csv")), row.names = FALSE)

# significance: one-way ANOVA + Tukey HSD post-hoc, WITHIN each experimental
# FAMILY (matching the MFI figures 3-5), at the WELL level (n = 6). A single
# ANOVA across all 21 conditions is inappropriate -- it dilutes the omnibus and
# over-corrects the contrast of interest. Each family shares the PBS control.
# Star = Tukey-adjusted P vs PBS within that condition's family.
conds_present <- unique(per_well$Condition)
FAMILIES <- list(
  ADE       = intersect(c("Neg_Cont", "Exo_Cont", "Exo_COV_Only", "Exo_TBI_Only", "Exo_COV_TBI"), conds_present),
  PosCtrl   = intersect(c("Neg_Cont", "LPS", "Zymosan", "Poly_AU", "FLA_BS", "ODN_DSL03"), conds_present),
  CU_CPT4a  = c("Neg_Cont", grep("^CU_CPT4a", conds_present, value = TRUE)),
  ODN_INH18 = c("Neg_Cont", grep("^ODN_INH18", conds_present, value = TRUE)))

run_family_tukey <- function(raw, metric_col, fam, reference = "Neg_Cont") {
  raw <- raw %>% filter(Condition %in% fam, is.finite(.data[[metric_col]]))
  raw$Condition <- factor(raw$Condition)
  if (nlevels(raw$Condition) < 3 || !reference %in% levels(raw$Condition)) return(tibble())
  fit <- aov(reformulate("Condition", response = metric_col), data = raw)
  aov_p <- summary(fit)[[1]][["Pr(>F)"]][1]
  tk <- as.data.frame(TukeyHSD(fit)$Condition); tk$pair <- rownames(tk)
  tk %>% filter(str_detect(pair, paste0("(^|-)", reference, "($|-)"))) %>%
    mutate(Condition = str_remove(pair, paste0("-?", reference, "-?")),
           p.adj = `p adj`, aov_p = aov_p) %>%
    filter(Condition != "") %>% select(Condition, p.adj, aov_p)
}

all_summaries <- list(); all_vs_control <- list()
for (ch in c("GFAP", "Vimentin")) {
  raw_ch <- per_well %>% filter(channel == ch)
  for (m in ALL_METRICS) {
    summ <- raw_ch %>% filter(is.finite(.data[[m]])) %>% group_by(Condition) %>%
      summarise(Mean = mean(.data[[m]]), SD = sd(.data[[m]]),
                SE = sd(.data[[m]])/sqrt(n()), n_wells = n(), .groups = "drop")
    all_summaries[[paste(ch, m)]] <- summ %>% mutate(channel = ch, metric = m)
    for (fname in names(FAMILIES)) {
      vc <- run_family_tukey(raw_ch, m, FAMILIES[[fname]])
      if (nrow(vc) > 0)
        all_vs_control[[paste(ch, m, fname)]] <- vc %>% mutate(channel = ch, metric = m, family = fname)
    }
  }
}
means <- bind_rows(all_summaries); vs_control <- bind_rows(all_vs_control)
write.csv(means, file.path(output_dir, paste0("supp_figure_6b_", TAG, "_HA_TLR_means.csv")), row.names = FALSE)
write.csv(vs_control, file.path(output_dir, paste0("supp_figure_6b_", TAG, "_HA_TLR_vs_control_tukey.csv")), row.names = FALSE)

ha_means <- means %>% filter(metric %in% HEATMAP_METRICS)
sig_lookup <- vs_control %>% filter(metric %in% HEATMAP_METRICS) %>%
  transmute(channel, metric, Condition, adj.p.value = p.adj)
sig_stars <- function(p) ifelse(is.na(p), "",
  ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))))

d <- ha_means %>% left_join(sig_lookup, by = c("channel", "metric", "Condition")) %>%
  mutate(sig_label = sig_stars(adj.p.value), metric_channel = paste(metric, channel, sep = "_"))

format_condition_label <- function(x) {
  x <- ifelse(x == "Neg_Cont", "PBS Control", x)
  x <- gsub("_COV_TBI$",  "-COVID-19 (+) mTBI (+) ADEs", x)
  x <- gsub("_COV_Only$", "-COVID-19 (+) mTBI (-) ADEs", x)
  x <- gsub("_TBI_Only$", "-COVID-19 (-) mTBI (+) ADEs", x)
  x <- gsub("_Cont$",     "-COVID-19 (-) mTBI (-) ADEs", x)
  x <- gsub("^Exo-", "", x); x <- gsub("_", "-", x); x
}
d <- d %>% mutate(Condition_label = format_condition_label(Condition))

wide_z <- d %>% group_by(metric_channel) %>% mutate(z = as.numeric(scale(Mean))) %>% ungroup() %>%
  select(Condition, metric_channel, z) %>% pivot_wider(names_from = metric_channel, values_from = z)
mat <- as.matrix(wide_z %>% select(-Condition)); rownames(mat) <- wide_z$Condition; mat[is.na(mat)] <- 0
dist_mat <- dist(mat, method = "euclidean"); hc <- hclust(dist_mat, method = "complete")
row_order <- rownames(mat)[seriation::get_order(seriation::seriate(dist_mat, method = "OLO", control = list(hclust = hc)))]
row_order_labels <- format_condition_label(row_order)

d <- d %>% mutate(
  Condition_label = factor(Condition_label, levels = row_order_labels),
  metric = factor(metric, levels = c("proc_thickness", "skeleton", "solidity", "circularity"),
                  labels = c("Process thickness", "Process length", "Solidity", "Circularity")),
  channel = factor(channel, levels = c("GFAP", "Vimentin")))

# Circularity is 0-1 -- needs 2 decimals; the tens/hundreds-scale metrics use 1.
fmt_val <- function(x, mname = NULL) {
  digits <- ifelse(!is.null(mname) && mname %in% c("Circularity", "Solidity"), 2, 1)
  r <- round(x, digits); ifelse(r == 0, sprintf("%.*f", digits, 0), sprintf("%.*f", digits, r))
}
d <- d %>% rowwise() %>% mutate(label = paste0(fmt_val(Mean, as.character(metric)), sig_label)) %>% ungroup()
# Color is normalized WITHIN each metric x channel (not pooled across both
# channels) so e.g. Vimentin's high-thickness outliers don't compress the GFAP
# column into a flat blue. Printed cell values are the real numbers; the two
# channel columns of a panel are each on their own min-max, so the colorbar is
# shown as a relative "lower -> higher (within marker)" scale.
d <- d %>% group_by(metric, channel) %>%
  mutate(norm = (Mean - min(Mean, na.rm = TRUE)) / (max(Mean, na.rm = TRUE) - min(Mean, na.rm = TRUE))) %>% ungroup()
metric_range <- d %>% group_by(metric) %>% summarise(mn = min(Mean, na.rm = TRUE), mx = max(Mean, na.rm = TRUE), .groups = "drop")

metric_colors <- c("Process thickness" = "#2E7D46", "Process length" = "#6B3FA0",
                    "Solidity" = "#1E6B6B", "Circularity" = "#B8860B")
metric_units  <- c("Process thickness" = "Process thickness (µm)", "Process length" = "Process length (µm)",
                   "Solidity" = "Solidity (0-1)", "Circularity" = "Circularity (0-1)")
# "reactive" direction: process thickness / solidity / circularity are HIGHER in
# reactive, but process length is higher in CONTROLS -- invert only that panel's
# fill so red consistently means "more reactive" across all four (printed values
# untouched).
INVERT_METRICS <- c("Process length")

# Highlight the TLR-agonist positive controls (red) and the COVID+/mTBI+ ADE
# condition itself (blue) in the row labels, to call out that COVID+/mTBI+
# clusters with the positive controls. Vector is in factor-level (= y-axis)
# order so ggplot maps each color to the right tick.
TREATED_CONDS <- c("ODN_INH18_COV_TBI", "CU_CPT4a_COV_TBI")  # COV/TBI + TLR inhibitor
treated_labels <- format_condition_label(TREATED_CONDS)
covtbi_label   <- format_condition_label("Exo_COV_TBI")
ylab_levels    <- levels(d$Condition_label)
ylab_colors    <- rep("black", length(ylab_levels))   # all row labels black

# y positions (discrete axis: level 1 = bottom) for the side annotations.
covtbi_pos <- match(covtbi_label, ylab_levels)
t_pos      <- sort(match(treated_labels, ylab_levels)); t_pos <- t_pos[!is.na(t_pos)]
# the transition set: COV/TBI (untreated, reactive) + its two TLR-inhibitor-
# treated versions (rescued, up in the control-like zone)
trans_pos  <- sort(c(covtbi_pos, t_pos))
# the two top controls to bracket
ctrl_pos   <- sort(match(format_condition_label(c("Neg_Cont", "Exo_Cont")), ylab_levels))
# the bottom reactive cluster (TLR-agonist positive controls + untreated COV/TBI)
POS_CONTROL_CONDS <- c("LPS", "Zymosan", "Poly_AU", "FLA_BS", "ODN_DSL03")
r_pos <- sort(match(format_condition_label(c(POS_CONTROL_CONDS, "Exo_COV_TBI")), ylab_levels))
r_pos <- r_pos[!is.na(r_pos)]

# All 4 heatmap panels share IDENTICAL geometry (same widths, same margins, no
# per-panel coordinate tricks) so every strip title lines up exactly over its
# two tile columns. The brackets live in a separate 5th "annotation" panel
# (below) that shares the same discrete y-scale, so its rows line up with the
# heatmap rows without having to distort any real panel's coordinate space.
build_panel <- function(mname, show_y) {
  sub <- d %>% filter(metric == mname); rng <- metric_range %>% filter(metric == mname)
  mid <- (rng$mn + rng$mx) / 2
  inverted <- mname %in% INVERT_METRICS
  sub <- sub %>% mutate(fill_val = if (inverted) 1 - norm else norm)
  # color is normalized within each marker column, so the bar is relative;
  # exact values are printed in every tile. Red end = "more reactive".
  if (inverted) {           # Process length: red = shorter = reactive
    break_labels <- c("longer", "", "shorter")
  } else {                  # red = higher = reactive
    break_labels <- c("lower", "", "higher")
  }
  g <- ggplot(sub, aes(x = channel, y = Condition_label, fill = fill_val)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = label), color = "black", fontface = "bold", size = 4.0) +
    scale_fill_gradientn(colors = c("#08306B", "#4292C6", "#F7F7F7", "#EF6548", "#99000D"),
      values = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1), breaks = c(0, 0.5, 1),
      labels = break_labels, name = metric_units[mname]) +
    guides(fill = guide_colorbar(title.position = "top", direction = "horizontal",
      barwidth = unit(3.4, "cm"), barheight = unit(0.35, "cm"))) +
    facet_wrap(~ metric) + labs(x = NULL, y = NULL) + theme_minimal(base_size = 16) +
    theme(strip.background = element_rect(fill = metric_colors[mname], color = NA),
      strip.text = element_text(color = "white", face = "bold", size = 13),
      panel.grid = element_blank(),
      axis.text.y = if (show_y) element_text(size = 12, face = "bold", color = ylab_colors) else element_blank(),
      axis.text.x = element_text(size = 15, face = "bold", color = "gray15"),
      legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9, face = "bold"),
      plot.margin = margin(10, 8, 10, 8))

  # outline the COV/TBI (untreated) + its two TLR-inhibitor-treated rows in
  # every panel, so the reader sees which cells the transition bracket links
  for (yp in trans_pos) {
    g <- g + annotate("rect", xmin = 0.5, xmax = 2.5, ymin = yp - 0.5, ymax = yp + 0.5,
                      fill = NA, color = "black", linewidth = 0.8)
  }
  g
}

# Annotation panel: same discrete y-scale (Condition_label, same level order)
# as the heatmap panels, plus an invisible fill legend of the SAME kind/size
# as the real ones (blank colors/labels) so patchwork reserves identical
# vertical space for it -- this is what keeps its rows aligned with the four
# real heatmap panels without any manual coordinate tweaking.
ann_df <- tibble(Condition_label = factor(ylab_levels, levels = ylab_levels), x = 0, v = 0)
build_annotation_panel <- function() {
  g <- ggplot(ann_df, aes(x = x, y = Condition_label, fill = v)) +
    geom_blank() +
    scale_x_continuous(limits = c(0, 4.2), expand = c(0, 0)) +
    scale_fill_gradientn(colors = c("white", "white"), limits = c(0, 1),
      breaks = c(0, 0.5, 1), labels = c("", "", ""), name = " ") +
    guides(fill = guide_colorbar(title.position = "top", direction = "horizontal",
      barwidth = unit(3.4, "cm"), barheight = unit(0.35, "cm"))) +
    facet_wrap(~ "") + labs(x = NULL, y = NULL) + theme_minimal(base_size = 16) +
    theme(strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(color = "white", face = "bold", size = 13),
      panel.grid = element_blank(), panel.background = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),
      legend.background = element_blank(), legend.text = element_blank(),
      legend.title = element_blank(),
      plot.margin = margin(10, 8, 10, 0)) +
    coord_cartesian(clip = "off")

  # bracket around the two top controls (PBS + COVID-/mTBI-)
  c_lo <- min(ctrl_pos) - 0.45; c_hi <- max(ctrl_pos) + 0.45
  g <- g +
    annotate("segment", x = 0.12, xend = 0.12, y = c_lo, yend = c_hi, linewidth = 0.9, color = "black") +
    annotate("segment", x = 0.05, xend = 0.12, y = c_lo, yend = c_lo, linewidth = 0.9, color = "black") +
    annotate("segment", x = 0.05, xend = 0.12, y = c_hi, yend = c_hi, linewidth = 0.9, color = "black") +
    annotate("text", x = 0.20, y = (c_lo + c_hi) / 2, hjust = 0, vjust = 0.5,
             size = 3.6, fontface = "bold", color = "black", label = "Controls")
  # ] bracket linking COV/TBI (untreated) to both TLR-inhibitor-treated rows
  t_lo <- min(trans_pos) - 0.45; t_hi <- max(trans_pos) + 0.45
  g <- g +
    annotate("segment", x = 0.12, xend = 0.12, y = t_lo, yend = t_hi, linewidth = 0.9, color = "black") +
    annotate("segment", x = 0.0, xend = 0.12, y = trans_pos[1], yend = trans_pos[1], linewidth = 0.9, color = "black") +
    annotate("segment", x = 0.0, xend = 0.12, y = trans_pos[2], yend = trans_pos[2], linewidth = 0.9, color = "black") +
    annotate("segment", x = 0.0, xend = 0.12, y = trans_pos[3], yend = trans_pos[3], linewidth = 0.9, color = "black") +
    annotate("text", x = 0.20, y = (t_lo + t_hi) / 2, hjust = 0, vjust = 0.5,
             size = 3.4, fontface = "bold", lineheight = 0.95, color = "black",
             label = "transition back\nto a control-\nlike state")
  # outer bracket: the bottom reactive cluster (COV/TBI groups with positive controls)
  r_lo <- min(r_pos) - 0.45; r_hi <- max(r_pos) + 0.45
  g <- g +
    annotate("segment", x = 1.55, xend = 1.55, y = r_lo, yend = r_hi, linewidth = 0.9, color = "black") +
    annotate("segment", x = 1.48, xend = 1.55, y = r_lo, yend = r_lo, linewidth = 0.9, color = "black") +
    annotate("segment", x = 1.48, xend = 1.55, y = r_hi, yend = r_hi, linewidth = 0.9, color = "black") +
    annotate("text", x = 1.63, y = (r_lo + r_hi) / 2, hjust = 0, vjust = 0.5,
             size = 3.4, fontface = "bold", lineheight = 0.95, color = "black",
             label = "Reactive —\nCOVID+/mTBI+\nADEs cluster\nwith the TLR-\nagonist positive\ncontrols")
  g
}

metric_levels <- levels(d$metric)
panels <- lapply(seq_along(metric_levels), function(i)
  build_panel(metric_levels[i], show_y = (i == 1)))
panels[[length(panels) + 1]] <- build_annotation_panel()
heatmap_plot <- wrap_plots(panels, nrow = 1, widths = c(1.18, 1, 1, 1, 0.85))
out_path <- file.path(tif_dir, paste0("supp_figure_6b_hypertrophy_heatmap_", TAG, ".pdf"))
ggsave(out_path, plot = heatmap_plot, width = 17.5, height = 16, dpi = 600, bg = "white")
message("Saved: ", out_path)
message("consolidation-axis redundancy (solidity/circularity/skeleton):"); print(redundancy)
