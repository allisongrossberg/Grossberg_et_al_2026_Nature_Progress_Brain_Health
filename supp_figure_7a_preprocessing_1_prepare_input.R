# =============================================================================
# Supplementary Figure 7a preprocessing 1 - Build MetaboAnalystR input table
#   from the raw HA metabolomics Global Report, restricted to the two
#   ADE-exposure groups: Cont Exo (COVID-19-/mTBI- control ADEs) vs
#   COV/TBI Exo (COVID-19+/mTBI+ ADEs).
# =============================================================================
#
# Description: Source data is the raw LC-MS peak intensity table
#   (untargeted global metabolomics, human primary astrocytes) from the
#   original vendor report; this script reshapes it into the
#   sample-by-metabolite matrix MetaboAnalystR expects, to compute the
#   underlying statistics for the glycolysis/gluconeogenesis claim
#   (2/3-phospho-D-glycerate and PEP increased in COV/TBI Exo vs Cont Exo
#   at 24h).
#
# Input: AG-68_Grossberg_Astrocytes_11June2024_Final Global Report.xlsx
#   ("Global Report" sheet) - see supp_figure_7a_input_data/ for a copy.
# Output: supp_figure_7a_input_data/supp_figure_7a_peak_table.csv (rows = samples,
#   col 1 = Sample, col 2 = Label/group, remaining cols = raw peak
#   intensities per named metabolite)
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")
input_dir <- "supp_figure_7/supp_figure_7_input_data"
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)

library(readxl)

src_path <- file.path(input_dir, "AG-68_Grossberg_Astrocytes_11June2024_Final Global Report.xlsx")
if (!file.exists(src_path)) stop("Need AG-68_Grossberg_Astrocytes_11June2024_Final Global Report.xlsx in ", input_dir)

gr <- read_excel(src_path, sheet = "Global Report", col_names = FALSE)

# Row 3 = group labels, row 4 = AG-ID sample codes, cols 8:37 = the 30
# biological samples (10 groups x 3 replicates); cols 38+ are QC
# (blanks/tech-mix/%CV), excluded. Rows 5:146 = the 142 named metabolites.
group_row <- as.character(unlist(gr[3, 8:37]))
sample_id <- as.character(unlist(gr[4, 8:37]))
compound  <- as.character(unlist(gr[5:146, 2]))
peak_mat  <- as.data.frame(gr[5:146, 8:37])
peak_mat  <- as.data.frame(lapply(peak_mat, as.numeric))
stopifnot(!anyNA(peak_mat))
rownames(peak_mat) <- compound
colnames(peak_mat) <- sample_id

keep_groups <- c("Cont Exo", "COV/TBI Exo")
keep <- group_row %in% keep_groups
stopifnot(sum(keep) == 6)  # 3 reps x 2 groups

samp_meta <- data.frame(
  Sample = sample_id[keep],
  # Control group listed second so MetaboAnalystR's fold change is
  # (COV/TBI Exo) / (Cont Exo), matching "compared to the control group".
  Label  = factor(group_row[keep], levels = c("COV/TBI Exo", "Cont Exo")),
  stringsAsFactors = FALSE
)

peak_sub <- peak_mat[, keep]
out <- cbind(samp_meta, t(peak_sub))
write.csv(out, file.path(input_dir, "supp_figure_7a_peak_table.csv"), row.names = FALSE)

cat("Wrote", nrow(out), "samples x", ncol(peak_sub), "metabolites\n")
print(table(samp_meta$Label))
