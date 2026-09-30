# =============================================================================
# Supplementary Figure 7a preprocessing 2 - Glycolysis fold-change/t-test reanalysis
#   (2/3-phospho-D-glycerate and phosphoenolpyruvate increased, COV/TBI Exo
#   vs Cont Exo, 24h HA metabolomics) using MetaboAnalystR, the R package
#   underlying the original MetaboAnalyst web-tool analysis.
# =============================================================================
#
# Description: Methods: median normalization, log10 transform, Pareto
#   scaling; fold change + t-test vs CONT; FC > 1.0 and raw p < 0.10.
#   Reproduces the original analyst's reported results almost exactly (9
#   significant metabolites, matching log2FC/p-values within rounding).
#   Student's (equal-variance) t-test is primary; Welch's (unequal-variance)
#   is kept as a reported robustness check.
#
# Input: supp_figure_7a_input_data/supp_figure_7a_peak_table.csv (built by
#   supp_figure_7a_preprocessing_1_prepare_input.R; 6 samples [3 Cont Exo,
#   3 COV/TBI Exo] x 142 metabolites, raw LC-MS peak intensities)
# Output: supp_figure_7a_output_data/supp_figure_7a_combined_results.csv (Student's,
#   primary), supp_figure_7a_output_data/supp_figure_7a_ttest_variant_comparison.csv
#   (Student vs Welch, target metabolites only),
#   supp_figure_7a_output_data/supp_figure_7a_published_number_check.csv (validation
#   against the analyst's reported numbers)
# =============================================================================

if (basename(getwd()) == "supp_figure_7") setwd("..")
input_dir  <- "supp_figure_7/supp_figure_7_input_data"
output_dir <- "supp_figure_7/supp_figure_7_output_data"
work_dir   <- file.path(output_dir, "metaboanalyst_workdir")
dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)

suppressMessages(library(MetaboAnalystR))

peak_table_path <- normalizePath(file.path(input_dir, "supp_figure_7a_peak_table.csv"))
target_metabolites <- c("2/3-Phospho-D-glycerate", "Phosphoenolpyruvate")

# Runs the full pipeline once for a given t-test variant and returns a
# combined fold-change + t-test result table.
run_pipeline <- function(equal.var) {
  old_wd <- getwd()
  setwd(work_dir)  # MetaboAnalystR reads/writes intermediate .qs files in cwd
  on.exit(setwd(old_wd))

  mSet <- InitDataObjects("pktable", "stat", FALSE, default.dpi = 72)
  mSet <- Read.TextData(mSet, peak_table_path, "rowu", "disc")
  mSet <- SanityCheckData(mSet)

  # Force the comparison direction explicitly: numerator = COV/TBI Exo,
  # denominator = Cont Exo (the control), matching "compared to the CONT
  # group" in the Methods -- independent of factor level ordering assigned
  # internally from the CSV.
  mSet$dataSet$comp.groups <- c("COV/TBI Exo", "Cont Exo")

  mSet <- Normalization(mSet, rowNorm = "MedianNorm", transNorm = "LogNorm",
                         scaleNorm = "ParetoNorm", ref = NULL, ratio = FALSE, ratioNum = 20)

  # Work around a state-passing quirk in this MetaboAnalystR build: when run
  # outside the web-server session (.on.public.web == FALSE), Normalization()
  # reads dataSet$prenorm.cls off a stale local copy of mSetObj and ends up
  # setting dataSet$cls to NULL, even though the normalized data itself is
  # correct. dataSet$orig.cls (set by SanityCheckData, untouched by
  # Normalization) still holds the correct per-sample group labels in the
  # same row order as dataSet$norm, so restore cls from it.
  stopifnot(is.null(mSet$dataSet$cls) || length(mSet$dataSet$cls) == nrow(mSet$dataSet$norm))
  mSet$dataSet$cls <- mSet$dataSet$orig.cls
  mSet$dataSet$comp.groups <- c("COV/TBI Exo", "Cont Exo")

  # fc.thresh = 1.0: the manuscript's own stated criterion (FC > 1.0), i.e.
  # no additional effect-size dead zone beyond the direction of change.
  mSet <- FC.Anal(mSet, fc.thresh = 1.0, cmp.type = 0, paired = FALSE, fc.method = "classical")

  mSet <- suppressMessages(suppressWarnings(Ttests.Anal(
    mSet, nonpar = FALSE, threshp = 0.10, paired = FALSE, equal.var = equal.var,
    pvalType = "raw", all_results = TRUE, tt.method = "classical"
  )))

  fc.all <- mSet$analSet$fc$fc.all
  fc.log <- mSet$analSet$fc$fc.log
  tt <- mSet$analSet$tt

  data.frame(
    metabolite  = names(fc.all),
    fold_change = as.numeric(fc.all),
    log2_fc     = as.numeric(fc.log),
    p_value     = as.numeric(tt$p.value[names(fc.all)]),
    stringsAsFactors = FALSE
  )
}

# Primary: Student's (equal-variance) t-test -- confirmed against the
# original analyst's published numbers (see header note).
student <- run_pipeline(equal.var = TRUE)
student$significant <- student$fold_change > 1.0 & student$p_value < 0.10
student <- student[order(student$p_value), ]

# Robustness check: Welch's (unequal variance) t-test
welch <- run_pipeline(equal.var = FALSE)
welch$significant <- welch$fold_change > 1.0 & welch$p_value < 0.10

write.csv(student, file.path(output_dir, "supp_figure_7a_combined_results.csv"), row.names = FALSE)

variant_compare <- merge(
  student[student$metabolite %in% target_metabolites, c("metabolite", "fold_change", "p_value", "significant")],
  welch[welch$metabolite %in% target_metabolites, c("metabolite", "p_value", "significant")],
  by = "metabolite", suffixes = c("_student", "_welch")
)
write.csv(variant_compare, file.path(output_dir, "supp_figure_7a_ttest_variant_comparison.csv"), row.names = FALSE)

# ---- Validation against the analyst's published numbers -------------------
published <- data.frame(
  metabolite = c("Phosphoenolpyruvate", "3-Methyleneoxindole", "acyl-C4-OH",
                 "2/3-Phospho-D-glycerate", "N-formyl kynurenine",
                 "butanoyl-l-carnitine (acyl-C4)", "GTP", "acyl-C18:2 (Linoleoyl-CoA)",
                 "Glutathione disulfide"),
  published_log2fc = c(1.09, 0.269, 0.889, 0.555, -0.831, 0.447, 0.247, 0.551, -0.966),
  published_p      = c(0.00802, 0.00839, 0.0392, 0.0610, 0.0606, 0.0699, 0.0999, 0.0954, 1.34e-4),
  stringsAsFactors = FALSE
)
check <- merge(published, student[, c("metabolite", "log2_fc", "p_value")], by = "metabolite")
check$log2fc_diff <- check$log2_fc - check$published_log2fc
check$p_diff <- check$p_value - check$published_p
check <- check[order(check$published_p), ]
write.csv(check, file.path(output_dir, "supp_figure_7a_published_number_check.csv"), row.names = FALSE)

cat("\n=== Validation: this script's Student's-t results vs the analyst's published numbers ===\n")
print(check)

cat("\n=== Target metabolites (Student's t-test, primary) ===\n")
print(student[student$metabolite %in% target_metabolites, ])

cat("\n=== Student vs Welch comparison (target metabolites) ===\n")
print(variant_compare)

cat("\n=== All significant under Student's t-test (FC > 1.0 & raw p < 0.10):", sum(student$significant), "of", nrow(student), "===\n")
print(student[student$significant, c("metabolite", "fold_change", "log2_fc", "p_value")])

cat("\n=== Rank of target metabolites by p-value (Student's), out of", nrow(student), "===\n")
print(which(student$metabolite %in% target_metabolites))

unlink(work_dir, recursive = TRUE)
