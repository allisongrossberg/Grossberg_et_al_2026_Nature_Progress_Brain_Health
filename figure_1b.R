# =============================================================================
# Figure 1b – Participant group distribution (pie chart)
# =============================================================================
#
# Description: Pie chart of participant counts by study group (four groups).
#
# Prerequisites: figure_1_output_data/COAST_Study_Data_Clean_Age_Groups_add_dates.csv
#
# Output: figure_1b_Group_Pie_Chart.pdf
# =============================================================================

cairo_pdf(tempfile(fileext = ".pdf"))  # Suppress default PDF device when running in R GUI
output_dir <- if (dir.exists("figure_1")) "figure_1/figure_1_output_data" else "figure_1_output_data"
tif_dir    <- if (dir.exists("figure_1")) "figure_1" else "."
COAST_Study_Data_Clean_Age_Groups_add_dates <- read.csv(
  file.path(output_dir, "COAST_Study_Data_Clean_Age_Groups_add_dates.csv"),
  stringsAsFactors = FALSE
)

library(ggplot2)

# -----------------------------------------------------------------------------
# Build pie chart
# ----------------------------------------------------------------------------- 

group_counts <- as.data.frame(table(COAST_Study_Data_Clean_Age_Groups_add_dates$qq_group))
names(group_counts) <- c("group", "participants")

df <- group_counts
df$percentage <- df$participants / sum(df$participants) * 100
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(df[, c("group", "participants", "percentage")], file.path(output_dir, "figure_1b_Group_Pie_Chart_data.csv"), row.names = FALSE)

# Create labels with percentage and n
df$label <- paste0(round(df$percentage, 1), "%\n(n = ", df$participants, ")")

# Create the pie chart
Group_Pie_Chart <- ggplot(df, aes(x = "", y = participants, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm"), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom", 
    legend.key.size = unit(1.5, "cm"),  # Increase legend key size
    legend.text = element_text(size = 20, face = "bold"),
    text = element_text(size = 14, face = "bold")
  ) +
  scale_fill_manual(values = c(
    "COVID-19 (+) mTBI (+)" = "#DD4726",
    "COVID-19 (-) mTBI (-)" = "#7583B7",
    "COVID-19 (+) mTBI (-)" = "#739D51",
    "COVID-19 (-) mTBI (+)" = "#D06FAC"
  )) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5),
            size = 8, fontface = "bold") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) 

ggsave(file.path(tif_dir, "figure_1b_Group_Pie_Chart.pdf"), Group_Pie_Chart, width = 15, height = 10, dpi = 600, device = grDevices::cairo_pdf)
ggsave(sub("\\.pdf$", ".tif", file.path(tif_dir, "figure_1b_Group_Pie_Chart.pdf")), Group_Pie_Chart, width = 15, height = 10, dpi = 600, device = "tiff")
