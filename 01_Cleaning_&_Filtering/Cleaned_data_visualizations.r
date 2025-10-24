# ===============================================================
# RAW DATA DESCRIPTIVE VISUALS WITH STATS
# ===============================================================


# ---------------------------------------------------------------
# SETUP
# ---------------------------------------------------------------

setwd("~/Library/CloudStorage/Box-Box/Proteomics_YvO_Master/Working")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
p_load(
  readxl, readr, dplyr, tidyr, ggplot2, ggpubr, stringr, purrr,
  RColorBrewer, scales
)

# ---------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------

cleaned_df <- read_csv("01_Cleaning_&_Filtering/01_cleaned.csv")

# Identify numeric sample columns (exclude metadata)
abundance_cols <- setdiff(names(cleaned_df), c("Gene", "N_seq", "uniprot_id", "Contaminant"))
dir.create("01_Cleaning_&_Filtering/Raw_Visualizations", showWarnings = FALSE)

# ===============================================================
# FIGURE 1 — TOTAL SIGNAL PER SAMPLE
# ===============================================================

total_signal <- cleaned_df %>%
  select(all_of(abundance_cols)) %>%
  summarise(across(everything(), ~sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "TotalSignal") %>%
  mutate(Group = case_when(
    grepl("^Old.*Pre", Sample)   ~ "Old_Pre",
    grepl("^Old.*Post", Sample)  ~ "Old_Post",
    grepl("^Young.*Pre", Sample) ~ "Young_Pre",
    grepl("^Young.*Post", Sample)~ "Young_Post",
    TRUE ~ "Other"
  ))


# Define pairwise comparisons of interest
my_comparisons <- list(
  c("Old_Post", "Old_Pre"),
  c("Old_Pre", "Young_Post"),
  c("Young_Pre", "Young_Post"),
  c("Old_Post", "Young_Post"),
  c("Old_Pre", "Young_Pre"),
  c("Old_Post", "Young_Pre")
)

# --- Boxplot with all statistical annotations
ggplot(total_signal, aes(x = Group, y = TotalSignal, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_compare_means(
    method = "anova",
    label.y = max(total_signal$TotalSignal) * 1.5,
    label.x = 0.65
  ) +  # Global ANOVA
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    size = 4
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Total Signal per Group (Raw Intensities)",
       x = "Experimental Group",
       y = "Total Signal (Sum of Intensities)") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
ggsave("01_Cleaning_&_Filtering/Raw_Visualizations/04_TotalSignal_per_Group.png",
       width = 7, height = 5, dpi = 300)

# ===============================================================
# FIGURE 2 — PROTEINS DETECTED PER SAMPLE
# ===============================================================

protein_detection <- cleaned_df %>%
  select(all_of(abundance_cols)) %>%
  summarise(across(everything(), ~sum(!is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "DetectedProteins") %>%
  mutate(Group = case_when(
    grepl("^Old.*Pre", Sample)   ~ "Old_Pre",
    grepl("^Old.*Post", Sample)  ~ "Old_Post",
    grepl("^Young.*Pre", Sample) ~ "Young_Pre",
    grepl("^Young.*Post", Sample)~ "Young_Post",
    TRUE ~ "Other"
  ))

# --- Define identical pairwise comparisons to Figure 1
my_comparisons <- list(
  c("Old_Post", "Old_Pre"),
  c("Old_Pre", "Young_Post"),
  c("Young_Pre", "Young_Post"),
  c("Old_Post", "Young_Post"),
  c("Old_Pre", "Young_Pre"),
  c("Old_Post", "Young_Pre")
)

# --- Boxplot with all statistical annotations
ggplot(protein_detection, aes(x = Group, y = DetectedProteins, fill = Group)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.6, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_compare_means(
    method = "anova",
    label.y = max(protein_detection$DetectedProteins) * 1.5,
    label.x = 0.65
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "t.test",
    label = "p.signif",
    size = 4
  ) +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Proteins Detected per Group",
       x = "Experimental Group",
       y = "Number of Detected Proteins") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

ggsave("01_Cleaning_&_Filtering/Raw_Visualizations/03_ProteinsDetected_per_Group.png",
       width = 7, height = 5, dpi = 300)


# ===============================================================
# FIGURE 3 — PROTEIN DETECTION SUMMARY BY GROUP
# ===============================================================

groups <- list(
  Old_Pre   = grep("^Old.*Pre", names(cleaned_df), value = TRUE),
  Old_Post  = grep("^Old.*Post", names(cleaned_df), value = TRUE),
  Young_Pre = grep("^Young.*Pre", names(cleaned_df), value = TRUE),
  Young_Post= grep("^Young.*Post", names(cleaned_df), value = TRUE)
)

detected_by_group <- sapply(groups, function(cols) {
  rowSums(!is.na(cleaned_df[, cols])) > 0
})

detection_summary <- data.frame(
  Group = names(groups),
  Detected = colSums(detected_by_group),
  Total = nrow(cleaned_df)
) %>%
  mutate(NotDetected = Total - Detected) %>%
  pivot_longer(cols = c("Detected", "NotDetected"),
               names_to = "Status", values_to = "Count")

ggplot(detection_summary, aes(x = Group, y = Count, fill = Status)) +
  geom_col(position = "stack", color = "black") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Protein Detection Summary by Group",
       x = "Experimental Group",
       y = "Protein Count",
       fill = "Status") +
  theme_minimal(base_size = 14)

ggsave("01_Cleaning_&_Filtering/Raw_Visualizations/02_ProteinDetection_by_Group.png",
       width = 7, height = 5, dpi = 300)


# ===============================================================
# FIGURE 4 — TOTAL SIGNAL PER SAMPLE (ORDERED)
# ===============================================================

ggplot(total_signal, aes(x = reorder(Sample, TotalSignal), y = TotalSignal, fill = Group)) +
  geom_col(color = "black") +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Total Signal per Sample (Ordered by Intensity)",
       x = "Sample", y = "Total Signal") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank())

ggsave("01_Cleaning_&_Filtering/Raw_Visualizations/01_TotalSignal_per_Sample.png",
       width = 7, height = 5, dpi = 300)


# ===============================================================
# PRINT SUMMARY
# ===============================================================

summary_table <- data.frame(
  Metric = c(
    "Total Proteins",
    "Total Samples",
    "Mean Proteins per Sample",
    "Mean Total Signal per Sample"
  ),
  Value = c(
    nrow(cleaned_df),
    length(abundance_cols),
    mean(protein_detection$DetectedProteins),
    mean(total_signal$TotalSignal)
  )
)

print(summary_table)

