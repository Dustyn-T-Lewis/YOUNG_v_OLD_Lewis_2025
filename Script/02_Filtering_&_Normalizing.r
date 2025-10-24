# ===============================================================
# INSTALL
# ===============================================================
setwd("~/Library/CloudStorage/Box-Box/Proteomics_YvO_Master/Working/Proteomics_YvO_Master")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

install.packages("devtools")
devtools::install_github("ByrumLab/proteoDA", 
                         dependencies = TRUE, 
                         build_vignettes = TRUE)
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
p_load(
  proteoDA, dplyr, stringr, ggplot2
)

# ===============================================================
# IMPORT & DATA SETUP
# ===============================================================
input_data <- read.csv("01_Cleaning_&_Filtering/01_cleaned.csv")

# Create a sample metadata dataframe from row names
sample_metadata <- data.frame(
  Sample = colnames(input_data[, 4:67])   # Sample names in columns 4–67
) |>
  mutate(
    AgeGroup  = str_extract(Sample, "^(Old|Young)"),
    Sample_ID = str_extract(Sample, "(?<=Old_|Young_)[A-Za-z0-9_]+(?=_Pre|_Post)"),
    Timepoint = str_extract(Sample, "(Pre|Post)"),
    Group = paste(AgeGroup, Timepoint, sep = "_")
  )
rownames(sample_metadata) = sample_metadata$Sample

annotation_data <- input_data |>
  select(Gene, uniprot_id, N_seq, Contaminant)

intensity_data <- input_data[,4:67]

# Create DAList
raw <- DAList(
  data       = intensity_data,
  annotation = annotation_data,
  metadata   = sample_metadata,
  design = NULL,
  eBayes_fit = NULL,
  results = NULL,
  tags = NULL)

# ===============================================================
# SAMPLE & PROTEIN FILTERING
# ===============================================================

                                  # Manual Group-Wise filtering function 
                                  filter_proteins_by_proportion_any <- function(DAList, min_prop = 0.66, grouping_column = "Group") {
                                      # Extract DAList data frames of interest
                                      dat <- DAList$data
                                      meta <- DAList$metadata
                                      groups <- meta[[grouping_column]]
                                      uniq_groups <- unique(groups)
                                      # Determine which proteins meet the threshold in at least (≥ 1 group)    
                                    keep <- apply(dat, 1, function(x) {
                                      # Copmute a '% present per group' = present/total (a fraction detection score)
                                      present_by_group <- sapply(uniq_groups, function(g) {
                                        samples <- which(groups == g)
                                        prop_present <- mean(!is.na(x[samples]))
                                        prop_present >= min_prop
                                      })
                                      any(present_by_group)  # *'any()' at least one group passes
                                    })
                                       # Subset DAList based on retained proteins   
                                    DAList$data <- dat[keep, ]
                                    DAList$annotation <- DAList$annotation[keep, ]
                                       # Filter Summary
                                    message("Filtered ", sum(!keep), "; kept ", sum(keep), " (2/3 Group-Wise Filter: ≥ ", min_prop*100, "% in ≥1 group).")
                                    return(DAList)
                                  }
# Filtering Function      
filtered <- raw |>
  filter_proteins_by_annotation(Contaminant != "+") |>
  zero_to_missing() |>
  filter_proteins_by_proportion_any(min_prop = 0.66, grouping_column = "Group")

# "Filtering by sample # = 2/3 largest group size (17)"      
#      filtered <- raw |>
#        filter_proteins_by_annotation(Contaminant != "+") |>
#        zero_to_missing() |>
#        filter_proteins_by_group(
#          min_reps = ceiling(0.66 * max(table(raw$metadata$Group))), # ≈ two-thirds rule
#          min_groups = 1,                                            # at least one group
#          grouping_column = "Group"
#        )
#      max(table(raw$metadata$Group)) * .66

# ===============================================================
# NORMALIZATION
# ===============================================================

write_norm_report(filtered,
                  grouping_column = "Group",
                  output_dir = "02_QC_report",
                  filename = "normalization.pdf",
                  overwrite = T,
                  suppress_zoom_legend = FALSE,
                  use_ggrastr = FALSE)

# methods = log2, median, mean, vsn, quantile, cycloess, rlr, gi
normalized <- normalize_data(filtered, 
                             norm_method = "cycloess")

# ===============================================================
# QC REPORT
# ===============================================================

write_qc_report(normalized,
                color_column = "Group",
                label_column = NULL,
                output_dir = "02_QC_report",
                filename = "QC_report.pdf",
                overwrite = T,
                top_proteins = 500,        # number of most variable proteins
                standardize = TRUE,
                pca_axes = c(1,2),          # first 2 PCs
                dist_metric = "euclidean",  # stats::dist for options
                clust_method = "complete",  # stats::hclust for options
                show_all_proteins = F)      # only proteins with missing data

df <- data.frame(
  Sample = rep(colnames(normalized$data), each = nrow(normalized$data)),
  Intensity = log2(as.vector(as.matrix(normalized$data))),
  Group = rep((normalized$metadata)$Group, each = nrow(normalized$data))
)

p <- ggplot(df, aes(Sample, Intensity, fill = Group)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "log2(Intensity)", title = "Sample Intensity Distributions") +
  coord_cartesian(ylim = c(3, 6))   

ggsave("02_QC_report/Manual_Violin_Plot.pdf", p, width = 10, height = 6)

# ===============================================================
# STATISTICAL MODEL 
# ===============================================================