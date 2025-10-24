# ===============================================================
# Setup
# ===============================================================

setwd("~/Library/CloudStorage/Box-Box/Proteomics_YvO_Master/Working")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)
p_load(readxl, dplyr, stringr)

# ===============================================================
# Load and clean raw data
# ===============================================================

# Read raw file (no column names)
raw <- read_excel("00_Inputs/00_raw.xlsx", col_names = FALSE)

# 1. Extract total spectra and column names
total_spectra <- raw[1,] %>% unlist() %>% as.numeric()
col_names <- raw[2,] %>% unlist() %>% as.character() %>%
  str_replace(".*OP-", "Old_OP_") %>%
  str_replace(".*O-", "Old_O_") %>%
  str_replace(".*YP-", "Young_YP_") %>%
  str_replace(".*Y-", "Young_Y_")

# 2. Clean dataframe and add names
cleaned_df <- raw[-c(1,2),]
colnames(cleaned_df) <- col_names

# 3. Keep only relevant columns
cleaned_df <- cleaned_df %>%
  select(where(~!all(is.na(.)))) %>%
  select(-First.Protein.Description) %>%
  rename(Gene = Genes, N_seq = N.Sequences)

# 4. Convert numeric columns
abundance_cols <- setdiff(names(cleaned_df), c("Gene", "N_seq"))
cleaned_df <- cleaned_df %>% mutate(across(all_of(abundance_cols), as.numeric))

# 5. Remove duplicates (keep highest N_seq)
cleaned_df <- cleaned_df %>%
  arrange(Gene, desc(N_seq)) %>%
  distinct(Gene, .keep_all = TRUE)

# ===============================================================
# Add Contaminant Column
# ===============================================================

contaminant_markers <- c(
  "HBA1","HBA2","HBB","HBG1","HBG2","HBM","HBD","HBZ",
  "ALB","TF","HP","HPX","PLG","APOA1","APOA2","APOB","APOC1","APOC2","APOC3","APOE",
  "FGA","FGB","FGG","AHSG","APOH","C3","C4A","C4B","C9","CRP","ORM1","ORM2",
  "SAA1","SAA2","SERPINA1","SERPINA3","SERPINC1","SERPIND1","SERPINF2","VTN","CLU",
  "SHBG","TTR","GC","IGHG1","IGHG2","IGHG3","IGHG4","IGHA1","IGHA2","IGHM","IGKC",
  "IGLC1","IGLC2","IGLC3","IGLC7","THBS1","ITGA2B","ITGB3","GP1BA","GP1BB","GP9",
  "PF4","PPBP","VWF","F13A1","F13B","KRT1","KRT2","KRT5","KRT9","KRT10","KRT14",
  "KRT16","KRT17","KRT18","KRT19","KRT20","FLG","IVL","TRY1","TRY2","PRSS1","PRSS2",
  "PRSS3","CSN1S1","CSN2","CSN3","LALBA","BSA","ALBU_BOVIN","CASA1","CASA2","CASA3"
)

cleaned_df <- cleaned_df %>%
  mutate(Contaminant = if_else(Gene %in% contaminant_markers, "+", "-"))

# ===============================================================
# Add UniProt ID Mapping
# ===============================================================

# 1. Load UniProt mapping file (must have 'Gene' and 'uniprot_id' columns)
id_map_df <- read_excel("00_Inputs/00_idmapping.xlsx", sheet = "Reviewed") %>%
  rename_with(~ str_replace_all(., " ", "_")) %>%
  distinct(Gene, .keep_all = TRUE)

# 2. Join to main dataframe
cleaned_df <- cleaned_df %>%
  left_join(id_map_df %>% select(Gene, uniprot_id), by = "Gene") %>%
  relocate(uniprot_id, .after = Gene)

# 3. Optional: Count unmapped genes
unmapped_count <- cleaned_df %>%
  filter(is.na(uniprot_id)) %>%
  nrow()

cat("Unmapped genes:", unmapped_count)
# ===============================================================
# Save Cleaned Data
# ===============================================================

if (!dir.exists("01_Cleaning_&_Filtering")) dir.create("01_Cleaning_&_Filtering")
write.csv(cleaned_df, "01_Cleaning_&_Filtering/01_cleaned.csv", row.names = FALSE)

