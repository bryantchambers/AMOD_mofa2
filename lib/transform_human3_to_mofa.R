library(tidyverse)

# 1. Load your stratified HUMAnN table (e.g., regrouped to KO)
# Replace 'path_to_file.tsv' with your actual filename
df <- read_tsv("path_to_file.tsv") %>%
  rename(Feature = 1)

# 2. Filter for stratified data and split the 'Feature' column
# We exclude the 'unclassified' or totals by looking for the pipe '|'
ko_species_matrix <- df %>%
  filter(str_detect(Feature, "\\|")) %>%
  separate(Feature, into = c("KO", "Species"), sep = "\\|")

# 3. Clean up names (Optional: removing 's__' prefix from MetaPhlAn names)
ko_species_matrix <- ko_species_matrix %>%
  mutate(Species = str_replace(Species, "s__", ""))

# 4. Reshape if you want a specific Sample/Species view for one KO
# Or keep as is for a "Long" format which is better for ggplot2
head(ko_species_matrix)

# 5. Example: Total abundance of a specific KO per species across samples
ko_summary <- ko_species_matrix %>%
  group_by(KO, Species) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop')

# 6. Convert to a Wide Matrix (KOs as rows, Species as columns) for one sample
# Replace 'Sample_A' with your actual sample column name
sample_matrix <- ko_summary %>%
  select(KO, Species, Sample_A) %>%
  pivot_wider(names_from = Species, values_from = Sample_A, values_fill = 0)
