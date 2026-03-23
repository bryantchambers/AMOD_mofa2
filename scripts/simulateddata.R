# ==============================================================================
# Generating Synthetic Multi-Omics Data for MOFA2
# ==============================================================================

# Load transformation library
# install.packages("compositions")
library(compositions)

set.seed(42) # For reproducibility

# 1. Setup Experimental Design
# ------------------------------------------------------------
n_samples <- 50
sample_names <- paste0("Sample_", 1:n_samples)

# Create a "Hidden Signal" (e.g., a Sediment Depth gradient)
# This is the "Latent Factor" we hope MOFA will rediscover.
depth_signal <- seq(-2, 2, length.out = n_samples) 

# 2. View 1: Taxonomy (mat_taxa)
# ------------------------------------------------------------
# We'll simulate 100 Genera. 
# Some will increase with depth, some will decrease, others are noise.
n_taxa <- 100
mat_taxa_raw <- matrix(rpois(n_samples * n_taxa, lambda = 10), 
                       nrow = n_taxa, ncol = n_samples)

# Inject signal into the first 20 taxa
for(i in 1:20) {
  mat_taxa_raw[i, ] <- rpois(n_samples, lambda = exp(2 + (depth_signal * runif(1, -1, 1))))
}

rownames(mat_taxa_raw) <- paste0("Genus_", 1:n_taxa)
colnames(mat_taxa_raw) <- sample_names

# CRITICAL STEP: CLR Transformation
# We add a small pseudocount (1) to handle zeros, then CLR.
mat_taxa <- t(clr(t(mat_taxa_raw + 1)))

# 3. View 2: Metabolism (mat_kegg)
# ------------------------------------------------------------
# 150 KEGG Orthologs (KOs). Let's assume these are continuous (e.g., TPM).
n_kegg <- 150
mat_kegg <- matrix(rnorm(n_samples * n_kegg, mean = 5, sd = 1), 
                   nrow = n_kegg, ncol = n_samples)

# Inject signal: KOs linked to the taxa above
for(i in 1:30) {
  mat_kegg[i, ] <- 5 + (depth_signal * rnorm(1, 2, 0.5)) + rnorm(n_samples, sd = 0.5)
}

rownames(mat_kegg) <- paste0("KO_", sprintf("%05d", 1:n_kegg))
colnames(mat_kegg) <- sample_names

# 4. View 3: Biosynthesis (mat_bgc)
# ------------------------------------------------------------
# 50 BGCs. We'll use binary presence/absence (Bernoulli distribution).
n_bgc <- 50
mat_bgc <- matrix(rbinom(n_samples * n_bgc, 1, prob = 0.2), 
                  nrow = n_bgc, ncol = n_samples)

# Inject signal: Certain BGCs appear only at specific depths
for(i in 1:10) {
  prob_signal <- 1 / (1 + exp(-(depth_signal * rnorm(1, 2, 0.5)))) # Logistic mapping
  mat_bgc[i, ] <- rbinom(n_samples, 1, prob = prob_signal)
}

rownames(mat_bgc) <- paste0("BGC_Cluster_", 1:n_bgc)
colnames(mat_bgc) <- sample_names

# 5. Combine into MOFA Data List
# ------------------------------------------------------------
data_list <- list(
  Taxonomy = mat_taxa,
  Metabolism = mat_kegg,
  Biosynthesis = mat_bgc
)

# Preview the data
str(data_list)