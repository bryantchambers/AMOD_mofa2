library(MOFA2)
library(dplyr)
library(tidyr)
library(igraph)
library(reticulate)

#reticulate::use_condaenv("amod_mofapy2")
reticulate::use_python("/maps/projects/caeg/people/gfx654/miniforge3/envs/amod_mofapy2/bin/python", required = TRUE)
# 2. Data Preparation (Conceptual)
# ------------------------------------------------------------------------------
# MOFA expects a list of matrices. 
# Rows = Features (Taxa, KOs, BGCs), Columns = Samples.
# IMPORTANT: Data must be normalized before this step!
# - Taxa: Agglomerated to a single rank (e.g., Genus) and CLR transformed.
# - KEGG: Log2(TPM + 1) transformed.
# - BGCs: Presence/Absence (1/0) or Log-normalized.

# Assuming you have three pre-processed matrices: mat_taxa, mat_kegg, mat_bgc
# We bundle them into a list. The names of the list become the "Views".
# data_list <- list(
#   Taxonomy = mat_taxa,
#   Metabolism = mat_kegg,
#   Biosynthesis = mat_bgc
# )

data_list_analysis <- list(
  Taxonomy = mat_taxa %>% unname() %>% as.matrix(),
  Metabolism = mat_kegg %>% unname() %>% as.matrix(),
  Biosynthesis = mat_bgc %>% unname() %>% as.matrix()
)

# 3. Initialize and Train the MOFA Model
# ------------------------------------------------------------------------------
# Create the MOFA object
MOFAobject <- create_mofa(data_list_analysis)

# Set Data Options: Tell MOFA what kind of statistical distributions to expect.
# We use "gaussian" for continuous/normalized data.
data_opts <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE # Scales views so one doesn't dominate the model

# Set Model Options: How many factors (hidden drivers) should it look for?
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8 # A good starting point. MOFA will drop inactive ones.

# Set Training Options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium" 

# Prepare and Train (This connects to Python's mofapy2 under the hood)
MOFAobject <- prepare_mofa(MOFAobject, data_opts, model_opts, train_opts)
MOFAobject <- run_mofa(MOFAobject, use_basilisk = FALSE)

# 4. Interpret and Select a Factor
# ------------------------------------------------------------------------------
# In a real workflow, you would correlate the factors with your metadata here.
# e.g., Correlate Factor scores with "Sediment Depth" or "Paleoecosystem State".
# Let's assume we found that Factor 3 strongly correlates with a key environmental shift.
target_factor <- "Factor3"

# 5. Extract Weights to Build the Network (Node & Edge Lists)
# ------------------------------------------------------------------------------
# We extract the weights (loadings) for all features across all views for Factor 3.
weights <- get_weights(MOFAobject, factors = target_factor, as.data.frame = TRUE)

# Filter for the most important features (the "signal").
# We take the absolute value of the weight (magnitude matters, sign dictates direction).
# Here, we keep features in the top 5% of weights for this specific factor.
threshold <- quantile(abs(weights$value), 0.95)
top_features <- weights %>% filter(abs(value) >= threshold)

# Create Node List
# We need to know what view each feature came from so Cytoscape can color them differently.
nodes <- top_features %>%
  select(id = feature, type = view, weight = value) %>%
  distinct()

# Create Edge List
# How do we connect them? In MOFA, if two features are highly weighted on the SAME factor,
# they are co-varying. We create edges between all top features in this factor.
# Note: This creates a fully connected sub-graph (clique) for this factor. 
# For multi-factor networks, you would calculate pairwise correlations of their weight profiles.
edges <- expand.grid(source = nodes$id, target = nodes$id) %>%
  filter(source != target) %>% # Remove self-loops
  # Optional: Remove duplicate undirected edges (A-B and B-A)
  mutate(combo = paste(pmin(source, target), pmax(source, target), sep="_")) %>%
  distinct(combo, .keep_all = TRUE) %>%
  select(source, target)

# Add an edge attribute (e.g., they belong to the Factor 3 interaction module)
edges$module <- target_factor

# 6. Build and Export for Cytoscape
# ------------------------------------------------------------------------------
# Create the igraph object
network_graph <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# Export to GraphML (Cytoscape's preferred format)
# In Cytoscape, simply go to File -> Import -> Network from File -> Select this GraphML.
write_graph(network_graph, "MOFA_Factor3_Network.graphml", format = "graphml")

# Alternatively, export clean CSVs for custom knowledge graph ingestion
write.csv(nodes, "MOFA_Node_List.csv", row.names = FALSE)
write.csv(edges, "MOFA_Edge_List.csv", row.names = FALSE)

print("Network construction complete. Ready for Cytoscape or Knowledge Graph integration.")