library(phyloseq)
# This function converts a phyloseq object into a list of matrices suitable for MOFA.
# It assumes the phyloseq object has an OTU table, a taxonomy table, and optionally a sample metadata table.
phyloseq_to_mofa <- function(physeq_obj, rank = "species") {
  if (!inherits(physeq_obj, "phyloseq")) {
    stop("Input must be a phyloseq object.")
  }

  # Agglomerate to Species level (or any desired rank) and CLR transform the OTU table.
  # This is critical for compositional data like microbiome counts.
  # You can change "Species" to "Genus" or any other rank depending on 
  # your data and research question.
  if (!is.null(tax_table(physeq_obj))) {
    physeq_obj_glom <- tax_glom(physeq_obj, taxrank = rank)
  } else {
    physeq_obj_glom <- physeq_obj  # Use original if no taxonomy for glomming
  }

  # 1. Extract the OTU table and convert to a matrix; throw an informative error if missing.
  otu_obj <- otu_table(physeq_obj_glom)
  if (is.null(otu_obj)) {
    stop("phyloseq object has no OTU table")
  }
  otu_matrix <- as(otu_obj, "matrix")
  if (taxa_are_rows(physeq_obj_glom)) {
    otu_matrix <- t(otu_matrix) # Transpose if taxa are rows
  }

  # 2. Extract the taxonomy table and convert to a matrix (if it exists)
  tax_matrix <- NULL
  tax_obj <- tryCatch(tax_table(physeq_obj_glom), error = function(e) NULL)
  if (!is.null(tax_obj)) {
    tax_matrix <- as(tax_obj, "matrix")
  }

  # 3. Extract the sample metadata and convert to a data frame (if it exists)
  sample_df <- NULL
  sample_obj <- tryCatch(sample_data(physeq_obj_glom), error = function(e) NULL)
  if (!is.null(sample_obj)) {
    sample_df <- as(sample_obj, "data.frame")
    # Keep row names in sync with sample names
    rownames(sample_df) <- sample_names(physeq_obj_glom)
  }

  # 4. Return a list of matrices/dataframes for MOFA
  return(list(
    OTU = otu_matrix,
    Taxonomy = tax_matrix,
    SampleData = sample_df
  ))
}
