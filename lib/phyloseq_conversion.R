library(phyloseq)
# This function converts a phyloseq object into a list of matrices suitable for MOFA.
# It assumes the phyloseq object has an OTU table, a taxonomy table, and optionally a sample metadata table.
phyloseq_to_mofa <- function(physeq) {
  # 1. Extract the OTU table and convert to a matrix
  otu_table <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) {
    otu_table <- t(otu_table) # Transpose if taxa are rows
  }

    # 2. Extract the taxonomy table and convert to a matrix (if it exists)
  if (!is.null(tax_table(physeq))) {
    tax_table <- as(tax_table(physeq), "matrix")
  } else {
    tax_table <- NULL
  }                 
    
    # 3. Extract the sample metadata and convert to a matrix (if it exists)             


    if (!is.null(sample_data(physeq))) {
    sample_data <- as(sample_data(physeq), "data.frame")
    } else {
    sample_data <- NULL
    }

    # 4. Return a list of matrices for MOFA
    return(list(
      OTU = otu_table,
      Taxonomy = tax_table,
      SampleData = sample_data
    ))
}
