#This function performs a Centered Log-Ratio (CLR) transformation on a matrix of compositional data, which is common in microbiome analyses. The CLR transformation helps to address the compositional nature of the data, making it more suitable for downstream analyses like MOFA.
mofa_prep_CLRtransform <- function(mat_taxa_raw) {
 
  # Add a small pseudocount to handle zeros (common in microbiome data)
     mat_taxa <- t(compositions::clr(t(mat_taxa_raw + 1)))
    return(mat_taxa)

}
