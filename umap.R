
umap_fun <- function(data.matrix, config = umap.defaults, seed = 1){
  # Set seed
  set.seed(seed) 
  # Run UMAP
  umap_out <- umap(data.matrix, config = config)
  # Extract coordinates
  umap_coords.df <- umap_out$layout %>%
    as.data.frame() %>%
    rownames_to_column(var = "id")
  # Return coordinates in dataframe
  return(umap_coords.df)
}
