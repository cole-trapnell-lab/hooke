#' adds umap coords to a data frame
#' 
add_umap_coords <- function(df, umap_centers) {
  
  from_df = left_join(tmp_df, umap_centers, by = c("from"= "cell_group")) %>% 
    rename_with(~gsub(pattern="_",replacement = "_from_", .x)) 
  
  to_df = left_join(tmp_df, umap_centers, by = c("to"= "cell_group")) %>% 
    rename_with(~gsub(pattern="_",replacement = "_to_", .x))
  
  return(full_join(from_df, to_df, by = c("from", "to")))
  
}
