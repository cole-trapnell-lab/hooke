#' adds umap coords to a data frame
#' 
add_umap_coords <- function(df, umap_centers) {
  
  from_df = left_join(df, umap_centers, by = c("from"= "cell_group")) %>% 
    rename_with(~gsub(pattern="_",replacement = "_from_", .x)) 
  
  to_df = left_join(df, umap_centers, by = c("to"= "cell_group")) %>% 
    rename_with(~gsub(pattern="_",replacement = "_to_", .x))
  
  return(full_join(from_df, to_df))
  
}

#' return the complementing edges
blacklist <- function(edges) {
  igraph::graph_from_data_frame(edges) %>% 
    igraph::complementer() %>% 
    igraph::as_data_frame()
}
