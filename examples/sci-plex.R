library(hooke)
library(monocle3)
library(devtools)
library(tidyr)
library(tibble)
library(dplyr)

devtools::load_all("~/OneDrive/UW/Trapnell/hooke/")
setwd("~/OneDrive/UW/Trapnell/hooke/examples/")


a549_cds <- readRDS("R_objects/A549_24hrs.RDS")
a549_coldata = a549_cds %>% colData %>% as.data.frame()
# recluster in umap space 

# res default 
a549_cds <- cluster_cells(a549_cds, reduction_method = "UMAP", resolution = 5e-5)

#dir.create("sci_plex")
# plot_cells(a549_cds) %>% ggsave(filename = "examples/sci_plex/sci_plex_cluster.png")



#k562_cds <- readRDS("R_objects/K562_24hrs.RDS")
# mcf7_cds <- readRDS("R_objects/MCF7_24hrs.RDS")

unique(a549_coldata$pathway_level_2) %>% sort()

a549_coldata %>% filter(pathway_level_2 == "Aurora kinase activity") %>% 
  pull(product_name) %>% unique()

hdac = a549_coldata %>% filter(pathway_level_2 == "Histone deacetylation")%>% 
  pull(catalog_number) %>% unique() 
# a549 cells ------------------------------------------------------------------

product_names = unique(colData(a549_cds)$product_name)
catalog_numbers = unique(colData(a549_cds)$catalog_number)

cat_to_prod = colData(a549_cds) %>% 
  as.data.frame() %>% 
  group_by(product_name, catalog_number) %>% 
  tally() %>% 
  select(-n)

drugs = c("Trametinib (GSK1120212)", "Epothilone A", "Epothilone B")

colData(a549_cds)$sample = NULL
colData(a549_cds)$log_dose = log10(colData(a549_cds)$dose+1)

# replace the -inf with 0
colData(a549_cds)$dose %>% unique()
colData(a549_cds)$log_dose %>% unique()
colData(a549_cds)$grouping = paste(colData(a549_cds)$catalog_number, 
                                   colData(a549_cds)$dose,
                                   colData(a549_cds)$replicate, 
                                   colData(a549_cds)$top_oligo_W, sep=".")

colData(a549_cds)$Cluster = clusters(a549_cds)
colData(a549_cds)$cluster = paste0("C", colData(a549_cds)$Cluster)

# colData(a549_cds) %>% 
#   as.data.frame() %>% 
#   filter(product_name == "Vehicle")
# # which column is embryo like
# colData(a549_cds) %>% 
#   as.data.frame() %>% 
#   group_by(product_name, dose,replicate, drug_dose, top_oligo_W) %>% 
#   tally() %>% 
#   filter(product_name =='Vehicle')
# colData(a549_cds) %>% 
#   as.data.frame() %>% 
  # filter(product_name == "Dacinostat (LAQ824)") 
  # filter(pathway_level_2 == "Histone deacetylation") %>% 
  # pull(product_name) %>% unique()
  # filter(grepl("cycle",pathway_level_1)) 
  
  
# drug dose is product_name * dose
a549_ccs = new_cell_count_set(a549_cds,
                         sample_group = "grouping",
                         cell_group = "cluster")


a549_ccs_coldata = colData(a549_ccs) %>% as.data.frame()
colData(a549_ccs)[["S0000"]] = 0

for (cn in catalog_numbers) {
  if (cn != "S0000") {
    a549_ccs_coldata = a549_ccs_coldata %>%
      mutate("{cn}":= ifelse(catalog_number==cn, log_dose, 0))
    colData(a549_ccs)[[cn]] = a549_ccs_coldata[[cn]]
  }
}

# make formula string
sub_catalog_numbers = catalog_numbers[catalog_numbers!="S0000"]
product_formula = paste(sub_catalog_numbers, collapse=" + ")
product_formula = paste0("~ ", product_formula)



# a549_ccm = new_cell_count_model(a549_ccs, model_formula_str = "~splines::ns(as.numeric(S1055), df=3) ", pseudocount=1)

# saveRDS(a549_ccm, "sci_plex/a549_ccm_spline.rds")
# saveRDS(a549_ccm, "sci_plex/a549_ccm.rds") # no spline

a549_ccm <- readRDS("sci_plex/a549_ccm.rds")


# recode the drugs
coef_results = coef(a549_ccm@best_model) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("catalog_number") %>% 
  left_join(cat_to_prod, by = "catalog_number") %>% 
  select(catalog_number, product_name,1:20)


# returns the catalog number for a given drug
return_catalog <- function(pn, cat_to_prod) {
  cat_to_prod %>% filter(product_name == pn) %>% pull(catalog_number)
}

# returns the product name for a given catalog number
return_drug_name <- function(cn, cat_to_prod) {
  cat_to_prod %>% filter(catalog_number == cn) %>% pull(product_name)
}


# make tibble with everything 0 except this col and value
make_tbl <- function(colname, value) {
  tbl = tibble(x = catalog_numbers, y=0)
  tbl_wide = spread(tbl, x, y)
  tbl_wide[[colname]] = value
  return(tbl_wide)
}

plot_series <- function(a549_ccm, drug_name, pvalue=1.0) {

  dose_0 = estimate_abundances(a549_ccm, make_tbl(drug_name, 0))
  dose_10 = estimate_abundances(a549_ccm, make_tbl(drug_name, log10(11)))
  dose_100 = estimate_abundances(a549_ccm, make_tbl(drug_name, log10(101)))
  dose_1000 = estimate_abundances(a549_ccm, make_tbl(drug_name, log10(1001)))
  dose_10000 = estimate_abundances(a549_ccm, make_tbl(drug_name, log10(10001)))

  cond_0_vs_10_tbl = compare_abundances(a549_ccm, dose_0, dose_10)
  cond_0_vs_100_tbl = compare_abundances(a549_ccm, dose_0, dose_100)
  cond_0_vs_1000_tbl = compare_abundances(a549_ccm, dose_0, dose_1000)
  cond_0_vs_10000_tbl = compare_abundances(a549_ccm, dose_0, dose_10000)
  
  # debug(plot_contrast)
  p = plot_contrast(a549_ccm, 
                cond_0_vs_10_tbl, 
                p_value_thresh = pvalue, 
                fc_limits = c(-3,3), 
                plot_labels = F) 
  ggsave(p, filename = paste0("sci_plex/",
           drug_name,"_0_10_pval_",pvalue,".png"), width = 7, height=7)
  
  p = plot_contrast(a549_ccm, 
                cond_0_vs_100_tbl, 
                p_value_thresh = pvalue, 
                fc_limits = c(-3,3), 
                plot_labels = F) 
  ggsave(p, filename = paste0("sci_plex/",
           drug_name,"_0_100_pval_", pvalue,".png"), width = 7, height=7)
  
  p = plot_contrast(a549_ccm, 
                cond_0_vs_1000_tbl, 
                p_value_thresh = pvalue, 
                fc_limits = c(-3,3), 
                plot_labels = F) 
  
  ggsave(p, filename = paste0("sci_plex/",
           drug_name,"_0_1000_pval_",pvalue,".png"), width = 7, height=7)
  
  # undebug(plot_contrast)
  p = plot_contrast(a549_ccm, 
                cond_0_vs_10000_tbl, 
                p_value_thresh = pvalue, 
                fc_limits = c(-3,3), 
                plot_labels = F)
  ggsave(p, filename = paste0("sci_plex/",
            drug_name,"_0_10000_pval_",pvalue,".png"), width = 7, height=7)
  

  # cond_0_vs_10_tbl %>% select(log_abund_x, log_abund_y, delta_log_abund, delta_p_value) %>% head()
  # cond_0_vs_100_tbl %>% select(log_abund_x, log_abund_y, delta_log_abund, delta_p_value) %>% head()
  # cond_0_vs_1000_tbl %>% select(log_abund_x, log_abund_y, delta_log_abund, delta_p_value) %>% head()
  # cond_0_vs_10000_tbl %>% select(log_abund_x, log_abund_y, delta_log_abund, delta_p_value) %>% head()
  
}


plot_10k <- function(drug_name, pvalue) {
  
  dose_0 = estimate_abundances(a549_ccm, make_tbl(drug_name, 0))
   dose_10000 = estimate_abundances(a549_ccm, make_tbl(drug_name, log10(10001)))
  cond_0_vs_10000_tbl = compare_abundances(a549_ccm, dose_0, dose_10000)
  
  p = plot_contrast(a549_ccm, 
                    cond_0_vs_10000_tbl, 
                    p_value_thresh = pvalue, 
                    fc_limits = c(-3,3), 
                    plot_labels = F)
  ggsave(p, filename = paste0("sci_plex/",
                              drug_name,"_hdac_0_10000_pval_",pvalue,".png"), width = 7, height=7)
}

for (d in hdac) {
  plot_10k(d, pvalue=0.05)
}


# Debug Vechicle -------------------------------------------------------------
# a549_tram_ccs = a549_ccs[,colData(a549_ccs)$product_name == "Trametinib (GSK1120212)"]
# a549_vehicle_ccs = a549_ccs[,colData(a549_ccs)$product_name == "Vehicle"]
# SingleCellExperiment::counts(a549_tram_ccs)
# SingleCellExperiment::counts(a549_vehicle_ccs)

plot_series(a549_ccm, "S2673", pvalue = 1)
plot_series(a549_ccm, "S1095", pvalue = 1)
plot_series(a549_ccm, "S1297", pvalue = 1)
plot_series(a549_ccm, "S1515", pvalue = 1)

plot_series(a549_ccm, "S2673", pvalue = 0.05)
plot_series(a549_ccm, "S1095", pvalue = 0.05)
plot_series(a549_ccm, "S1297", pvalue = 0.05)
plot_series(a549_ccm, "S1515", pvalue = 0.05)

plot_series(a549_ccm, "S1096", pvalue = 0.05)



plot_series(a549_ccm, "S1451", pvalue = 0.05)
plot_series(a549_ccm, "S1529", pvalue = 0.05)

return_catalog("Aurora A Inhibitor I",cat_to_prod) #"S1451"
return_catalog("Hesperadin",cat_to_prod) #"S1529"



# S2673# 0, 10, 100, 1000, 10000
return_catalog("Trametinib (GSK1120212)", cat_to_prod) #S2673

trametinib_dose_0 = estimate_abundances(a549_ccm, make_tbl("S2673", 0))
trametinib_dose_10 = estimate_abundances(a549_ccm, make_tbl("S2673", log10(11)))
trametinib_dose_100 = estimate_abundances(a549_ccm, make_tbl("S2673", log10(101)))
trametinib_dose_1000 = estimate_abundances(a549_ccm, make_tbl("S2673", log10(1001)))
trametinib_dose_10000 = estimate_abundances(a549_ccm, make_tbl("S2673", log10(10001)))
trametinib_cond_0_vs_10000_tbl = compare_abundances(a549_ccm, trametinib_dose_0, trametinib_dose_10000)

trametinib_dose_0 %>% select(log_abund, log_abund_sd)
trametinib_dose_10000 %>% select(log_abund, log_abund_sd)
trametinib_cond_0_vs_10000_tbl %>% select(S2673_x, S2673_y,
                                          log_abund_x, log_abund_sd_x, 
                                          log_abund_y, log_abund_sd_y, 
                                          delta_log_abund, delta_p_value)



plot_contrast(a549_ccm, 
              trametinib_cond_0_vs_10000_tbl, 
              p_value_thresh = 0.05) %>% 
  ggsave(filename = "sci_plex/trametinib_0_10000_pval05.png", width = 10, height=7)


my_plot_cells(a549_ccm, cond_b_vs_a_tbl = cond_0_vs_10_tbl)
my_plot_cells(a549_ccm, cond_b_vs_a_tbl = cond_10_vs_100_tbl)

return_catalog("Dacinostat (LAQ824)", cat_to_prod) #S1095




return_catalog("Epothilone A", cat_to_prod) #S1297
epoth_dose_0 = estimate_abundances(a549_ccm, make_tbl("S1297", 0))
epoth_dose_10000 = estimate_abundances(a549_ccm, make_tbl("S1297", log10(10001)))

epoth_cond_0_vs_10000_tbl = compare_abundances(a549_ccm, epoth_dose_0, epoth_dose_10000)

plot_contrast(a549_ccm, 
              epoth_cond_0_vs_10000_tbl, 
              p_value_thresh = 0.05) %>% 
  ggsave(filename = "sci_plex/epothilone_0_10000_pval05.png", width = 10, height=7)




return_catalog("Quisinostat (JNJ-26481585) 2HCl", cat_to_prod) #S1096
quis_dose_0 = estimate_abundances(a549_ccm, make_tbl("S1096", 0))
quis_dose_10 = estimate_abundances(a549_ccm, make_tbl("S1096", 1))
quis_dose_100 = estimate_abundances(a549_ccm, make_tbl("S1096", 2))
quis_dose_1000 = estimate_abundances(a549_ccm, make_tbl("S1096", 3))
quis_dose_10000 = estimate_abundances(a549_ccm, make_tbl("S1096", log10(10001)))

quis_cond_0_vs_10_tbl = compare_abundances(a549_ccm, quis_dose_0, quis_dose_10)
quis_cond_10_vs_100_tbl = compare_abundances(a549_ccm, quis_dose_10, quis_dose_100)
quis_cond_100_vs_1000_tbl = compare_abundances(a549_ccm, quis_dose_100, quis_dose_1000)
quis_cond_0_vs_10000_tbl = compare_abundances(a549_ccm, quis_dose_0, quis_dose_10000)

plot_contrast(a549_ccm, quis_cond_0_vs_10000_tbl, p_value_thresh = 0.05) %>% 
  ggsave(filename = "sci_plex/quisinostat_0_10000_pval05.png", width = 10, height=7)


plot_contrast(a549_ccm, quis_cond_0_vs_10_tbl, p_value_thresh = 1.0)
plot_contrast(a549_ccm, quis_cond_10_vs_100_tbl, p_value_thresh = 1.0)
plot_contrast(a549_ccm, quis_cond_100_vs_1000_tbl, p_value_thresh = 1.0)

return_catalog("Pracinostat (SB939)", cat_to_prod)
pracinostat_dose_0 = estimate_abundances(a549_ccm, make_tbl("S1515", 0))
pracinostat_dose_10 = estimate_abundances(a549_ccm, make_tbl("S1515", log10(10)))
pracinostat_dose_100 = estimate_abundances(a549_ccm, make_tbl("S1515", log10(100)))
pracinostat_dose_1000 = estimate_abundances(a549_ccm, make_tbl("S1515", log10(1000)))
pracinostat_dose_10000 = estimate_abundances(a549_ccm, make_tbl("S1515", log10(10000)))

pracinostat_cond_0_vs_10_tbl = compare_abundances(a549_ccm, pracinostat_dose_0, pracinostat_dose_10)
pracinostat_cond_0_vs_100_tbl = compare_abundances(a549_ccm, pracinostat_dose_0, pracinostat_dose_100)
pracinostat_cond_0_vs_1000_tbl = compare_abundances(a549_ccm, pracinostat_dose_0, pracinostat_dose_1000)
pracinostat_cond_0_vs_10000_tbl = compare_abundances(a549_ccm, pracinostat_dose_0, pracinostat_dose_10000)

plot_contrast(a549_ccm, pracinostat_cond_0_vs_10000_tbl, p_value_thresh = 0.05) %>% 
  ggsave(filename = "sci_plex/pracinostat_0_10000_pval05.png", width = 10, height=7)

# plot_contrast(a549_ccm, pracinostat_cond_0_vs_10_tbl, p_value_thresh = 1.0) %>% 
#   ggsave(filename = "sci_plex/pracinostat_0_10_pval_1.png", width = 10, height=7)
# plot_contrast(a549_ccm, pracinostat_cond_10_vs_100_tbl, p_value_thresh = 1.0) %>% 
#   ggsave(filename = "sci_plex/pracinostat_10_100_pval_1.png", width = 10, height=7)
# plot_contrast(a549_ccm, pracinostat_cond_100_vs_1000_tbl, p_value_thresh = 1.0) %>% 
#   ggsave(filename = "sci_plex/pracinostat_100_1000_pval_1.png", width = 10, height=7)
# plot_contrast(a549_ccm, pracinostat_cond_1000_vs_10000_tbl, p_value_thresh = 1.0) %>% 
#   ggsave(filename = "sci_plex/pracinostat_1000_10000_pval_1.png", width = 10, height=7)

plot_contrast(a549_ccm, pracinostat_cond_0_vs_10000_tbl, p_value_thresh = 0.01) 

corr_edges = collect_pln_graph_edges(a549_ccm,
                                     pracinostat_cond_0_vs_10000_tbl,
                                    log_abundance_thresh=-5)

corr_edges = corr_edges #%>% filter(from_delta_p_value <0.05, to_delta_p_value <0.05)
neg_edges = corr_edges %>% 
  filter(edge_type!="undirected", pcor < 0) %>% 
  select(from,to)
i=1

green_edges = lapply(1:nrow(neg_edges), function(i) { 
  source = neg_edges[i,]$from
  target = neg_edges[i,]$to
get_pos_pcor_edges(centroids(a549_ccs),
                   corr_edges,
                   alpha = 1,
                   beta = 1,
                   gamma = 1,
                   sum_weights=T) %>%
  calc_shortest_path(source, target) %>% select(from,to,weight)
})

a549_green_edge_df = do.call(rbind, green_edges) %>% add_umap_coords(centroids(a549_ccs))

a549_green_edge_df = green_edges[[1]]%>% add_umap_coords(centroids(a549_ccs))

bp = return_baseplot(a549_ccm, quis_cond_0_vs_10000_tbl)

bp  + geom_segment(data = a549_green_edge_df,
                   aes(x = umap_to_1,
                       y = umap_to_2,
                       xend=umap_from_1,
                       yend = umap_from_2),
                   color="black" , size=1.7) #+
  # geom_segment(data = a549_green_edge_df,
  #              aes(x = umap_from_1,
  #                  y = umap_from_1,
  #                  xend=(umap_to_1+umap_from_1)/2,
  #                  yend = (umap_to_2+umap_from_2)/2),
  #              color="black", size=1.7,
  #              linejoin='mitre',
  #              arrow = arrow(type="closed", angle=30, length=unit(1, "mm")))

# FC of sci-plex cell states --------------------------------------------------



# Heatmap of dose dependent effects (e.g. PLN coefficients) of each drug on each cell state, 
# showing that effect coefficients cluster according to compound class. 

# ------------------------------------------------------------------------------

pred_abund_mat = lapply(product_names, function(p) {
  dose = estimate_abundances(a549_ccm, tibble(product_name = p, log_dose = log10(1001)))
  dose$log_abund
}) %>% bind_cols()

names(pred_abund_mat) = product_names


pred_abund_mat = cbind(dose_0$log_abund, 
                       dose_d$log_abund, 
                       quisinostat_dose_10$log_abund)
# colnames(pred_abund_mat) = c("0", "10","100", "1000", "10000")

pheatmap::pheatmap(t(pred_abund_mat))

# -----------------------------------------------------------------------------

log_res = coef_df[["log_dose"]]
names(log_res) = 1:19

my_plot_cells(a549_ccm, residuals = log_res)

vehicle = coef_df[["product_nameVehicle"]]
names(vehicle) = 1:19
my_plot_cells(a549_ccm, residuals = vehicle)


dacinostat = coef_df[["product_nameDacinostat (LAQ824)"]]
names(dacinostat) = 1:19
my_plot_cells(a549_ccm, residuals = dacinostat)

quisinostat = coef_df[["product_nameQuisinostat (JNJ-26481585) 2HCl"]]
names(quisinostat) = 1:19
my_plot_cells(a549_ccm, residuals = quisinostat)

data.frame(vehicle, dacinostat, quisinostat) %>% t() %>% pheatmap::pheatmap()

# ------------------

# try to recreate supplementary figure 6 ---------------------------------------

# rows are drugs by group
# columns are clusters numbers clustered 
# log2 fold enrichment over vehicle

coef_df = as.data.frame(coef(a549_ccm@best_model))

std_error_df = standard_error(a549_ccm@best_model)
intercept = coef_df %>% select("(Intercept)") %>%rownames_to_column("cluster") 

std_error_t = std_error_df%>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("catalog_number") %>% 
  mutate("catalog_number" = gsub("as.numeric", "",catalog_number)) %>% 
  mutate("catalog_number" = gsub("[[:punct:]]", "",catalog_number)) %>% 
  left_join(cat_to_prod, by = "catalog_number") %>% 
  filter(!catalog_number %in% c("Intercept")) %>% 
  column_to_rownames("product_name") %>% 
  select(-catalog_number) %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column("cluster") %>% 
  pivot_longer(-cluster)%>% rename("se"=value)


coef_df_t = coef_df %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("catalog_number") %>% 
  mutate("catalog_number" = gsub("as.numeric", "",catalog_number)) %>% 
  mutate("catalog_number" = gsub("[[:punct:]]", "",catalog_number)) %>% 
  left_join(cat_to_prod, by = "catalog_number") %>% 
  filter(!catalog_number %in% c("Intercept")) %>% 
  column_to_rownames("product_name") %>% 
  select(-catalog_number) 

# are each drug term significant are they adding anything to the model

# which ones are not significantly different than 0 
# calculate std error using wald test


adj_heatmap = coef_df_t %>%
  t() %>% as.data.frame %>% 
  rownames_to_column("cluster") %>% 
  pivot_longer(-cluster) %>%
  left_join(intercept, by="cluster") %>% 
  left_join(std_error_t, by=c("cluster",'name')) %>% 
  mutate(pvalue = pnorm(abs(value), sd = se, lower.tail=FALSE)) %>% 
  mutate(value = ifelse(pvalue <= 0.05, value, 0)) %>% 
  select(-c(`(Intercept)`, se, pvalue)) %>% 
  pivot_wider(names_from = "name", values_from = "value") %>% 
  column_to_rownames("cluster")

  
mat = as.data.frame(t(adj_heatmap))

pathway_df = colData(a549_cds) %>% 
    as.data.frame() %>%
  group_by(product_name, pathway_level_1) %>% tally() %>% select(-n) %>% 
  column_to_rownames("product_name")


pathway2_df = colData(a549_cds) %>% 
  as.data.frame() %>%
  group_by(product_name, pathway_level_1, pathway_level_2) %>% tally() %>% select(-n) %>% 
  column_to_rownames("product_name")



proliferation_index_df = colData(a549_cds) %>% 
  as.data.frame() %>%
  group_by(cluster) %>% summarize(PI=mean(proliferation_index)) %>% 
  column_to_rownames("cluster")

pathway_df = pathway_df %>% arrange(pathway_level_1)
mat = mat[rownames(pathway_df),]

res = pheatmap::pheatmap(mat, cluster_rows = F, cluster_cols = T, 
                         clustering_method = "ward.D2",
                         annotation_row = pathway_df, 
                         annotation_col = proliferation_index_df, 
                         treeheight_row = 0, treeheight_col = 0)
res
#switch 


row_order = rownames(mat)[res$tree_row$order]
col_order = colnames(mat)[res$tree_col$order]
  
mat = mat %>% 
    rownames_to_column("product_name") %>% 
    pivot_longer(-product_name)

mat$product_name <- factor(mat$product_name, levels=rownames(pathway_df))
mat$name <- factor(mat$name, levels=col_order)

ggplot(mat, aes(x=name,y = product_name, fill=value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#3234a9", mid = "white", 
                       high= "darkred", na.value="black", name="") 
  

# subset to just epigenetic regulation 
er_df = pathway_df %>% filter(pathway_level_1 == "Epigenetic regulation")
sub_mat = mat[rownames(er_df),]

path_df = pathway2_df%>% 
  filter(pathway_level_1 == "Epigenetic regulation") 
  

sub_res = pheatmap::pheatmap(sub_mat, cluster_rows = T, cluster_cols = T, 
                         clustering_method = "ward.D2",
                         annotation_row = path_df, 
                         annotation_col = proliferation_index_df, 
                         treeheight_row = 0, treeheight_col = 0)
sub_res
