library(ggplot2)
library(tibble)
library(dplyr)
library(monocle3)
library(hooke)
library(msigdbr)

#pap_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/macrophages/PAP/monos.al.cds_2021-11-01.RDS")

# Peak matrix with LSI looks a lot better than ArchR gene scores in terms of UMAP structure
atac_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/sci-plex_ATAC_RNA/ATAC/sciPlexATAC_peaks_cds.RDS")
#atac_cds = preprocess_cds(atac_cds, method="LSI")
#atac_cds = reduce_dimension(atac_cds, preprocess_method="LSI")
#atac_cds = cluster_cells(atac_cds)
#plot_cells(atac_cds, color_cells_by="treatment") + facet_wrap(~treatment)

#atac_cds = atac_cds[!is.na(rowData(atac_cds)$gene_short_name) & rowData(atac_cds)$gene_short_name != "NA",] # NA is actually coded as a character string!!
atac_cds = preprocess_cds(atac_cds, method="LSI")
atac_cds = reduce_dimension(atac_cds, preprocess_method="LSI")
atac_cds = cluster_cells(atac_cds)
plot_cells(atac_cds, color_cells_by="treatment") + facet_wrap(~treatment)

agg_atac_mat = monocle3:::my.aggregate.Matrix(counts(atac_cds), rowData(atac_cds) %>% as.data.frame %>% dplyr::pull(nearestGene))

#atac_cds = preprocess_cds(atac_cds, method="LSI")
#atac_cds = reduce_dimension(atac_cds, preprocess_method="LSI")
#atac_cds = cluster_cells(atac_cds)

#atac_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/sci-plex_ATAC_RNA/ATAC/sciPlexATAC_genescore_cds.RDS")
#atac_cds = preprocess_cds(atac_cds)
#atac_cds = reduce_dimension(atac_cds)
#atac_cds = cluster_cells(atac_cds)

rna_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/sci-plex_ATAC_RNA/RNA/sciPlex2_cds_NB5processed.RDS")
rna_cds = preprocess_cds(rna_cds)
rna_cds = reduce_dimension(rna_cds)
rna_cds = cluster_cells(rna_cds)

#cluster_gene_fit = fit_models(rna_cds, "~cluster", cores=4)
cluster_markers = top_markers(rna_cds, genes_to_test_per_group = 1000, marker_sig_test=TRUE, cores=8)
sig_cluster_markers = cluster_markers %>% filter(marker_test_q_value < 0.05)
###

#atac_cds = atac_cds[!duplicated(rowData(atac_cds)$gene_short_name),]
#atac_cds = preprocess_cds(atac_cds, method="LSI")
#atac_cds = reduce_dimension(atac_cds, preprocess_method="LSI")
#plot_cells(atac_cds)


atac_genescore_mat = agg_atac_mat
#row.names(atac_genescore_mat) = rowData(atac_cds)$gene_short_name

atac_cell_metadata = colData(atac_cds) %>% as.data.frame
atac_cell_metadata$assay = "ATAC"
atac_cell_metadata$treatment[atac_cell_metadata$dose == 0] = "vehicle"

atac_cell_metadata$grouping = paste(atac_cell_metadata$assay,
                                    atac_cell_metadata$treatment,
                                    atac_cell_metadata$Relative_dose,
                                    atac_cell_metadata$replicate,
                                    atac_cell_metadata$well_oligo, sep=".")
atac_cell_metadata = atac_cell_metadata %>% dplyr::select(assay, treatment, dose, replicate, well_oligo, grouping)
atac_gene_metadata = data.frame(row.names=row.names(atac_genescore_mat), gene_short_name=row.names(atac_genescore_mat))
atac_cds = new_cell_data_set(atac_genescore_mat, cell_metadata = atac_cell_metadata, gene_metadata=atac_gene_metadata)

atac_cds = preprocess_cds(atac_cds, method="LSI")
atac_cds = reduce_dimension(atac_cds, preprocess_method="LSI")
atac_cds = cluster_cells(atac_cds)
plot_cells(atac_cds, color_cells_by="treatment") + facet_wrap(~treatment)

subset_atac_cds = atac_cds[intersect(sig_cluster_markers$gene_id, row.names(atac_cds)),]
subset_atac_cds = estimate_size_factors(subset_atac_cds)
subset_atac_cds = preprocess_cds(subset_atac_cds, method="LSI")
subset_atac_cds = reduce_dimension(subset_atac_cds, preprocess_method="LSI")
subset_atac_cds = cluster_cells(subset_atac_cds)
plot_cells(subset_atac_cds, color_cells_by="treatment") + facet_wrap(~treatment)


#atac_gene_metadata = rowData(atac_cds) %>% as.data.frame
#row.names(atac_gene_metadata) = rowData(atac_cds)$gene_short_name
#atac_gene_metadata = atac_gene_metadata[row.names(atac_genescore_mat),]
#atac_gene_metadata = data.frame(row.names=row.names(atac_genescore_mat), gene_short_name=row.names(atac_genescore_mat))
#atac_cds = new_cell_data_set(atac_genescore_mat, cell_metadata = atac_cell_metadata, gene_metadata=atac_gene_metadata)

#atac_genes = rowData(rna_cds)[cluster_markers$gene_id,]$gene_short_name %>% unique
#atac_genes = intersect(atac_genes, row.names(atac_cds))
#de_genes_atac_cds = atac_cds[atac_genes,]

#de_genes_atac_cds = preprocess_cds(de_genes_atac_cds)
#de_genes_atac_cds = reduce_dimension(de_genes_atac_cds)
#plot_cells(de_genes_atac_cds, color_cells_by="dose") + facet_wrap(~treatment)

#size_factors(rna_cds) = 1
# rna_expr_mat = monocle3::aggregate_gene_expression(rna_cds,
#                                                    rowData(rna_cds) %>% as.data.frame %>% dplyr::select(id, gene_short_name),
#                                                    norm_method="size_only",
#                                                    pseudocount=0,
#                                                    scale_agg_values=FALSE)

rna_cds = rna_cds[!duplicated(rowData(rna_cds)$gene_short_name),]

rna_expr_mat = exprs(rna_cds)
rna_cell_metadata = colData(rna_cds) %>% as.data.frame
rna_cell_metadata$treatment = rna_cell_metadata$treatmentRNA_broad
rna_cell_metadata$assay = "RNA"

rna_cell_metadata$grouping = paste(rna_cell_metadata$assay,
                                   rna_cell_metadata$treatment,
                                   rna_cell_metadata$dose,
                                   rna_cell_metadata$replicate,
                                   rna_cell_metadata$well_oligo, sep=".")
rna_cell_metadata = rna_cell_metadata %>% dplyr::select(assay, treatment, dose, replicate, well_oligo, grouping)

rna_gene_metadata = rowData(rna_cds) %>% as.data.frame
row.names(rna_expr_mat) = rna_gene_metadata$gene_short_name
row.names(rna_gene_metadata) = rna_gene_metadata$gene_short_name
rna_cds = new_cell_data_set(rna_expr_mat, cell_metadata = rna_cell_metadata, gene_metadata=rna_gene_metadata)


atacAndRNA_cds = combine_cds(list(atac_cds, rna_cds), sample_col_name="experiment", keep_all_genes=FALSE)

atacAndRNA_cds = atacAndRNA_cds[intersect(cluster_markers$gene_id, row.names(atacAndRNA_cds)),]

colData(atacAndRNA_cds)$perturb_group = interaction(colData(atacAndRNA_cds)$assay,
                                                    colData(atacAndRNA_cds)$treatment,
                                                    colData(atacAndRNA_cds)$dose)
#perturb_counts = table(colData(atacAndRNA_cds)$perturb_group)
#atacAndRNA_cds = atacAndRNA_cds[,colData(atacAndRNA_cds)$perturb_group %in% names(perturb_counts)[which(perturb_counts > 50)]]
#colData(atacAndRNA_cds)$perturb_group = droplevels(colData(atacAndRNA_cds)$perturb_group)

colData(atacAndRNA_cds)$log_dose = log10(colData(atacAndRNA_cds)$dose + 0.1)

atacAndRNA_cds_coldata = colData(atacAndRNA_cds) %>% as.data.frame
for (cn in unique(atacAndRNA_cds_coldata$treatment)) {
  if (cn != "vehicle") {
    atacAndRNA_cds_coldata = atacAndRNA_cds_coldata %>%
      mutate("{cn}":= ifelse(treatment==cn, log_dose, 0))
    #colData(atacAndRNA_cds)[[cn]] = a549_ccs_coldata[[cn]]
  }
}
colData(atacAndRNA_cds)$SAHA = atacAndRNA_cds_coldata$SAHA
colData(atacAndRNA_cds)$Dex = atacAndRNA_cds_coldata$Dex
colData(atacAndRNA_cds)$BMS = atacAndRNA_cds_coldata$BMS
colData(atacAndRNA_cds)$Nutlin = atacAndRNA_cds_coldata$Nutlin


atacAndRNA_cds = estimate_size_factors(atacAndRNA_cds)


atacAndRNA_cds = preprocess_cds(atacAndRNA_cds, method="LSI", num_dim=50)
atacAndRNA_cds = align_cds(atacAndRNA_cds,
                           preprocess_method="LSI",
                           alignment_group="assay"
                           )

atacAndRNA_cds = reduce_dimension(atacAndRNA_cds, preprocess_method="Aligned")

plot_cells(atacAndRNA_cds, color_cells_by="assay") + facet_wrap(~treatment)

#colData(atacAndRNA_cds)$log_dose = log(colData(atacAndRNA_cds)$dose + 0.1)
# Harmony-based:
# Try using harmony between PC and UMAP levels of reduction and see if that works better?
harmonyPCs = harmony::HarmonyMatrix(reducedDims(atacAndRNA_cds)$LSI,
                                    as.data.frame(colData(atacAndRNA_cds)),
                                    #vars_use = c("assay", "treatment", "dose"),
                                    #vars_use = c("assay", "Dex", "BMS", "SAHA", "Nutlin"),
                                    #vars_use = c("assay", "treatment"),
                                    vars_use = c("assay"),
                                    theta=0,
                                    do_pca = FALSE,
                                    #theta=c(4, 4),
                                    #lambda=c(0.5, 1),
                                    verbose=TRUE)

harmonyCDS = atacAndRNA_cds
reducedDims(harmonyCDS)$PCA = harmonyPCs

harmonyCDS = reduce_dimension(harmonyCDS, preprocess_method="PCA", reduction_method="UMAP")
plot_cells(harmonyCDS, color_cells_by="assay") + facet_wrap(~treatment)

# drug dose is product_name * dose
a549_ccs = new_cell_count_set(harmonyCDS,
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








