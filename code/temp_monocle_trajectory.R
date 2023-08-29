pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "Seurat", "CytoTRACE", "reticulate", "monocle3")
for (pkg in pkgs){
  if(pkg == "reticulate"){
    #suppressPackageStartupMessages(library(pkg, character.only = T))
    #use_python('/opt/simplehpc/miniconda3/bin/python')    ## another version of python has already
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = T))
  }
}






reticulate::virtualenv_create(envname = 'python_environment', 
                              python= '/opt/simplehpc/miniconda3/bin/python')
reticulate::virtualenv_install("python_environment", 
                               packages=c('pandas','catboost'))
reticulate::use_virtualenv("python_environment",required = TRUE)


mef2d <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/mef2d.rds")

sct_assay <- GetAssayData(mef2d)
rna_assay <- GetAssayData(mef2d, assay = "RNA", slot = "data")
meta_mef2d <- mef2d@meta.data
# check annotation 
table(meta_mef2d[, c("orig.ident", "cell_types2")])
# 069: mostly pre-B, some tumor1, matureB and pre-proB, little NK-T
# M2-4: mostly tumor cells, some NK-T, and pre-B, erythroid and matureB
# E1: mostly tumor cells, many NK-T, some matureB and preB

# calculate contribution?

mef2d_sub <- subset(mef2d, `cell_types2` %notin% c("NK_T", "erythroid_cell", "myeloid"))

cell_color <- c("HSC_MPP" = "indianred1", "CLP" = "indianred3", "pre_proB" = "darkseagreen1",
                 "proB" = "lightgreen", "preB" = "darkolivegreen1",  "preB_I" = "yellowgreen",
                 "preB_II" = "chartreuse4", "imatureB" = "#3CB371", "matureB" = "darkgreen",
                 "erythroid_cell" = "lightpink2", "NK_T" = "#9575aa", "myeloid" = "deeppink4", "tumor1" = "grey",
                 "tumor_cell1" = "yellow4", "tumor_cell2" = "aquamarine2", "tumor_cell3" = "lightseagreen")

samples <- unique(meta_mef2d$orig.ident)
sample_color <- c("#ee5d5a", "#9f2f6a", "#34956c", "#6b6498","#b29bc9", "#71a3a2", 
                  "#c88978",  "#a15217", "#ce662f", "#cc9b32", "#53096a", "#6569c2", "#66676c")[1:length(samples)]
names(sample_color) <- samples


psub1 <- DimPlot(mef2d_sub, group.by="cell_types2", reduction='umap', label = T) + NoLegend() + 
  scale_color_manual(values = cell_color)

ggsave(psub1, filename = "/cluster/home/yjliu_jh/projects/temp/umap_allsamples.pdf", width = 5, height = 4)

psubsep1 <- DimPlot(mef2d_sub, split.by="cell_types2", reduction='umap', label = F) + NoLegend() + 
  scale_color_manual(values = cell_color)

ggsave(psubsep1, filename = "/cluster/home/yjliu_jh/projects/temp/umap_all_sep.pdf", width = 12, height = 3)




# CytoTRACE plot ----------------------------------------------------------

cyt_res_mef2d_sub <- CytoTRACE(as.matrix(mef2d_sub[["RNA"]]@counts), ncores = 12)

# plot the CytoTRACE results

 pdf("/cluster/home/yjliu_jh/projects/temp/plotCytoTRACE_umap_sub.pdf", width = 12, height = 4)
 plotCytoTRACE(cyt_res_mef2d_sub, colors = NULL, emb = mef2d_sub@reductions$umap@cell.embeddings, 
               phenotype = cell_type##, outputDir = "./tumor_b_cell_umap_"
               )
 dev.off()

cell_type <- as.character(mef2d_sub@meta.data$cell_types2)
orig_ident <- as.character(mef2d_sub@meta.data$orig.ident)
names(cell_type) <-rownames(mef2d_sub@meta.data)
names(orig_ident) <- rownames(mef2d_sub@meta.data)

meta_sub <- mef2d_sub@meta.data
emb <- as.matrix(mef2d_sub@reductions$umap@cell.embeddings)
mat <- t(cyt_res_mef2d_sub$exprMatrix)[rownames(emb), ]
cyto <- cyt_res_mef2d_sub$CytoTRACE[rownames(emb)]

cyto_plot_datfm <- data.frame(emb, CytoTRACE_score = cyto, HDAC9 = mat[, "HDAC9"], 
                              MEF2D = mat[, "MEF2D"], meta_sub)


temp_color <- RColorBrewer::brewer.pal(11, "Spectral")
temp_color[6] <- "gold"
rbPal <- colorRampPalette(temp_color)

p_cyt <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
  ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
                                                guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
                                                breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
  theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
        axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
        legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) # + ggtitle("CytoTRACE Score")
p_cyt_celltypes <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
  ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
                                                guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
                                                breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
  theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
        axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
        legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) + facet_wrap(~cell_types, nrow = 1)
p_cyt_sample <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
  ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
                                                guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
                                                breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
  theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
        axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
        legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) + facet_wrap(~`orig.ident`, nrow = 1)

p_cyt_celltypes_sample <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE_score)) + geom_point(size = .5) + 
  ggpubr::theme_pubr() + scale_colour_gradientn(name = "CytoTRACE score", colours = rev(rbPal(50)),
                                                guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"),
                                                breaks = seq(0, 1, 0.2), labels = c("0.0 (More diff.)", 0.2, 0.4, 0.6, 0.8, "1.0 (Less diff.)")) +
  theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
        plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
        axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
        legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) + facet_wrap(~`orig.ident` + cell_types)
ggsave("cytotrace_split_celltype_sample_idv.pdf", (p_cyt_celltypes/p_cyt_sample + plot_layout(guides = "collect")), width = 11, height = 7)

ggsave("cytotrace_split_celltype_sample.pdf", p_cyt_celltypes_sample, width = 11, height = 11)

# p_celltype <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = cell_types)) + geom_point(size = .5) + 
#   ggpubr::theme_pubr() + scale_colour_manual(values = celltypes_color[c("Tumor cell1", "Tumor cell2", "Tumor cell3", "B cell")]) +
#   theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
#         plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
#         axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
#         legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) + ggtitle("Cell types")
# 
# p_sample <- ggplot(cyto_plot_datfm, aes(x = UMAP_1, y = UMAP_2, color = orig.ident)) + geom_point(size = .5) + 
#   ggpubr::theme_pubr() + scale_colour_manual(values = sample_color) +
#   theme(legend.text = ggplot2::element_text(size = 6), legend.title = ggplot2::element_text(size = 7), 
#         plot.title = ggplot2::element_text(size = 10, hjust = 0.5), axis.title.x = ggplot2::element_text(size = 8), 
#         axis.title.y = ggplot2::element_text(size = 7), axis.text = ggplot2::element_text(size = 6), 
#         legend.position = "right", plot.margin = ggplot2::unit(c(0.5, 1, 0.5, 1), "cm")) + ggtitle("Samples")


# monocle3 trajectory analysis --------------------------------------------


data <- GetAssayData(mef2d_sub, assay = "RNA", slot = 'counts')
# meta
fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
# built CDS
cds <- new_cell_data_set(data, cell_metadata = meta_sub, gene_metadata = fd)

##Pre-process the data
#Monocle principal components
#preprocess_cds equivalent of seurat NormalizeData+ScaleData+RunPCA
# cds <- cds %>% preprocess_cds(., num_dim = 50) %>% 
#   reduce_dimension(., preprocess_method = "PCA", umap.n_neighbors = 30) %>% 
#   cluster_cells(., k = 20)

# reduce_dimension, umap.fast_sgd, umap.min.dist may largely affect branches distribution
cds_test <- cds %>% preprocess_cds(., num_dim = 50) %>% 
  reduce_dimension(., preprocess_method = "PCA", umap.n_neighbors = 50, umap.fast_sgd = F, umap.min_dist = 1.5) %>% 
  cluster_cells(., k = 30, partition_qval = .05) %>% learn_graph(close_loop = T)

# cds.embed <- cds@int_colData$reducedDims$UMAP
# int.embed <- Embeddings(sel_sce, reduction = "umap")
# int.embed <- int.embed[rownames(cds.embed),]
# cds@int_colData$reducedDims$UMAP <- int.embed
# p1 <- plot_cells(cds, color_cells_by = "cell_types", label_groups_by_cluster=FALSE, label_leaves=F,
#                  label_branch_points=F, group_label_size=4, cell_size=0.7) + theme(legend.position="none") + 
#   scale_color_manual(values = celltypes_color) + ggtitle("Original root")
# 
# cds <- cds %>% learn_graph(., close_loop = T) %>% order_cells()
# p2 <- plot_cells(cds, color_cells_by = "cell_types", label_groups_by_cluster=FALSE, label_leaves=F,
#            label_branch_points=F, group_label_size=4, cell_size=0.7) + theme(legend.position="none") + 
#   scale_color_manual(values = celltypes_color) + ggtitle("Reset root -- grouped by cell types")
# plot_cells(cds_test, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=F,
#            label_branch_points=F, group_label_size=4, cell_size=0.7) + 
#   theme(legend.position="none") + ggtitle("Reset root -- pseudotime")
# 

cds_test.embed <- cds_test@int_colData$reducedDims$UMAP
int.embed <- Embeddings(mef2d_sub, reduction = "umap")
int.embed <- int.embed[rownames(cds_test.embed),]
cds_test@int_colData$reducedDims$UMAP <- int.embed
# p1 <- plot_cells(cds_test, color_cells_by = "cell_types", label_groups_by_cluster=FALSE, label_leaves=F,
#                  label_branch_points=F, group_label_size=4, cell_size=0.7) + theme(legend.position="none") + 
#   scale_color_manual(values = celltypes_color) + ggtitle("Original root")


readr::write_rds(cds_test, "/cluster/home/yjliu_jh/projects/temp/cds_test.rds")
cds_test <- readr::read_rds("/cluster/home/yjliu_jh/projects/temp/cds_test.rds")

cds_test_reset <- cds_test %>% learn_graph(., close_loop = T) %>% order_cells()
readr::write_rds(cds_test_reset, "/cluster/home/yjliu_jh/projects/temp/cds_test_reset.rds")
cds_test_reset <- readr::read_rds("/cluster/home/yjliu_jh/projects/temp/cds_test_reset.rds")

p2 <- plot_cells(cds_test_reset, color_cells_by = "cell_types2", label_groups_by_cluster=FALSE, label_leaves=F,
                 label_branch_points=F, group_label_size=4, cell_size=0.7) + theme(legend.position="none") + 
  scale_color_manual(values = cell_color)
p3 <- plot_cells(cds_test_reset, color_cells_by = "pseudotime", label_groups_by_cluster=FALSE, label_leaves=F,
                 label_branch_points=F, group_label_size=4, cell_size=0.7) + 
  theme(legend.position="right", axis.text = element_text(size = 7), axis.title = element_text(size = 8), legend.text = element_text(size = 6)
        , legend.title = element_text(size = 8))# + 
# scale_colour_gradientn(name = "Pseudotime", colours = colorRampPalette(c("purple","yellow"))(50),
#                                  guide = ggplot2::guide_colourbar(ticks.colour = "black", ticks.linewidth = 1, frame.colour = "black"))


psub1 <- DimPlot(mef2d_sub, group.by="cell_types2", reduction='umap', label = T, label.size = 2) + 
  NoLegend() + scale_color_manual(values = cell_color)

p_traj <- psub1 + p_cyt + p3 + plot_layout(nrow = 1)
ggsave("/cluster/home/yjliu_jh/projects/temp/trajectory_cytotrace_monocle3_0109.pdf", p_traj, width = 10, height = 3)

# cds <- preprocess_cds(cds, num_dim = 100)
# cds <- reduce_dimension(cds)
# cds <- cluster_cells(cds, resolution=1e-5)

ciliated_cds_pr_test_res <- graph_test(cds_test_reset, neighbor_graph="principal_graph", cores=12)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds_test_reset[pr_deg_ids,], resolution=1e-2)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_test_reset)), 
                                cell_types = cds_test_reset@colData$cell_types)
agg_mat <- aggregate_gene_expression(cds_test_reset, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
















colors_cell <- c("Stellate" = "#F2ADBE",
                 "TFF1_Ductal" =  "#E8AB56", "MUC6_Ductal" = "#e8d056",
                 "MUC4_Ductal" = "#e88756", "KRT17_Ductal" = "#f26e50",
                 "KRT14_Ductal_cycling" = "#B4604E",
                 "KRT14_Ductal" = "#c7897b", "Ductal_acinar" = "#dfaea2",
                 "Acinar_ductal" = "#D6BDA1","Acinar" = "#edc9af")

#Extract data, phenotype data, and feature data from the SeuratObject
data <- ipmn_ductal@assays$RNA@counts

pd <- new('AnnotatedDataFrame', data = ipmn_ductal@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
ipmn_ductal_mono <- newCellDataSet(data,
                                   phenoData = pd,
                                   featureData = fd,
                                   lowerDetectionLimit = 0.5,
                                   expressionFamily = negbinomial.size())

ipmn_ductal_mono <- ipmn_ductal_mono %>% estimateSizeFactors() %>% estimateDispersions()

ipmn_ductal_mono <- detectGenes(ipmn_ductal_mono, min_expr = 0.1)
print(head(fData(ipmn_ductal_mono))) #17702
expressed_genes <- row.names(subset(fData(ipmn_ductal_mono),
                                    num_cells_expressed >= 10))#16978


diff_test_res <- differentialGeneTest(ipmn_ductal_mono[expressed_genes,],
                                      fullModelFormulaStr="~cell_types",cores=1)
head(diff_test_res)
deg <- subset(diff_test_res, qval < 0.001)
deg <- deg[order(deg$qval,decreasing=F),]

write_csv(diff_test_res, "ipmn_ductal_mono_diff_test_res.csv")
write_csv(deg, "ipmn_ductal_mono_DEG.csv")
readr::write_rds(diff_test_res, "diff_test_res.rds")
readr::write_rds(deg, "deg.rds")

#####
#choose genes that define a cell's progress
ordering_genes <- row.names(deg[1:2000,])

ipmn_ductal_mono_deg <- setOrderingFilter(ipmn_ductal_mono, ordering_genes)

pdf("ipmn_ductal_monocle_plot_ordering_genes.pdf")
plot_ordering_genes(ipmn_ductal_mono_deg)
dev.off()

#reduce data dimensionality
ipmn_ductal_mono_deg <- reduceDimension(ipmn_ductal_mono_deg, max_components = 2,
                                        method = 'DDRTree')

readr::write_rds(ipmn_ductal_mono_deg, "ipmn_ductal_mono_deg.rds")

#order cells along the trajectory
ipmn_ductal_mono_deg <- orderCells(ipmn_ductal_mono_deg)

pdf("plot_cell_clusters_monocle_pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(ipmn_ductal_mono_deg,color_by="Pseudotime", size=1,show_backbone=TRUE)
dev.off()

pdf("plot_cell_clusters_monocle_celltype.pdf",width = 6,height = 5)
plot_cell_trajectory(ipmn_ductal_mono_deg,color_by="cell_types", size=1,show_backbone=TRUE) +
  scale_color_manual(values = colors_cell) + theme(legend.position = "right")
dev.off()

pdf("plot_cell_clusters_monocle_celltype_faceted.pdf",width = 15,height = 8)
plot_cell_trajectory(ipmn_ductal_mono_deg, color_by = "cell_types") +
  facet_wrap("~cell_types", nrow = 2) +
  scale_color_manual(values = colors_cell) + theme(legend.position = "right")
dev.off()

pdf("plot_cell_clusters_monocle_celltype_faceted_sample.pdf",width = 15,height = 4)
plot_cell_trajectory(ipmn_ductal_mono_deg, color_by = "cell_types") +
  facet_wrap("~orig.ident", nrow = 1) +
  scale_color_manual(values = colors_cell) + theme(legend.position = "right")
dev.off()


pdf("plot_cell_clusters_monocle_celltype_tree.pdf",width = 10,height = 5)
p1 <- plot_cell_trajectory(ipmn_ductal_mono_deg, color_by = "cell_types") +
  scale_color_manual(values = colors_cell) +
  theme(legend.position = "none", panel.border = element_blank())
p2 <- plot_complex_cell_trajectory(ipmn_ductal_mono_deg, color_by = "cell_types") +
  scale_color_manual(values = colors_cell) +
  theme(legend.title = element_blank())
print(p1|p2)
dev.off()

df <- pData(ipmn_ductal_mono_deg)
p <- ggplot(df,aes(Pseudotime,colour = cell_types, fill = cell_types)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) + theme_classic2() +
  scale_color_manual(values = colors_cell) +
  scale_fill_manual(values = colors_cell)
ggsave("monocle_celltype_density.pdf", p, height = 4)

df <- pData(ipmn_ductal_mono_deg)
smplec_col <- c("cjl" = "#6EBEE6", "gyz" = "#F2ADBE", "wxf" = "#BE96BE", "wzl" = "#E8AB56", "wzl_n" = "#B4DCB4")
p <- ggplot(df,aes(Pseudotime,colour = orig.ident, fill = orig.ident)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) + theme_classic2() +
  scale_color_manual(values = smplec_col) +
  scale_fill_manual(values = smplec_col)
ggsave("monocle_sample_density.pdf", p, height = 4)

readr::write_rds(ipmn_ductal_mono_deg, "ipmn_ductal_mono_deg.rds")

diff_time_res <- differentialGeneTest(ipmn_ductal_mono_deg[ordering_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_time_res <- diff_time_res[,c(5,2,3,4,1,6,7)]
readr::write_csv(diff_time_res, "diff_time_res.csv")
diff_time_genes <- diff_time_res %>% pull(gene_short_name) %>% as.character()

pdf("monocle_pseudotime_heatmap.pdf", height = 8)
p <- plot_pseudotime_heatmap(ipmn_ductal_mono_deg[diff_time_genes,],
                             num_clusters = 5,
                             cores = 1,
                             show_rownames = T,
                             return_heatmap = T)
print(p)
dev.off()
readr::write_rds(p,"monocle_pseudotime_heatmap_obj.rds")
#p$tree_row %>% str()
clusters_tree <- cutree(p$tree_row, k = 5)
clusters_tree_df <- data.frame(clusters_tree)
clusters_tree_df[,1] <- as.character(clusters_tree_df[,1])
colnames(clusters_tree_df) <- "cluster"
clusters_tree_df <- clusters_tree_df %>% as_tibble(rownames = "gene")
#table(clusters_tree_df)
readr::write_csv(clusters_tree_df, "clusters_tree_df.csv")

pdf("plot_cell_clusters_monocle_state.pdf",width = 5,height = 4)
plot_cell_trajectory(ipmn_ductal_mono_deg, color_by = "State") +
  scale_color_manual(values = show_me_the_colors("single_cell2")) + theme(legend.position = "right")
dev.off()

ids <- bitr(clusters_tree_df$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
clusters_tree_df <- merge(clusters_tree_df, ids, by.x='gene', by.y='SYMBOL')
degs_clus <- split(clusters_tree_df$ENTREZID, clusters_tree_df$gene_cluster)

## GO
gos <- compareCluster(degs_clus,
                      fun = "enrichGO",
                      OrgDb = "org.Hs.eg.db",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05
)
p <- dotplot(gos, showCategory = 20)
p <- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5,
  hjust = 0.5,
  size = 6),
  axis.text.y = element_text(
    angle = 15,
    size = 6),) + scale_size(range = c(0, 3))
ggsave(glue::glue("pseudotime_deg_k5_go.pdf"),
       p, width = 6, height = 8)
readr::write_rds(gos,"pseudotime_deg_k5_gos.rds")
readr::write_csv(as.data.frame(gos), glue::glue("pseudotime_deg_k5_gos.csv"))

## KEGG
keggs <- compareCluster(degs_clus,
                        fun = "enrichKEGG",
                        organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(keggs)
p <- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5,
  hjust = 0.5,
  size = 6),
  axis.text.y = element_text(
    angle = 15,
    size = 6),) + scale_size(range = c(0, 3))
ggsave(glue::glue("pseudotime_deg_k5_kegg.pdf"),
       p, width = 5, height = 5)
readr::write_rds(keggs, glue::glue("pseudotime_deg_k5_keggs.rds"))
readr::write_csv(as.data.frame(keggs), glue::glue("pseudotime_deg_k5_keggs.csv"))


# branch ------------------------------------------------------------------
#trace('buildBranchCellDataSet', edit = T, where = asNamespace("monocle"))
BEAM_res <- BEAM(ipmn_ductal_mono_deg[ordering_genes,], branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write_csv(BEAM_res, "BEAM_res.csv")
readr::write_rds(BEAM_res, "BEAM_res.rds")

pdf("plot_cell_clusters_monocle_state_heatmap.pdf",width = 5,height = 7)
p <- plot_genes_branched_heatmap(ipmn_ductal_mono_deg[row.names(subset(BEAM_res,
                                                                       qval < 1e-4)),],
                                 branch_point = 1,
                                 num_clusters = 4,
                                 cores = 1,
                                 use_gene_short_name = T,
                                 show_rownames = T,
                                 return_heatmap = T)
dev.off()

readr::write_rds(p,"plot_monocle_state_brunch_heatmap_obj.rds")

#p$tree_row %>% str()
clusters_tree <- cutree(p$ph_res$tree_row, k = 4)
clusters_tree_df <- data.frame(clusters_tree)
clusters_tree_df[,1] <- as.character(clusters_tree_df[,1])
colnames(clusters_tree_df) <- "cluster"
clusters_tree_df <- clusters_tree_df %>% as_tibble(rownames = "gene")
#table(clusters_tree_df)
readr::write_csv(clusters_tree_df, "brunch1_clusters_tree_df.csv")

ids <- bitr(clusters_tree_df$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
clusters_tree_df <- merge(clusters_tree_df, ids, by.x='gene', by.y='SYMBOL')
degs_clus <- split(clusters_tree_df$ENTREZID, clusters_tree_df$cluster)

## GO
gos <- compareCluster(degs_clus,
                      fun = "enrichGO",
                      OrgDb = "org.Hs.eg.db",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05
)
p <- dotplot(gos, showCategory = 20)
p <- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5,
  hjust = 0.5,
  size = 6),
  axis.text.y = element_text(
    angle = 15,
    size = 6),) + scale_size(range = c(0, 3))
ggsave(glue::glue("branch1_deg_k4_go.pdf"),
       p, width = 6, height = 8)
readr::write_rds(gos,"branch1_deg_k4_gos.rds")
readr::write_csv(as.data.frame(gos), glue::glue("branch1_deg_k4_gos.csv"))

## KEGG
keggs <- compareCluster(degs_clus,
                        fun = "enrichKEGG",
                        organism = "hsa", pvalueCutoff = 0.05
)
p <- dotplot(keggs)
p <- p + theme(axis.text.x = element_text(
  angle = 45,
  vjust = 0.5,
  hjust = 0.5,
  size = 6),
  axis.text.y = element_text(
    angle = 15,
    size = 6),) + scale_size(range = c(0, 3))
ggsave(glue::glue("branch1_deg_k4_kegg.pdf"),
       p, width = 5, height = 5)
readr::write_rds(keggs, glue::glue("branch1_deg_k4_keggs.rds"))
readr::write_csv(as.data.frame(keggs), glue::glue("branch1_deg_k4_keggs.csv"))
