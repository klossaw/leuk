
# - after sample_clustering.R


# 
rarb_temp <- fusions_rar[fusions_rar$gene3 %in% "RARB", ]$sample_id
rarg_more <- setdiff(fusions_rar[fusions_rar$gene3 %in% "RARG", ]$sample_id, rarg_samples)


# set annotation for first filtering
meta_plot$call <- "AML"
meta_plot$call[rownames(meta_plot) %in% pml_rara1] <- "APL"
meta_plot$age <- as.numeric(meta_plot$age)

# set colors
col_fun <- circlize::colorRamp2(c(0, 10, 100), c("#FFEEEE", "#FFBBBB", "#FF0000"))
anno_col <- c("RARG" = "yellow", "RARB" = "#FF8534", "TTMV-RARA" = "#E25385",
              "pml-others" = "#63E278", "RARA-others" = "lightgreen",
              "APL" = "#8864B6", "APL-like" = "steelblue", "AML" = "#96B2F3")
anno_def <- ggsci::pal_nejm("default")(8) 

# set annotation
ha2 <- HeatmapAnnotation(
  df = as.data.frame(meta_plot[, c(1:5, 8)]),
  show_legend = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  col = list(age = col_fun, 
             call = anno_col,
             gender = setNames(anno_def[1:2], c("female", "male")))
)

# plot using the same seed
set.seed(123)
h <- ComplexHeatmap::Heatmap(dat_mat[1:1000, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = F, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = T,
                             column_title = "aml_heatmap",
                             column_names_side = "top",
                             column_km = 26, border = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha2)

ht <- draw(h)
pdf("/cluster/home/yjliu_jh/projects/leukemia/output/cluster_big_exp1000_2.pdf", width = 16, height = 9)
ht
dev.off()


# check for samples
cl <- map_int(column_order(ht), length)
ci <- NULL
for (i in 1:length(cl)) {
  ci <- c(ci, paste0(rep("G", cl[i]), i))
}
ordered_samples <- data.frame(samples = colnames(dat_mat)[unlist(column_order(ht))], 
                              subgroup = ci)
meta_1 <- meta_plot
meta_1$samples <- rownames(meta_1)
ordered_samples <- left_join(ordered_samples, meta_1[, c("samples", "anno", "call")])

#table(ordered_samples[, c("subgroup", "anno")])
outlier_apl_samples <- ordered_samples$samples[ordered_samples$call %in% "APL" & ordered_samples$subgroup %notin% paste0("G", 22:25)]
aapl_samples <- ordered_samples$samples[ordered_samples$call %in% "AML" & ordered_samples$subgroup %in% paste0("G", 22:25)]
aapl_rarg <- intersect(aapl_samples, rarg_samples)


# update annotation
meta_plot$call <- "AML"
meta_plot$call[rownames(meta_plot) %in% aapl_samples] <- "APL-like"
meta_plot$call[rownames(meta_plot) %in% aapl_rarg] <- "RARG"
meta_plot$call[rownames(meta_plot) %in% ttmv_samples] <- "TTMV-RARA"



# 
small_dat <- dat[, rownames(meta_plot)[meta_plot$datasets %in% c("suzhou", "beataml") | meta_plot$ident != "AML"]]
vadj <- matrixStats::rowVars(as.matrix(small_dat))
dat_exp_h <- cbind(small_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_small <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_small <- meta_plot[colnames(dat_mat_small), ]


smaller_dat <- dat[, rownames(meta_plot)[meta_plot$call %notin% c("AML")]]
vadj <- matrixStats::rowVars(as.matrix(smaller_dat))
dat_exp_h <- cbind(smaller_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_smaller <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_smaller <- meta_plot[colnames(dat_mat_smaller), ]




ha3 <- HeatmapAnnotation(
  df = as.data.frame(meta_small[, c(1:5, 8)]),
  show_legend = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  col = list(age = col_fun, 
             call = anno_col,
             gender = setNames(anno_def[1:2], c("female", "male")))
)

set.seed(123)
h <- ComplexHeatmap::Heatmap(dat_mat_small[1:1000, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = F, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = T,
                             column_title = "aml_heatmap",
                             column_names_side = "top",
                             column_km = 12, border = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha3)


ht <- draw(h)
pdf("/cluster/home/yjliu_jh/projects/leukemia/output/cluster_small_exp1000_1.pdf", width = 16, height = 9)
ht
dev.off()




ha4 <- HeatmapAnnotation(
  df = as.data.frame(meta_smaller[, c(1:5, 8)]),
  show_legend = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  col = list(age = col_fun, 
             call = anno_col,
             gender = setNames(anno_def[1:2], c("female", "male")))
)

set.seed(123)
h <- ComplexHeatmap::Heatmap(dat_mat_smaller[1:1500, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = F, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = T,
                             column_title = "aml_heatmap",
                             column_names_side = "top",
                             column_km = 9, border = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha4)


ht <- draw(h)
pdf("/cluster/home/yjliu_jh/projects/leukemia/output/cluster_smaller_exp1500_1.pdf", width = 16, height = 9)
ht
dev.off()



smallest_dat <- dat[, rownames(meta_plot)[meta_plot$call %notin% c("AML", "APL")]]
vadj <- matrixStats::rowVars(as.matrix(smallest_dat))
dat_exp_h <- cbind(smallest_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_smallest <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_smallest <- meta_plot[colnames(dat_mat_smallest), ]



ha5 <- HeatmapAnnotation(
  df = as.data.frame(meta_smallest[, c(1:5, 8)]),
  show_legend = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
  col = list(age = col_fun, 
             call = anno_col,
             gender = setNames(anno_def[1:2], c("female", "male")))
)

set.seed(123)
h <- ComplexHeatmap::Heatmap(dat_mat_smallest[1:1000, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = T,
                             show_column_names = F, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = T,
                             column_title = "aml_heatmap",
                             column_names_side = "top",
                             column_km = 5, border = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha5)


ht <- draw(h)
pdf("/cluster/home/yjliu_jh/projects/leukemia/output/cluster_smallest_exp1000_1.pdf", width = 16, height = 9)
ht
dev.off()









