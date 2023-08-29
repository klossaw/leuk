pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2", "ggthemes", "EnsDb.Hsapiens.v86",
          "jhuanglabRNAseq", "clusterProfiler", "DOSE", "org.Hs.eg.db", "rstatix")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


new_samples <- c(paste0("B00", 1:7), "B009", "B010")

setwd("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/exp")
what1 <- jhuanglabRNAseq::getexpSalmon("leu_j", "apll", "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/salmon")

new_exp <- read_csv("datExpBeforeSVA.csv")

mut_apll_new <- mut_int[, colnames(mut_int) %in% c(mut_cols, "tool"), with = F]
temp <- mut_apll_new[, toString(tool), by = setdiff(names(mut_apll_new), "tool")]
colnames(temp)[5] <- "tools"
mutations_new <- rbind(mutations_new, as.data.frame(temp)[, mut_cols])



# combine exp value
new_exp <- new_exp[, -1] %>% as.data.frame() %>% na.omit() %>%
  remove_rownames() %>% column_to_rownames(var = "gene_name")
new_exp <- as.matrix(new_exp[rownames(dat), ])
dat_new <- cbind(dat, new_exp)

small_dat <- dat[, rownames(meta_plot)[meta_plot$datasets %in% c("suzhou", "beataml") | meta_plot$ident != "AML" | 
                                         rownames(meta_plot) %in% c(apll_samples, apll_rar_samples)]]
small_dat <- cbind(small_dat, new_exp)
vadj <- matrixStats::rowVars(as.matrix(small_dat))
dat_exp_h <- cbind(small_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
mat_small_new <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_small <- meta_plot[colnames(dat_mat_small), ]
meta_small$pheno <- "AML"
meta_small$pheno[rownames(meta_small) %in% pml_rara1] <- "APL"
meta_small$pheno[rownames(meta_small) %in% apll_samples] <- "APL-like"
meta_small$pheno[rownames(meta_small) %in% apll_rar_samples] <- "RARA-others"
meta_add <- data.frame(gender = NA_character_, age = NA_integer_, rna_type = "unknown",
                       disease_type = "AML", datasets = "aplsz_new", ident = "AML",
                       call = "AML", pheno = "APL-like", sample_id = new_samples)
meta_add <- meta_add %>% remove_rownames() %>% column_to_rownames(var = "sample_id")
meta_small_new <- rbind(meta_small, meta_add)



ms2 <- as.data.frame(meta_small_new)[, c("gender", "age", "rna_type", "disease_type", "datasets", "call", "pheno")]
ms2$sample_id <- rownames(ms2)
ms2$age <- as.numeric(ms2$age)
meta_plot3x <- enrich_anno(ms2, mutall_list, mut_genes)
meta_plot4x <- enrich_anno(meta_plot3x, itd_list, itd_genes)
meta_plot4x <- enrich_anno(meta_plot4x, fusion_list, fus_genes[-7])
meta_plot5 <- enrich_anno2(meta_plot4x, mut2_list, gene2, anno_col = "mutation")

ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot5[, -c(3, 4, 8)]),
  show_legend = c(TRUE, TRUE, FALSE, T, T, rep(F, length(mut_genes) + length(itd_genes) + length(fus_genes[-7])), T, T),
  col = c(list(age = col_fun, 
               call = anno_col,
               pheno = anno_col,
               gender = setNames(anno_def[1:2], c("female", "male"))),
          auto_anno(meta_plot5, c(mut_genes, itd_genes, fus_genes#, asxl1_hs
          ), c("magenta", "ivory")),
          auto_anno(meta_plot5, c(gene2, "sub_groups", "subgroups_230712"), all_colors))
)

# output plots


readr::write_rds(dat_new, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/dat_apll_new.rds")
readr::write_rds(mat_small_new, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/exp_apll_new.rds")
readr::write_rds(meta_small_new, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/meta_apll_new.rds")
readr::write_rds(ha_x, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/ha_apll_new.rds")
readr::write_rds(small_dat, "/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/exp/small_dat.rds")


quick_heatmap(small_dat, ha = ha_x, width = 25, height = 15, outdir = "plots", 
              column_split = 12, top_var_percent = 0.95)



smaller_dat <- dat_new[, rownames(meta_small_new)[meta_small_new$dataset %notin% c("beataml")]]
vadj <- matrixStats::rowVars(as.matrix(smaller_dat))
dat_exp_h <- cbind(smaller_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_smaller <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta5_smaller <- meta_plot5[colnames(dat_mat_smaller), ]

ha_sm <- HeatmapAnnotation(
  df = as.data.frame(meta_plot5[, -c(3, 4, 8)]),
  show_legend = c(TRUE, TRUE, FALSE, T, T, rep(F, length(mut_genes) + length(itd_genes) + length(fus_genes[-7])), T, T),
  col = c(list(age = col_fun, 
               call = anno_col,
               pheno = anno_col,
               gender = setNames(anno_def[1:2], c("female", "male"))),
          auto_anno(meta_plot5, c(mut_genes, itd_genes, fus_genes#, asxl1_hs
          ), c("magenta", "ivory")),
          auto_anno(meta_plot5, c(gene2, "sub_groups", "subgroups_230712"), all_colors))
)










set.seed(123)
h <- ComplexHeatmap::Heatmap(mat_small_new[1:1000, ], 
                             name = " ",
                             show_row_names = FALSE, 
                             use_raster = F,
                             show_column_names = F, 
                             clustering_method_columns = 'ward.D2',
                             clustering_method_rows = 'ward.D2',
                             show_row_dend = F,
                             show_column_dend = T,
                             column_title = "apll_heatmap",
                             column_names_side = "top",
                             column_km = 12, border = F,
                             heatmap_legend_param = list(
                               legend_direction = "horizontal", 
                               legend_width = unit(3, "cm")),
                             top_annotation = ha_x)

ht <- draw(h)
pdf("/cluster/home/yjliu_jh/projects/leukemia/output/apll_new_test01.pdf", width = 25, height = 15)
ht
dev.off()



# npm1 short duplications
mut_apll <- mutations_new[mutations_new$sample_id %in% rownames(meta_plot5)[meta_plot5$pheno %in% "APL-like"], ]
mut_apll_npm1 <- mut_apll[mut_apll$Hugo_Symbol %in% "NPM1", ]
mut_apl_npm1 <- mutations_new[mutations_new$sample_id %in% rownames(meta_plot5)[meta_plot5$pheno %in% "APL"] & mutations_new$Hugo_Symbol %in% "NPM1", ]




("/cluster/home/yjliu_jh/projects/leu_j/analysis/apll/human/rnaseq/fusion/starfusion/{sample}/star-fusion.fusion_predictions.tsv")







