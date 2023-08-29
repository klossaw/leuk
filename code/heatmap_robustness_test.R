pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# after other scripts
# sth else

# 
#/cluster/home/jhuang/projects/leukemia/data/aplsz/human/rnaseq/rename.R

big_path <- "/cluster/home/jhuang/projects/leukemia/data/suzhou/human/rnaseq"
objs_big <- file.info(list.files(big_path, full.names = T))
objs_big <- objs_big[!grepl("bak", rownames(objs_big)), ]
objs_big <- objs_big[!grepl("dup", rownames(objs_big)), ]
rownames(objs_big) <- sub(big_path, "", rownames(objs_big))


hzb <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/data/suzhou/样本汇总表.xlsx")
ztt <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/data/suzhou/张彤彤/张彤彤112例RNA-SEQ.xls")
colnames(hzb)[13] <- "curator"
colnames(ztt)[3] <- "sample_id"



small_path <- "/cluster/home/jhuang/projects/leukemia/data/aplsz/human/rnaseq"













# load data matrix and annotation
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
cld <- colData(rds) |> as.data.frame() %>% replace(., is.na(.), "")
fil <- cld$disease_type  %in% c("AML")
se <- rds[,fil]
metadata <- colData(se)
fpkm <- assay(se, "tpm_remove_rnatype") %>% as.data.frame() %>% dplyr::select(metadata$sample_id) # exp matrix
dat <- fpkm |> as.data.frame() |> tibble::rownames_to_column("gene_name") |>
  filter_black_heatmap_genes() |> tibble::remove_rownames() |> tibble::column_to_rownames("gene_name")


# apl_anno1
apl_anno <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/112RNA-SEQ.xls")
colnames(apl_anno) <- c("id", "sample_id", "name", "age", "gender", "ident", "fab_subtype")


# apl_anno2


# ttmv_data
ttmv_slim <- unique(ttmv_data[, c("gene", "dataset", "sample_id")])
ttmv_tmp <- ttmv_slim[grepl("FHIT", ttmv_slim$gene), ]
ttmv_freq <- table(ttmv_slim$gene)
ttmv_freq <- sort(ttmv_freq, decreasing = T)
ttmv_freq_all <- table(ttmv_data$gene)
ttmv_freq_all <- sort(ttmv_freq_all, decreasing = T)
candidates <- c("HNRNPL", names(ttmv_freq)[1:20][!grepl(".1", names(ttmv_freq)[1:20])])

# ttmv candidate list
sample_ttmv_list <- list()
sample_ttmv_list[["HNRNPL"]] <- unique(ttmv_data$sample_id[ttmv_data$gene %in% "AC008982.1/HNRNPL"])
for (i in 2:length(candidates)) {
  sample_ttmv_list[[candidates[i]]] <- unique(ttmv_data$sample_id[ttmv_data$gene %in% candidates[i]])
}


# itd data
itd1 <- readr::read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/leukemia/itd.csv")


# metadata for plot
meta_plot2 <- as.data.frame(metadata)[, c("gender", "age", "rna_type", "disease_type", "datasets")]
meta_plot2$sample_id <- rownames(meta_plot2)
meta_plot2$age <- as.numeric(meta_plot2$age)








# ======= set functions ========

all_colors <- c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", 
                "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
                "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3",
                "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB",
                "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", "#008B8B",
                "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4",
                "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347",
                "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887")
col_fun <- circlize::colorRamp2(c(0, 10, 100), c("#FFEEEE", "#FFBBBB", "#FF0000"))
anno_col <- c("RARG" = "yellow", "RARB" = "#FF8534", "TTMV-RARA" = "#E25385",
              "pml-others" = "#63E278", "RARA-others" = "lightgreen",
               "APL" = "#8864B6", "APL-like" = "steelblue", "AML" = "#96B2F3")
anno_def <- ggsci::pal_nejm("default")(8)

# intakes a annotation data, a sample containing data, and columns to be added to annotation
enrich_anno <- function(anno, data, extra, sample_col = "sample_id") {
  for (col in extra){
    anno[[col]] <- ifelse(anno[[sample_col]] %in% data[[col]], "exist", "non-exist")
  }
  anno
}


auto_anno <- function(data, cols, colors) {
  anno_list <- list()
  for (col in cols) {
    elements <- sort(unique(data[[col]]))
    anno_list[[col]] <- setNames(colors[1:length(elements)], elements)
  }
  anno_list
}


meta_plot3 <- enrich_anno(meta_plot2, sample_ttmv_list, candidates)


test_heat <- function(data, anno, path, k = 26, top = 1000, title = "aml_heatmap",
                      seed = 2023, show_col = F, show_row = F, raster = T, bd = F) {
  set.seed(seed)
  h <- ComplexHeatmap::Heatmap(data[1:top, ], 
                               name = " ",
                               show_row_names = show_row, 
                               use_raster = raster,
                               show_column_names = show_col, 
                               clustering_method_columns = 'ward.D2',
                               clustering_method_rows = 'ward.D2',
                               show_row_dend = F,
                               show_column_dend = T,
                               column_title = "aml_heatmap",
                               column_names_side = "top",
                               column_km = k, border = bd,
                               heatmap_legend_param = list(
                                 legend_direction = "horizontal", 
                                 legend_width = unit(3, "cm")),
                               top_annotation = anno)
  
  ht <- draw(h)
  pdf(path, width = 16, height = 9)
  print(ht)
  dev.off()
}



# ======== function block end =======


# set annotation
ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot3[, c("gender", "age", "datasets", candidates)]),
  show_legend = c(TRUE, TRUE, FALSE, rep(T, length(candidates))),
  col = append(list(age = col_fun, 
             gender = setNames(anno_def[1:2], c("female", "male"))),
             auto_anno(meta_plot3, candidates, c("magenta", "ivory")))
)

# dat
vadj <- matrixStats::rowVars(as.matrix(dat))
dat_exp_h <- cbind(dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 

# plot using the quick function
out_path <- "/cluster/home/yjliu_jh/projects/temp/testxxxxheatmap.pdf"
test_heat(dat_mat, ha_x, out_path)














# sth else


library(ngstk)
setwd("~/projects/leukemia/data/aplsz/human/rnaseq")
format_filenames(files_dir = "~/projects/leukemia/data/aplsz/human/rnaseq",
                 pattern = "*.gz",do.rename = T,
                 replace = list(old = c("_L00\\d+","-G1D2J17XKF1-H004K984_S\\d+","_001.fastq.gz","1.fq.gz", "2.fq.gz","RRR","RR"), 
                                new = c("","",".fq.gz","R1.fq.gz", "R2.fq.gz","R","R")))




# ====== real robustness function ======




robust_test <- function(, times = 10) {
  
  
  
  
}


# loop before function:
# 1: check randomness of functions
# 2: make sure sample annotations



col_fun <- circlize::colorRamp2(c(0, 10, 100), c("#FFEEEE", "#FFBBBB", "#FF0000"))
anno_col <- c("RARG" = "yellow", "RARB" = "#FF8534", "TTMV-RARA" = "#E25385",
              "pml-others" = "#63E278", "RARA-others" = "lightgreen",
              "APL" = "#8864B6", "APL-like" = "steelblue", "AML" = "#96B2F3")
anno_def <- ggsci::pal_nejm("default")(8)

meta_plot$age <- as.numeric(meta_plot$age)

small_dat <- dat[, rownames(meta_plot)[meta_plot$datasets %in% c("suzhou", "beataml") | meta_plot$ident != "AML"]]
vadj <- matrixStats::rowVars(as.matrix(small_dat))
dat_exp_h <- cbind(small_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_small <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_small <- meta_plot[colnames(dat_mat_small), ]


ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot3[, c("gender", "age", "datasets", candidates)]),
  show_legend = c(TRUE, TRUE, FALSE, rep(T, length(candidates))),
  col = append(list(age = col_fun, 
                    gender = setNames(anno_def[1:2], c("female", "male"))),
               auto_anno(meta_plot3, candidates, c("magenta", "ivory")))
)




param <- data.frame(split = rep(11:15, each = 2), topvar = c(0.95, 0.98))


for (i in 10) {
  quick_heatmap(dat_mat_small, ha = ha_small, height = 9, outdir = paste0("aml/new", i), 
                column_split = 13, top_var_percent = 0.98)
}





auto_selection <- 







# select strange genes from heatmap
heatmap_genes <- readxl::read_excel("/cluster/home/yjliu_jh/aml/new1/sampleinfo_0.95.xlsx", sheet = 2)
od_samples <- readxl::read_excel("/cluster/home/yjliu_jh/aml/new1/sampleinfo_0.95.xlsx", sheet = 1)


# finish sample selection and make a ppt

# make a dav update ppt



hx <- ComplexHeatmap::Heatmap(dat_mat_small[heatmap_genes$gene_name[2300:2800], od_samples$sample_id],
                              top_annotation = ha_small, cluster_rows = F, cluster_columns = F)

ht <- draw(hx)
pdf("/cluster/home/yjliu_jh/projects/temp/testxsssxxxheatmape4.pdf", width = 16, height = 9)
print(ht)
dev.off()



















function (dat, top_var_percent = 0.9, outdir = "tall", ha = NULL, 
          heatmap_title = "leukemia", use_raster = F, column_split = 10, 
          height = 6, width = 16, cluster_columns = TRUE) 
{
  dat_var <- apply(dat, 1, var)
  checkdir(outdir)
  message(glue("We are processing percent: {top_var_percent}"))
  index <- dat_var >= quantile(dat_var, probs = top_var_percent)
  exp_plot <- data.matrix(dat[index, ])
  hm <- jhHeatmap(as.matrix(exp_plot), glue("{outdir}/jhHeatmap_{top_var_percent}.pdf"), 
                  width = width, height = height, Colv = cluster_columns)
  pdf(glue("{outdir}/heatmap_{top_var_percent}.pdf"), height = height, 
      width = width)
  mat_scaled <- t(hm$carpet)[row.names(exp_plot), colnames(exp_plot)]
  ht <- ComplexHeatmap::Heatmap(mat_scaled, name = heatmap_title, 
                                use_raster = use_raster, col = circlize::colorRamp2(c(-1.8, 
                                                                                      0, 1.8), c("blue", "white", "red")), row_order = rev(hm$rowInd), 
                                column_order = hm$colInd, cluster_rows = rev(hm$rowDendrogram), show_row_names = FALSE, 
                                show_column_names = FALSE, top_annotation = ha, column_split = column_split, 
                                column_gap = unit(0.3, "mm"), column_title = NULL)
  ht <- ComplexHeatmap::draw(ht, annotation_legend_side = "left", 
                             heatmap_legend_side = "right")
  dev.off()
  cl <- map_int(column_order(ht), length)
  ci <- NULL
  for (i in 1:length(cl)) {
    ci <- c(ci, paste0(rep("G", cl[i]), i))
  }
  output_info <- data.frame(sample_id = colnames(mat_scaled)[unlist(column_order(ht))], 
                            sub_groups = ci)
  heatmapeGenes <- data.frame(gene_name = rownames(mat_scaled)[unlist(row_order(ht))])
  list2excel(list(sampleinfo = output_info, heatmapeGenes = heatmapeGenes), 
             glue("{outdir}/sampleinfo_{top_var_percent}.xlsx"))
}





