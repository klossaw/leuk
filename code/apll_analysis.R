pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# get apll annotations from clinical
apll_s <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/apll.xlsx", col_names = F)
colnames(apll_s) <- c("patient_name", "sample_seq_id")
apll_s$sample_seq_id <- sub("-", "_", apll_s$sample_seq_id)
apll_s$id_short <- sub("\\_.*", "", apll_s$sample_seq_id)
apll_s$id_short[grepl("^[[:digit:]]+", apll_s$id_short)] <- paste0("s", apll_s$sample_seq_id[grepl("^[[:digit:]]+", apll_s$id_short)])
apll_s$id_short[apll_s$id_short %in% "s202213312"] <- apll_s$sample_seq_id[apll_s$id_short %in% "s202213312"]
# match data to our sample_id
toj_1 <- cld[cld$datasets %in% c("aplsz"), c("patient_id", "sample_id")]
toj_1$patient_id <- sub("aplsz@", "", toj_1$sample_id)
toj_1$id_short <- sub("\\_.*", "", toj_1$patient_id)
toj_1$id_short[toj_1$id_short %in% "s202213312"] <- toj_1$patient_id[toj_1$id_short %in% "s202213312"]
# get apl-like samples by phenotype
apll_s2 <- left_join(apll_s, toj_1)
apll_s2$sample_seq_id[is.na(apll_s2$sample_id)]
apll_samples <- na.omit(apll_s2$sample_id)

# get other samples
# ! no need to read because we already had larger annotations 
# apll_m <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/apll.xlsx", sheet = 2, col_names = F)
# apll_m <- apll_m[1:8, 1:3]
# colnames(apll_m) <- c("patient_name", "sample_seq_id", "fusion")
# apll_m$id_short <- sub("\\-.*", "", apll_m$sample_seq_id)
# apll_m2 <- left_join(apll_m, toj_1)
# apll_m2$sample_seq_id[is.na(apll_m2$sample_id)]
# apll_rar_samples1 <- na.omit(apll_m2$sample_id)
# apll_rar_samples2 <- rownames(meta_small)[grepl("rarg", rownames(meta_small))]
# apll_rar_samples3 <- fusions_rar[!grepl("PML", fusions_rar$fusion) & fusions_rar$dataset %in% "aplsz", ]$sample_id
# apll_rar_samples <- unique(c(apll_rar_samples1, apll_rar_samples2, apll_rar_samples3))

# combine with newly added samples
new_samples <- c(paste0("B00", 1:7), "B009", "B010")
apll_samples <- c(apll_samples, new_samples)
# output
readr::write_rds(apll_samples, "/cluster/home/yjliu_jh/projects/leukemia/data/anno/apll_samples_clinical.rds")


# modify meta_plot from other scripts





small_dat <- dat[, rownames(meta_plot)[meta_plot$datasets %in% c("suzhou", "beataml") | meta_plot$ident != "AML" | 
                                         rownames(meta_plot) %in% c(apll_samples, apll_rar_samples)]]
vadj <- matrixStats::rowVars(as.matrix(small_dat))
dat_exp_h <- cbind(small_dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat_small <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 
meta_small <- meta_plot[colnames(dat_mat_small), ]


meta_small$pheno <- "AML"
meta_small$pheno[rownames(meta_small) %in% pml_rara1] <- "APL"
meta_small$pheno[rownames(meta_small) %in% apll_samples] <- "APL-like"
meta_small$pheno[rownames(meta_small) %in% apll_rar_samples] <- "RARA-others"












# continue add annotations of other RARA  also RARG and RARB 






# set annotation
ha_sm <- HeatmapAnnotation(
  df = as.data.frame(meta_small[, c("gender", "age", "datasets", "pheno", "call", "anno")]),
  show_legend = c(TRUE, TRUE, FALSE, T, T, T),
  col = list(age = col_fun, 
             gender = setNames(anno_def[1:2], c("female", "male")),
             pheno = anno_col, call = anno_col, anno = anno_col)
)



# plot using the quick function
out_path <- "/cluster/home/yjliu_jh/projects/temp/testapllheatmap2.pdf"
test_heat(dat_mat_small, ha_sm, out_path)



metasmall2 <- meta_small[meta_small$pheno %notin% "RARA-others", ]
dat_mat_small2 <- dat_mat_small[, colnames(dat_mat_small) %in% rownames(metasmall2)]
ha_sm2 <- HeatmapAnnotation(
  df = as.data.frame(metasmall2[, c("gender", "age", "datasets", "pheno", "call", "anno")]),
  show_legend = c(TRUE, TRUE, FALSE, T, T, T),
  col = list(age = col_fun, 
             gender = setNames(anno_def[1:2], c("female", "male")),
             pheno = anno_col, call = anno_col, anno = anno_col)
)



metasmall3 <- metasmall2[metasmall2$dataset %in% c("suzhou", "aplsz"), ]
dat_mat_small3 <- dat_mat_small2[, colnames(dat_mat_small2) %in% rownames(metasmall3)]
ha_sm3 <- HeatmapAnnotation(
  df = as.data.frame(metasmall3[, c("gender", "age", "datasets", "pheno", "call", "anno")]),
  show_legend = c(TRUE, TRUE, FALSE, T, T, T),
  col = list(age = col_fun, 
             gender = setNames(anno_def[1:2], c("female", "male")),
             pheno = anno_col, call = anno_col, anno = anno_col)
)
out_path <- "/cluster/home/yjliu_jh/projects/temp/testapllheatmap4.pdf"
apl_hts <- test_heat2(dat_mat_small3, ha_sm3, out_path, k = 7)


cl <- map_int(column_order(apl_hts), length)
ci <- NULL
for (i in 1:length(cl)) {
  ci <- c(ci, paste0(rep("G", cl[i]), i))
}
output_info <- data.frame(sample_id = colnames(dat_mat_small3)[unlist(column_order(apl_hts))], 
                          sub_groups = ci)
hzb <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/data/suzhou/样本汇总表.xlsx")
colnames(hzb)[13] <- "curator"
hzb3 <- as.data.frame(hzb[, 1:3])
output_info$SampleID <- sub(".*\\_", "", output_info$sample_id)
output_info <- left_join(output_info, hzb3)
output_info$SampleID[grepl("aplsz", output_info$sample_id)] <- output_info$sample_id[grepl("aplsz", output_info$sample_id)]

list2excel(list(sampleinfo = output_info[, -1]), glue("/cluster/home/yjliu_jh/projects/temp/clustering_small_group.xlsx"))





# plot using the quick function
out_path <- "/cluster/home/yjliu_jh/projects/temp/testapllheatmap3.pdf"
apl_ht <- test_heat2(dat_mat_small2, ha_sm2, out_path, k = 18)













# get count
leu_path <- "/cluster/home/jhuang/projects/leukemia/analysis"
datasets <- setdiff(list.dirs(leu_path, F, F), c("diagnosis", "bak", "meta", "phs000922")) ## unique(cld$datasets)
hugo_anno <- readr::read_delim("/cluster/home/yjliu_jh/projects/leukemia/data/hgnc_complete_set_2023-07-01.txt",
                               col_types = cols(intermediate_filament_db = col_character()))
hugo_anno <- hugo_anno[, c("symbol", "locus_group", "ensembl_gene_id", "entrez_id")]
colnames(hugo_anno)[3] <- "gene_id"
# build the count matrix
counts_list <- list()
for (dataset in datasets){
  counts <- read_csv(glue::glue("{leu_path}/{dataset}/human/rnaseq/exp/tables/{dataset}_human_counts.csv"))
  counts_fil <- counts[, -2] %>% left_join(hugo_anno) %>% na.omit() %>% 
    dplyr::select(-c(gene_id, locus_group, entrez_id)) %>% as.data.frame() %>%
    remove_rownames() %>% column_to_rownames(var = "symbol") %>% round() %>% as.matrix()
  colnames(counts_fil) <- paste0(dataset, "@", colnames(counts_fil))
  counts_list[[dataset]] <- counts_fil
}
all_counts <- as.data.frame(counts_list)
all_counts <- do.call(cbind, counts_list)
readr::write_rds(all_counts, "/cluster/home/yjliu_jh/projects/leukemia/data/counts8608.rds")
# filter based on different criteria 
sub_counts <- all_counts[, rownames(meta_small)]





# function to order and add columns to the results
add_col <- function(data){
  data$gene_name <- rownames(data)
  data$absFC <- abs(data$log2FoldChange)
  data <- data[order(data$pvalue), ]
  data2 <- left_join(as.data.frame(data), hugo_anno[, c("symbol", "entrez_id")], by = c("gene_name" = "symbol"))
  data$entrez_id <- data2$entrez_id
  data
}

# check functions
check_ego_gene <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego <- ego[order(ego$p.adjust), ]
  ego
}

check_ego_gene2 <- function(data){
  ego <- enrichGO(gene = data, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  ego
}




# auto-de functions (no need for QC in current dataset)
auto_de <- function(data, anno, compare){
  dds <- DESeq(DESeqDataSetFromMatrix(data, anno, ~ group))
  as.data.frame(add_col(results(dds, c("group", compare))))
}

meta_small$group <- meta_small$pheno
test <- auto_de(sub_counts, meta_small, c("APL-like", "APL"))
test2 <- auto_de(sub_counts, meta_small, c("RARA-others", "APL"))

apll_up <- rownames(test[test$log2FoldChange > 1, ])[1:200]
apll_down <- rownames(test[test$log2FoldChange < -1, ])[1:200]
rara_up <- rownames(test2[test2$log2FoldChange > 1, ])[1:200]
rara_down <- rownames(test2[test2$log2FoldChange < -1, ])[1:200]


temp <- dotplot(enrichGO(gene = rara_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
         ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2))

pdf(glue("/cluster/home/yjliu_jh/projects/temp/rara_vs_apl_up.pdf"), height = 8, width = 6)
print(temp)
dev.off()



