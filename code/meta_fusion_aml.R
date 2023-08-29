pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "leukemia"
dataset <- "meta"
species <- "human"
workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq/exp/tables") |> checkdir()
setwd(workdir)


# load total leukemia dataset
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
# tpm_raw tpm_remove_rnatype tpm_remove_cohort tpm
cld <- colData(rds) |> as.data.frame() %>% replace(., is.na(.), "")
fil <- cld$remove == "no" &  cld$disease_type  %in% c("AML") & cld$cell_lines == ""
fil <- cld$disease_type  %in% c("AML")
fusion_info <- metadata(rds)$fusion_info$fusion_groups
se <- rds[, fil]
metadata <- colData(se)


# load sampleinfo 
sinfo <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/docs/meta/sampleinfo/sampleinfo_meta.xlsx", guess_max = 10000)


# load clustering results
clusterings <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
clusterings <- as.data.frame(clusterings)
colnames(clusterings)[2] <- "rna_group"

# load fusion
fusions_all <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/fusions_all.rds")
fusions_new <- fusions_all  ## change to new sample_id
fusions_new$sample_id <- paste0(fusions_new$dataset, "@", fusions_new$sample_id)

# separate the samples into groups and perform functions

meta_plot <- as.data.frame(metadata)
meta_plot$fusion_genes_rna[meta_plot$rna_group %in% c("G1", "G2") & meta_plot$fusion_genes_rna == ""] <- "RARA-others"
meta_plot$fusion_plot <- meta_plot$fusion_genes_rna

meta_part <- left_join(meta_plot, clusterings)
#meta_part <- meta_part[meta_part$sample_id %in% clusterings$sample_id[clusterings$rna_group %in% paste0("G", 15:26)], ]


# 1 top fusions appeared
find_top <- function(sample_list, fusion_dat, top = 20){
  fusions <- unique(fusion_dat[fusion_dat$sample_id %in% sample_list, c(2, 7)])
  temp <- as.data.frame(sort(table(fusions$fusion_anno), decreasing = T)[1:top])
  temp$prop <- temp$Freq / length(sample_list)
  temp
}


# 2 whitelist fusions  and corresponding proportion
fusions_new$fusion_anno <- fusions_new$fusion
fusions_new$fusion_anno[grepl("KMT2A", fusions_new$fusion_anno)] <- "KMT2A-r"
fusions_new$fusion_anno[grepl("NUP98", fusions_new$fusion_anno)] <- "NUP98-r"
fusions_new$fusion_anno[grepl("MECOM", fusions_new$fusion_anno)] <- "MECOM-r"
fusions_new$fusion_anno[grepl("DUX4", fusions_new$fusion_anno)] <- "DUX4-r"
fusions_new$fusion_anno[grepl("PAX5", fusions_new$fusion_anno)] <- "PAX5-r"
fusions_new$fusion_anno[grepl("MLLT10", fusions_new$fusion_anno)] <- "MLLT10-r"
fusions_new$fusion_anno[grepl("CRLF2", fusions_new$fusion_anno)] <- "CRLF2-r"
fusions_new$fusion_anno[grepl("LMO1", fusions_new$fusion_anno)] <- "LMO1-r"
fusions_new$fusion_anno[grepl("LMO2", fusions_new$fusion_anno)] <- "LMO2-r"
# may add more new alterations... maybe some canonical in other leukemias?



# select likely gene fusions
# groups <- unique(meta_part$rna_group)
# f_list <- list()
# for(i in 1:length(groups)){
#   f_list[[groups[i]]] <- find_top(meta_part$sample_id[meta_part$rna_group %in% groups[i]], fusions_new, 20)
# }
# f_top <- bind_rows(f_list, .id = "rna_group")
# f_top_wide <- f_top[, -3] %>% pivot_wider(names_from = "Var1", values_from = "prop", values_fill = NA)
# f_top_filtered <- f_top_wide[, c(1, which(colSums(is.na(f_top_wide)) >= 13))]
# ftf <- f_top_filtered %>% pivot_longer(-rna_group) %>% na.omit() %>% as.data.frame()
# colnames(ftf)[2:3] <- c("fusion", "prop")
# tj <- as.data.frame(table(meta_part$rna_group))
# colnames(tj) <- c("rna_group", "sample_count")
# ftf <- left_join(ftf, tj)
# ftf %>% df2excel("/cluster/home/yjliu_jh/projects/temp/aml_g1-26.xlsx")
# f_list %>% df2excel("/cluster/home/yjliu_jh/projects/temp/g1-26_all.xlsx")  ## may further investigate
# ftf2 <- ftf[ftf$prop > 0.85, ]






# alter annotation
# g15_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% "G15"], 
#                               fusions_new$sample_id[fusions_new$fusion_anno %in% "CBFB-MYH11"])
# 
# 
# meta_plot$fusion_plot[meta_plot$sample_id %in% g15_anno_samples] <- "aml_CBFB_MYH11"
# 
# g17_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% "G17"], 
#                              fusions_new$sample_id[fusions_new$fusion_anno %in% "NUP98-r"])
# 
# meta_plot$fusion_plot[meta_plot$sample_id %in% g17_anno_samples] <- "aml_NUP98_rearrangement"
# 
# g19_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% "G19"], 
#                               fusions_new$sample_id[fusions_new$fusion_anno %in% "CBFA2T3-GLIS2"])
# 
# meta_plot$fusion_plot[meta_plot$sample_id %in% g19_anno_samples] <- "aml_CBFA2T3-GLIS2"


g8_9_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% c("G8", "G9")], 
                                fusions_new$sample_id[fusions_new$fusion_anno %in% "KMT2A-r"])
meta_plot$fusion_plot[meta_plot$sample_id %in% g8_9_anno_samples] <- "aml_KMT2A-r"


g16_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% c("G16")], 
                           fusions_new$sample_id[fusions_new$fusion_anno %in% "PCAT18-KCTD1"])
meta_plot$fusion_plot[meta_plot$sample_id %in% g16_anno_samples] <- "aml_PCAT18-KCTD1"






# update to rds
indices <- match(meta_plot$sample_id, cld$sample_id)
cld$fusion_genes_rna[indices] <- meta_plot$fusion_plot
colData(rds)$fusion_genes_rna <- cld$fusion_genes_rna
readr::write_rds(rds, rds_fn)



# output xlsx
metadatap <- as.data.frame(colData(rds)) %>% dplyr::select("sample_id", "datasets", "age", "gender", 
                                                           starts_with("subgroups_"), starts_with("subtype_"), starts_with("PMID"), "fusion_genes", "mutations", "disease_type", 
                                                           "rna_type", "vital_status", "os_time", "fusion_genes_rna")

metadatap %>% df2excel("/cluster/home/jhuang/projects/leukemia/docs/meta/sampleinfo/sampleinfo_meta.xlsx")







# "NPM1", "CEBPA" mutation ...











