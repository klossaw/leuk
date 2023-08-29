pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "SummarizedExperiment", "jhuanglabRNAseq", "ComplexHeatmap")
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
meta_part <- left_join(as.data.frame(metadata), clusterings)

# read prepared fusion data
fusions_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/fusions_new.rds")

# update rds to remove history bad assignments
test <- unique(meta_part[, c("fusion_genes_rna", "rna_group")])
test <- test[test$fusion_genes_rna != "", ]

keep_groups <- c("G1", "G2", "G8", "G9", "G15", "G16", "G17", "G19")  ## keep before-annotated samples
other_groups <- unique(setdiff(test$rna_group, keep_groups))
metadata$fusion_genes_rna[meta_part$rna_group %in% other_groups] <- ""

# annotate 13 14 clusters as Boqiang Hu did not integrate his results to rds
g13_14_anno_samples <- intersect(meta_part$sample_id[meta_part$rna_group %in% c("G13", "G14")], 
                               fusions_new$sample_id[fusions_new$fusion_anno %in% c("RUNX1-RUNX1T1", "RUNX1T1-RUNX1")])
metadata$fusion_genes_rna[metadata$sample_id %in% g13_14_anno_samples] <- "aml_RUNX1_RUNX1T1"


# 
remaining_groups <- setdiff(unique(meta_part$rna_group), c(keep_groups, "G13", "G14"))





# 





# check for known fusions?

fusions_new_aml <- inner_join(fusions_new, meta_part[, c("sample_id", "rna_group")])
table(fusions_new_aml$rna_group[fusions_new_aml$fusion %in% "DEK-NUP214"])  ## may have data loss
table(fusions_new_aml$rna_group[fusions_new_aml$fusion %in% "ZNF292-PNRC1"])

table(fusions_new_aml$rna_group[fusions_new_aml$fusion %in% "RBM15-MRTFA"])

# first check for samples without canonical fusions assigned to that group







# identical(meta_plot$sample_id, meta_part$sample_id)
# table(meta_plot$fusion_plot[meta_part$rna_group %in% c("G16")])

g13_14_left_samples <- metadata$sample_id[metadata$fusion_genes_rna %notin% "aml_RUNX1_RUNX1T1" & meta_part$rna_group %in% c("G13", "G14")]
temp_13_14 <- fusions_new[fusions_new$sample_id %in% g13_14_left_samples, ]
g8_9_left_samples <- metadata$sample_id[metadata$fusion_genes_rna %notin% "aml_KMT2A-r" & meta_part$rna_group %in% c("G8", "G9")]
temp_8_9 <- fusions_new[fusions_new$sample_id %in% g8_9_left_samples, ]

temp_8_9 <- temp_8_9[temp_8_9$fusion %notin% blacklist_fusions, ]
g_8_9 <- temp_8_9 %>% group_by(sample_id, gene5) %>% summarise(count5 = n())
g_8_9 <- g_8_9 %>% group_by(gene5) %>% summarise(count = n())
g_8_9 <- g_8_9[order(g_8_9$count, decreasing = T), ]


blacklist_fusions <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/blacklist_fusions.rds")

head(sort(table(temp_13_14$fusion[temp_13_14$fusion %notin% blacklist_fusions]), decreasing = T))  ## no fusion patterns
head(sort(table(temp_8_9$fusion[temp_8_9$fusion %notin% blacklist_fusions]), decreasing = T))  ## some non-exclusive fusion patterns

head(sort(table(temp_8_9$gene5[temp_8_9$fusion %notin% blacklist_fusions]), decreasing = T))  ## MALAT1 500!
head(sort(table(temp_13_14$gene5[temp_13_14$fusion %notin% blacklist_fusions]), decreasing = T))  ## max OAZ1 3 B2M 2

# table 





fusion_fc <- read.delim("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/fusion/tables/leukemia_fusioncatcher.tsv")


