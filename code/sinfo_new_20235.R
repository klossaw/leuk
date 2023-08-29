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



# read in metadata
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
meta <- rds %>% colData() %>% as.data.frame()



# read in newest sampleinfo
sinfo1 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_leu_jh3.xlsx", guess_max = 10000)
sinfo1 <- as.data.frame(sinfo1)
sinfo1$sid <- paste0(sinfo1$datasets, "@", sinfo1$sample_id)


ids <- intersect(meta$sample_id[meta$age == "" | meta$gender == ""], sinfo1$sid)
ind1 <- which(meta$sample_id %in% ids)
ind1x <- match(ids, sinfo1$sid)

meta[ind1, "gender"] <- tolower(sinfo1$gender[ind1x])
meta[ind1, "age"] <- sinfo1$age[ind1x]
to_join <- unique(sinfo1[, c("sid", "vital_status", "os_time")])
colnames(to_join)[1] <- "sample_id"
meta2 <- left_join(meta, to_join)

colData(rds)$age <- meta2$age
colData(rds)$gender <- meta2$gender
colData(rds)$vital_status <- meta2$vital_status
colData(rds)$os_time <- meta2$os_time


# add new fusion_genes_rna called by rnaseq
fusion_rara <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/pml_rara.rds")
meta2 <- left_join(meta, fusion_rara)
colData(rds)$fusion_genes_rna <- meta2$fusion_genes_rna


# output
readr::write_rds(rds, rds_fn)
