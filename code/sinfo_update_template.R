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
sinfo <- readxl::read_excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_leu_jh3.xlsx", guess_max = 10000)




# add new columns
# the 'fusion' part may include
# 1 fusion annotations from original data sources (maybe mixed and needs to be categorized/standardized)
# 2 fusion annotations called by our analysis pipeline
# 3 other fusion annotations from methods including clustering, etc.
# therefore, for convenience, use separate columns for these fusions and make adjustments later






# fill blanks for existing columns










