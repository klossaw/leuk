# load packages
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "SummarizedExperiment")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# read data
rds_fn <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
meta_rds <- metadata(rds)
cld <- colData(rds)

# select incorrect-labelled samples
rarg_samples <- setdiff(rownames(cld)[grepl("rarg", rownames(cld))], rownames(cld)[grepl("rarg@rarg", rownames(cld))])
rarg_ind <- match(rarg_samples, rownames(cld))

# manually change coldata / metadata / assay of affected samples
# metadata should be updated if still needed!
rownames(colData(rds))[rarg_ind] <- sub("rarg", "aplsz", rownames(colData(rds))[rarg_ind])
colData(rds)[rarg_ind, "datasets"] <- "aplsz"
colData(rds)[rarg_ind, "sample_id"] <- sub("rarg", "aplsz", colData(rds)[rarg_ind, "sample_id"])
colnames(assay(rds))[rarg_ind] <- sub("rarg", "aplsz", colnames(assay(rds))[rarg_ind])
readr::write_rds(rds, rds_fn)


