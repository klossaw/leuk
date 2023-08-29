pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


# 1 list all files
mut_path <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/mutations"
files = list.files(path = mut_path, pattern = "gz$", recursive = T, full.names = F)

# 2 generate all correct file names
suffix <- c("_fil.HaplotypeCaller.maf.gz", "_fil.snp.varscan2.maf.gz",
            "_fil.speedseq.maf.gz", ".indel.varscan2.maf.gz", ".rnaindel.maf.gz")
temp <- str_split(files, "\\/") %>% map_chr(~ paste0(.[c(1:3, 3)], collapse = "/")) %>% unique()
temp <- temp[grepl("maf/", temp)]
files_all <- paste0(rep(temp, each = 5), suffix)

# 3 locate and delete extra files with "fil_" pattern
to_remove <- setdiff(files, files_all)[grepl("fil_", setdiff(files, files_all))]
to_remove <- paste0(mut_path, "/", to_remove)
file.remove(to_remove)

