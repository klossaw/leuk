pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# 1 get all file names under the folder
mut_path <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/mutations"
files = list.files(path = mut_path, pattern = "gz$", recursive = T, full.names = F)
suffix <- c("_fil.HaplotypeCaller.maf.gz", "_fil.snp.varscan2.maf.gz",
            "_fil.speedseq.maf.gz", ".indel.varscan2.maf.gz", ".rnaindel.maf.gz")
temp <- str_split(files, "\\/") %>% map_chr(~ paste0(.[c(1:3, 3)], collapse = "/")) %>% unique()
temp <- temp[grepl("maf/", temp)]
files_all <- paste0(rep(temp, each = 5), suffix)
not_generated <- sub(".maf.gz", "", setdiff(files_all, files))
curr_ng_list <- str_split(basename(not_generated), "\\.", n = 2, simplify = T)
clist1 <- data.frame(datasets = sub("\\/.*", "", not_generated), 
                     patient_id = curr_ng_list[, 1], tool = curr_ng_list[, 2])
# modify special data names
special_ind <- grepl("XX.", not_generated)
clist1[special_ind, 2] <- paste0(clist1[special_ind, 2], ".", clist1[special_ind, 3])
clist1[special_ind, 2] <- sub("^(.*)[.].*", "\\1", clist1[special_ind, 2])
clist1[special_ind, 3] <- gsub("^.{0,2}", "", clist1[special_ind, 3])
list2excel(list(not_completed = clist1[!grepl("_fil", clist1$patient_id), ], 
                fil_not_completed = clist1[grepl("_fil", clist1$patient_id), ]),
           "/cluster/home/yjliu_jh/projects/leukemia/output/curr0716_1.xlsx")


# 2 get all file names with fastq
leu_data_path <- "/cluster/home/jhuang/projects/leukemia/data"
leu_datasets <- setdiff(list.dirs(leu_data_path, F, F), c("diagnosis", "bak", "meta", "cell_lines")) 
data_dirs <- paste0(leu_data_path, "/", leu_datasets, "/human/rnaseq")
data_dirs <- data_dirs[file.exists(data_dirs)]
pts <- "_1.fastq.gz|_2.fastq.gz|_R1.fq.gz|_R2.fq.gz|_2.fq.gz|_1.fq.gz|.fq.gz"
datainfo <- function (dir) {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[8], 
             patient_id = sub(pts, "", list.files(dir, pattern = ".fq.gz|.fastq.gz")))
}
curr_all_data <- unique(bind_rows(lapply(data_dirs, datainfo)))
curr_all_data_processed <- as.data.frame(str_split(temp, "\\/", simplify = T))[, c(1, 3)]
colnames(curr_all_data_processed) <- c("datasets", "patient_id")
remaining_data <- anti_join(curr_all_data, curr_all_data_processed)
list2excel(list(not_started = remaining_data),
           "/cluster/home/yjliu_jh/projects/leukemia/output/curr0716_2.xlsx")

# 3 delete extra files
to_remove <- setdiff(files, files_all)[grepl("fil_", setdiff(files, files_all))]
to_remove <- paste0(mut_path, "/", to_remove)
file.remove(to_remove)


