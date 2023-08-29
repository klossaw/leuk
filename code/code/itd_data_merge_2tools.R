pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq",
          "jsonlite")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

# ======== itd data from filt3r ================
leu_ana_path <- "~/projects/leukemia/analysis"
leu_datasets <- base::setdiff(list.dirs(leu_ana_path, F, F), c("diagnosis", "bak", "meta", "cell_lines"))
data_dirs <- paste0(leu_ana_path, "/", leu_datasets, "/human/rnaseq/splicing/filt3r/", leu_datasets)
data_dirs <- data_dirs[file.exists(data_dirs)]
datainfo <- function (dir) {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[5],
             patient_id = sub(".json", "", list.files(dir, pattern = ".json")))
}
curr_all_data <- unique(bind_rows(lapply(data_dirs, datainfo)))

test_list <- list()
for (i in 1:nrow(curr_all_data)) {
  test3 <- readLines(glue("{leu_ana_path}/{curr_all_data[i, 1]}/human/rnaseq/splicing/filt3r/{curr_all_data[i, 1]}/{curr_all_data[i, 2]}.json"), warn = F)
  test2 <- sub("NaN|nan", 999, test3) %>% fromJSON()
  test1 <- test2$details
  test1[["event_reads"]] <- test2$nb_reads_in_events
  test1[["percentage"]] <- test2$percentage
  test1[["sample_id"]] <- glue("{curr_all_data[i, 1]}@{curr_all_data[i, 2]}")
  test1[["dataset"]] <- curr_all_data[i, 1]
  test_list[[i]] <- test1
}

flt3_itd_all <- bind_rows(test_list)
flt3_itd_fil <- flt3_itd_all[flt3_itd_all$is_duplication == "TRUE", ]

output_dir <- "~/projects/leukemia/analysis/meta/human/rnaseq/splicing/filt3r/tables" ## change output dir here
readr::write_tsv(flt3_itd_fil, glue("{output_dir}/filt3r.tsv"))



# ======= find samples without filt3r results ==========
leu_data_path <- "~/projects/leukemia/data"
data_dirs <- paste0(leu_data_path, "/", leu_datasets, "/human/rnaseq")
data_dirs <- data_dirs[file.exists(data_dirs)]
pts <- "_1.fastq.gz|_2.fastq.gz|_R1.fq.gz|_R2.fq.gz|_2.fq.gz|_1.fq.gz|.fq.gz"
datainfo <- function (dir) {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[5],
             patient_id = sub(pts, "", list.files(dir, pattern = ".fq.gz|.fastq.gz")))
}
curr_all_data_fq <- unique(bind_rows(lapply(data_dirs, datainfo)))
anti_join(curr_all_data_fq, curr_all_data)



# ======== itd data from hamlet_itd ================
data_dirs2 <- paste0(leu_ana_path, "/", leu_datasets, "/human/rnaseq/splicing/hamlet_itd/", leu_datasets)
data_dirs2 <- data_dirs2[file.exists(data_dirs2)]
datainfo_new <- function (dir, x = "_FLT3.csv") {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[8],
             path = list.files(list.dirs(dir, full.names = T, recursive = F), pattern = x, full.names = T))
}
curr_all_flt3 <- unique(bind_rows(lapply(data_dirs2, datainfo_new)))
curr_all_kmt2a <- unique(bind_rows(lapply(data_dirs2, datainfo_new, x = "_KMT2A.csv")))
all_flt3 <- bind_rows(lapply(curr_all_flt3$path, read_delim, col_types = "ddddddddcd"))
all_kmt2a <- bind_rows(lapply(curr_all_kmt2a$path, read_delim, col_types = "ddddddddcd"))


output_dir <- "~/projects/leukemia/analysis/meta/human/rnaseq/splicing/hamlet_itd/tables" 
readr::write_tsv(all_flt3, glue("{output_dir}/flt3_hamlet.tsv"))
readr::write_tsv(all_kmt2a, glue("{output_dir}/kmt2a_hamlet.tsv"))


# ======= find samples without hamlet_itd results ==========
curr_all_flt3$patient_id <- str_split(curr_all_flt3$path, pattern = "/", simplify = T)[, 14]
curr_all_kmt2a$patient_id <- str_split(curr_all_kmt2a$path, pattern = "/", simplify = T)[, 14]
anti_join(curr_all_data_fq, curr_all_flt3[, c(1, 3)])
anti_join(curr_all_data_fq, curr_all_kmt2a[, c(1, 3)])