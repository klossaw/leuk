pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq",
          "jsonlite")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}



# ========== ttmv data ==============
ttmv_path <- "/cluster/home/yliang_jh/projects/mRNA/TTMV"
datasets <- dir_ls(ttmv_path, type = "directory") %>% path_file()
datasets <- datasets[!datasets %in% c(".snakemake")]
ttmv_list <- list()
for (dataset in datasets){
  csv <- glue("{ttmv_path}/{dataset}/{dataset}.summary.csv")
  if (file.exists(csv)){
    ttmv_list[[dataset]] <- read_csv(csv, col_types = "cddddddccccc")  ## prevent empty csv default as character
  }
}
ttmv_data <- bind_rows(ttmv_list) 
ttmv_data$sample_id <- paste0(ttmv_data$dataset, "@", ttmv_data$sample)
readr::write_rds(ttmv_data, "/cluster/home/yjliu_jh/projects/leukemia/data/ttmv_data.rds")



# ========= itd data =================
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
cld <- colData(rds) %>% as.data.frame()
cld$sample_name <- sub(".*@", "", cld$sample_id)
itd1 <- readr::read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/leukemia/itd.csv")
itd1 <- left_join(itd1, cld[, c("sample_name", "sample_id")])
itd1$itd[itd1$itd %in% "flt3"] <- "FLT3-ITD"
itd1$itd[itd1$itd %in% "kmt2a"] <- "KMT2A-PTD"
itd1 <- as.data.frame(itd1)
readr::write_rds(itd1, "/cluster/home/yjliu_jh/projects/leukemia/data/itd_data.rds")



# ======== itd data 2 ================
leu_ana_path <- "/cluster/home/jhuang/projects/leukemia/analysis"
leu_datasets <- setdiff(list.dirs(leu_ana_path, F, F), c("diagnosis", "bak", "meta", "cell_lines")) 
data_dirs <- paste0(leu_ana_path, "/", leu_datasets, "/human/rnaseq/splicing/filt3r/", leu_datasets)
data_dirs <- data_dirs[file.exists(data_dirs)]
datainfo <- function (dir) {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[8], 
             patient_id = sub(".json", "", list.files(dir, pattern = ".json")))
}
curr_all_data <- unique(bind_rows(lapply(data_dirs, datainfo)))

test_list <- list()
for (i in 1:nrow(curr_all_data)) {
  test3 <- readLines(glue("{leu_ana_path}/{curr_all_data[i, 1]}/human/rnaseq/splicing/filt3r/{curr_all_data[i, 1]}/{curr_all_data[i, 2]}.json"), warn = F)
  test2 <- sub("NaN|nan", 999, test3) %>% fromJSON() 
  test1 <- do.call(cbind, test2[c(6, 5, 7)])
  test1[["sample_id"]] <- glue("{curr_all_data[i, 1]}@{curr_all_data[i, 2]}")
  test_list[[i]] <- test1
}

flt3_itd_all <- bind_rows(test_list)
flt3_itd_fil <- flt3_itd_all[flt3_itd_all$is_duplication == "TRUE", ]
#readr::write_rds(flt3_itd_fil, "/cluster/home/yjliu_jh/projects/leukemia/data/itd_data2.rds")
readr::write_tsv(flt3_itd_fil, "/cluster/home/yjliu_jh/projects/leukemia/data/itd_data2.tsv")



# ======= find samples without filt3r results ==========
leu_data_path <- "/cluster/home/jhuang/projects/leukemia/data"
data_dirs <- paste0(leu_data_path, "/", leu_datasets, "/human/rnaseq")
data_dirs <- data_dirs[file.exists(data_dirs)]
pts <- "_1.fastq.gz|_2.fastq.gz|_R1.fq.gz|_R2.fq.gz|_2.fq.gz|_1.fq.gz|.fq.gz"
datainfo <- function (dir) {
  data.frame(datasets = str_split(dir, pattern = "/", simplify = T)[8], 
             patient_id = sub(pts, "", list.files(dir, pattern = ".fq.gz|.fastq.gz")))
}
curr_all_data_fq <- unique(bind_rows(lapply(data_dirs, datainfo)))
anti_join(curr_all_data_fq, curr_all_data)
