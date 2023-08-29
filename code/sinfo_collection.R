
leu <- readr::read_rds("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds")
sampleinfo_r <- as.data.frame(colData(leu))


info_path <- "/cluster/home/yjliu_jh/projects/leu_j/data/leukemia_anno_new.xlsx"
collected_info_sheets <- readxl::excel_sheets(info_path)
collected_info <- lapply(collected_info_sheets, function(X) readxl::read_excel(info_path, sheet = X))
names(collected_info) <- collected_info_sheets


coln_info <- list()
dtst <- character()
for (i in 1:length(collected_info)) {
  coln_info[[i]] <- colnames(collected_info[[i]])
  dtst <- c(dtst, rep(collected_info_sheets[i], ncol(collected_info[[i]])))
}

col_match <- data.frame(dataset = dtst, col = unlist(coln_info))
colns <- unique(col_match$col)




# manually alter that xlsx before doing this step! 
# collect os columns and change months to days
survival_cols <- colns[grep("vital", colns, ignore.case = T)]
survival_cols2 <- colns[grep("outcome", colns, ignore.case = T)]
survival_cols3 <- colns[grep("os", colns, ignore.case = T)]
survival_cols3 <- survival_cols3[c(1, 4, 8, 14, 18, 19)]
survival_cols4 <- colns[grep("surv", colns, ignore.case = T)]




scols <- c(survival_cols, survival_cols2, survival_cols3, survival_cols4)
sdata <- col_match[col_match$col %in% scols, ]
sdata <- sdata[-c(2, 5, 6, 9, 10), ]
sdata <- sdata[c(1, 4), ]
sdata$os <- c("new os", "os_days")

temp_all <- data.frame(matrix(NA, 0, 3))
colnames(temp_all) <- c("sample_id", "vital_status", "os_time")

for (i in 1:nrow(sdata)){
  temp <- data.frame(sample_id = collected_info[[sdata$dataset[i]]]$sample_id,
                     vital_status = collected_info[[sdata$dataset[i]]][[sdata$col[i]]],
                     os_time = collected_info[[sdata$dataset[i]]][[sdata$os[i]]])
  temp_all <- rbind(temp_all, temp)
}



all2_sur <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_ClinicalData_Phase_II_Discovery_20221108.xlsx")
all1_sur <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_ClinicalData_Phase_I_20221108.xlsx")
all3_sur <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_ClinicalData_Phase_III_20221108.xlsx")
all2_sur <- as.data.frame(all2_sur[, c("TARGET USI", "Vital Status", "Overall Survival Time in Days")])
all1_sur <- as.data.frame(all1_sur[, c("TARGET USI", "Vital Status", "Overall Survival Time in Days")])
all3_sur <- as.data.frame(all3_sur[, c("TARGET USI", "Vital Status", "Overall Survival Time in Days")])

all_sur <- rbind(all1_sur, all2_sur)
all_sur <- rbind(all_sur, all3_sur)

all2_cl <- read.delim("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_mRNA-seq_Phase2_20220110.sdrf.txt")
all3_cl <- read.delim("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_mRNA-seq_Phase3_20181025.sdrf.txt")
all1_cl <- read.delim("/cluster/home/yjliu_jh/projects/leukemia/data/anno/TARGET_ALL_mRNA-seq_Phase1_20170609.sdrf.txt")

all2_cl <- unique(all2_cl[, c("Source.Name", "Comment.SRA_RUN.", "Comment.SRA_RUN..1")])
all3_cl <- unique(all3_cl[, c("Source.Name", "Comment.SRA_RUN.", "Comment.SRA_RUN..1")])
all1_cl <- unique(all1_cl[, c("Source.Name", "Comment.SRA_RUN.", "Comment.SRA_RUN..1")])

all_cl <- rbind(all1_cl, all2_cl)
all_cl <- rbind(all_cl, all3_cl)

colnames(all_cl) <- c("target", "srr1", "srr2")
all_cl$srr1[all_cl$srr1 == ""] <- NA
all_cl$srr2[all_cl$srr2 == ""] <- NA
all_cl <- all_cl %>% mutate(sample_id = coalesce(srr1, srr2))


colnames(all_sur) <- c("target", "vital_status", "os_time")
all_cl <- all_cl[, c(1, 4)]

p22 <- read.delim("/cluster/home/yjliu_jh/projects/leukemia/data/anno/clinical_all_p2.tsv")
all2_sur_extra <- p22[, c("case_submitter_id", "vital_status", "days_to_last_follow_up")]
colnames(all2_sur_extra) <- c("target", "vital_status", "os_time")
all_sur <- rbind(all_sur, all2_sur_extra)

xall <- unique(inner_join(all_cl, all_sur))
xall <- na.omit(xall)
xall$os_time <- as.numeric(xall$os_time)
xall <- na.omit(xall)
temp_all2 <- rbind(temp_all, xall[, 2:4])
temp_all2$vital_status[temp_all2$vital_status %in% "D"] <- "Dead"
temp_all2$vital_status[temp_all2$vital_status %in% "A"] <- "Alive"

sinfo_fil <- inner_join(sinfo_all, unique(temp_all2))
mat_fil <- mat[, colnames(mat) %in% sinfo_fil$sample_id]

readr::write_rds(sinfo_fil, "/cluster/home/yjliu_jh/projects/leukemia/data/ball_sampleinfo_with_survival.rds")
readr::write_rds(mat_fil, "/cluster/home/yjliu_jh/projects/leukemia/data/ball_expression_matrix_with_survival.rds")



# get mef2d and tcf3 fusion samples, compare with other fusion samples, calculate differences


sinfo_sub <- sinfo_fil[, c("sample_id", "vital_status", "os_time")]









reslist <- list()
# loop
for (i in 1:length(hcc_samples)) {
  curr_sample <- hcc_samples[i]
  diff_res <- readr::read_rds(glue::glue("/cluster/home/yjliu_jh/projects/hcc/data/ggdiff_075_{curr_sample}.rds"))
  res <- bind_rows(diff_res, .id = "gene1")
  res <- res[order(res$abs_diff, decreasing = T), ]
  colnames(res)[4] <- "gene2"
  join1 <- counts_new[counts_new$samples %in% curr_sample, c("gene", "count05")]
  colnames(join1) <- c("gene1", "count1")
  join2 <- counts_new[counts_new$samples %in% curr_sample, c("gene", "count05", "count067", "count075")]
  colnames(join2) <- c("gene2", "count2_05", "count2_067", "count2_075")
  res <- left_join(res, join1)
  res <- left_join(res, join2)
  res[, 8:10] <- abs(res[, 8:10] - res[, 7])
  res$mincount <- with(res, pmin(count2_05, count2_067, count2_075))
  # some loose thresholding procedures
  res_fil <- res[res$abs_diff > 0.2 | res$abs_diff < -0.1, ]
  reslist[[i]] <- res_fil
}

resall <- bind_rows(reslist)
resall <- resall[, c(1:7, 11)]
readr::write_rds(resall, "/cluster/home/yjliu_jh/projects/hcc/output/dist/integrated/gene_gene/ggdiff_075.rds")













