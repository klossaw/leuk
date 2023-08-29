










ball_sinfo <- readxl::read_excel("/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/BALL/all_MEF2D/fusion_tsv1/sampleinfo_0.95.xlsx")
ball_sinfo <- as.data.frame(ball_sinfo)
colnames(ball_sinfo)[2] <- "rna_group"
fil <- cld$disease_type %in% c("B-ALL")
metadata_ball <- as.data.frame(cld[fil, ])
meta_part_ball <- left_join(metadata_ball, ball_sinfo)



tall_sinfo <- readxl::read_excel("/cluster/home/yjliu_jh/projects/temp/sampleinfo_0.8.xlsx")
tall_sinfo <- as.data.frame(tall_sinfo)
colnames(tall_sinfo)[2] <- "rna_group"
fil <- cld$disease_type %in% c("T-ALL")
metadata_tall <- as.data.frame(cld[fil, ])
meta_part_tall <- left_join(metadata_tall, tall_sinfo)



groups <- unique(meta_part_tall$rna_group)
f_list <- list()
for(i in 1:length(groups)){
  f_list[[groups[i]]] <- find_top(meta_part_tall$sample_id[meta_part_tall$rna_group %in% groups[i]], fusions_new, 20)
}
f_top <- bind_rows(f_list, .id = "rna_group")
f_top_wide <- f_top[, -3] %>% pivot_wider(names_from = "Var1", values_from = "prop", values_fill = NA)
f_top_filtered <- f_top_wide[, c(1, which(colSums(is.na(f_top_wide)) >= 13))]
ftf <- f_top_filtered %>% pivot_longer(-rna_group) %>% na.omit() %>% as.data.frame()
colnames(ftf)[2:3] <- c("fusion", "prop")
tj <- as.data.frame(table(meta_part_tall$rna_group))
colnames(tj) <- c("rna_group", "sample_count")
ftf <- left_join(ftf, tj)
ftf_tall %>% df2excel("/cluster/home/yjliu_jh/projects/temp/tall-groups.xlsx")




# first build dds  then loop to calculate contrast
counts_aml <- round(as.matrix(assay(se)))
counts_aml[counts_aml < 0] <- 0
groupings <- data.frame(sample_id = colnames(counts_aml))
groupings <- left_join(groupings, clusterings)

# also to find annotated fusions/mutations group from sampleinfo

all_groups <- unique(groupings$rna_group)
diff_gene_res <- list()

for (i in 1:length(all_groups)){
  curr_group <- as.character(all_groups[i])
  groupings_x <- groupings
  groupings_x$rna_group <- ifelse(groupings_x$rna_group %in% curr_group, curr_group, "other_groups")
  dds1 <- DESeqDataSetFromMatrix(countData = counts_aml,
                                 colData = groupings_x,
                                 design = ~ rna_group)
  keep <- rowSums(counts(dds1) >= 2) > ceiling(0.2 * ncol(counts_aml))
  dds1 <- dds1[keep,]
  keep2 <- rowSums(counts(dds1) == 0) < ceiling(0.5 * ncol(counts_aml))
  dds1 <- dds1[keep2,]
  suppressMessages(dds <- DESeq(dds1, betaPrior = TRUE)) 
  diff_gene_res[[curr_group]] <- results(dds, contrast = list(curr_group, "other_groups"))
}
readr::write_rds(diff_gene_res, "/cluster/home/yjliu_jh/projects/temp/aml_groups_de.rds")





out_dir <- "/cluster/home/yjliu_jh/projects/hyperion/output/data"
sparse_score_all <- readr::read_rds(glue::glue("{out_dir}/sparse/score_all_march22.rds"))

