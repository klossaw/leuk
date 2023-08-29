pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "SummarizedExperiment", "jhuanglabRNAseq")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "leukemia"
workdir <- glue("~/projects/{project}/data")
setwd(workdir)

# load mutations called by three tools
all_mutations <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/all_mutations.rds")
# load current clusterings for B-ALL
path0 <- "/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/BALL/all_MEF2D/230702/sampleinfo_0.95.xlsx"
path1 <- "/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx"
path2 <- "/cluster/home/ztao_jh/projects/leukemia/analysis/meta/human/figures/heatmap/tall/step1/select_heatmap/sampleinfo_0.9.xlsx"
clu_ball <- readxl::read_excel(path0)  ## may change the path to other groups
clu_ball <- as.data.frame(clu_ball)
colnames(clu_ball)[2] <- "rna_group"
ball_mut <- all_mutations[all_mutations$sample_id %in% clu_ball$sample_id, ]

# === 1 top mutations appeared - gene level ===
# function
find_top <- function(sample_list, mut_dat, col = "Hugo_Symbol", top = 20){
  mutations <- unique(mut_dat[mut_dat$sample_id %in% sample_list, c(1, 2)])
  temp <- as.data.frame(sort(table(mutations[[col]]), decreasing = T)[1:top])
  temp$prop <- temp$Freq / length(sample_list)
  temp
}

# pre-filtering, remove mutations widely across the genome 
gene_freq_ball <- sort(table(unique(ball_mut[, 1:2])$Hugo_Symbol), decreasing = T)
gfb_fil <- gene_freq_ball[gene_freq_ball < 0.67 * nrow(clu_ball)]
ball_mut_fil <- ball_mut[ball_mut$Hugo_Symbol %in% names(gfb_fil), ]
# - use dplyr::filter(nchar(tool) > 3) if [tool > 1] needed 
# select top genes for each cluster
groups <- unique(clu_ball$rna_group)
f_list <- list()
for(i in 1:length(groups)){
  top_genes <- find_top(clu_ball$sample_id[clu_ball$rna_group %in% groups[i]], ball_mut_fil, top = 200)
  f_list[[groups[i]]] <- top_genes[top_genes$prop >= 0.75, ]
}
# get candidate genes 
f_top <- bind_rows(f_list, .id = "rna_group")
f_top_wide <- f_top[, -3] %>% pivot_wider(names_from = "Var1", values_from = "prop", values_fill = NA)
f_top_filtered <- f_top_wide[, c(1, which(colSums(is.na(f_top_wide)) >= (nrow(f_top_wide) - 2)))]
ftf <- f_top_filtered %>% pivot_longer(-rna_group) %>% na.omit() %>% as.data.frame()
colnames(ftf)[2:3] <- c("mutation", "prop")

# enhance specificity
ball_mut_freq <- unique(ball_mut_fil[, 1:2]) %>% group_by(Hugo_Symbol) %>%
  summarise(freq = n() / nrow(clu_ball)) %>% as.data.frame()
colnames(ball_mut_freq)[1] <- "mutation"
ftf <- left_join(ftf, ball_mut_freq)
tj <- as.data.frame(table(clu_ball$rna_group))
colnames(tj) <- c("rna_group", "sample_count")
#tj$freq_group <- tj$sample_count / nrow(clu_ball)
ftf <- left_join(ftf, tj)
n_ball <- nrow(clu_ball)
ftf <- ftf %>% mutate(diff = prop - (n_ball * freq - sample_count * prop) / (n_ball - sample_count))
ftf_fil <- ftf[ftf$diff > 0.6, ]
ftf_fil %>% df2excel("/cluster/home/yjliu_jh/projects/temp/ball_specific_mutation.xlsx")


# === 2 top mutations appeared - locus level ===

# select top genes for each cluster

ball_mut2 <- ball_mut
ball_mut2$mut <- paste0(ball_mut2$Hugo_Symbol, "_", ball_mut2$HGVSp_Short)
ball_mut2 <- ball_mut2[, c(1, 7, 6)]

gene_freq_ball <- sort(table(unique(ball_mut2[, 1:2])$mut), decreasing = T)
gfb_fil <- gene_freq_ball[gene_freq_ball < 0.5 * nrow(clu_ball)]
ball_mut_fil2 <- ball_mut2[ball_mut2$mut %in% names(gfb_fil), ]

f_list <- list()
for(i in 1:length(groups)){
  top_genes <- find_top(clu_ball$sample_id[clu_ball$rna_group %in% groups[i]], ball_mut_fil2, 
                        col = "mut", top = 200)
  f_list[[groups[i]]] <- top_genes
}
# get candidate genes 
f_top <- bind_rows(f_list, .id = "rna_group")
f_top_wide <- f_top[, -3] %>% pivot_wider(names_from = "Var1", values_from = "prop", values_fill = NA)
f_top_filtered <- f_top_wide[, c(1, which(colSums(is.na(f_top_wide)) >= (nrow(f_top_wide) - 2)))]
ftf <- f_top_filtered %>% pivot_longer(-rna_group) %>% na.omit() %>% as.data.frame()
colnames(ftf)[2:3] <- c("mutation", "prop")

# enhance specificity
ball_mut_freq <- unique(ball_mut_fil2[, 1:2]) %>% group_by(mut) %>%
  summarise(freq = n() / nrow(clu_ball)) %>% as.data.frame()
colnames(ball_mut_freq)[1] <- "mutation"
ftf <- left_join(ftf, ball_mut_freq)
ftf <- left_join(ftf, tj)
ftf <- ftf %>% mutate(diff = prop - (n_ball * freq - sample_count * prop) / (n_ball - sample_count))
ftf_fil <- ftf[ftf$diff > 0.5, ]
ftf_fil %>% df2excel("/cluster/home/yjliu_jh/projects/temp/ball_specific_mutation_locus.xlsx")







#cmaf <- distinct(cmaf)
#cmaf$software[cmaf$software %in% "speedseq"] <- "SS"
#temp <- as.data.table(cmaf)[, toString(software), by = eval(names(cmaf)[1:5])]








# 



