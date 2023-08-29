pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "SummarizedExperiment", "jhuanglabRNAseq", "ComplexHeatmap")
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

fpkm <- assay(se, "tpm_remove_rnatype") %>% as.data.frame() %>% dplyr::select(metadata$sample_id) # exp matrix
dat <- fpkm |> as.data.frame() |> tibble::rownames_to_column("gene_name") |>
  filter_black_heatmap_genes() |> tibble::remove_rownames() |> tibble::column_to_rownames("gene_name")


allcolour = c("#DC143C", "#0000FF", "#20B2AA", "#FFA500",
                       "#9370DB", "#98FB98", "#F08080", "#1E90FF", "#7CFC00",
                       "#FFFF00", "#808000", "#FF00FF", "#FA8072", "#7B68EE",
                       "#9400D3", "#800080", "#A0522D", "#D2B48C", "#D2691E",
                       "#87CEEB", "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD",
                       "#008B8B", "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A",
                       "#4682B4", "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE",
                       "#FF6347", "#6A5ACD", "#9932CC", "#8B008B", "#8B4513",
                       "#DEB887")

anno_col <- ggsci::pal_nejm("default")(8)
meta_plot <- metadata
meta_plot$fusion_plot <- metadata$fusion_genes_rna

ha <-  HeatmapAnnotation(subgroups = meta_plot$subgroups,
                         subgroups_hjy = meta_plot$subgroups_hjy,
                         fusion_type = meta_plot$fusion_plot,
                         rna_type = meta_plot$rna_type,
                         gender = meta_plot$gender,
                         col = list(rna_type = c("total_RNA" = anno_col[1],
                                                 "mRNA" = anno_col[2]),
                                    gender = c("male" = anno_col[1],
                                               "female" = anno_col[2]),
                                    disease_type = setNames(allcolour[1:length(unique(meta_plot$disease_type))],
                                                            unique(meta_plot$disease_type))),
                         show_legend = rep(FALSE, 2),
                         simple_anno_size = unit(0.2, "cm"),
                         annotation_name_gp = gpar(fontsize = 6))

quick_heatmap(dat, ha = ha, outdir = "aml/step1", column_split = 26, top_var_percent = 0.95)




sampleinfo <- readxl::read_excel("aml/step1/sampleinfo_0.9.xlsx", sheet = 1)
idx <- match(sampleinfo$sample_id, rds$sample_id)
rds$subgroups_hjy[idx] <- glue("aml_{sampleinfo$sub_groups}")
#unique(rds$subgroups_hjy)
readr::write_rds(rds, rds_fn)






# ====== temp g1 g2 mef2d =====

sampleinfo <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx", sheet = 1)
g1_2_samples <- sampleinfo$sample_id[sampleinfo$sub_groups %in% c("G1", "G2")]
temp2 <- fusions_new[fusions_new$sample_id %in% g1_2_samples, ]
fus12mef <- temp2[grepl("MEF2D", temp2$gene5), ]


xxe <- c("BCL9", "SS18", "FOXJ2", "CSF1R", "DAZAP1", "STAT6", "HNRNPUL1", "HNRNPH1", "HNRNPM")

temp3 <- fus12mef[fus12mef$gene3 %in% xxe, ]

temp1 <- fusions_new[grepl("MEF2D", fusions_new$gene5), ]
temp11 <- temp1[temp1$gene3 %in% xxe, ]

meta_plot$fusion_plot[meta_plot$sample_id %in% temp11$sample_id] <- "aml_MEF2D"





# ===== add whitelist fusions anyway

find_top2 <- function(sample_list, fusion_dat, top = 20){
  fusions <- unique(fusion_dat[fusion_dat$sample_id %in% sample_list, c(2, 7)])
  t0 <-  table(fusions$fusion_anno)
  whitelist <- c("CBYB-MYH11", "RUNX1-RUNX1T1", "DEK-NUP214", "BCR-ABL1",
                 "KMT2A-r", "MECOM-r", "NUP98-r", "PAX5-r", "DUX4-r", "CRLF2-r", 
                 "MLLT10-r", "LMO1-r", "LMO2-r")
  temp <- as.data.frame(sort(t0, decreasing = T)[1:top])
  temp2 <- t0[]
  temp$prop <- temp$Freq / length(sample_list)
  temp
}





# ========

count_all <- numeric()
for (i in 1:length(hcc_samples)){
  hcc_temp <- get(paste0(hcc_samples[i], "_two"))
  count_all[i] <- nrow(hcc_temp@images$slice1@coordinates)
}
count_all <- data.frame(samples = hcc_samples, spotcount = count_all)

mgdiff <- left_join(mg_diff, count_all)
mgdiff$abs_count <- pmin(mgdiff$spotcount - mgdiff$count05, mgdiff$count05)
mgdiff <- as.data.frame(mgdiff)
mgsub <- mgdiff[mgdiff$abs_count > 100 & mgdiff$count >= 3, ]
mgsub2 <- left_join(mgsub[, c("mz", "gene")], mgdiff)






