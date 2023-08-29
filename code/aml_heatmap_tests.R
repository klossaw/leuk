# ======== first half - check identifiers/sample_id =========




















# ======== second half - heatmap test with mutation data ========
itd1 <- readr::read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/leukemia/itd.csv")
all_mutations <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/all_mutations.rds")
fusions_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/fusions_new.rds")


all_mut_fil <- all_mutations[all_mutations$tool %in% "HC, VS, SS", ]

mut_genes <- c("KRAS", "WT1", "CEBPA", "NPM1", "DNMT3A", "CYB5R2", "ADAMTS5", "CEBPD", "COL6A5", "COBL", "SLC22A1")
mutfil_list <- list()
for (i in 1:length(mut_genes)) {
  mutfil_list[[mut_genes[i]]] <- unique(all_mut_fil$sample_id[all_mut_fil$Hugo_Symbol %in% mut_genes[i]])
}

mutall_list <- list()
for (i in 1:length(mut_genes)) {
  mutall_list[[mut_genes[i]]] <- unique(all_mutations$sample_id[all_mutations$Hugo_Symbol %in% mut_genes[i]])
}

fus_genes <- c("PML-RARA", "KMT2A-r", "RUNX1-RUNX1T1", "CBFB-MYH11", "PCAT18-KCTD1", "NUP98-r", "CBFA2T3-GLIS2")
fusion_list <- list()
for (i in 1:length(fus_genes)) {
  fusion_list[[fus_genes[i]]] <- unique(fusions_new$sample_id[fusions_new$fusion_anno %in% fus_genes[i]])
}






itd_list <- list()
temp_genes <- c("flt3", "kmt2a")
for (i in 1:length(temp_genes)) {
  pids <- unique(itd1$sample_name[itd1$itd %in% temp_genes[i]])
  itd_list[[temp_genes[i]]] <- cld$sample_id[cld$patient_id %in% pids]
}

meta_plot3x <- enrich_anno(meta_plot2, mutall_list, mut_genes)
meta_plot4x <- enrich_anno(meta_plot3x, itd_list, temp_genes)
meta_plot4x <- enrich_anno(meta_plot4x, fusion_list, fus_genes)
meta_plot4x <- left_join(meta_plot4x, clusterings[, 1:2])

groups <- unique(meta_plot4x$rna_group)

ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot4x[, c("gender", "age", "datasets", mut_genes, temp_genes, fus_genes, "rna_group")]),
  show_legend = c(TRUE, TRUE, FALSE, rep(F, length(mut_genes) + length(temp_genes) + length(fus_genes)), F),
  col = append(append(list(age = col_fun, 
                    gender = setNames(anno_def[1:2], c("female", "male")),
                    rna_group = setNames(all_colors[1:length(groups)], groups)), 
               auto_anno(meta_plot4x, mut_genes, c("magenta", "ivory"))), 
               append(auto_anno(meta_plot4x, temp_genes, c("magenta", "ivory")),
                      auto_anno(meta_plot4x, fus_genes, c("magenta", "ivory"))))
)

out_path <- "/cluster/home/yjliu_jh/projects/temp/testxsssxxxheatmap2.pdf"
test_heat(dat_mat, ha_x, out_path)


quick_heatmap(dat, ha = ha_x, height = 9, outdir = "aml/new", column_split = 26, top_var_percent = 0.95)












# ======= extra - functions ==========














