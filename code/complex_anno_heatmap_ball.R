pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", 
          "jhuanglabRNAseq", "ComplexHeatmap", "SummarizedExperiment", "matrixStats")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
#project <- "project"
#dataset <- "dataset"
#species <- "human"
#workdir <- glue("~/projects/{project}/analysis/{dataset}/{species}/rnaseq")
#setwd(workdir)


# ======== initial data load ============

# load data: mutation, fusion, ttmv, itd
mutations_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/temp/all_mutations_0716.rds")
fusions_new <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/fusions_new.rds")
ttmv_data <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/ttmv_data.rds")
itd_data <- readr::read_rds("/cluster/home/yjliu_jh/projects/leukemia/data/itd_data.rds")
# ball specific fusions
fusions_new$fusion_anno[fusions_new$gene5 %in% "MEF2D"] <- "MEF2D-r"
fusions_new$fusion_anno[fusions_new$gene3 %in% "HLF"] <- "HLF-r"
fusions_new$fusion_anno[grepl("DUX4", fusions_new$fusion)] <- "DUX4-r"
fusions_new$fusion_anno[fusions_new$gene3 %in% "ZNF384"] <- "ZNF384-r"

# load clusterings info
all_clustering_anno <- read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/need_to_add_info_230712.csv")
aml_anno <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
ball_anno <- readxl::read_excel("/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/BALL/all_MEF2D/230710/sampleinfo_0.95.xlsx")
tall_anno <- readxl::read_excel("/cluster/home/ztao_jh/projects/leukemia/analysis/meta/human/figures/heatmap/tall/step1/select_heatmap22lab/sampleinfo_0.9.xlsx")
tall_anno <- left_join(tall_anno, all_clustering_anno[, c("sample_id", "subgroups_230712")])
ball_anno <- left_join(ball_anno, all_clustering_anno[, c("sample_id", "subgroups_230712")])
aml_anno <- left_join(aml_anno, all_clustering_anno[, c("sample_id", "subgroups_230712")])

# load data: expression profile
metadata <- read_csv("/cluster/home/ylxie_jh/projects/leukemia/analysis/meta/jhuang/human/figures/heatmap/allBALL_fusion_mutation_metadata0712.csv")
#BALL total
seeds <- 2023
rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
cld <- colData(rds)
se <- rds[, metadata$sample_id]
fpkm <- assay(se, "tpm_remove_rnatype") %>% as.data.frame() %>% dplyr::select(metadata$sample_id) # exp matrix
dat <- fpkm %>% as.data.frame() %>% tibble::rownames_to_column("gene_name") %>%
  filter_black_heatmap_genes() %>% tibble::remove_rownames() %>% tibble::column_to_rownames("gene_name")
metadata <- metadata %>% as.data.frame() %>% `rownames<-`(.$sample_id) %>% .[colnames(dat), ]


# ========= functions ============

# add new columns to annotation data.frame (exist or not)
enrich_anno <- function(anno, data, extra = NULL, sample_col = "sample_id") {
  if (is.null(extra)) { extra <- names(data) }
  for (col in extra){
    anno[[col]] <- ifelse(anno[[sample_col]] %in% data[[col]], "exist", "non-exist")
  }
  anno
}

# set colors for columns
auto_anno <- function(data, cols, colors) {
  anno_list <- list()
  for (col in cols) {
    elements <- sort(unique(data[[col]]))
    anno_list[[col]] <- setNames(colors[1:length(elements)], elements)
  }
  anno_list
}

# add new columns to annotation data.frame (based on annotation col)
enrich_anno2 <- function(anno, data, extra, sample_col = "sample_id", anno_col = "anno") {
  for (col in extra){
    new_col <- left_join(anno, data[[col]][, c(sample_col, anno_col)])
    anno[[col]] <- new_col[[anno_col]]
  }
  anno
}

# add separated columns to annotation data.frame (based on single annotation col with potential duplicates)
enrich_anno3 <- function(anno, data, sample_col = "sample_id", anno_col = "anno") {
  spl <- split(data[[sample_col]], f = data[[anno_col]])
  enrich_anno(anno, spl)
}



# ========== extra functions designed for different data ===========

# select hotspots for mutations
mut_hotspot <- function (gene, mut_data, col = "HGVSp_Short", hotspot = NULL, k = 10, sum = T) {
  mut_gene <- mut_data[mut_data$Hugo_Symbol %in% gene, ] %>% as.data.frame()
  if (is.null(hotspot)) {
    hotspot <- names(sort(table(mut_gene[[col]]), decreasing = T))[1:k]
  } 
  mut_gene_fil <- mut_gene[mut_gene[[col]] %in% hotspot, ] 
  if (sum == T) {
    mut_gene_fil <- mut_gene_fil[, c("sample_id", col)] %>% unique() %>% group_by(sample_id) %>% 
      summarise(mutation = paste(.data[[col]], collapse = ", "))
  } else { mut_gene_fil$mutation <- paste0(mut_gene_fil$Hugo_Symbol, "_", mut_gene_fil[[col]]) }
  mut_gene_fil
}


# some others when needed?
# separate itd types
# fuzzy-end or not
# ...


# ========== plot color settings =============

# set colors
all_colors <- c("#DC143C", "#0000FF", "#20B2AA", "#FFA500", "#9370DB", 
                "#98FB98", "#F08080", "#1E90FF", "#7CFC00", "#FFFF00",
                "#808000", "#FF00FF", "#FA8072", "#7B68EE", "#9400D3",
                "#800080", "#A0522D", "#D2B48C", "#D2691E", "#87CEEB",
                "#40E0D0", "#5F9EA0", "#FF1493", "#0000CD", "#008B8B",
                "#FFE4B5", "#8A2BE2", "#228B22", "#E9967A", "#4682B4",
                "#32CD32", "#F0E68C", "#FFFFE0", "#EE82EE", "#FF6347",
                "#6A5ACD", "#9932CC", "#8B008B", "#8B4513", "#DEB887")
                         

col_fun <- circlize::colorRamp2(c(0, 10, 100), c("#FFEEEE", "#FFBBBB", "#FF0000"))
anno_col <- c("RARG" = "yellow", "RARB" = "#FF8534", "TTMV-RARA" = "#E25385",
              "pml-others" = "#63E278", "RARA-others" = "lightgreen",
              "APL" = "#8864B6", "APL-like" = "steelblue", "AML" = "#96B2F3")
anno_def <- ggsci::pal_nejm("default")(8)


# ========== plot annotation settings ===========

# set gene list (change here)
mut_genes <- c("KRAS", "WT1", "DNMT3A", "NUTM1", "PAX5")
fus_genes <- c("KMT2A-r", "BCR-ABL1", "ETV6-RUNX1", "TCF3-PBX1", "HLF-r", "DUX4-r",
               "MEF2D-r", "BCL-MYC", "IGH-CEBPE", "ZNF384-r")
itd_genes <- c("FLT3-ITD", "KMT2A-PTD")

# list construction
mutall_list <- list()
for (i in 1:length(mut_genes)) {
  mutall_list[[mut_genes[i]]] <- unique(mutations_new$sample_id[mutations_new$Hugo_Symbol %in% mut_genes[i]])
}

fusion_list <- list()
for (i in 1:length(fus_genes)) {
  fusion_list[[fus_genes[i]]] <- unique(fusions_new$sample_id[fusions_new$fusion_anno %in% fus_genes[i]])
}

itd_list <- list()
for (i in 1:length(itd_genes)) {
  itd_list[[itd_genes[i]]] <- unique(itd_data$sample_id[itd_data$itd %in% itd_genes[i]])
}

# ttmv list
ttmv_slim <- unique(ttmv_data[, c("gene", "dataset", "sample_id")])
ttmv_freq <- sort(table(ttmv_slim$gene), decreasing = T)
candidates <- c("HNRNPL", names(ttmv_freq)[1:20][!grepl(".1", names(ttmv_freq)[1:20])])
sample_ttmv_list <- list()
sample_ttmv_list[["HNRNPL"]] <- unique(ttmv_data$sample_id[ttmv_data$gene %in% "AC008982.1/HNRNPL"])
for (i in 2:length(candidates)) {
  sample_ttmv_list[[candidates[i]]] <- unique(ttmv_data$sample_id[ttmv_data$gene %in% candidates[i]])
}

# ============ mutation hotspot settings =========
# set mutation (hotspots)
# according to the proteinpaint tool:
idh1_hotspots <- c("p.R132H", "p.R132C", "p.R132G", "p.R132S", "p.R132L", "p.R132V",
                   "p.R132I", "p.P33S", "p.R100Q")
idh2_hotspots <- c("p.R140Q", "p.R172K", "p.R172S", "p.R172M", "p.R172G", "p.R172W",
                   "p.R140W", "p.R140L", "p.R140G")
asxl1_hotspots <- c("p.G645fs", "", "p.R634_I641fs", "p.R693*", "p.L19_E1splice_region",
                    "p.Y591fs", "p.Y591*", "p.K85R")

idh1_mut <- mut_hotspot("IDH1", mutations_new, hotspot = idh1_hotspots)
idh2_mut <- mut_hotspot("IDH2", mutations_new, hotspot = idh2_hotspots)
asxl1_mut_sep <- mut_hotspot("ASXL1", mutations_new, sum = F)
pax5_mut <- mut_hotspot("PAX5", mutations_new, hotspot = "p.P80R")
ikzf1_mut <- mut_hotspot("IKZF1", mutations_new, hotspot = "p.N159Y")
zeb2_mut <- mut_hotspot("ZEB2", mutations_new, hotspot = "p.H1038R")

mut2_list <- list(idh1_mut, idh2_mut, pax5_mut, ikzf1_mut, zeb2_mut)
gene2 <- c("IDH1", "IDH2", "PAX5", "IKZF1", "ZEB2")
names(mut2_list) <- gene2


# ========== final annotation data ============
# set annotation
metadata <- colData(se)
meta_plot2 <- as.data.frame(metadata)[, c("gender", "age", "rna_type", "disease_type", "datasets")]
meta_plot2$sample_id <- rownames(meta_plot2)
meta_plot2$age <- as.numeric(meta_plot2$age)
meta_plot3x <- enrich_anno(meta_plot2, mutall_list, mut_genes)
meta_plot4x <- enrich_anno(meta_plot3x, itd_list, itd_genes)
meta_plot4x <- enrich_anno(meta_plot4x, fusion_list, fus_genes)
meta_plot4x <- left_join(meta_plot4x, ball_anno)  ## switch to other subtype if changes were made above (line:35)
meta_plot5 <- enrich_anno2(meta_plot4x, mut2_list, gene2, anno_col = "mutation")
#meta_plot5 <- enrich_anno3(meta_plot5, asxl1_mut_sep, anno_col = "mutation")
asxl1_hs <- unique(asxl1_mut_sep$mutation)

ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot5[, -c(3, 4, 6)]),
  show_legend = c(TRUE, TRUE, FALSE, rep(F, length(mut_genes) + length(itd_genes) + length(fus_genes)),
                  F, T, rep(T, length(gene2)) #rep(F, length(asxl1_hs))),
  ),
  col = c(list(age = col_fun, 
               gender = setNames(anno_def[1:2], c("female", "male"))),
          auto_anno(meta_plot5, c(mut_genes, itd_genes, fus_genes#, asxl1_hs
                                  ), c("magenta", "ivory")),
          auto_anno(meta_plot5, c(gene2, "sub_groups", "subgroups_230712"), all_colors))
)

# output plots
quick_heatmap(dat, ha = ha_x, width = 28, height = 16, outdir = "leu/ball/new", column_split = 20, top_var_percent = 0.95)
























# aml code
clusterings <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
quick_heatmap(dat, ha = ha_x, height = 9, outdir = "aml/new", column_split = 26, top_var_percent = 0.95)

# ball code


ball_samples <- readxl::read_excel("BALL/all_MEF2D/230710/sampleinfo_0.95.xlsx") %>%
  dplyr::rename("rna_group" = "sub_groups")


# tall code
quick_heatmap(dat, ha = ha2, outdir = glue("~/projects/{project}/analysis/{dataset}/{species}/figures/heatmap/tall/step1/select_heatmap22lab"),
              column_split = 27, top_var_percent = 0.9, height = 20)
clustersheeet <- glue("~/projects/{project}/analysis/{dataset}/{species}/figures/heatmap/tall/step1/select_heatmap22/sampleinfo_0.9.xlsx") %>%
  readxl::read_excel(sheet = 1)





