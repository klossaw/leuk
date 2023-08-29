pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "optparse",
          "vroom", "jhtools", "glue", "ggsci", "tidyverse", "dplyr", "Seurat",
          "ggridges", "viridis", "circlize", "ComplexHeatmap", "matrixStats")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


liver <- readr::read_rds("/cluster/home/jhuang/projects/liver/analysis/HRA001748/human/tenx/rds/seu_A124.counts_anno.rds")




"/cluster/home/yjliu_jh/projects/leukemia/output/aml_samples_anno.xlsx"

allx_list <- append(mutall_list, append(fusion_list, itd_list))


colnames(metainfo)[3] <- "label"
metainfo$label[metainfo$label %in% c("aml_NUP98_rearrangement", "aml_NUP98")] <- "aml_NUP98-r"
metainfo$label[metainfo$rna_group %in% c("G3", "G4", "G10", "G12") & metainfo$sample_id %in% allx_list[["NPM1"]]] <- "aml_NPM1"
metainfo$label[metainfo$rna_group %in% c("G3", "G4") & metainfo$sample_id %in% allx_list[["flt3"]]] <- paste0(
  metainfo$label[metainfo$rna_group %in% c("G3", "G4") & metainfo$sample_id %in% allx_list[["flt3"]]], " FLT3-ITD")
metainfo$label[metainfo$label %in% "NA FLT3-ITD"] <- "FLT3-ITD"
metainfo$label[metainfo$rna_group %in% c("G11", "G17") & metainfo$sample_id %in% allx_list[["NUP98"]]] <- "aml_NUP98-r"
metainfo$label[metainfo$rna_group %in% c("G18") & metainfo$sample_id %in% allx_list[["CEBPA"]]] <- "aml_CEBPA"
metainfo$label[metainfo$rna_group %in% c("G17") & metainfo$sample_id %in% allx_list[["CEBPD"]]] <- paste0(
  metainfo$label[metainfo$rna_group %in% c("G17") & metainfo$sample_id %in% allx_list[["CEBPD"]]], " CEBPD")
metainfo$label[metainfo$label %in% "NA CEBPD"] <- "CEBPD mutation"
metainfo$label[metainfo$rna_group %in% c("G20") & metainfo$sample_id %in% allx_list[["COL6A5"]]] <- "COL6A5 mutation"


list2excel(list(sample_anno = data.frame(metainfo)), "/cluster/home/yjliu_jh/projects/leukemia/output/aml_samples_anno_curr.xlsx")











# zt code

pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes",
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "circlize", "ComplexHeatmap", "SummarizedExperiment", "jhuanglabRNAseq",
          "viridis","ggrepel","limma")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "leukemia"
dataset <- "meta"
species <- "human"
workdir <- glue("~/projects/{project}/data")
setwd(workdir)
seeds <- 2023

pth1 <- "/cluster/home/yjliu_jh/projects/leukemia/analysis/"
pth2 <- "/human/rnaseq/mutations/"
datasets <- list.files("~/projects/leukemia/analysis", full.names = F)
pth_dataset <- paste0(pth1, datasets, pth2, datasets,"/maf")
datasets <- datasets[file.exists(pth_dataset)]
all_path <- list.files(pth_dataset, full.names = T, pattern = "fil", recursive = T)

safe_read_tsv <- safely(read_tsv)
safe_read <- function(pth_i){
  tmpp <- str_split(pth_i, pattern = "/")
  dataset <- tmpp[[1]][8]
  sample <- tmpp[[1]][14]
  tools <- str_extract(tmpp[[1]][15], pattern = "(?<=fil.).+(?=.maf.gz)")
  tmpcsv <- safe_read_tsv(pth_i, show_col_types = F, progress = F, col_types = rep("c", times = 200) %>% paste0(collapse = ""))
  if(is.null(tmpcsv$error)){
    out <- tmpcsv$result
    out$dataset <- dataset
    out$sample <- sample
    out$tools <- tools
    out <- out %>% dplyr::select(dataset, sample, tools, everything())
  }else{
    out <- NULL
  }
  out
}

samplesheet <- tibble(
  full_name = all_path
) %>%
  tidyr::separate_wider_delim(delim = "/", cols = "full_name", names = 1:15 %>% as.character(), too_many = "drop") %>%
  dplyr::select(8, 14, 15)
colnames(samplesheet) <- c("datasets", " samples", "tools")
samplesheet$tools <- str_extract(samplesheet$tools, pattern = "(?<=fil.).+(?=.maf.gz)")
samplesheet$path <- all_path
all_list <- parallel::mclapply(all_path, safe_read, mc.cores = 50L)
all_frame <- all_list %>% list_rbind()


sel_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "HGVSp_Short",
              "Variant_Classification", "Variant_Type", "Consequence", "PUBMED", "FILTER",
              "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
              "SIFT", "PolyPhen", "IMPACT", "CLIN_SIG", "Existing_variation", "Protein_position")
all_fil <- all_frame[, c("dataset", "sample", "tools", sel_cols)]

all_fil_small <- all_fil[, c("dataset", "sample", "tools", "Hugo_Symbol", "HGVSc", "HGVSp_Short")]
all_fil_small$tools[all_fil_small$tools %in% "HaplotypeCaller"] <- "HC"
all_fil_small$tools[all_fil_small$tools %in% "snp.varscan2"] <- "VS"
all_fil_small$tools[all_fil_small$tools %in% "speedseq"] <- "SS"
all_fil_small <- all_fil_small %>% dplyr::distinct()
afs <- as.data.table(all_fil_small)[, toString(tools), by = setdiff(names(all_fil_small), "tools")]
colnames(afs)[6] <- "tools"
afs$sample_id <- paste0(afs$dataset, "@", afs$sample)
afs <- afs[, -c(1:2)]
readr::write_rds(afs, "/cluster/home/yjliu_jh/projects/temp/all_mutations_0716.rds")





fil <- cld$disease_type %in% c("AML")
se <- rds[, fil]
metadata <- colData(se)
fpkm <- assay(se, "tpm_remove_rnatype") %>% as.data.frame() %>% dplyr::select(metadata$sample_id) # exp matrix
dat <- fpkm |> as.data.frame() |> tibble::rownames_to_column("gene_name") |>
  filter_black_heatmap_genes() |> tibble::remove_rownames() |> tibble::column_to_rownames("gene_name")

meta_plot2 <- as.data.frame(metadata)[, c("gender", "age", "rna_type", "disease_type", "datasets")]
meta_plot2$sample_id <- rownames(meta_plot2)
meta_plot2$age <- as.numeric(meta_plot2$age)




# according to the proteinpaint tool:
idh1_hotspots <- c("p.R132H", "p.R132C", "p.R132G", "p.R132S", "p.R132L", "p.R132V",
                   "p.R132I", "p.P33S", "p.R100Q")
idh2_hotspots <- c("p.R140Q", "p.R172K", "p.R172S", "p.R172M", "p.R172G", "p.R172W",
                   "p.R140W", "p.R140L", "p.R140G")
asxl1_hotspots <- c("p.G645fs", "", "p.R634_I641fs", "p.R693*", "p.L19_E1splice_region",
                    "p.Y591fs", "p.Y591*", "p.K85R")


idh1_mut <- afs[afs$Hugo_Symbol %in% "IDH1" & afs$HGVSp_Short %in% idh1_hotspots, ]
idh2_mut <- afs[afs$Hugo_Symbol %in% "IDH2" & afs$HGVSp_Short %in% idh2_hotspots, ]
asxl1_mut <- afs[afs$Hugo_Symbol %in% "ASXL1" & afs$HGVSp_Short %in% asxl1_hotspots, ]



mut3_list <- list()
mut3_list[["IDH1"]] <- unique(idh1_mut$sample_id)
mut3_list[["IDH2"]] <- unique(idh2_mut$sample_id)
mut3_list[["ASXL1"]] <- unique(asxl1_mut$sample_id)
gene3 <- c("IDH1", "IDH2", "ASXL1")
meta_plot3 <- enrich_anno(meta_plot2, mut3_list, gene3)










clusterings <- readxl::read_excel("/cluster/home/jhuang/projects/leukemia/analysis/meta/human/figures/heatmap/aml/step1/sampleinfo_0.95.xlsx")
clusterings <- as.data.frame(clusterings)
colnames(clusterings)[2] <- "rna_group"
meta_plot3 <- left_join(meta_plot3, clusterings[, 1:2])
groups <- unique(meta_plot3$rna_group)

# set annotation   ## functions see heatmap_robustness_test.R
ha_x <- HeatmapAnnotation(
  df = as.data.frame(meta_plot3[, c("gender", "age", "datasets", gene3, "rna_group")]),
  show_legend = c(TRUE, TRUE, FALSE, rep(T, length(gene3))),
  col = append(list(age = col_fun, 
                    gender = setNames(anno_def[1:2], c("female", "male")),
                    rna_group = setNames(all_colors[1:length(groups)], groups)),
               auto_anno(meta_plot3, gene3, c("magenta", "ivory")))
)

# dat
vadj <- matrixStats::rowVars(as.matrix(dat))
dat_exp_h <- cbind(dat, vadj)
var_adj <- as.data.frame(dat_exp_h) %>% dplyr::arrange(desc(vadj))
dat_mat <- var_adj %>% dplyr::select(-last_col()) %>% t() %>% scale() %>% t() %>% 
  as.matrix() 

# plot using the quick function
out_path <- "/cluster/home/yjliu_jh/projects/temp/test3genesheatmap.pdf"
test_heat(dat_mat, ha_x, out_path)


quick_heatmap(dat, ha = ha_x, height = 9, outdir = "/cluster/home/yjliu_jh/projects/temp/htm", column_split = 26, top_var_percent = 0.95)



