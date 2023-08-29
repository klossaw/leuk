

# === code from zt ===
library(infercnv)

# ========================= infercnv =============================



mef2d <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merge/mef2d.rds")

sce_tumor <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/merged_4_tumors.rds")
sce_hd <- readr::read_rds("/cluster/home/yjliu_jh/projects/mef2d/data/healthy_donors.rds")
Idents(sce_tumor) <- "cell_types_broad"
sce_tumor_sub <- subset(sce_tumor, `cell_types_broad` %notin% c("NK_T", "erythroid_cell", "myeloid"))
Idents(sce_hd) <- "cell_types"
sce_hd_sub <- subset(sce_hd, `cell_types` %notin% c("NK_T", "erythroid_cell", "myeloid"))
rm(list = c("sce_tumor", "sce_hd"))
gc()
sce_hd_sub <- sce_hd_sub[, sample(colnames(sce_hd_sub), 8000)]

sce_all <- readr::read_rds("/cluster/home/jhuang/projects/mef2d/analysis/zwn/human/tenx/seurat/merged/all_merged.rds")
new_cell_id <- str_sub(colnames(sce_all), 1 + str_locate(colnames(sce_all), ":")[, 1] / 2)
sce_all <- RenameCells(sce_all, new.names = new_cell_id)

mef_id <- colnames(sce_tumor_sub)
hd_id <- colnames(sce_hd_sub)
sce_all_sub <- sce_all[, c(mef_id, hd_id)]
rm(list = c("sce_all"))
gc()


case <- unique(sce_tumor_sub$orig.ident)
control <- unique(sce_hd_sub$orig.ident)

sce_all_sub@meta.data$celltype <- "healthy"
sce_all_sub$celltype[sce_all_sub$orig.ident %in% case] <- sce_tumor_sub$cell_types_broad


raw_counts_matrixcounts <- as.matrix(sce_all_sub@assays$RNA@counts)


gencode42 <- rtracklayer::import("/cluster/home/yjliu_jh/share/gencode.v42.annotation.gtf")
total <- as.data.frame(gencode42)
total <- total[, c("gene_name", "seqnames", "start", "end")]
total <- total[!duplicated(total$gene_name), ]
rownames(total) <- total$gene_name

odgenes <- intersect(rownames(total), rownames(raw_counts_matrixcounts))
raw_counts_matrixcounts <- raw_counts_matrixcounts[odgenes, ]
total <- total[match(odgenes, rownames(total)), c("seqnames", "start", "end")]

write.table(total, file = '/cluster/home/yjliu_jh/share/ref_gencodev42_slim_all.txt',
            row.names = T, col.names = F, quote = F, sep = "\t")



data <- data.frame(V1 = names(sce_all_sub@active.ident), V2 = sce_all_sub@meta.data$celltype)

write.table(data, file = '/cluster/home/yjliu_jh/share/tempanno_all_cells.txt',
            row.names = F, col.names = F, quote = F, sep = '\t')

# 开始构建cnv对象
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrixcounts,
                                    annotations_file = "/cluster/home/yjliu_jh/share/tempanno_all_cells.txt",
                                    delim = "\t",
                                    gene_order_file = "/cluster/home/yjliu_jh/share/ref_gencodev42_slim_all.txt",
                                    ref_group_names = c('proB'))

readr::write_rds(infercnv_obj, "/cluster/home/yjliu_jh/share/infercnv_obj_all.rds")



infercnv_obj <- readr::read_rds("/cluster/home/yjliu_jh/share/infercnv_obj_all.rds")
options(scipen = 100)
infercnv_obj2 = infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir = "/cluster/home/yjliu_jh/projects/mef2d/output/drug/cnv",  
                             denoise = T, 
                             HMM = T,
                             cluster_by_groups = T,
                             tumor_subcluster_partition_method = "random_trees",
                             output_format = "pdf",
                             num_threads = 50) 

readr::write_rds(infercnv_obj2, "/cluster/home/yjliu_jh/share/infercnv_obj.rds")


# 手动读入个cnv的结果看看
run.final.infercnv <- readRDS('/cluster/home/ztao_jh/projects/mef2d/analysis/zwn/human/output_dir_Bcell/run.final.infercnv_obj')
# run.final.infercnv@count.data |> density() |> plot()
run.final.infercnv@expr.data |> density() |> plot()

# Heatmap(run.final.infercnv@expr.data,
#         cluster_rows = F,cluster_columns = F)

# 提取cnv值计算平均cnv
expr.data <- run.final.infercnv@expr.data
celldata <- data.frame(
  cells = colnames(expr.data),
  means = apply(expr.data,2,mean)
) 



