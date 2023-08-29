pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr",
          "lubridate", "ComplexHeatmap", "maftools", "DESeq2")
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}

rds_fn <- "~/projects/leukemia/analysis/meta/human/rnaseq/exp/tables/leukemia.rds"
rds <- readr::read_rds(rds_fn)
cld <- colData(rds) |> as.data.frame() 
samples <- cld[, c("patient_id", "dataset")]

mafloc <- "/cluster/home/jhuang/projects/leukemia/analysis"




# ====== subset and merge MAFs for each tool ======

# dv <- list()
# for (i in 1:length(samples)){
#   tryCatch({
#     temp <- read.maf(glue("{mafloc}/{samples[i]}/{samples[i]}.DeepVariant.maf"))
#     temp_3 <- temp@data[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
#                             "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
#                             "Variant_Classification", "Variant_Type", "Consequence", "PUBMED", "FILTER",
#                             "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
#                             "SIFT", "PolyPhen", "IMPACT", "CLIN_SIG", "Existing_variation", "Protein_position")]
#     colnames(temp_3)[22] <- "VEP_Impact"
#     dv[[samples[i]]] <- temp_3
#   }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# }
# readr::write_rds(dv, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/DeepVariant.rds")

sel_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
              "Variant_Classification", "Variant_Type", "Consequence", "PUBMED", "FILTER",
              "HGVSc", "EXON", "INTRON", "t_depth", "t_ref_count", "t_alt_count", "BIOTYPE",
              "SIFT", "PolyPhen", "IMPACT", "CLIN_SIG", "Existing_variation", "Protein_position")


speedseq <- list()
for(i in 1:nrow(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i, 2]}/human/rnaseq/mutations/{samples[i, 2]}/maf/{samples[i, 1]}/{samples[i, 1]}.speedseq.maf.gz"))
  temp_3 <- temp@data[, ..sel_cols]
  colnames(temp_3)[22] <- "VEP_Impact"
  speedseq[[samples[i, 1]]] <- temp_3
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
readr::write_rds(speedseq, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/speedseq.rds")

hc <- list()
for(i in 1:nrow(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i, 2]}/human/rnaseq/mutations/{samples[i, 2]}/maf/{samples[i, 1]}/{samples[i, 1]}.HaplotypeCaller.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    hc[[samples[i, 1]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
readr::write_rds(hc, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/HaplotypeCaller.rds")

ri <- list()
for(i in 1:nrow(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i, 2]}/human/rnaseq/mutations/{samples[i, 2]}/maf/{samples[i, 1]}/{samples[i, 1]}.rnaindel.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    ri[[samples[i, 1]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
readr::write_rds(ri, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/rnaindel.rds")

vs2s <- list()
for(i in 1:nrow(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i, 2]}/human/rnaseq/mutations/{samples[i, 2]}/maf/{samples[i, 1]}/{samples[i, 1]}.snp.varscan2.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    vs2s[[samples[i, 1]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
readr::write_rds(vs2s, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/snp.varscan2.rds")

vs2i <- list()
for(i in 1:nrow(samples)){
  tryCatch({
    temp <- read.maf(glue("{mafloc}/{samples[i, 2]}/human/rnaseq/mutations/{samples[i, 2]}/maf/{samples[i, 1]}/{samples[i, 1]}.indel.varscan2.maf.gz"))
    temp_3 <- temp@data[, ..sel_cols]
    colnames(temp_3)[22] <- "VEP_Impact"
    vs2i[[samples[i, 1]]] <- temp_3
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
readr::write_rds(vs2i, "/cluster/home/yjliu_jh/projects/leu_j/data/2023/indel.varscan2.rds")




# ====== further filtering process ======

mut_ri <- dplyr::bind_rows(ri, .id = 'Samples')
mut_vs2s <- dplyr::bind_rows(vs2s, .id = 'Samples')
mut_vs2i <- dplyr::bind_rows(vs2i, .id = 'Samples')
mut_ri$tool <- "RI"
mut_vs2s$tool <- "VS"
mut_vs2i$tool <- "VS"


hc <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/HaplotypeCaller.rds")
mut_hc <- dplyr::bind_rows(hc, .id = 'Samples')
dv <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/DeepVariant.rds")
mut_dv <- dplyr::bind_rows(dv, .id = 'Samples')
speedseq <- readr::read_rds("/cluster/home/yjliu_jh/projects/leu_j/data/speedseq.rds")
mut_ss <- dplyr::bind_rows(speedseq, .id = 'Samples')

mut_hc$tool <- "HC"
mut_ss$tool <- "SS"
mut_dv$tool <- "DV"


mut_int <- bind_rows(mut_hc, mut_ss, mut_dv, mut_ri, mut_vs2i, mut_vs2s)

system.time(temp <- as.data.table(mut_int)[, toString(tool), by = setdiff(names(mut_int), "tool")])
readr::write_rds(temp, "/cluster/home/yjliu_jh/projects/leu_j/data/merged.rds")


system.time(temp <- as.data.table(mut_int)[, toString(tool), by = eval(names(mut_int)[1:8])])
system.time(tempx <- left_join(temp, mut_int))
tempx <- as.data.frame(tempx)[, !(names(tempx) %in% "tool")]   ## -which() or ! not working for data.table
colnames(tempx)[colnames(df) %in% "V1"] <- "Tool"
readr::write_rds(tempx, "/cluster/home/yjliu_jh/projects/leu_j/data/merged_v2.rds")



temp1 <- tempx
temp1 <- temp1[temp1$FILTER %in% c("PASS", "."), ]

# subset based on sift and polyphen? 
temp1$SIFT <- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", temp1$SIFT, perl = T)
temp1$PolyPhen <- gsub("(?<=\\()[^()]*(?=\\))(*SKIP)(*F)|.", "", temp1$PolyPhen, perl = T)
temp1$SIFT <- as.numeric(temp1$SIFT)
temp1$PolyPhen <- as.numeric(temp1$PolyPhen)


rna_edit <- readr::read_delim("/cluster/home/yjliu_jh/projects/leu_j/data/TABLE1_hg38.txt")
#rna_edit$edit <- paste0(rna_edit$Ref, "to", rna_edit$Ed)

editdata <- rna_edit[, c(1:7)]
colnames(editdata)[1:3] <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
temp1 <- anti_join(temp1, editdata)





head(sort(table(unique(temp1[temp1$Hugo_Symbol %in% vitalgenes, 1:2])$Hugo_Symbol), decreasing = T), 14)



ind <- (temp1$SIFT <= 0.05 | temp1$PolyPhen >= 0.85 | temp1$VEP_Impact %in% "HIGH")
ind[is.na(ind)] <- FALSE
temp2 <- temp1[ind, ]

# readr::write_rds(temp2, "/cluster/home/yjliu_jh/projects/leu_j/data/tempmut.rds")
sinfo2 <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/sampleinfo_20220809.xlsx")
groupings <- sinfo2[, 1:3]
colnames(groupings)[1] <- "Samples"

temp3 <- left_join(temp2, groupings)


quick_check <- function(dat, group){
  head(sort(table(unique(dat[dat$Hugo_Symbol %in% vitalgenes & dat$subgroups %in% group, 1:2])$Hugo_Symbol), decreasing = T), 20)
}
quick_check(temp3, "G7")


# G1 characteristic CEBPA    G2 AML-ETO   G3 PML-RARA   G4 NPM1  G5 NPM1 KMT2A with FLT3-ITD
# G6 not obvious (NO)   G7 NPM1 KMT2A  similar to G5?   G8 CBFB-MYH11   G9 NO    G10 NO maybe BCR-ABL1
#
#



temp2 <- temp %>% dplyr::filter(nchar(V1) > 3)


hamlet_leu <- read_csv("/cluster/home/yliang_jh/projects/leukemia_linxiangjie/hamlet/snv-indels_variants_all.csv")
hamlet_leu_fil <- hamlet_leu[is.na(hamlet_leu$filters), ]
hamlet_leu_fil <- hamlet_leu_fil[hamlet_leu_fil$IMPACT %notin% c("MODIFIER"), ]

hamlet_leu_fil <- left_join(hamlet_leu_fil, sinfo2[, 1:3])
readr::write_rds(hamlet_leu_fil, "/cluster/home/yjliu_jh/projects/leu_j/data/mut_all_hamlet.rds")

ham_fusions



head(sort(table(unique(hamlet_leu_fil[hamlet_leu_fil$subgroups %in% "G5", 1:2])$gene_symbol), decreasing = T), 20)



#

tempxx <- temp1[, 1:14] %>% group_by_at(1:10) %>% dplyr::filter(n() > 1)







ss3 <- temp2 %>% group_by(Samples, Hugo_Symbol, HGVSc) %>% summarise(count = n())
ss3$count <- 1  ## multiple transcripts of same gene? 
ss3 <- ss3 %>% group_by(Hugo_Symbol, HGVSc) %>% summarise(tcount = sum(count)) %>% 
  group_by(Hugo_Symbol) %>% mutate(allcount = sum(tcount), freq = tcount/allcount) 

flags <- read.delim("/cluster/home/yjliu_jh/projects/leu_j/data/FLAGS.txt")
ss4 <- ss3[ss3$freq < 0.8, ]
ss4 <- ss4[ss4$tcount < 100, ]
ss4 <- ss4 %>% group_by(Hugo_Symbol) %>% summarise(newcount = sum(tcount))
ss4 <- ss4[ss4$Hugo_Symbol %notin% unique(flags$FLAGS), ]

ssx <- ss4[order(ss4$newcount, decreasing = T), ]


cols_6 <- as.character(readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 
                                          6, skip = 3, n_max = 1, col_names = FALSE))
leu_drivers <- readxl::read_excel("/cluster/home/yjliu_jh/projects/leu_j/data/ref_ng.xlsx", 6, 
                                  skip = 5, col_names = cols_6)

head(match(leu_drivers$Gene, ssx$Hugo_Symbol), 180)



head(rna_edit[rna_edit$Gene.wgEncodeGencodeBasicV34 %in% "MECOM", ])


head(hamlet_dav[, c("")])

editdata <- rna_edit[, c(1:7)]
colnames(editdata)[1:3] <- c("CHROM", "POS", "REF")

temph <- inner_join(hamlet_dav_fil[, 1:10], editdata)
head(as.data.frame(temph))
temph$flag <- mapply(grepl, temph$Ed, temph$alleles)
temph <- temph[temph$flag == TRUE, 1:10]
tempa <- anti_join(hamlet_dav_fil, temph)
tempa1 <- tempa[tempa$IMPACT %notin% c("MODIFIER"), ]

vitalgenes <- c("MECOM", "TET2", "DNMT3A", "IDH2", "RUNX1", "TP53", "NRAS", 
                "FLT3", "ASXL1", "CEBPA", "KIT", "NPM1", "WT1", "IDH1")

hamlet_res <- readr::read_rds("/cluster/home/yliang_jh/projects/mRNA/hamlet/leukemia/CRA000588/CRA000588.summary.rds")
hamlet_dav <- read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/hematology_dav/snv-indels_variants_all.csv")
hamlet_dav_fil <- hamlet_dav[is.na(hamlet_dav$filters), ]
head(sort(table(unique(hamlet_dav_fil[, 1:2])$gene_symbol), decreasing = T), 20)

fka <- read_csv("/cluster/home/yliang_jh/projects/mRNA/hamlet/test/Frankenstein/snv-indels/Frankenstein.variants_all.csv")



r <- unclass(lsf.str(envir = asNamespace("jhuanglabRNAseq"), all = T))


#r <- r[-grep("\\[", r)]
#r <- r[-grep("<-", r)]

# create functions in the Global Env. with the same name
# for(name in r) eval(parse(text=paste0(name, '<-jhuanglabRNAseq:::', name)))

dummydata <- data.frame(a = c(1, 3, 5, 1, 3, 5), 
                        b = c(1, 2, 3, 4, 2, 4),
                        c = c("A", "A", "A", "B", "C", "C"))
dummydata %>% group_by_at(1:2) %>% filter(n() > 1)




