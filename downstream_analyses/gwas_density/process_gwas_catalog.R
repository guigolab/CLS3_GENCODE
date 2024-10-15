library(dplyr)
library(optparse)
library(stringr)
library(tidyr)

file = read.delim("gwas.catalog.extHap.NAs.noInt.tsv", header = T, row.names = NULL, sep = "\t", quote = "")

file_merged = file[ file$MERGED == 1 & file$SNP_ID_CURRENT != "NA", ]
file_merged$SNP_ID_CURRENT = paste0("rs", file_merged$SNP_ID_CURRENT)
file_merged = file_merged[, c("CHR_ID","CHR_POS", "SNP_ID_CURRENT", "CONTEXT", "INTERGENIC", "MAPPED_TRAIT", "MAPPED_TRAIT_URI")]
colnames(file_merged)[3] = "SNPS"
dim(file_merged)

file_not_merged = file[ ! rownames(file) %in% rownames(file_merged), ]
file_not_merged = file_not_merged[, c("CHR_ID","CHR_POS", "SNPS", "CONTEXT", "INTERGENIC", "MAPPED_TRAIT", "MAPPED_TRAIT_URI")]
dim(file_not_merged)

data = rbind(file_merged, file_not_merged)
data$CHR_ID = paste0("chr", data$CHR_ID)
data = data[data$CHR_ID != "chrNA",]
data = unique(data)

print("Writing on file")
write.table(data, "../data/gwas_condensed_catalog.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
saveRDS(data,"../output/gwas_condensed_catalog.rds")

print("Get clustering")
clusters = read.delim("gwas.catalog.semantic_clustering.th0.8.bed7", header = F, row.names = NULL, sep = "\t", quote = "")
colnames(clusters) = c("CHR_ID","CHR_POSmin1","CHR_POS","SNPS","score","strand","efo","group","cluster")
data = merge(data[, c("CHR_ID","CHR_POS","SNPS")], clusters[,c("CHR_ID","CHR_POS","SNPS","cluster")], by = c("CHR_ID","CHR_POS","SNPS"), all.x = T)
write.table(unique(data), "gwas_coordinates_resolvedName_catalog.tsv", col.names = T, row.names = F, quote = F, sep = "\t")
saveRDS(unique(data),"gwas_coordinates_resolvedName_catalog.rds")

bedformat = data.frame("CHR_ID"=data$CHR_ID,"CHR_Start"=data$CHR_POS-1,"CHR_Stop"=data$CHR_POS,"SNP_ID_CURRENT"=data$SNPS,"TRAIT"=data$cluster )
write.table(unique(bedformat), ".gwas_coordinates_resolvedName_catalog.bed", col.names = F, row.names = F, quote = F, sep = "\t")
