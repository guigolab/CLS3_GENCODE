library(dplyr)
library(optparse)
library(stringr)
library(tidyr)

file = read.delim("gwas.catalog.efo.bed7", header = F, row.names = NULL, sep = "\t", quote = "")
colnames(file) = c("chr","start","end","rsid","score","strand","trait","EFO.URI")

print("Split catalog by term ID")
file = as.data.frame(file %>% mutate(EFO.URI = str_split(EFO.URI, ",")) %>% unnest(EFO.URI))
file$EFO.URI = sapply(str_split(file$EFO.URI,"/"),tail,1)

mapping = read.delim("trait_mappings", header = T, row.names = NULL, sep = "\t", quote = "")
mapping$EFO.URI = sapply(str_split(mapping$EFO.URI, "/"), tail, 1)
mapping = unique(mapping[, c("EFO.URI","EFO.term")])

file = merge(file, mapping, by = "EFO.URI", all.x = T)
file$trait = ifelse( !is.na(file$EFO.term), file$EFO.term, file$trait)
file = unique(file[, -c(1, ncol(file))])

write.table(file, "gwas.catalog.mappedefo.bed7", col.names = T, row.names = F, quote = F, sep = "\t")
