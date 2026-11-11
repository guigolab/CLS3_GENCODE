packages = c("optparse", "Seurat", "stringr")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Sample list", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default=NULL, help="Tissue tag", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

taus = as.data.frame(matrix( ncol = 3, nrow = 0))

mappings = paste0("../stats/", opt$tissue, ".mappings.tsv")
if ( ! file.exists(mappings) ) { print(paste0("Skip: ", opt$tissue)); exit; }
cluster_map = read.table(mappings, header = T, row.names = NULL)

for ( i in readLines(opt$input) )
{
    print(i)
    path_novel = paste0("./matrices_rds/", i, "_v47Novel.rds")
    path_old = paste0("./matrices_rds/", i, "_v47woNovel.rds")

    if ( ! file.exists(path_novel) | ! file.exists(path_old) | ! i %in% cluster_map$sample ) { print(paste0("Skip: ", i)); next; }

    current_map = cluster_map[cluster_map$sample == i, ]
    current_map$old_cluster = sapply(str_split(current_map$old_cluster, "[_]"), "[", 2)
    current_map$new_cluster = sapply(str_split(current_map$new_cluster, "[_]"), "[", 2)
    min_cluster = min(length(current_map$new_cluster), length(current_map$old_cluster))
    if ( min_cluster < 3 ) { print(paste0("Skip: ", i, ". Few samples.")); }
    
    matrix_novel = readRDS(path_novel)
    matrix_woNovel = readRDS(path_old)

    clusters_novel = as.data.frame(matrix_novel$seurat_clusters)
    clusters_novel = split(rownames(clusters_novel), clusters_novel[,1])

    clusters_woNovel = as.data.frame(matrix_woNovel$seurat_clusters)
    clusters_woNovel = split(rownames(clusters_woNovel), clusters_woNovel[,1])

    concordant = 0
    for ( x in names(clusters_novel) )
    {
        match = current_map[current_map$new_cluster == x, "old_cluster"]
        concordant = concordant + length( intersect( unlist(clusters_woNovel[match]), unlist(clusters_novel[x]) ) )   
    }

    tau = (concordant - (ncol(matrix_novel) - concordant)) / ncol(matrix_novel)
    taus = rbind(taus, c(i, tau, min_cluster))
}    

colnames(taus) = c("sample", "tau", "min_cluster")
write.table(taus, paste0("./stats/tau.",opt$tissue,".tsv"), quote = F, row.names = F, sep = "\t")
