packages = c("data.table", "optparse", "Seurat", "ggplot2", "stringr", "ggpubr", "dplyr", "pheatmap", "grid", "gridExtra")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Sample list", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default=NULL, help="Tissue tag", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

pdf(paste0("./output/Jaccard_",opt$tissue,".pdf"), width = 15)
for ( sample in readLines(opt$input) )
{
    path = paste0("./matrices_rds/", sample, "_v47Novel.rds")
    path_old = paste0("./matrices_rds/", sample, "_v47woNovel.rds")

    if ( ! file.exists(path) | ! file.exists(path_old)) { print(paste0("Skip: ", sample)); next; }

    matrix = readRDS(path)
    matrix_old = readRDS(path_old)

    # Extract cluster labels
    tmp = as.data.frame(matrix$seurat_clusters)
    new_cluster_labels = data.frame(cell = rownames(tmp), cluster_new = tmp[,1])
    new_cluster_list = split(rownames(tmp), tmp[,1])

    tmp = as.data.frame(matrix_old$seurat_clusters)
    old_cluster_labels = data.frame(cell = rownames(tmp), cluster_old = tmp[,1])
    old_cluster_list = split(rownames(tmp), tmp[,1])

    ## ROWS - With Novel, COLUMNS - Without novel
    jaccard = as.data.frame(matrix(nrow = length(names(new_cluster_list)), ncol = length(names(old_cluster_list))))
    colnames(jaccard) = names(old_cluster_list)
    rownames(jaccard) = names(new_cluster_list)

    intersection = as.data.frame(matrix(nrow = length(names(new_cluster_list)), ncol = length(names(old_cluster_list))))
    colnames(intersection) = names(old_cluster_list)
    rownames(intersection) = names(new_cluster_list)
    size_new = as.data.frame(matrix(ncol = 2, nrow = nrow(intersection)))
    colnames(size_new) = c("Cluster", "Size")
    rownames(size_new) = rownames(intersection)

    size_old = as.data.frame(matrix(ncol = 2, nrow = ncol(intersection)))
    colnames(size_old) = c("Cluster", "Size")
    rownames(size_old) = colnames(intersection)

    for ( c_old in names(old_cluster_list) )
    {
        for ( c_new in names(new_cluster_list) )
        {
            a = unlist(new_cluster_list[c_new])
            b = unlist(old_cluster_list[c_old])
            jaccard[c_new, c_old] = round( length(intersect(a, b)) / length(union(a,b)), 2 )
            intersection[c_new, c_old] = round( length(intersect(a, b)), 0 )

            size_new[c_new, ] = c(c_new, length(a))
            size_old[c_old, ] = c(c_old, length(b))
        }
    }

    colnames(jaccard) = paste0("Without_",colnames(jaccard))
    rownames(jaccard) = paste0("Novel_",rownames(jaccard))
    colnames(intersection) = paste0("Without_",colnames(intersection))
    rownames(intersection) = paste0("Novel_",rownames(intersection))

    write.csv(jaccard, paste0("./stats/",sample,"_jaccard.csv"))
    write.csv(intersection, paste0("./stats/",sample,"_intersection.csv"))

    size_old$tag = "OLD"
    size_new$tag = "NEW"
    sizes = rbind(size_new, size_old)
    write.table(sizes, paste0("./stats/",sample,"_clusterSize.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

    if ( length(names(old_cluster_list)) < 3 | length(names(new_cluster_list)) < 3 ) { next; }
    j = pheatmap(jaccard, display_numbers = T, cluster_rows = FALSE, cluster_cols = FALSE, main = "Jaccard Index", silent = TRUE, legend = F)
    i = pheatmap(intersection, display_numbers = round(intersection, 0), cluster_rows = FALSE, cluster_cols = FALSE, main = "Intersection Size", silent = TRUE, legend = F)

    size_new$Cluster = factor(size_new$Cluster, levels = seq(0,nrow(size_new)))
    s1 = ggplot(size_new, aes( x = Cluster, y = as.numeric(Size))) + geom_bar(stat = "identity", fill = "lightblue" ) + ggtitle(" Clustering solution - With Novel lncRNAs") + ylab("#Cells") + geom_text(aes(label = Size), vjust = -0.3, size = 3)

    size_old$Cluster = factor(size_old$Cluster, levels = seq(0,nrow(size_old)))
    s2 = ggplot(size_old, aes( x = Cluster, y = as.numeric(Size))) + geom_bar(stat = "identity", fill = "lightgray" ) + ggtitle(" Clustering solution - Without Novel lncRNAs") + ylab("#Cells") + geom_text(aes(label = Size), vjust = -0.3, size = 3)

    grid.arrange(j[[4]], i[[4]], ggarrange(s1, s2, ncol = 1, nrow = 2), ncol = 3, nrow = 1, top = textGrob(sample))
    
}
dev.off()



cl_maps = as.data.frame(matrix(ncol = 4, nrow = 0))

for ( sample in readLines(opt$input) )
{
    path=paste0("./stats/",sample,"_jaccard.csv")
    if ( ! file.exists(path) ) { print(paste0("Skip: ", sample)); next; }

    jaccard = read.csv(path, header = T, row.names = 1)
    jaccard$new_cluster = rownames(jaccard)
    jaccard = melt(jaccard, value.name = "jaccard_index", id = "new_cluster", variable.name = "old_cluster")

    jaccard = jaccard[order(jaccard$jaccard_index, decreasing = T), ]
    jaccard$new_cluster = as.character(jaccard$new_cluster)
    jaccard$old_cluster = as.character(jaccard$old_cluster)

    while ( nrow(jaccard) > 0 )
    {
        new = jaccard[1, "new_cluster"]
        old = jaccard[1, "old_cluster"]
        cl_maps = rbind(cl_maps, c(new, old, jaccard[1, "jaccard_index"], sample))

        jaccard = jaccard[jaccard$new_cluster != new, ]
    }
    
}
colnames(cl_maps) = c("new_cluster", "old_cluster", "jaccard_index", "sample")
write.table(cl_maps, paste0("./stats/", opt$tissue, ".mappings.tsv"), row.names = F, col.names = T, quote = F, sep = "\t")
