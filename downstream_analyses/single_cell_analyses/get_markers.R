packages = c("data.table", "optparse", "Seurat", "Matrix", "stringr", "dplyr", "ggplot2", "presto")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
    make_option(c("-i", "--input_path"), type="character", default=NULL, help="Path to MTX Matrix", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default="NA", help="Tissue", metavar="character"),
    make_option(c("-f", "--tau_filter"), type="logical", default=TRUE, help="Filter for tau; default T", metavar="character"),
    make_option(c("-x", "--tau_threshold"), type="numeric", default=0.7, help="TAU threshold", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

tissue = opt$tissue
samples = readLines(opt$input_path)

#Gene Sets
genenames = data.frame(gene_id = sapply(str_split(readLines(paste0("./data/", samples[0], "/kb_out_nac/counts_unfiltered/cells_x_genes.genes.txt")), "[.]"), head, 1),
                     gene_name = sapply(str_split(readLines(paste0("./data/", samples[0], "/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt")), "[.]"), head, 1))

cls = readLines("cls.loci.disjoint.ids")
cls = sapply(str_split(cls, "[.]"),"[",1)

lncrna = readLines("lncRNA.loci.v27.disjoint.ids")
lncrna = sapply(str_split(lncrna, "[.]"),"[",1)

pcg = readLines("protein_coding.loci.v47.disjoint.ids")
pcg = sapply(str_split(pcg, "[.]"),"[",1)

lncrna = genenames[ genenames$gene_id %in% lncrna, "gene_name"]
cls = genenames[ genenames$gene_id %in% cls, "gene_name"]
pcg = genenames[ genenames$gene_id %in% pcg, "gene_name"]

if(opt$tau_filter) { 
    tau = read.table(paste0("./cluster_mapping/tau_analysis/stats/tau.", tissue ,".tsv"), header = T, row.names = NULL) 
    pass = tau[tau$tau >= 0.7 & tau$min_cluster > 2, "sample"]
    discarded = length(samples)-length(pass)

    print(paste0("Low TAU or 2 clusters: ", discarded ))
    samples = samples[samples %in% pass]
}

cluster_info = read.table("../postprocess/allclusters_flag.tsv", header = T, sep = "\t", row.names = NULL) 
clusters_info$ids = paste(clusters_info$tissue, clusters_info$sample, clusters_info$cluster, sep = "_")
bad_clusters = clusters_info[clusters_info$mitochondrial.flag | clusters_info$malat.flag | clusters_info$intronic.flag, "ids"]

#Run for each sample pf the same tissue
for (sample in samples) {

    input = paste0("./matrices_rds/", sample, "_v47Novel.rds")
    path = paste0("./data/", sample)

    if ( ! file.exists(input) ) { print(paste0("Skipping ", sample)); next; }

    tryCatch({
        #Get matrix
        matrixv47 = readRDS(input)
        clusters = paste(tissue, sample, unique(Idents(matrixv47)), sep = "_")
        clusters = clusters[! clusters %in% bad_clusters]
        
        if ( length(clusters) <= 2) { print(paste0("Skip - few cluster: ", sample)); next; }
        markers = FindAllMarkers(matrixv47, logfc.threshold = 0.5, min.pct = 0.2, return.thresh = 0.05) 
        markers$cluster = paste(tissue, sample, markers$cluster, sep = "_")
        #Remove flagged clusters
        markers = markers[markers$p_val_adj < 0.05 & ! markers$cluster %in% bad_clusters, ]
        markers$biotype = ifelse(markers$gene %in% cls, "CLS", ifelse(markers$gene %in% lncrna, "lncRNAv27", ifelse(markers$gene %in% pcg, "protein-coding", "Other")))
        write.table(markers, paste0("./output/", sample, "_markers.tsv"), quote = F, sep = "\t", row.names = F)

    }, error = function(e) {
        cat("An error occurred:\n", e$message, "\n")
    })

}
