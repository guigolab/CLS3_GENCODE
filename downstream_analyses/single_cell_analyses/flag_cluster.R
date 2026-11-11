packages = c("data.table", "optparse", "Seurat", "Matrix", "stringr", "dplyr", "DoubletFinder","ggplot2", "future", "grid", "gridExtra")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
    make_option(c("-i", "--input_path"), type="character", default=NULL, help="Path to MTX Matrix", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default="NA", help="Tissue", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

novel_lncrnas = readLines("./intergenic_genes.names.txt")
lncrnas = readLines("./lncRNAs_v47.noTEC.names.txt")
proteins = readLines("./pc_v47.names.txt")

if (!grepl(".tsv", opt$input_path)) { samples = c(opt$input_path) } else { samples = readLines(opt$input_path) }
tissue = opt$tissue

path_input = "/./matrices_rds/"
path_output = "./matrices_rds/tag_clusters/"

stats = as.data.frame(matrix(ncol = 13, nrow = 0))
colnames(stats) = c("sample","tissue","cluster","size","mitochondrial.median","mitochondrial.flag","malat.median","malat.flag","intronic.median","intronic.flag","cls.expressed","lnc.expressed","pc.expressed")

for (sample in samples) {
    print(sample)

    if ( ! file.exists(paste0(path_input, sample, "_v47Novel.rds")) ) 
    { 
        print("Missing file.")
        stats = rbind(stats, data.frame(sample = sample, tissue = tissue, cluster = NA, size = NA, mitochondrial.median = NA, mitochondrial.flag = NA,
                                        malat.median = NA, malat.flag = NA, intronic.median = NA, intronic.flag = NA, cls.expressed = NA, lnc.expressed = NA, pc.expressed = NA))
        next;
    }

    matrixv47 = readRDS(paste0(path_input, sample, "_v47Novel.rds"))
    intronic_cell = readRDS(paste0("./intronic_reads/", sample, "_nascent.bycell.rds"))

    malat_filter = ifelse( tissue == "blood", 0, as.numeric(quantile(matrixv47[["SCT"]]$data["MALAT1", ], 0.05)))
    mitoc_filter = as.numeric(quantile(matrixv47@meta.data$percent.mitochondrial, 0.95))
      intronic_filter = as.numeric(quantile(unlist(intronic_cell), 0.95))

    malat_avgs = mitoc_avgs = intrn_avgs = list()
    sizes = lnc = pc = cls = list()
    clusters = sort(unique(Idents(matrixv47)))

    for (cx in clusters)
    {
        cls_cells = WhichCells(matrixv47, idents = cx)

        malat1_expr = FetchData(matrixv47, vars = "MALAT1", cells = cls_cells)
        malat_avgs[[cx]] = median(malat1_expr$MALAT1, na.rm = T)

        mito_expr = FetchData(matrixv47, vars = "percent.mitochondrial", cells = cls_cells)
        mitoc_avgs[[cx]] = median(mito_expr$percent.mitochondrial, na.rm = T)
 
        sizes[[cx]] = length(cls_cells)

        intr_expr = intronic_cell[cls_cells]
        intrn_avgs[[cx]] = median(unlist(intr_expr), na.rm = T)

        counts = GetAssayData(matrixv47, layer = "counts")
        counts = counts[, cls_cells]
        counts = counts[ rowSums(counts >= 1) >= 0.05*ncol(counts) , ]  #Expressed in 5% of the cells
        genes = rownames(counts)

        lnc[[cx]] = sum(genes %in% lncrnas)
        cls[[cx]] = sum(genes %in% novel_lncrnas)
        pc[[cx]] =  sum(genes %in% proteins)
    }

    outliers_mlt = names(malat_avgs[malat_avgs < malat_filter])
    outliers_mtc = names(mitoc_avgs[mitoc_avgs > mitoc_filter])
    outliers_ntr = names(intrn_avgs[intrn_avgs > intronic_filter])

    tmp = data.frame( sample = sample, tissue = tissue, cluster = clusters, size = unlist(sizes), mitochondrial.median = unlist(mitoc_avgs), 
               mitochondrial.flag = ifelse(as.character(clusters) %in% outliers_mtc, T, F),
               malat.median = unlist(malat_avgs), malat.flag = ifelse(as.character(clusters) %in% outliers_mlt, T, F), intronic.median = unlist(intrn_avgs), 
               intronic.flag = ifelse(as.character(clusters) %in% outliers_ntr, T, F), cls.expressed = unlist(cls), lnc.expressed = unlist(lnc), pc.expressed = unlist(pc))
    stats = rbind(stats, tmp)
}

write.table(stats, paste0(path_output, tissue, ".clusterFlags.tsv"), row.names = F, quote = F, sep = "\t")
