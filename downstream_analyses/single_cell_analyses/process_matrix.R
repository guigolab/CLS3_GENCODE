#IMPORT PACKAGES
packages = c("data.table", "optparse", "Seurat", "Matrix", "stringr", "dplyr", "DoubletFinder","ggplot2", "future", "grid", "gridExtra")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
    make_option(c("-i", "--input_path"), type="character", default=NULL, help="Path to MTX Matrix", metavar="character"),
    make_option(c("-t", "--tissue"), type="character", default="NA", help="Tissue", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
set.seed(1234)

#SET ENVIRONMENT
plan("multisession", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) #Gives error, so going to 10GiB.

#SET INPUTS
if (!grepl(".tsv", opt$input_path)) { samples = c(opt$input_path) } else { samples = readLines(opt$input_path) }
tissue = opt$tissue

#Run for each sample of the same tissue
for (sample in samples) {
    print(sample)
    path = paste0("./cca_human_v47/human/data/", sample)

    tryCatch({
        #Get expression matrix
        expression_matrix = readMM(file = paste0(path, "/kb_out_nac/counts_filtered/cells_x_genes.total.mtx"))
        expression_matrix = Matrix(t(expression_matrix), sparse = TRUE)
        colnames(expression_matrix) = readLines(paste0(path, "/kb_out_nac/counts_filtered/cells_x_genes.barcodes.txt"))
        rownames(expression_matrix) = sapply(str_split(readLines(paste0(path, "/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt")), "[.]"), head, 1)
        matrixv47 = CreateSeuratObject(counts = expression_matrix)

        #Get simulated expression matrix without lncRNAs
        expression_matrix_nolnc = readMM(file = paste0(path, "/kb_out_nolnc/cells_x_genes.total.mtx"))
        expression_matrix_nolnc = Matrix(t(expression_matrix_nolnc), sparse = TRUE)
        colnames(expression_matrix_nolnc) = readLines(paste0(path, "/kb_out_nolnc/cells_x_genes.barcodes.txt"))
        rownames(expression_matrix_nolnc) = sapply(str_split(readLines(paste0(path, "/kb_out_nolnc/cells_x_genes.genes.names.txt")), "[.]"), head, 1)
        matrixv47_nolnc = CreateSeuratObject(counts = expression_matrix_nolnc)

        #Apply QC: Mitochrondrial percentage
        matrixv47 = PercentageFeatureSet(matrixv47, pattern = "^MT-", col.name = "percent.mitochondrial")
        matrixv47_nolnc = PercentageFeatureSet(matrixv47_nolnc, pattern = "^MT-", col.name = "percent.mitochondrial")

        ## INVESTIGATION OF LOW-QUALITY CELLS ######################################################
        test = matrixv47
        filter_feature = median(test@meta.data$nFeature_RNA, na.rm = TRUE) + 3*mad(test@meta.data$nFeature_RNA, na.rm = TRUE)
        test@meta.data$raw_malat1 = as.vector(test[["RNA"]]$counts["MALAT1", ])

        filter_mitochondrial = median(test@meta.data$percent.mitochondrial, na.rm = TRUE) + 3*mad(test@meta.data$percent.mitochondrial, na.rm = TRUE)
        filter_mitochondrial = ifelse(filter_mitochondrial < 5, 5, ifelse(filter_mitochondrial > 80, 80, filter_mitochondrial))

        test = subset(test, subset = nFeature_RNA >= 200 & nFeature_RNA <= filter_feature & percent.mitochondrial <= filter_mitochondrial)

        #Default SCTrasform, top 3000 most variable features
        test = SCTransform(test, vars.to.regress = "percent.mitochondrial", verbose = F, conserve.memory = TRUE)
        test@meta.data$normal_malat1 = as.vector(test[["SCT"]]$data["MALAT1", ])

        malat_filter = ifelse( tissue == "blood", 0, median(test[["SCT"]]$data["MALAT1", ], na.rm = TRUE) - 5*mad(test[["SCT"]]$data["MALAT1", ], na.rm = TRUE))
        test = subset(test, subset = normal_malat1 >= malat_filter)

        test = RunPCA(test, verbose = FALSE)
        test = FindNeighbors(test, dims = 1:30, verbose = FALSE)
        test = FindClusters(test, verbose = FALSE)
        test = RunUMAP(test, dims = 1:30, verbose = FALSE)

        # pK Identification
        sweep.res.test = paramSweep(test, PCs = 1:30, sct = T)
        sweep.stats_test = summarizeSweep(sweep.res.test, GT = FALSE)
        bcmvn = find.pK(sweep.stats_test)
        # Optimal pK is the max of the bimodality coefficient (BCmvn) distribution
        optimal.pk = bcmvn %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
        optimal.pk = as.numeric(as.character(optimal.pk[[1]]))

        # Homotypic Doublet Proportion Estimate
        homotypic.prop = modelHomotypic(test@meta.data$ClusteringResults)
        nExp_poi = round(optimal.pk*nrow(test@meta.data)) 
        nExp_poi.adj = round(nExp_poi*(1-homotypic.prop))

        # Run DoubletFinder with varying classification stringencies
        test = doubletFinder(test, PCs = 1:10, pN = 0.25, pK = optimal.pk, nExp = nExp_poi.adj, sct = T)
        colnames(test@meta.data)[grepl('DF.classifications.*', colnames(test@meta.data))] = "doublet_finder"

        test = subset(test, subset = doublet_finder == "Singlet")

        pass_barcodes = colnames(test)
        ################################################################

        ## FILTERED MATRIX; GO WITH ANALYSES ####################################################################
        matrixv47 = subset(matrixv47, cells = pass_barcodes)
        matrixv47_nolnc = subset(matrixv47_nolnc, cells = pass_barcodes)

        #Select most variable features as requested: Default SCTrasform, top 3000 most variable features
        matrixv47 = SCTransform(matrixv47, vars.to.regress = "percent.mitochondrial", verbose = F, conserve.memory = TRUE)
        matrixv47_nolnc = SCTransform(matrixv47_nolnc, vars.to.regress = "percent.mitochondrial", verbose = F, conserve.memory = TRUE)

        variable_genes = data.frame(residual_variance = matrixv47@assays$SCT@SCTModel.list$counts@feature.attributes$residual_variance)
        rownames(variable_genes) = rownames(matrixv47)
        variable_genes = variable_genes[order(variable_genes$residual_variance, decreasing = T), , drop = F]
        variable_genes$rank = seq(1:nrow(variable_genes))
        write.table(variable_genes, paste0("./matrices_rds/", sample, "_variableFeatures.tsv"), row.names = T, col.names = T, quote = F, sep ="\t")

        #Apply QC: MALAT1 Expression
        matrixv47@meta.data$normal_malat1 = as.vector(matrixv47[["SCT"]]@data["MALAT1", ])
        matrixv47_nolnc@meta.data$normal_malat1 = as.vector(matrixv47_nolnc[["SCT"]]@data["MALAT1", ])

        # These are now standard steps in the Seurat workflow for visualization and clustering
        matrixv47 = RunPCA(matrixv47, verbose = FALSE)
        matrixv47 = FindNeighbors(matrixv47, dims = 1:30, verbose = FALSE)
        matrixv47 = FindClusters(matrixv47, verbose = FALSE)
        matrixv47 = RunUMAP(matrixv47, dims = 1:30, verbose = FALSE)

        matrixv47_nolnc = RunPCA(matrixv47_nolnc, verbose = FALSE)
        matrixv47_nolnc = FindNeighbors(matrixv47_nolnc, dims = 1:30, verbose = FALSE)
        matrixv47_nolnc = FindClusters(matrixv47_nolnc, verbose = FALSE)
        matrixv47_nolnc = RunUMAP(matrixv47_nolnc, dims = 1:30, verbose = FALSE)

        saveRDS(matrixv47, paste0("./matrices_rds/", sample, "_v47Novel.rds"))
        saveRDS(matrixv47_nolnc, paste0("./matrices_rds/", sample, "_v47woNovel.rds"))

    }, error = function(e) {
        print(e$message)
    }, finally = { 
        dev.off() 
    })
}
