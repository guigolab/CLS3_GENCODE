packages = c("data.table", "optparse", "Seurat", "Matrix", "stringr", "dplyr")
invisible(lapply(packages, library, character.only = TRUE))
option_list = list(
    make_option(c("-i", "--input_path"), type="character", default=NULL, help="Path to MTX Matrix", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
set.seed(1234)

plan("multisession", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)

#SET INPUTS
if (!grepl(".tsv", opt$input_path)) { samples = c(opt$input_path) } else { samples = readLines(opt$input_path) }

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
        expression_matrix[is.na(expression_matrix)] = 0
        expression_matrix = colSums(expression_matrix)
        
        expression_matrix_nascent = readMM(file = paste0(path, "/kb_out_nac/counts_filtered/cells_x_genes.nascent.mtx"))
        expression_matrix_nascent = Matrix(t(expression_matrix_nascent), sparse = TRUE)
        colnames(expression_matrix_nascent) = readLines(paste0(path, "/kb_out_nac/counts_filtered/cells_x_genes.barcodes.txt"))
        rownames(expression_matrix_nascent) = sapply(str_split(readLines(paste0(path, "/kb_out_nac/counts_unfiltered/cells_x_genes.genes.names.txt")), "[.]"), head, 1)
        expression_matrix_nascent[is.na(expression_matrix_nascent)] = 0
        expression_matrix_nascent = colSums(expression_matrix_nascent)
        expression_matrix_nascent = expression_matrix_nascent[names(expression_matrix)]
        
        nascents = round(expression_matrix_nascent/expression_matrix*100, 2)
        saveRDS(nascents, paste0("./intronic_reads/", sample, "_nascent.rds"))

    }, error = function(e) {
        print(e$message)
    })
}
