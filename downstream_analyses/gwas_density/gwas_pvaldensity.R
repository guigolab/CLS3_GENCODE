###LOAD PACKAGE
packages = c("dplyr", "optparse", "ggplot2", "rstatix", "ggpubr", "ggbreak", "patchwork", "tidyverse")
invisible(lapply(packages, library, character.only = TRUE))

option_list = list(
  make_option(c("-t", "--test"), type="character", default=NULL, help="BED file from bedtools intersect", metavar="character"),  
  make_option(c("-c", "--reference"), type="character", default=NULL, help="BED file from bedtools intersect for reference", metavar="character"),
  make_option(c("-r", "--random"), type="character", default=NULL, help="BED file from bedtools intersect for randomised set", metavar="character"),
  make_option(c("-p", "--control"), type="character", default=NULL, help="BED file from bedtools intersect for reference", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="", help="Title", metavar="character"),
  make_option(c("-m", "--gene_mapping"), type="character", default="/users/rg/tperteghella/Gencode/Analysis/GWAS_hit/data/mapping_ic_ID.tsv", help="Title", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

format = c("chrcls","startcls","endcls","gene_id","chr","start","pos","trait","pvalue","assigned","overlap")

if ( ! file.exists(paste0("densities_gwas",opt$name,".rds")) )
{
    print("Collect data")
    mapping = read.table(opt$gene_mapping)
    colnames(mapping) = c("gene_id","ics")

    test = read.delim(opt$test, header = F, row.names = NULL, quote = "", sep = "\t")
    colnames(test) = c("chrcls","startcls","endcls","ics","chr","start","pos","trait","pvalue","assigned","overlap")
    test = merge(test, mapping, by = "ics")
    test = test[, format]

    randomised = read.delim(opt$random, header = F, row.names = NULL, quote = "", sep = "\t")
    colnames(randomised) = format  #No need to merge, already gene level

    reference = read.delim(opt$reference, header = F, row.names = NULL, quote = "", sep = "\t")
    colnames(reference) = format  #No need to merge, already gene level

    control = read.delim(opt$control, header = F, row.names = NULL, quote = "", sep = "\t")
    colnames(control) = format   #No need to merge, already gene level

    print("Shape dataframes")
  
    #change ICS to GENE_ID to go back to Intron Chain level
    void = as.data.frame(test[test$pos == -1, ] %>% group_by(gene_id) %>% summarize( hits = 0, pruned = 0, startcls = min(startcls), endcls = max(endcls) ))
    test = as.data.frame(test[test$pos != -1, ] %>% group_by(gene_id) %>% summarize( hits = n_distinct(chr,pos), pruned = n_distinct(assigned), startcls = min(startcls), endcls = max(endcls) ))
    test = rbind(test,void) %>% group_by(gene_id) %>% summarize( hits = sum(hits), pruned = sum(pruned), startcls = min(startcls), endcls = max(endcls) )
    test$tag = rep("CLS_Intergenic_v27", nrow(test))
    test$length = test$endcls-test$startcls
    test$density = (test$hits/test$length)*1000
    test$densitypruned = (test$pruned/test$length)*1000

    void = as.data.frame(randomised[randomised$pos == -1, ] %>% group_by(gene_id) %>% summarize( hits = 0, pruned = 0, startcls = min(startcls), endcls = max(endcls) ))
    randomised = as.data.frame(randomised[randomised$pos != -1, ] %>% group_by(gene_id) %>% summarize( hits =  n_distinct(chr,pos), pruned = n_distinct(assigned), startcls = min(startcls), endcls = max(endcls) ))
    randomised = rbind(randomised,void) %>% group_by(gene_id) %>% summarize( hits = sum(hits), pruned = sum(pruned), startcls = min(startcls), endcls = max(endcls) )
    randomised$tag = rep("Decoy", nrow(randomised))    
    randomised$length = randomised$endcls-randomised$startcls
    randomised$density = (randomised$hits/randomised$length)*1000
    randomised$densitypruned = (randomised$pruned/randomised$length)*1000

    void = as.data.frame(reference[reference$pos == -1, ] %>% group_by(gene_id) %>% summarize( hits = 0, pruned = 0, startcls = min(startcls), endcls = max(endcls) ))
    reference = as.data.frame(reference[reference$pos != -1, ] %>% group_by(gene_id) %>% summarize( hits =  n_distinct(chr,pos), pruned = n_distinct(assigned), startcls = min(startcls), endcls = max(endcls) ))
    reference = rbind(reference,void) %>% group_by(gene_id) %>% summarize( hits = sum(hits), pruned = sum(pruned), startcls = min(startcls), endcls = max(endcls) )
    reference$tag = rep("lncRNA_v27", nrow(reference))    
    reference$length = reference$endcls-reference$startcls
    reference$density = (reference$hits/reference$length)*1000
    reference$densitypruned = (reference$pruned/reference$length)*1000

    void = as.data.frame(control[control$pos == -1, ] %>% group_by(gene_id) %>% summarize( hits = 0, pruned = 0, startcls = min(startcls), endcls = max(endcls) ))
    control = as.data.frame(control[control$pos != -1, ] %>% group_by(gene_id) %>% summarize( hits =  n_distinct(chr,pos), pruned = n_distinct(assigned), startcls = min(startcls), endcls = max(endcls) ))
    control = rbind(control,void) %>% group_by(gene_id) %>% summarize( hits = sum(hits), pruned = sum(pruned), startcls = min(startcls), endcls = max(endcls) )
    control$tag = rep("coding_v27", nrow(control))
    control$length = control$endcls-control$startcls
    control$density = (control$hits/control$length)*1000
    control$densitypruned = (control$pruned/control$length)*1000

    final = do.call(rbind, list(test, randomised, reference, control))
    saveRDS(final, paste0("densities_gwas",opt$name,".rds"))
} else {
  final = readRDS(paste0("densities_gwas",opt$name,".rds"))
}

final$tag = factor(final$tag, levels = c("lncRNA_v27","CLS_Intergenic_v27","Decoy","coding_v27"))

#with zeroes
Ttest = final %>% pairwise_t_test(as.formula(densitypruned ~ tag), paired = FALSE)

bxp = ggplot(final) + geom_violin(aes(x = tag, y = densitypruned, fill = tag), outlier.shape = NA) + ggtitle(paste("GWAS density - Pruned catalog", opt$name)) +
  scale_y_continuous(limits = quantile(final$densitypruned, c(0, 0.96))) +
  theme(text = element_text(size = 25), legend.position = "none") + stat_summary(fun=mean, mapping = aes( x = tag, y = density ), geom="point", shape=20, size=8, color="red", fill="red") +
  scale_fill_manual(values=c("#5B4A8D","#BC96E6","#D8B4E2","lightblue")) + xlab("") + scale_x_discrete(guide = guide_axis(n.dodge=2)) + 
  theme(axis.text.x.top = element_blank())

#Add statistical test p-values
Ttest = Ttest %>% add_xy_position(x = "tag")
Ttest$y.position = c(0.15,0.2,0.25,0.3,0.35,0.4)#-0.12
bxp = bxp + stat_pvalue_manual(Ttest, label = "p.adj", step.increase = 0)
