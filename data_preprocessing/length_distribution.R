library(ggplot2)
library(data.table)
library(scales)
library(wesanderson)
library(dplyr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
cbPalette=c("")

print(paste0("Reading: ",args[1]))
colnames(plot) = c('id', 'seqTech', 'cap', 'sizeFrac', 'sample', 'length')
plot = plot[, -1]

plot$seqTech=ifelse(grepl("ont", plot$seqTech), "ONT", "PacBio")
plot$seqTech=factor(plot$seqTech, levels=c("ONT","PacBio"))

plot$cap = ifelse(plot$cap == "Hv3" | plot$cap == "Mv2", "post-capture", "pre-capture")
plot$species = ifelse(grepl("H", plot$cap), "Human", "Mouse")
plot$sample = str_sub(plot$sample, 1, str_length(plot$sample)-6)

datSumm = plot %>% group_by(seqTech, sizeFrac, cap, sample) %>% summarise(n=n(), med=median(length))

summaryStats = transform(datSumm, LabelN = paste0('N=', comma(n)), LabelM = paste0( 'Median=', comma(med)))
geom_textSize = 25
lineSize = 2

pdf(gsub(".tsv", ".pdf", args[1]), bg = "white", width = 100, height = 150)
print(
   ggplot(data=plot, aes(x=length)) + geom_histogram(binwidth=100, position='identity') + scale_x_continuous(limits=c(0,5000), labels=comma) + ylab("Count") + theme_bw() + ylim(c(0,3500)) +
       scale_fill_manual(cbPalette) + xlab("") + theme(axis.text.x = element_text(size = 70, colour = "black", vjust=0.5), 
       axis.text.y = element_text(size = 70,  colour = "black"), legend.text = element_text(size=0), plot.title = element_text(size = 0), legend.title =element_text(size=0, color="white"), 
       axis.title.x = element_text(size = 0), axis.title.y = element_text(size = 0), legend.text.align=0.5, strip.text.x=element_text(size=80), strip.text.y=element_text(size=80)) + 
       theme(legend.position="top") + facet_grid( seqTech + sample ~ cap, scales="free_y") + geom_vline(data = summaryStats, aes(xintercept=med), color='#ff0055', linetype='solid', size=lineSize) + 
       geom_text(data = summaryStats, aes(label = LabelN, x = Inf, y = Inf), hjust=1, vjust=1.5, size=geom_textSize, fontface = 'bold') +
       geom_text(data = summaryStats, aes(label = LabelM, x = Inf, y = Inf), hjust=1, vjust=3.5, size=geom_textSize, fontface = 'bold', color='#ff0055') + scale_y_continuous(labels = scales::comma)
)
dev.off()
