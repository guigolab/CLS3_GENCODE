if (!require("ggplot2")) install.packages("ggplot2")

s<-c("Hv3","Mv2")

for(spec in s)
{	options(scipen = 999)
	infile<-paste("stats/",spec,".MergedTargetRegions.proportionDetected_AllTissues",sep="")
	data_all <- read.table(infile,header = TRUE)
	data <- subset(data_all, (regionType=="detectedNovel" | regionType=="detectedNonNovel") & class!="uncaptured")
	print(data)
	data$class=factor(data$class, levels=c("NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers","allCatalogs"))
	data$regionType=factor(data$regionType, levels=c("detectedNonNovel","detectedNovel"))
	#data$spec=factor(data$spec, levels=c("Hv3"),labels=c("Human"))
	colors_labels <- c(rep("darkmagenta",6), rep("forestgreen",4), rep("lightgoldenrod3",2),rep("black",1))
	library(ggplot2)
	cbPalette=c("#DCF2F1","#d6e751") ##DCF2F1 "#B0578D")

	if(spec == "Hv3")
        {       coord=coord_cartesian(ylim = c(0, 80))
        }else{
                coord=coord_cartesian(ylim = c(0, 85))
        }

	ggplot(data,aes(y=propDetectedRegions,x=class,fill=regionType)) +
	geom_bar(position="stack", stat="identity",width=0.8)+  xlab(" ")+
	scale_fill_manual(values=cbPalette)+
	ylab ("% regions detected")+
	coord+
	theme_bw(base_size=12,base_family="Helvetica")+theme(legend.title=element_blank(), axis.text.x=element_text(face="bold",angle = 45, colour = colors_labels, vjust = 1, hjust=1)) + 

	geom_text(aes(label = paste(data$detectedRegions,"\n\n"),colour="#catalog regions",fontface="bold"),size=2,position = position_stack(0.6))+
	geom_text(aes(label = paste(data$detectedICs),colour="#transcript models",fontface="bold"), size=2,position = position_stack(0.6))
	#+ geom_text(aes(label = paste0((data$detectedRegions,colour="red"),"\n(",data$detectedLoci,")\n(",data$detectedTMs,")")), size=2,position = position_stack(0.6)) #+facet_wrap(~ref, ncol=2)
	outfile<-paste("plots/",spec,".MergedTargetRegions.proportionDetected.vennbar.tiff",sep="")
	ggsave(outfile, device = 'tiff', width = 7, height = 4)
}
