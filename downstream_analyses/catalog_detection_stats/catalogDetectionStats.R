options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("ggplot2", "gridExtra", "grid", "tiff")
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!(pkg %in% installed)) {
    install.packages(pkg)
        }
      }
library(ggplot2)
library(gridExtra)
library(grid)
library(tiff)

s<-c("Hv3","Mv2")
options(scipen = 999)
#for(spec in s)
#{	
	#infile<-paste("stats/",spec,".MergedTargetRegions.proportionDetected_AllTissues",sep="")
	infile<-paste("stats/MergedTargetRegions.proportionDetected_AllTissues",sep="")
        data_all <- read.table(infile,header = TRUE)
	data <- subset(data_all, (regionType=="detectedNovel" | regionType=="detectedNonNovel") & class!="uncaptured")
	print(data)
	data$class=factor(data$class, levels=c("NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers","allCatalogs"))
	data$regionType=factor(data$regionType, levels=c("detectedNonNovel","detectedNovel"))
	#data$spec=factor(data$spec, levels=c("Hv3"),labels=c("Human"))
	colors_labels <- c(rep("darkmagenta",6), rep("forestgreen",4), rep("lightgoldenrod3",2),rep("black",1))
	cbPalette=c("#cee3f0","#eaf759") #"#DCF2F1","#d6e751"

	#if(spec == "Hv3")
        #{       coord=coord_cartesian(ylim = c(0, 80))
        #}else{
        #        coord=coord_cartesian(ylim = c(0, 85))
        #}
	coord=coord_cartesian(ylim = c(0, 85))

	ggplot(data,aes(y=propDetectedRegions,x=class,fill=regionType)) +
	geom_bar(position="stack", stat="identity",width=0.8)+  xlab(" ")+
	scale_fill_manual(values=cbPalette)+
	ylab ("% regions detected")+
	coord+
	theme_bw(base_size=5,base_family="Helvetica")+
	theme(legend.position="bottom",legend.title=element_blank(), axis.text.x=element_text(face="bold", family="Helvetica", size=5, angle = 45, colour = colors_labels, vjust = 1, hjust=1)) + 
	geom_text(aes(label = paste(data$detectedRegions,"\n\n"),colour="#catalog regions",fontface="bold"),size=1.67,family = "Helvetica",position = position_stack(0.6))+
	geom_text(aes(label = paste(data$detectedICs),colour="#transcript models",fontface="bold"), size=1.67,family = "Helvetica",position = position_stack(0.6))+
  	scale_color_manual(values = c("#catalog regions" = "#f55858",
                                "#transcript models" = "#577af7")) +
	guides(color = guide_legend(title = "Region Type"))+
	

#+ geom_text(aes(label = paste0((data$detectedRegions,colour="red"),"\n(",data$detectedLoci,")\n(",data$detectedTMs,")")), size=2,position = position_stack(0.6)) 
	facet_wrap(~spec, nrow=2)
	#outfile<-paste("plots/",spec,".MergedTargetRegions.proportionDetected.vennbar.tiff",sep="")
	outfile<-paste("plots/MergedTargetRegions.proportionDetected.vennbar.tiff",sep="") 
	ggsave(outfile, width = 130, height = 115, units = "mm", dpi = 1200)
#}

img1 <- rasterGrob(as.raster(readTIFF(paste0("plots/MergedTargetRegions.proportionDetected.vennbar.tiff"))), interpolate = TRUE) 

        final_width = 130
        final_height = 115

        outfile <- paste("plots/MergedTargetRegions.proportionDetected.vennbar.pdf",sep="")
        pdf(outfile, width=final_width, height=final_height)
        grid.arrange(img1, ncol = 1)
        dev.off()
