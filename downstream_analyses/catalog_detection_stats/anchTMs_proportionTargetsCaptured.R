#################correct format out add targetS col in bash script before proceeding
###class\tregions\ttotal\tprop
options(scipen = 999)

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


#s<-c("Hv3","Mv2")
#type=c("repeated","actual")
type=c("actual")
#for(spec in s)
#{	
	for(typ in type)
	{	if(typ == "repeated")
		{	infile<-paste("stats/MergedTargetRegions.proportionAllTogether",sep="")
			outfile<-paste("plots/MergedTargetRegions.proportionAllTogether.vennbar.tiff",sep="")
		}else{	infile<-paste("stats/MergedTargetRegions.proportionAllTogetherACTUAL",sep="")
			outfile<-paste("plots/MergedTargetRegions.proportionAllTogetherACTUAL.vennbar.tiff",sep="")
		}
        	data <- read.table(infile,header = TRUE)
		data$class=factor(data$class, levels=c("pre-capture","common","post-capture","allDetected"), labels=c("pre-capture","common","post-capture","all detected"))

		library(ggplot2)

		ggplot(data,aes(fill="#F0AFA5",y=prop,x=class)) +
		geom_bar(position=position_dodge(), stat="identity")+  xlab(" ") +#, width=0.5)+
		ylab ("%target regions")+
		coord_cartesian(ylim = c(0, 50))+
		theme_bw(base_size=5,base_family="Helvetica")+	
		theme(legend.position = "none", legend.title=element_blank(), axis.text.x=element_text(face="bold",family="Helvetica", size=5,angle = 45, vjust = 1, hjust=1))+
		geom_text(aes(label = after_stat(round(y, digits=2))),vjust=-1,size=1.67,fontface="bold",family = "Helvetica") +
		facet_wrap(~spec, ncol=1)
                #theme(panel.spacing = unit(1, "cm"))
		ggsave(outfile, width = 50, height = 115, units = "mm", dpi = 1200)
	}



       img1 <- rasterGrob(as.raster(readTIFF(paste0("plots/MergedTargetRegions.proportionAllTogetherACTUAL.vennbar.tiff"))), interpolate = TRUE)
       # img2 <- rasterGrob(as.raster(readTIFF(paste0("plots/Mv2.MergedTargetRegions.proportionAllTogetherACTUAL.vennbar.tiff"))), interpolate = TRUE)

        final_width = 50
        final_height = 115

        outfile <- paste("plots/MergedTargetRegions.proportionACTUAL.pdf",sep="")
        pdf(outfile, width=final_width, height=final_height)
        grid.arrange(img1, ncol = 1)
        dev.off()

