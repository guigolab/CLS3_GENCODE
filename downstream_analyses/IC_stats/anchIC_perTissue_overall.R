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

s <- c("Hv3","Mv2")
tec <- c("overall","pacBioSII","ont")
nov <- c("ALL","N")

for(spec in s)
{ for(novelty in nov)
  { infile<-paste("stats/",spec,"_refined_overall_",novelty,"_anchIC_capstats",sep="")
    data <- read.table(infile,header = TRUE)
    data$capture=factor(data$capture, levels=c("preCap","common","postCap"), c("pre-capture","common","post-capture"))
    category_totals <- aggregate(stats ~ categ, data, sum) 
    cbPalette=c("#7fc7a2", "#357955", "#9bd5f4")
        ggplot(data,aes(fill=capture,y=stats,x=categ)) +
        geom_bar(position="stack", stat="identity", width=1)+
        geom_text(data = category_totals,aes(x = categ, y = stats, label = format(stats, big.mark = ",")),vjust = -0.5,size = 1.67,fontface="bold",family = "Helvetica", color = "gray30", inherit.aes = FALSE) +
        labs(x = NULL, y = NULL) +
	scale_fill_manual(values=cbPalette)+
        coord_cartesian(ylim = c(0,category_totals$stats))+
        theme_bw(base_size=5,base_family="Helvetica")+
        theme(legend.position = "none", axis.text.x = element_blank(),
                                        axis.text = element_blank(),      # remove y-axis text
                                        axis.ticks = element_blank(),       # remove axis ticks
                                        panel.grid = element_blank(),       # remove gridlines
                                        panel.border = element_blank(),     # remove plot border
                                        axis.title = element_blank()        # remove axis titles
                )
        outfile<-paste("plots/",spec,"_",novelty,"_anchIC.vennbar.pdf",sep="")
        ggsave(outfile, width = 20, height = 45, units = "mm")
  }
}


for(spec in s)
{for(novelty in nov)
{for(tech in tec)
{	infile<-paste("stats/",spec,"_refined_",tech,"_",novelty,"_anchIC_stats",sep="")
        data <- read.table(infile,header = TRUE)
	data$capture=factor(data$capture, levels=c("pre-capture","common","post-capture"))
	if(spec == "Hv3")
	{	data$tissue=factor(data$tissue, levels=c("EmbryoBrain","AdultBrain","EmbryoHeart","AdultHeart","EmbryoLiver","AdultLiver","EmbryoiPSC","AdultWBlood","AdultTestis","PlacentaPlacenta","AdultCpoolA","AdultTpoolA"),labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","iPSC","WBlood","Testis","Placenta","CPool","TPool"))
		breaks_vec <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 17, 19)
	}
	else{
		data$tissue=factor(data$tissue, levels=c("EmbryoBrain","AdultBrain","EmbryoHeart","AdultHeart","EmbryoLiver","AdultLiver","EmbryoESC","AdultWBlood","AdultTestis","AdultTpoolA"),labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","EmbESC","WBlood","Testis","TPool"))
		breaks_vec <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15)
	}
	manual_breaks <- scale_x_continuous(
        		breaks = breaks_vec, #c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 17, 19, 21),
                        labels = levels(data$tissue),
                        expand = c(0, 1),
                        minor_breaks = NULL
                        )
        if(novelty=="ALL")
	{	lim = 160000
	}
	else{			
		lim = 30000
	}      
	library(ggplot2)
	cbPalette=c("#ffe24d","#ffa500","#ec3b7b","#bfbfbf")
	cbPalette=c("#7fc7a2", "#357955", "#9bd5f4", "#ff9cb5", "#ff4f7c", "#ff7a4f")
	ggplot(data,aes(fill=capture,y=stats,x=categ2)) +
	geom_bar(position="stack", stat="identity",color="gray22", linewidth = 0.2, width=1)+ xlab(" ") +
	scale_fill_manual(values=cbPalette)+
	ylab ("# transcript models")+
#coord_cartesian(ylim = c(0, 80000))+
	coord_cartesian(ylim = c(0, lim))+
	theme_bw(base_size=5,base_family="Helvetica")+theme(legend.position = "none", axis.text.x = element_text(face="bold",family="Helvetica", size=6,angle=45, vjust=1, hjust=1))+
	manual_breaks
	#guides(fill = "none")
	outfile<-paste("plots/",spec,"_perTissue_",tech,"_",novelty,"_anchIC.vennbar.pdf",sep="")
	ggsave(outfile, width =50, height = 45, units = "mm")
}	
}
}
# Save as a PDF
#for(novelty in nov)
#{	img1 <- rasterGrob(as.raster(readTIFF(paste0("plots/Hv3_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE)
#	img2 <- rasterGrob(as.raster(readTIFF(paste0("plots/Hv3_perTissue_overall_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE)
#        img3 <- rasterGrob(as.raster(readTIFF(paste0("plots/Hv3_perTissue_pacBioSII_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE) 
#	img4 <- rasterGrob(as.raster(readTIFF(paste0("plots/Hv3_perTissue_ont_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE) 
#	img5 <- rasterGrob(as.raster(readTIFF(paste0("plots/Mv2_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE)
#        img6 <- rasterGrob(as.raster(readTIFF(paste0("plots/Mv2_perTissue_overall_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE)
#	img7 <- rasterGrob(as.raster(readTIFF(paste0("plots/Mv2_perTissue_pacBioSII_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE#)
#	img8 <- rasterGrob(as.raster(readTIFF(paste0("plots/Mv2_perTissue_ont_",novelty,"_anchIC.vennbar.tiff"))), interpolate = TRUE)
#
#	final_width = 180
#        final_height = 90

#        outfile <- paste("plots/",novelty,"_anchIC.pdf",sep="")
#        pdf(outfile, width=final_width, height=final_height)
#        grid.arrange(img1, img2, img3, img4, img5, img6, img7, img8, ncol = 4, nrow = 2)
#        dev.off()
#}
