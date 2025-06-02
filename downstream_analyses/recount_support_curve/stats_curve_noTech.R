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

s<-c("Hv3","Mv2")
for(spec in s)
{	#data <- read.table(file="stats/filterThreshold_vs_supportedTMs.txt",header = TRUE)
	infile<-paste("stats/",spec,"_filterThreshold_vs_supportedTMs_3classes.txt",sep="")
	data_all <- read.table(infile,header = TRUE)
	data <- subset(data_all,tech=="all_ICs"|tech=="intergenic_ICs"|tech=="novel_ICs")
	data$tech=factor(data$tech, levels=c("all_ICs","intergenic_ICs","novel_ICs"))
	LINES <- c("all_ICs"="longdash","novel_ICs"="longdash","intergenic_ICs"="longdash")

	cbPalette=c("all_ICs"="#000000","novel_ICs"="#ED9455","intergenic_ICs"="#bf62a6")
	library(ggplot2)

	ggplot(data,aes(y=percsupportedICs,x=filterThreshold,group=tech,linetype=factor(tech))) +
	geom_line(aes(color=factor(tech)))+
#geom_line()+
	scale_color_manual(values=cbPalette)+
	scale_linetype_manual(values = LINES)+
	labs(linetype="noveltyStatus",color="noveltyStatus")+
#guides(color = guide_legend(title = "Users By guides"))+
#geom_line()+
#scale_linetype_manual(values=c("normal","normal","normal","normal","twodash","twodash","twodash","twodash","dotted","dotted","dotted","dotted"))+

#geom_point(aes(color=tech))+
	coord_cartesian(ylim = c(0, 100))+ 
	ylab ("% transcripts supported") +
	xlab ("# recount reads support")+
	theme_bw(base_size=5,base_family="Helvetica")+
        theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size = 5, family = "Helvetica"), legend.key.size = unit(0.3, "cm"), axis.text.x=element_text(face="bold",family="Helvetica", size=5, vjust = 1, hjust=1))

#ggsave('recountSupportFILTERS.bar.tiff', device = 'tiff', width = 6, height = 4) 
	outfile<-paste("plots/",spec,"_recountSupportFILTERS_3classesCurve_noTech",".tiff",sep="")
	ggsave(outfile, width = 80, height = 80, units = "mm", dpi = 1200)
}

	img1 <- rasterGrob(as.raster(readTIFF("plots/Hv3_recountSupportFILTERS_3classesCurve_noTech.tiff")), interpolate = TRUE)
        img2 <- rasterGrob(as.raster(readTIFF("plots/Mv2_recountSupportFILTERS_3classesCurve_noTech.tiff")), interpolate = TRUE)

        final_width = 160
        final_height = 80

        outfile <- paste("plots/recountSupportFILTERS_3classesCurve_noTech.pdf",sep="") 
        pdf(outfile, width=final_width, height=final_height)
        grid.arrange(img1, img2, ncol = 2, nrow=1)
        dev.off()
