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

type<-c("NEWTSSs","decoyTSSs","proteincodingTSSs","lncRNATSSs")

for (t in type)
{	options(scipen = 999)
	if(t == "NEWTSSs")
        {       data_all <- read.table(file="stats/completenessSupport_novelTSSs.txt",header = TRUE)
         #       xlabel="#novel TSSs"
                max_coord=100000
		title = "Novel CLS\n80,284"
		ylab="# TSSs"
	}else if(t == "decoyTSSs")
	{	data_all <- read.table(file="stats/completenessSupport_decoyTSSs.txt",header = TRUE)
           #     xlabel="#decoy TSSs"
                max_coord=65000
		title = "Decoy models\n60,143"
		ylab=NULL
	}else if(t == "proteincodingTSSs")
        {       data_all <- read.table(file="stats/completenessSupport_proteincodingTSSs.txt",header = TRUE)
          #      xlabel="#protein coding TSSs"
                max_coord=100000
		title = "protein-coding (v27)\n83,455"
		ylab=NULL
        }else if(t == "lncRNATSSs")
        {       data_all <- read.table(file="stats/completenessSupport_lncRNATSSs.txt",header = TRUE)
         #       xlabel="#lncRNA TSSs"
                max_coord=15000
		title = "lncRNAs (v27)\n12,361"
		ylab=NULL
        }

	print(data_all)
	data <- subset(data_all, proCapfilter=="MTE5")
	data$categ=factor(data$categ, levels=c("5_support","no_support"))
	data$spec=factor(data$spec, levels=c("Hv3"),labels=c("Human"))
	data$dataset=factor(data$dataset, levels=c("cage","proCap","CageProCap"))

	library(ggplot2)
	cbPalette=c("#e5b3e5","#a6a6a6")
	ggplot(data,aes(fill=categ,y=count,x=spec)) +
	geom_bar(position="stack", stat="identity",width=0.9)+  xlab(" ")+# , width=0.5+
	scale_fill_manual(values=cbPalette)+
	labs(title = title) +
	ylab (ylab)+
	coord_cartesian(ylim = c(0, max_coord))+
	theme_bw(base_size=5,base_family="Helvetica")+
	theme(plot.title = element_text(hjust = 0.5, size = 7, family = "Helvetica"),
		legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size = 5, family = "Helvetica"), legend.key.size = unit(0.3, "cm"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())+
	geom_text(aes(label = round(data$perc,digits=2)),position=position_stack(0.5), size=1.67,fontface="bold",family = "Helvetica")+
	facet_wrap(~dataset, ncol=3)
	if(t == "NEWTSSs")
        {       ggsave('plots/Hv3_completeness_novelTSSs_5p.tiff', width = 45, height = 45, units = "mm", dpi = 1200)
        }else if(t == "decoyTSSs")
        {       ggsave('plots/Hv3_completeness_decoyTSSs_5p.tiff', width = 45, height = 45, units = "mm", dpi = 1200)
        }else if(t == "proteincodingTSSs")
        {       ggsave('plots/Hv3_completeness_proteincodingTSSs_5p.tiff', width = 45, height = 45, units = "mm", dpi = 1200)
        }else if(t == "lncRNATSSs")
        {       ggsave('plots/Hv3_completeness_lncRNATSSs_5p.tiff', width = 45, height = 45, units = "mm", dpi = 1200)
	}
}

img1 <- rasterGrob(as.raster(readTIFF("plots/Hv3_completeness_novelTSSs_5p.tiff")), interpolate = TRUE)
img2 <- rasterGrob(as.raster(readTIFF("plots/Hv3_completeness_lncRNATSSs_5p.tiff")), interpolate = TRUE)
img3 <- rasterGrob(as.raster(readTIFF("plots/Hv3_completeness_proteincodingTSSs_5p.tiff")), interpolate = TRUE)
img4 <- rasterGrob(as.raster(readTIFF("plots/Hv3_completeness_decoyTSSs_5p.tiff")), interpolate = TRUE)

final_width = 180
final_height = 45

outfile <- paste("plots/Hv3_completeness_TSSs_5p.pdf",sep="")
pdf(outfile, width=final_width, height=final_height)
grid.arrange(img1, img2, img3, img4, ncol = 4)
dev.off()
