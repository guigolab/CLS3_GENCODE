options(scipen = 999)
library(ggplot2)
library(gridExtra)
library(grid)
library(tiff)

spec1 <- c("H","M")
tec <- c("pacBio","ont")
ttype <- c("perTissue","collated")

for(plot in ttype)
	{for(val in spec1)
	{	for(T in tec)
		{	print(val)
			print(T)
			data_all <- read.table(file="stats/Enrichment_phase3",header = TRUE)
			if(plot == "perTissue")
			{	data <- subset(data_all, spec==paste("CapTrap_",val,sep="") & tech==T & EnrichCateg=="OTRnonERCC")
				print(data)
				data$tech=factor(data$tech,levels=c("pacBio","ont"),labels=c("PacBio","ONT"))
				data$cap=factor(data$cap,levels=c("pre","post"),labels=c("pre-capture","post-capture"))
				if(val == "H")
				{	data$tissue <- factor(data$tissue, levels=c("EmbBrain01Rep1","Brain03Rep1","EmbHeart01Rep1","Heart02Rep1","EmbLiver01Rep1","Liver01Rep1","iPSC01Rep1","WBlood01Rep1","Testis01Rep1","Placenta01Rep1","TpoolA01Rep1","CpoolA01Rep1","allTissues"), labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","iPSC","WBlood","Testis","Placenta","Tpool","CPool","all tissues"), ordered=TRUE)	
#					colors_labels<-c("#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#331D2C","#331D2C","#331D2C","#331D2C","#090580")
					breaks_vec <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 17, 19, 21)
				}
				if(val == "M")
				{	data$tissue <- factor(data$tissue, levels=c("EmbBrain01Rep1","Brain01Rep1","EmbHeart01Rep1","Heart01Rep1","EmbLiver01Rep1","Liver01Rep1","EmbSC01Rep1","WBlood01Rep1","Testis01Rep1","TpoolA01Rep1","allTissues"), labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","EmbSC","WBlood","Testis","Tpool","all tissues"), ordered=TRUE)
#					colors_labels<-c("#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#331D2C","#331D2C","#090580")
					breaks_vec <- c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 17)
				}
				manual_breaks <- scale_x_continuous(
                                                        breaks = breaks_vec, #c(1, 2, 4, 5, 7, 8, 10, 11, 13, 15, 17, 19, 21),
                                                        labels = levels(data$tissue),
                                                        expand = c(0, 1),
                                                        minor_breaks = NULL
                                                        )
                                myplot <- ggplot(
                                           data,aes(y=OTRperc,x=category,fill=cap)
                                           )
                                text_layer <- geom_text(
                                                aes(x = category, y = max(OTRperc) + 10, label = sprintf("%sx", ifelse(Enrichment < 10, sprintf("%.1f", round(Enrichment, 1)), round(Enrichment)))),
                                                fontface = "bold",
                                                size = 1.67,
                                                angle = 45,
						family = "Helvetica"
                                                )
			}
			if(plot == "collated")
			{
				data <- subset(data_all, spec==paste("CapTrap_",val,sep="") & tech==T & EnrichCateg=="OTRnonERCC" & tissue=="allTissues")
				print(data)
				data$cap=factor(data$cap,levels=c("pre","post"))	
				if (T == "ont") {Trevised <- "ONT"} else {Trevised <- "PacBio"}
				data$tissue <- factor(data$tissue,labels=c(Trevised))
				#colors_labels<-c("#331D2C")
				manual_breaks <- NULL
				myplot <- ggplot(
                                           data,aes(y=OTRperc,x=tissue,fill=cap)
                                           )
				text_layer <- geom_text( 
						aes(x = data$tissue, y = max(OTRperc) + 5, label = sprintf("%sx", ifelse(Enrichment < 10, sprintf("%.1f", round(Enrichment, 1)), round(Enrichment)))),
                                                fontface = "bold",
                                                size = 1.67,
                                                angle = 45,
						family = "Helvetica"
                                                )
			}	
	cbPalette=c("#7fc7a2", "#9bd5f4")
	
	plot_width <- ifelse(plot == "perTissue", 85, ifelse(plot == "collated", 45))
        plot_height <- ifelse(plot == "perTissue", 45, ifelse(plot == "collated", 40))

	if (val == "H") {species <- "Human"} else {species <- "Mouse"}
	plot_title <- paste(species, T, "OTR", plot)

	p <- myplot + #ggplot(data,aes(y=OTRperc,x=category,fill=cap)) +
	geom_bar(position="dodge", stat="identity",width=0.9)+  xlab(" ")+
	scale_fill_manual(values=cbPalette)+
	ylab ("% reads on target")+
	coord_cartesian(ylim = c(0, max(data$OTRperc)+15))+
	theme_bw(base_size=5,base_family="Helvetica")+theme(legend.position = "none", axis.text.x=element_text(face="bold",family="Helvetica", size=7,angle = 45,vjust = 1, hjust=1))+
	manual_breaks+
	text_layer +
	geom_text(aes(label = sprintf("%s%%", ifelse(OTRperc < 10, sprintf("%.1f", round(OTRperc, 1)), round(OTRperc)))),position=position_dodge(width=0.9), vjust=-0.25, size=1.67, family = "Helvetica", hjust=0.5)

	#ggtitle(plot_title)

	outfile <- paste("tmp.",val,"_",T,"_",plot,"_OTR.pdf",sep="")
	ggsave(outfile, plot=p, width = plot_width, height = plot_height, units = "mm")
	outfile <- paste("tmp.",val,"_",T,"_",plot,"_OTR.tiff",sep="")
        ggsave(outfile, plot=p, width = plot_width, height = plot_height, units = "mm", dpi = 1200)
	}}
}

for(plot in ttype)
{	img1 <- rasterGrob(as.raster(readTIFF(paste0("tmp.H_pacBio_",plot,"_OTR.tiff"))), interpolate = TRUE)
	img2 <- rasterGrob(as.raster(readTIFF(paste0("tmp.H_ont_",plot,"_OTR.tiff"))), interpolate = TRUE)
	img3 <- rasterGrob(as.raster(readTIFF(paste0("tmp.M_pacBio_",plot,"_OTR.tiff"))), interpolate = TRUE)
	img4 <- rasterGrob(as.raster(readTIFF(paste0("tmp.M_ont_",plot,"_OTR.tiff"))), interpolate = TRUE)

# Save as a PDF
	if(plot == "perTissue")
	{	final_width = 180 
		final_height = 80
	}
	if(plot == "collated")
	{	final_width = 90
		final_height = 80
	}
	outfile <- paste("OTR_",plot,".pdf",sep="")
	pdf(outfile, width=final_width, height=final_height)
	grid.arrange(img1, img2, img3, img4, ncol = 2, nrow = 2)
	dev.off()
}
