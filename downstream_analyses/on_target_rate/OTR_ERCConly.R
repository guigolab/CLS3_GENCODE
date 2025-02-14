options(scipen = 999)
library(ggplot2)
library(qpdf) 

pdf_list <- c()  

spec1 <- c("H","M")
tec <- c("pacBio","ont")
ttype <- c("perTissue","collated")

for(plot in ttype)
	{for(val in spec1)
	{	for(T in tec)
		{	print(val)
			print(T)
			data_all <- read.table(file="Enrichment_phase3test",header = TRUE)
			if(plot == "perTissue")
			{	data <- subset(data_all, spec==paste("CapTrap_",val,sep="") & tech==T & EnrichCateg=="OTRonlyERCC")
				print(data)
				data$tech=factor(data$tech,levels=c("pacBio","ont"),labels=c("pacBio","ONT"))
				data$cap=factor(data$cap,levels=c("pre","post"))
				if(val == "H")
				{	data$tissue <- factor(data$tissue, levels=c("EmbBrain01Rep1","Brain03Rep1","EmbHeart01Rep1","Heart02Rep1","EmbLiver01Rep1","Liver01Rep1","iPSC01Rep1","WBlood01Rep1","Testis01Rep1","Placenta01Rep1","TpoolA01Rep1","CpoolA01Rep1","allTissues"), labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","iPSC","WBlood","Testis","Placenta","Tpool","CPool","all tissues"), ordered=TRUE)	
#					colors_labels<-c("#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#331D2C","#331D2C","#331D2C","#331D2C","#090580")
				}
				if(val == "M")
				{	data$tissue <- factor(data$tissue, levels=c("EmbBrain01Rep1","Brain01Rep1","EmbHeart01Rep1","Heart01Rep1","EmbLiver01Rep1","Liver01Rep1","EmbSC01Rep1","WBlood01Rep1","Testis01Rep1","TpoolA01Rep1","allTissues"), labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","EmbSC","WBlood","Testis","Tpool","all tissues"), ordered=TRUE)
#					colors_labels<-c("#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#A78295","#331D2C","#331D2C","#331D2C","#090580")
				}
			}
			if(plot == "collated")
			{	data <- subset(data_all, spec==paste("CapTrap_",val,sep="") & tech==T & EnrichCateg=="OTRonlyERCC" & tissue=="allTissues")
				print(data)
				data$cap=factor(data$cap,levels=c("pre","post"))
				data$tissue <- factor(data$tissue,labels=c(T))
#				colors_labels<-c("#331D2C")
			}

	cbPalette=c("#7fc7a2", "#9bd5f4")
	
	plot_width <- ifelse(plot == "perTissue", 7, ifelse(plot == "collated", 3, 6))
        plot_height <- ifelse(plot == "perTissue", 4, ifelse(plot == "collated", 3, 5))
        pdf_filename <- paste0("tmp.plot_", plot, "_", val, "_", T, "-ERCConly.pdf")
        pdf_list <- c(pdf_list, pdf_filename)
        pdf(file = pdf_filename, width = plot_width, height = plot_height)

	plot_title <- paste(val, T, "OTR", plot, "- ERCC only")

	p <- ggplot(data,aes(y=OTRperc,x=tissue,fill=cap)) +
	geom_bar(position="dodge", stat="identity",width=0.9)+  xlab(" ")+
	scale_fill_manual(values=cbPalette)+
	ylab ("% reads on target")+
	coord_cartesian(ylim = c(0, max(data$OTRperc)+10))+
	theme_bw(base_size=12,base_family="Helvetica")+theme(legend.title=element_blank(), axis.text.x=element_text(face="bold",angle = 90,vjust = 1, hjust=1))+
	geom_text(aes(tissue, max(OTRperc)+5, label=paste("x",round(Enrichment,1))),fontface="bold",size=3.5,angle=45)+
	geom_text(aes(label=paste(round(OTRperc,1),"%",sep="")), position=position_dodge(width=0.9), vjust=-0.25, size=2.5, hjust=0.5)+
        ggtitle(plot_title)

        print(p)

#	outfile <- paste(val,"_",T,"_",plot,"_OTRonlyERCC.tiff",sep="")
#	if(plot == "perTissue")
#	{	
#	ggsave(outfile, device = 'tiff', width = 7, height = 4)
#	}
#	if(plot == "collated")
#	{	ggsave(outfile, device = 'tiff', width = 3, height = 3)
#	}
	dev.off()
	}}
}

# Merge all PDFs into a single file
if (length(pdf_list) > 1) {
    pdf_combine(input = pdf_list, output = "OTR-ERCConly_plots.pdf")
    cat("Merged PDF saved as: OTR-ERCConly_plots.pdf\n")
} else {
    cat("Only one PDF generated, no need to merge.\n")
}
