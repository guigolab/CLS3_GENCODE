options(scipen = 999)

options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("ggplot2", "gridExtra", "grid", "tiff", "ggrepel")
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
library(ggrepel)

#s<-c("Hv3","Mv2")
t<-c("p") #,"c")
#for (sp in s)
#{      	
for (type in t)
{	options(scipen = 999)
	data_all <- read.table("stats/targs_biotypes",header = TRUE)
	#data <- subset(data_all, class != "allCatalogs" & spec == sp & type == "detected")
	data <- subset(data_all, class != "allCatalogs" & type == "detected")
        #print(sp)
	#print(data$perct)
	#if(spec=="Hv3")
	#{	coord=30000		
	#}else{	coord=25000
	#}
	data$categ=factor(data$categ, levels=c("nonIntersecting","lncRNA","protein_coding","pseudogene","Others"), labels=c("non-exonic","lncRNA","protein_coding","pseudogene","others"))
	data$class=factor(data$class, levels=c("NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers","allCatalogs"))


	#data$y_value <- ifelse(type == "p", data$perc, data$count)
	if(type=="p")
	{	data$y_value = data$perct
		data$lab = data$count
		thresh<-5 #200 
		v<-""
		ylabel<-"% target regions detected"
	}#else{	data$y_value = data$count
	#	data$lab = data$perct
	#	thresh<-5
	#	v<-"%"
	#	ylabel<-"# target regions"
	#}
	print(thresh)
	
		cbPalette=c("#f3afd4","#01a101","#373790","#FF33FF","#B2BEB5") #"#8C3333","#FAE392","#00DFA2","#FCE9F1")
		ggplot(data, aes(fill = categ, y = y_value, x = class)) +	
		geom_bar(position="stack", stat="identity")+  xlab(" ") +#, width=0.5)+
		scale_fill_manual(values=cbPalette)+
		ylab (ylabel)+
		theme_bw(base_size=5,base_family="Helvetica")+
		theme(legend.position = "bottom", legend.title=element_blank(), legend.text = element_text(size = 5, family = "Helvetica"), legend.key.size = unit(0.3, "cm"), axis.text.x=element_text(face="bold",family="Helvetica", size=5,angle = 45, vjust = 1, hjust=1))+	
		geom_text(aes(label = ifelse(y_value > thresh, paste0(round(lab, digits = 1), v), "")), size=1.67,color="white",fontface="bold",family = "Helvetica",position = position_stack(0.5))+
        	facet_wrap(~spec, ncol=2) +
		theme(panel.spacing = unit(1, "cm"))
      		
		outfile <- paste("plots/targs_biotypes_",type,".tiff",sep="")
		ggsave(outfile, width = 170, height = 65, units = "mm", dpi = 1200)	
}
#}

        img1 <- rasterGrob(as.raster(readTIFF("plots/targs_biotypes_p.tiff")), interpolate = TRUE)

        final_width = 170
        final_height = 65

        outfile <- paste("plots/targs_biotypes_p.pdf",sep="")
        pdf(outfile, width=final_width, height=final_height)
        grid.arrange(img1)#, img2, img3, img4, img5, img6, img7, img8, ncol = 4, nrow = 2)
        dev.off()

