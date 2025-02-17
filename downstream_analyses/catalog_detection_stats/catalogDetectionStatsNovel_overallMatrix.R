#t<-c("propDetectedRegions")
ov <- c("AllCatalogs","AllTissues","uncaptured")
t<-c("detectedICs","propDetectedRegions")
s<-c("Hv3","Mv2")
library(ggplot2)

for(cat in ov)
{   for(spec in s)
{	for(typ in t)
	{	options(scipen = 999)
		infile<-paste("stats/",spec,".refined_MergedTargetRegions.proportionDetected.Novel_",cat,sep="")
		data_all <- read.table(infile,header = TRUE)
		data <- subset(data_all, regionType=="detectedNovel")	
		print(cat)	 
		if(spec == "Hv3")
		{      
			if(cat == "AllCatalogs" || cat == "uncaptured")
			{	data$tissue=factor(data$tissue, levels=c("SIDHBrE","SIDHBrA","SIDHHeE","SIDHHeA","SIDHLiE","SIDHLiA","SIDHWbE","SIDHWbA","SIDHTeA","SIDHPlP","SIDHTpA","SIDHCpA"),labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","iPSC","WBlood","Testis","Placenta","TPool","CPool"))		
				if(cat == "uncaptured")
				{	
					data$class=factor(data$class, levels=c("uncaptured"),labels=c("uncaptured"))	
					uplim=20000
				}
				else{	#AllCatalogs
					data$class=factor(data$class, levels=c("allCatalogs"),labels=c("allCatalogs"))
					uplim=40000
				}
			
				colors_labels <- c(rep("black",1))
			}else
			{	data$tissue=factor(data$tissue, levels=c("allTissues"),labels=c("allTissues"))
				data$class=factor(data$class, levels=c("allCatalogs","NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers"))
				uplim=75000
				colors_labels <- c(rep("darkmagenta",6), rep("forestgreen",4), rep("lightgoldenrod3",2))		
			}
			#uplim=22000
			if(typ == "propDetectedRegions")
			{	uplim=65
			}
		}else
		{	if(cat == "AllCatalogs" || cat == "uncaptured")
                        {	data$tissue=factor(data$tissue, levels=c("SIDMBrE","SIDMBrA","SIDMHeE","SIDMHeA","SIDMLiE","SIDMLiA","SIDMWbE","SIDMWbA","SIDMTeA","SIDMTpA"),labels=c("EmbBrain","Brain","EmbHeart","Heart","EmbLiver","Liver","EmbSC","WBlood","Testis","TPool"))	
				colors_labels <- c(rep("black",1))
				#uplim=17000
				if(cat == "uncaptured")
				{	data$class=factor(data$class, levels=c("uncaptured"),labels=c("uncaptured"))
					uplim=30000
				}
				else{
					data$class=factor(data$class, levels=c("allCatalogs"),labels=c("allCatalogs"))
					uplim=36000
				}
			}else
			{	data$tissue=factor(data$tissue, levels=c("allTissues"),labels=c("allTissues"))
				data$class=factor(data$class, levels=c("allCatalogs","NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers"))
				colors_labels <- c(rep("darkmagenta",6), rep("forestgreen",4), rep("lightgoldenrod3",2))
				uplim=75000
			}
			
			if(typ == "propDetectedRegions")
			{       uplim=65
			}
		}
		#data$class=factor(data$class, levels=c("NONCODE", "miTranscriptome", "fantomCat", "refSeq", "gencodeLncRna", "bigTranscriptome", "CMfinderCRSs", "phyloCSF", "GWAScatalog", "UCE", "fantomEnhancers", "VISTAenhancers"))  ##scRNA excluded
		#colors_labels <- c(rep("darkmagenta",6), rep("forestgreen",4), rep("lightgoldenrod3",2))
		if(typ == "detectedICs")
		{
			outfile<-paste("plots/",spec,".refined_MergedTargetRegions.NovelICsDetected_",cat,".marix.tiff",sep="")
			p<-ggplot(data,aes(x=tissue,y=class,fill=detectedICs))
			text<-geom_text(aes(label = detectedICs),size=2)
			matrix_palette = c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
			title=ggtitle("# novel transcripts detected")
		}else
		{	
			outfile<-paste("plots/",spec,".refined_MergedTargetRegions.NovelpercDetectedRegions_",cat,".marix.tiff",sep="")
			p<-ggplot(data,aes(x=tissue,y=class,fill=propDetectedRegions))
			text<-geom_text(size=2,aes(label = round(propDetectedRegions,digits=1)))
			matrix_palette = c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
			title=ggtitle("% target regions detected in novel transcripts")
		}
		print(data)		
		
	#	ggplot(data,aes(x=tissue,y=class,fill=))+
		p+
		geom_tile()+
		scale_fill_gradientn(colours=matrix_palette,limits=c(0,uplim))+
		ylab ("target catalog")+
		theme_bw(base_size=12,base_family="Helvetica")+
		theme(legend.title=element_blank(), axis.text.y=element_text(face="bold", colour = colors_labels, vjust = 1, hjust=1), axis.text.x=element_text(face="bold",angle=45,vjust = 1, hjust=1))+
		text+
		title
		if(cat == "AllCatalogs" || cat == "uncaptured")
		{	ggsave(outfile, device = 'tiff', width = 5, height = 1.5)
		}else
		{	ggsave(outfile, device = 'tiff', width = 4, height = 4)
		}
	}
   }
}
