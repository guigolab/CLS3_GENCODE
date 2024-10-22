library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(ComplexUpset)
library(tidyverse)
library(stringr)
library(scales)

show_hide_scale = scale_color_manual(values=c('show'='black', 'hide'='transparent'), guide='none')

upset_basic = function(masterTable, nameplot, title, w, h, s, x)
{
  #Format the table
  masterTable$tech = substr(masterTable$sampleID,8,8)
  masterTable$ONT = ifelse(masterTable$tech == "O", TRUE, FALSE)
  masterTable$PacBio = ifelse(masterTable$tech == "P", TRUE, FALSE)
  masterTable$PreCapture = ifelse(masterTable$capture == "PreCapture" | masterTable$capture == "Common", TRUE, FALSE)
  masterTable$PostCapture = ifelse(masterTable$capture == "PostCapture" | masterTable$capture == "Common", TRUE, FALSE)
  masterTable$adult = ifelse(masterTable$stage == "A", TRUE, FALSE)
  masterTable$embryo = ifelse(masterTable$stage == "E" | masterTable$stage == "P", TRUE, FALSE)
  masterTable = masterTable %>% group_by(anchTM) %>% summarize( ONT = any(ONT), PacBio = any(PacBio), adult = any(adult), embryo = any(embryo), 
                                                                PreCapture = any(PreCapture), PostCapture = any(PostCapture), 
                                                                tissue = as.character(length(unique(tissue))))
  tot = nrow(masterTable)
  masterTable = masterTable %>% rename("pre-capture" = "PreCapture", "post-capture" = "PostCapture") #Rename as decided
  cond = c("ONT","PacBio","pre-capture","post-capture","adult","embryo")
  
  stage_tech_capture = upset(masterTable, cond, name = title, sort_sets = FALSE,
         set_sizes=( upset_set_size() + geom_text(aes(label=paste0(round(after_stat(count)/tot*100,0),"%")), size = 8, hjust=1.1, stat='count') + 
                       theme( element_text(size = 13), axis.text.x = element_text(angle=90, size = 12), axis.title.x = element_blank()) + 
                       scale_y_reverse(labels = scales::comma) + expand_limits(y=s)),
         base_annotations=list(
           'Intersection size'=intersection_size(
             text_colors=c(on_background='black', on_bar='black'),
             text_mapping=aes(
               label=paste0(!!upset_text_percentage(),'\n','(',comma(!!get_size_mode('exclusive_intersection')),')'),
               y=ifelse(!!get_size_mode('exclusive_intersection') > x, !!get_size_mode('exclusive_intersection') - 0.2*x, !!get_size_mode('exclusive_intersection')),
               colour='on_background', size = 6), color = 'black') + ylim(c(0,x)) + theme(legend.position = "none", axis.title.y = element_text(size = 20), axis.text.y = element_text(size = 15))),
         annotations = list(
           Design = list(
             aes = aes(x = intersection, fill = tissue),
             geom = list(
               geom_bar(stat = 'count', position = 'fill', na.rm = TRUE),
               theme(text = element_text(size = 20), axis.text.x = element_text(size = 0)),
               ylab("# Tissues of origin"),
               geom_text(
                 aes(label = !!aes_percentage(relative_to = 'intersection'), 
                     color=ifelse(as.numeric(sapply(str_split(!!aes_percentage(relative_to = 'intersection'),"%"), head,1)) >= 5, 'show', 'hide')), size = 6,
                 stat = 'count',
                 position = position_fill(vjust = .5)), 
               show_hide_scale,
               scale_fill_manual("", values = c("1"="#c7522a","2"="#e5c185","3"="#f0daa5","4"="#ded39e","5"="#b8cdab","6"="#74a892","7"="#008585","8"="#004343"))
             ))), 
         #Set queries for assignig colors to groups in the intersection matrix and barplot
         queries = list(
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'adult'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'adult'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'pre-capture', 'adult'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'post-capture', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'adult'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'post-capture', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'post-capture', 'adult'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'post-capture', 'adult'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'post-capture', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'adult'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'adult'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'post-capture', 'adult'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'post-capture', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'adult', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'adult', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'pre-capture', 'adult'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'post-capture', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'adult'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'post-capture', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'post-capture', 'adult'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'pre-capture', 'post-capture', 'adult'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'post-capture', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'adult'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'pre-capture', 'adult'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'post-capture', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'post-capture', 'adult', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'pre-capture', 'post-capture', 'adult', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'PacBio', 'post-capture', 'adult', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT', 'post-capture', 'adult', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'post-capture', 'adult', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT','pre-capture', 'post-capture', 'adult', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio','pre-capture', 'post-capture', 'adult', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT','PacBio','pre-capture','adult', 'embryo'),
             color = 'purple',
             fill = 'purple',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT','pre-capture', 'adult', 'embryo'),
             color = 'darkblue',
             fill = 'darkblue',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio', 'pre-capture', 'adult', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('PacBio','post-capture', 'adult', 'embryo'),
             color = 'darkred',
             fill = 'darkred',
             only_components = c('intersections_matrix')
           ),
           upset_query(
             intersect = c('ONT','post-capture', 'adult', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT','PacBio','post-capture', 'adult'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT','pre-capture','post-capture', 'adult', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('ONT','pre-capture', 'adult', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio','pre-capture', 'adult', 'embryo'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio','post-capture', 'adult', 'embryo'),
             color = 'lightblue',
             fill = 'lightblue',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio','post-capture','pre-capture', 'adult', 'embryo'),
             color = 'darkgreen',
             fill = 'darkgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio','ONT','pre-capture', 'adult'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           ),
           upset_query(
             intersect = c('PacBio','pre-capture', 'adult'),
             color = 'lightgreen',
             fill = 'lightgreen',
             only_components = c('Intersection size')
           )
         )
    ) + theme(
    text = element_text(size = 13),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
    )

    ggsave(nameplot, plot = stage_tech_capture, width = w, height = h, device='tiff', dpi=700)
    print(nameplot)
}

## HUMAN
ucs = fread("Hv3_splicedmastertable_refined.gtf.gz", header = F)
colnames(ucs) = c("anchTM","capture","sampleID","tissue","stage")

ics = fread("Hv3_unsplicedmastertable_refined.gtf.gz", header = F)
colnames(ics) = c("anchTM","capture","sampleID","tissue","stage")

chains = rbind(ics, ucs)
upset_basic(chains, "human_upset_plot_techcapdstage.allIC.png", "CLS Transcripts", 23, 12, 550000, 120000)
 
## MOUSE
ucs = fread("Mv2_splicedmastertable_refined.gtf.gz", header = F)
colnames(ucs) = c("anchTM","capture","sampleID","tissue","stage")

ics = fread("Mv2_unsplicedmastertable_refined.gtf.gz", header = F)
colnames(ics) = c("anchTM","capture","sampleID","tissue","stage")

chains = rbind(ics, ucs)
upset_basic(chains, "mouse_upset_plot_techcapdstage.allIC.png", "CLS Transcripts", 23, 12, 500000, 110000)
quit()
