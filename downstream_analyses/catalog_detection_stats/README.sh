####Target regions detected per target class/catalog (barplot)
#---Plots the % target regions detected overall, further differentiating the novel and known regions, in addition to displaying the number of catalog regions and transcripts (CLS3 transcripts/intron chains)

bash rcd_catalogDetectionStats.sh

#final plot -> https://docs.google.com/presentation/d/1S4pvtA8Lpl8P_fMwfMA2GnKOKp8mIoZaAQWCFaMuGiw/edit#slide=id.g30a87e2963f_0_0

####Novel transcripts detected per target class per tissue (matrices)
#---Plots novel transcripts (CLS3 transcripts/intron chains) as well as the proportion of target regions that help detect novel transcripts, thanks to the target catalog

bash rcd_catalogDetection_perTissue.sh
#final plots: https://docs.google.com/presentation/d/1S4pvtA8Lpl8P_fMwfMA2GnKOKp8mIoZaAQWCFaMuGiw/edit#slide=id.g2ffe14a0607_0_472 
#https://docs.google.com/presentation/d/1S4pvtA8Lpl8P_fMwfMA2GnKOKp8mIoZaAQWCFaMuGiw/edit#slide=id.g2fa50d3fcc5_1_16
