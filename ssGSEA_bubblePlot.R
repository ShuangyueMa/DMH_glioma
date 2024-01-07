install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci) 

######The figure shows the ssGSEA score in different cell types across NMF subtypes
mypal3=c("#003399","#CC0000","#CC6633")
ssGSEA_bubblePlot<-ggplot(data=ssGSEA_bubble,aes(Comparation,cellType,size=level,color=sig))+
    geom_point(alpha=.75,shape=16)+
    scale_size("-log10_Pval",limits=c(0,85),range = c(4.5,7.05))+
    labs(title="Differences in cell infiltration level(ssGSEA score)", x="", y="")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.position="right", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),        
           )+scale_color_manual(values=mypal3[c(1,3)])

ssGSEA_bubblePlot<-ssGSEA_bubblePlot+guides(color = guide_legend(override.aes = list(size = 10)))
ssGSEA_bubblePlot
save(ssGSEA_bubblePlot,file="ssGSEA_bubblePlot.Rdata")

