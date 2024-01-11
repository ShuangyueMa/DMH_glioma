install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggbreak) 

######The figure shows the correlation between β-values of 8 survival-related CpGct and gene expressions
high_riskPlot<-ggplot(wai,aes(site,gene))+
    geom_point(data=wai,aes(site,gene,size=q,color=Rho),
               shape=16,show.legend=T)+
    scale_color_gradient(low="#638EBB",high="#CC3333")+
    labs(size="-log10 (q-value)",color="spearman's correlation")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        axis.text.x=element_text(angle=45,hjust = 0.95),
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"),
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
    )
high_riskPlot

low_riskPlot<-ggplot(wai_low,aes(site,gene))+
    geom_point(data=wai_low,aes(site,gene,size=q,color=Rho),
               shape=16,show.legend=T)+
    scale_color_gradient(low="#638EBB",high="#CC3333")+
    labs(size="-log10 (q-value)",color="spearman's correlation")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        axis.text.x=element_text(angle=45,hjust = 0.95),
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"),
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
    )
low_riskPlot
high_riskPlot

HriskP<-high_riskPlot+coord_flip()+theme(legend.position = "none")
LriskP<-low_riskPlot+coord_flip()+theme(legend.position = "none",axis.text.x = element_blank())
library(cowplot)
plot_grid(LriskP,HriskP,ncol=1,align="hv",rel_heights = c(1,1))
