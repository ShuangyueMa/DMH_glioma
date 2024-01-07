#####show colors
install.packages("colourpicker")
library(colourpicker)

######The figure shows the correlation degree in different intermediate DNA methylation levels
col3=c("#FF8C00", "#EE2C2C", "#0000CD")
library(ggplot2)
#gbmloli
GBM_Lollipop<-ggplot(data=gbmloli,aes(x=cut,y=meanV))+
    geom_errorbar(aes(ymin=minV,ymax=maxV),width=0,color=rep(col3[3],9),size=1.35)+
    geom_point(size=3.5, color=rep(col3[1],9))+
    theme_classic()+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none",
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
    )+
    geom_hline(yintercept =0.5,linetype="dashed",size=0.7,color="gray40")+
    geom_hline(yintercept =1.5,linetype="dashed",size=0.7,color="gray40")

#lggloli
LGG_Lollipop<-ggplot(data=lggloli,aes(x=cut,y=meanV))+
    geom_errorbar(aes(ymin=minV,ymax=maxV),width=0,color=rep(col3[3],9),size=1.35)+
    geom_point(size=3.5, color=rep(col3[2],9))+
    theme_classic()+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none",
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        )+
    geom_hline(yintercept =20,linetype="dashed",size=0.7,color="gray40")+
    geom_hline(yintercept =45,linetype="dashed",size=0.7,color="gray40")

library(cowplot)
plot_grid(GBM_Lollipop,LGG_Lollipop,ncol=1,rel_heights = c(1,1),align="hv")
save(lggloli,gbmloli,LGG_Lollipop,GBM_Lollipop,file="GBM_LGG_LollipopPlot.Rdata")
