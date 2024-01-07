install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 

######The figure shows the correlation between intermediate DNA methylation and tumor immune microenvironment
mypal3=c("#FF8C00", "#EE2C2C", "#0000CD")

#Correlation between PIM score and ESTIMATE score
cor_tu1<-cor.test(DF$PIM,DF$ESTIMATE.combined.score,method="spearman")
PIM_ESti<-ggplot(data = DF, mapping = aes(x = PIM,
                                     y = ESTIMATE.combined.score))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=1.2,alpha=0.4)+
    scale_fill_distiller(palette = "YlOrBr")+
    xlab("PIM")+
    ylab("ESTIMATE score")+
    theme_few()+
    theme(legend.position = "none")+
    xlim(0.08,0.26)+
    ylim(-1669.96,7626.79)+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu1$estimate,2),"p=",formatC(cor_tu1$p.value,2),sep=" "),family="Arail",
                 parse=F,x=0.19,y=3945.00,color="black",size=5.5,fontface="italic")

PIM_ESti
#Correlation between PIM score and immune score
cor_tu1<-cor.test(DF$PIM,DF$ESTIMATE.stromal.score,method="spearman")
PIM_immune<-ggplot(data = DF, mapping = aes(x = PIM,
                                     y = ESTIMATE.immune.score))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=1.2,alpha=0.4)+
        scale_fill_distiller(palette = "YlOrBr")+
    xlab("PIM")+
    ylab("immune cell infiltration level")+
    theme_few()+
    theme(legend.position = "none")+
    xlim(0.08,0.26)+
    ylim(-1122.40,4285.21)+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme(
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu1$estimate,2),"p=",formatC(cor_tu1$p.value,2),sep=" "),family="Arail",
                 parse=F,x=0.19,y=4285.21,color="black",size=5.5,fontface="italic")
PIM_immune
#Correlation between PIM score and stromal score
cor_tu1<-cor.test(DF$PIM,DF$ESTIMATE.stromal.score,method="spearman")
PIM_stromal<-ggplot(data = DF, mapping = aes(x = PIM,
                                     y = ESTIMATE.stromal.score))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=1.2,alpha=0.4)+
        scale_fill_distiller(palette = "YlOrBr")+
    xlab("PIM")+
    ylab("stromal cell infiltration level")+
    theme_few()+
    theme(legend.position = "none")+
    xlim(0.08,0.26)+
    ylim(-617.88,3945.00)+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        #plot.margin = unit(c(1.5,1.5,1.5,1.5),"cm"),
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu1$estimate,2),"p=",formatC(cor_tu1$p.value,2),sep=" "),family="Arail",
                 parse=F,x=0.19,y=3945.00,color="black",size=5.5,fontface="italic")

PIM_stromal

save(PIM_stromal,PIM_immune,PIM_ESti,file="PIM_tumorImmuneMicro_Correlation.Rdata")
