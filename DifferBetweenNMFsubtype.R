install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 

######The figure shows the PIM differences between NMF subtypes
mypal3=c("#003399","#CC3333","#CC6633")
compaired<-list(c("Cluster1","Cluster2"),c("Cluster2","Cluster3"),c("Cluster1","Cluster3"))

#boxplot shows the PIM differences between NMF subtypes
PIM_box<-ggplot(PIM_NMF,aes(NMF_cluster,PIM,col=NMF_cluster))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="PIM")+
    scale_x_discrete(limits=c("Cluster1","Cluster2","Cluster3"))+
    scale_color_manual(values = mypal3)+
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
    )+geom_signif(comparisons=compaired,step_increase = 0.1,map_signif_level = T,y_position=max(ppdata$PIM),textsize=5.5,parse=F,color="black",tip_length = 0)
PIM_box

######The figure shows the immune score differences between NMF subtypes
ymax3=max(boxImmune$ESTIMATE.immune.score)
ymin3=min(boxImmune$ESTIMATE.immune.score)
ESTIMATE.immune.score_plot<-ggplot(boxImmune,aes(NMF_cluster,ESTIMATE.immune.score,col=NMF_cluster))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    scale_color_manual(values=mypal3)+
    labs(x="",y="immune cell infiltration level")+
    scale_x_discrete(limits=c("Cluster1","Cluster2","Cluster3"))+
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
    )+geom_signif(comparisons=compaired,step_increase = 0.1,map_signif_level = T,y_position=ymax3,textsize=5.5,parse=F,color="black",tip_length = 0)

ESTIMATE.immune.score_plot

######The figure shows the stromal score differences between NMF subtypes
ymax2=max(boxStromal$ESTIMATE.stromal.score)
ymin2=min(boxStromal$ESTIMATE.stromal.score)
ESTIMATE.stromal.score_plot<-ggplot(boxStromal,aes(NMF_cluster,ESTIMATE.stromal.score,col=NMF_cluster))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="stromal cell infiltration level")+
    scale_color_manual(values=mypal3)+
    scale_x_discrete(limits=c("Cluster1","Cluster2","Cluster3"))+
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
        axis.line = element_line(color = "black", size =1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+geom_signif(comparisons=compaired,step_increase = 0.1,map_signif_level = T,y_position=ymax2,textsize=5.5,parse=F,color="black",tip_length = 0)

ESTIMATE.stromal.score_plot

#####The figure shows the ESTIMATE score differences between NMF subtypes
ymax5=max(boxEsti$ESTIMATE.combined.score)
ymin5=min(boxEsti$ESTIMATE.combined.score)
ESTIMATE.combined.score_plot<-ggplot(boxEsti,aes(NMF_cluster,ESTIMATE.combined.score,col=NMF_cluster))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="tumor purify")+
    scale_color_manual(values=mypal3)+
    ylim(0,10000)+
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
    )+geom_signif(comparisons=compaired,step_increase = 0.1,map_signif_level = T,y_position=10000,textsize=5.5,parse=F,color="black",tip_length = 0)

ESTIMATE.combined.score_plot


######The figure shows the DNA methylation differences between NMF subtypes
avemethplot=ggplot(aveMethdata,aes(NMF_cluster,meanMethLevel,col=NMF_cluster))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    scale_color_manual(values=mypal3)+
    labs(x="",y="Methylation Level")+
    scale_x_discrete(limits=c("Cluster1","Cluster2","Cluster3"))+
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
    )+geom_signif(comparisons=compaired,step_increase = 0.08,map_signif_level = T,y_position=max(aveMethdata$meanMethLevel),textsize=5.5,parse=F,color="black",tip_length = 0)

avemethplot

save(PIM_NMF,PIM_box,PIM_immune,ESTIMATE.immune.score_plot,PIM_stromal,ESTIMATE.stromal.score_plot,
     PIM_ESti,ESTIMATE.combined.score_plot,aveMethdata,avemethplot,file="DifferBetweenNMFsubtype_boxplot.Rdata")
