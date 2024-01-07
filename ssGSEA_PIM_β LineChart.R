install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 

######The figure shows the CV differences of PIM score, DNA methylation levels, and ssGSEA score between NMF subtypes
mypal3=c("#003399","#CC3333","#CC6633")
CV_lineChart<-ggplot(data =CVnmf,aes(x=PIM_group,y=score,group = Score_type,color=NMF_Cluster))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    xlab("")+
    ylab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+scale_color_manual(values=mypal3)
CV_lineChart

######The figure shows the Entropy differences of PIM score, DNA methylation levels, and ssGSEA score between NMF subtypes
Entropy_lineChart<-ggplot(data =CVnmf,aes(x=PIM_group,y=score,group = Score_type,color=NMF_Cluster))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    xlab("")+
    ylab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+scale_color_manual(values=mypal3)
Entropy_lineChart

######The figure shows the Gini differences of PIM score, DNA methylation levels, and ssGSEA score between NMF subtypes
Gini_lineChart<-ggplot(data =CVnmf,aes(x=PIM_group,y=score,group = Score_type,color=NMF_Cluster))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    xlab("")+
    ylab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+scale_color_manual(values=mypal3)
Gini_lineChart

######The figure shows the PIM distribution between NMF subtypes
PIM_lineChart<-ggplot(data =guog,aes(x=PIM_group,y=Number,group = NMF_Cluster,color=NMF_Cluster))+
    geom_point(size=2)+
    geom_line(size=1.2)+
    xlab("")+
    ylab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="none", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+scale_color_manual(values=mypal3)
PIM_lineChart

save(CV_lineChart,Entropy_lineChart,Gini_lineChart,PIM_lineChart,file="LineChart.Rdata")
