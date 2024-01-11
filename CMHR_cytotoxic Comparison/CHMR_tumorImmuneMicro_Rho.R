install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggbreak,ggsignif) 

######The figure shows the correlation between CHMR score and tumor immune microenvironment 
#boxplot
cyt_box<-ggplot(shengt,aes(fenzu,cyt,col=fenzu))+
    stat_boxplot(geom="errorbar",width=0.3,size=1.3,position = position_dodge(0.6),color=mypal3[c(3,1)])+ 
    geom_boxplot(size=1.3,width=0.65,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    scale_color_manual(values=mypal3[c(3,1)])+
    labs(x="",y="")+
    scale_x_discrete(limits=c("Low","High"))+
    ylim(0,max(shengt$cyt)+0.025)+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme(
        legend.position="none", 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"),
        axis.title.x=element_text(family="Arail",size=16,face="bold"),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),    
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+geom_signif(comparisons=compaired,step_increase = 0,map_signif_level = T,y_position=max(shengt$cyt)-0.5,textsize=6.5,parse=F,tip_length = 0,color="black")
cyt_box

#scatterplot
cor_tuT<-cor.test(as.numeric(shengt$CMHR),as.numeric(shengt$cyt),method="spearman")
CMHR_cytCor<-ggplot(data = shengt, mapping = aes(x = CMHR,
                                                 y = cyt))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=0.5,alpha=0.2)+
    scale_fill_distiller(palette = "YlOrBr")+
    ylim(0,2)+
    xlim(-2.5,3)+
    xlab(NULL)+
    ylab(NULL)+
    theme_few()+
    theme(legend.position = "none")+
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
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tuT$estimate,2),"p=",formatC(cor_tuT$p.value,2),sep=" "),family="Arail",
               parse=F,x=2,y=3,color="black",size=5.5,fontface="italic")

CMHR_cytCor

cor_tu1<-cor.test(as.numeric(CHMRscatter1[,1]),as.numeric(CHMRscatter1[,2]),method="spearman")
ScatterCHMR1<-ggplot(data =CHMRscatter1, mapping = aes(x =CMHR,y = CHMRscatter1[,2]))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=0.5,alpha=0.2)+
    scale_fill_distiller(palette = "YlOrBr")+
    xlab(NULL)+
    ylab(NULL)+
    xlim(-2,3.5)+
    theme_few()+
    theme(legend.position = "none")+
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
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu1$estimate,2),"p=",formatC(cor_tu1$p.value,3),sep=" "),family="Arail",
               parse=F,color="black",x=-1,y=3000,size=5.5,fontface="italic")
ScatterCHMR1

cor_tu2<-cor.test(as.numeric(CHMRscatter2[,1]),as.numeric(CHMRscatter2[,2]),method="spearman")
ScatterCHMR2<-ggplot(data =CHMRscatter2, mapping = aes(x =CMHR,y = CHMRscatter2[,2]))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=0.5,alpha=0.2)+
    scale_fill_distiller(palette = "YlOrBr")+
    xlab(NULL)+
    xlim(-2,3.5)+
    ylab(NULL)+
    theme_few()+
    theme(legend.position = "none")+
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
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu2$estimate,2),"p=",formatC(cor_tu2$p.value,3),sep=" "),family="Arail",
               parse=F,color="black",x=-1,y=3000,size=5.5,fontface="italic")
ScatterCHMR2

cor_tu3<-cor.test(as.numeric(CHMRscatter3[,1]),as.numeric(CHMRscatter3[,2]),method="spearman")
ScatterCHMR3<-ggplot(data =CHMRscatter3, mapping = aes(x =CMHR,y = CHMRscatter3[,2]))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=0.5,alpha=0.2)+
    scale_fill_distiller(palette = "YlOrBr")+
    xlab(NULL)+
    xlim(-2,3.5)+
    ylab(NULL)+
    theme_few()+
    theme(legend.position = "none")+
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
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu3$estimate,2),"p=",formatC(cor_tu3$p.value,3),sep=" "),family="Arail",
               parse=F,color="black",x=-1,y=3000,size=5.5,fontface="italic")
ScatterCHMR3

cor_tu4<-cor.test(as.numeric(CHMRscatter4[,1]),as.numeric(CHMRscatter4[,2]),method="spearman")
ScatterCHMR4<-ggplot(data =CHMRscatter4, mapping = aes(x =CMHR,y = CHMRscatter4[,2]))+
    stat_density_2d(geom = "polygon", 
                    aes(alpha = ..level..,fill = ..level..))+
    geom_point(size=0.5,alpha=0.2)+
    xlim(-2,3.5)+
    scale_fill_distiller(palette = "YlOrBr")+
    xlab(NULL)+
    ylab(NULL)+
    theme_few()+
    theme(legend.position = "none")+
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
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )+annotate("text",label=paste("R=",round(cor_tu4$estimate,2),"p=",formatC(cor_tu4$p.value,3),sep=" "),family="Arail",
               parse=F,color="black",x=-1,y=0.6,size=5.5,fontface="italic")
ScatterCHMR4

library(cowplot)
plot_grid(cyt_box,CMHR_cytCor+theme(axis.text.y = element_blank()),ScatterCHMR1+theme(axis.text.y = element_blank())
          ,ScatterCHMR2+theme(axis.text.y = element_blank()),ScatterCHMR3+theme(axis.text.y = element_blank()),
          ScatterCHMR4+theme(axis.text.y = element_blank()),ncol=3,align="hv")
