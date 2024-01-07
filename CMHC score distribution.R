install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 
mypal3=c("#003399","#CC3333","#CC6633")

######The figure shows β-value and CMHC cell score for six immune cells in different genome position
CMHC_sitebeta_plot<-ggplot(data=CMHC_sitebeta,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCmono")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+scale_color_manual(values=c(mypal3))
CMHC_sitebeta_plot

CMHCDistrimono_plot<-ggplot(data=CMHCDistrimono,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCmono")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+scale_color_manual(values=c(mypal3))
CMHCDistrimono_plot


CMHCDistriCD8_plot<-ggplot(data=CMHCDistriCD8,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCCD8")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2)
    )+scale_color_manual(values=c(mypal3))
CMHCDistriCD8_plot

CMHCDistriCD4_plot<-ggplot(data=CMHCDistriCD4,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCCD4")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        
    )+scale_color_manual(values=c(mypal3))
CMHCDistriCD4_plot


CMHCDistriB_plot<-ggplot(data=CMHCDistriB,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCB")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+scale_color_manual(values=c(mypal3))
CMHCDistriB_plot

CMHCDistriB_plot<-ggplot(data=CMHCDistriB,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCB")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+scale_color_manual(values=c(mypal3))
CMHCDistriB_plot

CMHCDistrigran_plot<-ggplot(data=CMHCDistrigran,aes(x=shunxu,y=mean,color=group))+
    geom_line(aes(group=group),size=1.5)+ 
    geom_point(shape=16,size=3,stroke=1)+
    ylab("CMHCB")+
    xlab(NULL)+
    theme_bw()+
    theme(
        legend.position="none",
        axis.text.x=element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+scale_color_manual(values=c(mypal3))
CMHCDistrigran_plot

library(cowplot)
plot_grid(CMHC_sitebeta_plot,CMHCDistriB_plot,CMHCDistriCD4_plot,CMHCDistriCD8_plot,CMHCDistrigran_plot,CMHCDistrimono_plot,
          ncol=1,rel_heights = c(1,1,1,1,1,1),align="hv")