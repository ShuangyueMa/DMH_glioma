#####show colors
install.packages("colourpicker")
library(colourpicker)
install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif)


######The figure shows the distribution of β values in different cell types and tumor bulk tissue
library(ggplot2)
#CD19 Bcell
Bcell_plot=ggplot()+geom_histogram(data = Bcell,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("B cell")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
Bcell_plot

#Tcell
Tcell_plot=ggplot()+geom_histogram(data = Tcell,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("T cell")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
Tcell_plot

#macrophage
macro_plot=ggplot()+geom_histogram(data = macro,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("Macrophage")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
macro_plot

#monocyte
mono_plot=ggplot()+geom_histogram(data = mono,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("Monocyte")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
mono_plot

#neutrophil
neutro_plot<-ggplot()+geom_histogram(data =neu,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("Neutrophil")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())

neutro_plot

#NK cell
NK_plot<-ggplot()+geom_histogram(data =NKcell,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("Natural Killer cell")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())#
NK_plot


##cell line
#U-87
u871_plot=ggplot()+geom_histogram(data = u871,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("U87")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
u871_plot

#U-251
U251=data.frame(data=as.numeric(beta.m[,1]))
u251_plot=ggplot()+geom_histogram(data =U251,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("U251")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top",
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
u251_plot


##tissue
tissue_plot<-ggplot()+geom_histogram(data = tumor_tissue,aes(x=data),stat = 'bin',bins =30,fill='gray',color='black',size=0.25)+
    geom_density(data=tumor_tissue,aes(x=data))+
    theme_bw()+
    theme(axis.title = element_blank())+
    scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1.0),labels=c("0","0.2","0.4","0.6","0.8","1.0"))+
    scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle("tumor tissue")+
    geom_vline(xintercept = c(0.2,0.6),linetype="dashed",color="gray25",size=1)+
    theme(
        plot.margin = unit(c(0.48,0.25,0.25,0.25),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.ticks=element_blank())
tissue_plot

library(cowplot)
library(ggpubr)
library(patchwork)
cell_9<-plot_grid(Bcell_plot,Tcell_plot,NK_plot,
                  macro_plot,mono_plot,neutro_plot,u871_plot,u251_plot,
                  tissue_plot,ncol=3,rel_widths = c(1,1,1),rel_heights = c(1,1,1),
                  align="hv")
cell_patch<-annotate_figure(cell_9,left = text_grob('Frequency',
                                                       size = 15, rot = 90),
                               bottom=text_grob('DNA methylation level (β)',
                                                size = 15, rot = 0))
save(Bcell,Tcell,macro,mono,neu,NKcell,u871,u251,tumor_tissue,
    cell_patch,Bcell_plot,Tcell_plot,NK_plot,macro_plot,mono_plot,neutro_plot,
     u871_plot,u251_plot,tissue_plot,file="InterMethyHistogram.Rdata")
