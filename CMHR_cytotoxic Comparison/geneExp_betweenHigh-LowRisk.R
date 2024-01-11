######The figure shows gene expression differences between high- and low-risk groups
install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggbreak,ggsignif) 

gene_boxlist<-list()
for(i in 1:25){
    box_160<-EXP[na.omit(match(pipeiindex[i],rownames(EXP))),na.omit(match(rownames(shengt),colnames(EXP)))]
    box_160=data.frame(gene=as.numeric(box_160),sample=colnames(box_160))
    box_160$group<-shengt$fenzu
    y.v<-max(box_160$gene)
    compaired<-list(c("Low","High"))
    jiyin_box<-ggplot(box_160,aes(group,gene,col=group))+
        stat_boxplot(geom="errorbar",width=0.3,size=1.3,position = position_dodge(0.6),color=mypal3[c(3,1)])+ 
        geom_boxplot(size=1.3,width=0.65,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
        geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
        scale_color_manual(values=mypal3[c(3,1)])+
        labs(x="",y="")+
        scale_x_discrete(limits=c("Low","High"))+
        ggtitle(pipeiindex[i])+
        ylim(0,y.v+0.025)+
        theme_classic()+theme(legend.position = 0)+
        theme_bw()+
        theme(
            legend.position="none", 
            axis.text.x = element_blank(),
            plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
            panel.border = element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )+geom_signif(comparisons=compaired,step_increase = 0,map_signif_level = T,y_position=y.v-1.5,textsize=6.5,parse=F,tip_length = 0,color="black")
    
    gene_boxlist[[i]]<-jiyin_box
}

library(cowplot)
plot_grid(gene_boxlist[[1]],gene_boxlist[[2]],gene_boxlist[[3]],gene_boxlist[[4]],gene_boxlist[[5]],
          gene_boxlist[[6]],gene_boxlist[[7]],gene_boxlist[[8]],gene_boxlist[[9]],gene_boxlist[[10]],
          gene_boxlist[[11]],gene_boxlist[[12]],gene_boxlist[[13]],gene_boxlist[[14]],gene_boxlist[[15]],
          gene_boxlist[[16]],gene_boxlist[[17]],gene_boxlist[[18]],gene_boxlist[[19]],gene_boxlist[[20]],
          gene_boxlist[[21]],gene_boxlist[[22]],gene_boxlist[[23]],gene_boxlist[[24]],gene_boxlist[[25]],
          ncol=5,nrow=5,algn="hv",rel_widths = c(1,1,1,1,1),rel_heights = c(1,1,1,1,1))


