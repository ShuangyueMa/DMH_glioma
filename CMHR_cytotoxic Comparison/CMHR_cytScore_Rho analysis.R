install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggbreak) 

######The figure shows the correlation between gene expression and CHMR score as well as cytotoxic score
mypal3=c("#003399","#CC3333","#CC6633")
geneCMHRCyt_bar<-ggplot(geneCHMRCyt, aes(Ge,score,fill=type)) +
    geom_col(position = position_dodge(width = 0),width = 0.8, size = 0.5, colour = 'black') + 
    theme_classic()+ 
    coord_flip()+
    xlab("")+
    ylab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    ylim(-1,1)+
    scale_fill_manual(values=c("#F6C782","#6777B2","#CC6633","#CC3333","#003399"))+
    scale_y_break(c(-0.2, 0.2), scales = 1.5)+
    scale_x_discrete(labels=unique(geneCHMRCyt[order(geneCHMRCyt$Ge),1]))+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(
        legend.position="none", 
        legend.title =element_blank(), 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.line.y = element_blank(),
        axis.ticks=element_blank(), 
    )
geneCMHRCyt_bar    
