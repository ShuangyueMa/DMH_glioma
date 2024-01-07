install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 
mypal3=c("#003399","#CC3333","#CC6633")

######The figure shows association between β values of Cell-specific CpG sites and gene expression
PALM=data.frame(cg1=as.numeric(TCGA[match(geneCor[i,1],rownames(TCGA)),]),
                ge=as.numeric(EXP[match(geneCor[i,2],rownames(EXP)),match(colnames(TCGA),colnames(EXP))]))
PALM_scatter<-ggplot(PALM,aes(cg1,ge))+
    geom_point(aes(color=cg1),size=4)+
    guides(color="none")+ 
    scale_color_gradient(low=mypal3[3],high=mypal3[2])+
    geom_smooth(method =lm,se=T,color=mypal3[1],size=1.2)+
    xlab(geneCor[i,1])+
    ylab(geneCor[i,2])+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme(
        plot.margin = unit(c(0.45,0.45,0.45,0.45),"cm"),
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
    )+annotate("text",label=paste("R=",geneCor[i,4],"p=",formatC(geneCor[i,5]),sep="  "),family="Arail",
               parse=F,x=max(data[,1]-0.5),y=max(data[,2]-0.5),color="black",size=5.5,fontface="italic")

PALM_scatter

SCAMP4_scatter<-ggplot(SCAMP4,aes(cg1,ge))+
    geom_point(aes(color=cg1),size=4)+
    guides(color="none")+ 
    scale_color_gradient(low=mypal3[3],high=mypal3[2])+
    geom_smooth(method =lm,se=T,color=mypal3[1],size=1.2)+
    xlab(geneCor[i,1])+
    ylab(geneCor[i,2])+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme(
        plot.margin = unit(c(0.45,0.45,0.45,0.45),"cm"),
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
    )+annotate("text",label=paste("R=",geneCor[i,4],"p=",formatC(geneCor[i,5]),sep="  "),family="Arail",
               parse=F,x=max(data[,1]-0.5),y=max(data[,2]-0.5),color="black",size=5.5,fontface="italic")

SCAMP4_scatter