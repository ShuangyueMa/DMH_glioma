install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,pROC,cowplot) 

######The figure compares the predictive performance for CHMR score, cytotosic score, GZMA expression, and PRF1 expression
library(pROC)
TCGA_AUCCI<-data.frame()
for(j in 1:4){
    roclist <- roc(yuce_TCGA[,j] ~ CMHR+cyt+GZMA+PRF1,data= yuce_TCGA,legacy.axes=T,smooth=T,plot=T,ci=T,auc=T)
    auc1<- round(auc(roclist[[1]]),2)
    auc2<- round(auc(roclist[[2]]),2)
    auc3<- round(auc(roclist[[3]]),2)
    auc4<-round(auc(roclist[[4]]),2)
    lab<- paste0("AUC_CHMR = ",auc1, "\n", "AUC_cytotoxic = ",auc2, "\n", "AUC_GZMA = ",auc3,"\n","AUC_PRF1 =",auc4)
    p1<-ggroc(roclist,legacy.axes=T,linetype=1,size=1.5)+
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),color="black", linetype="dashed",size=0.5)+
        scale_colour_manual(values = c( mypal3,"#666699"))+
        theme_classic()+
        theme_bw()+
        theme(
            legend.title = element_blank(),
            legend.text = element_text(colour="black",family="Arail",size=14,face="bold"),
            legend.position="none",
            axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
            axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
            axis.title.y=element_text(family="Arail",size=16,face="bold"), 
            axis.title.x=element_text(family="Arail",size=16,face="bold"), 
            plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
            panel.border=element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
        )+annotate( "text",x = 0.75,y = 0.2,
                    label=lab,
                    size= 3)+ylab(NULL)+xlab(NULL)
    TCGA_ROC_plot[[j]]<-p1
    auc_ci<-data.frame()
    for(i in 1:4){
        auc_ci<-rbind(auc_ci,as.vector(as.numeric(roclist[[i]]$ci)))
    }
    colnames(auc_ci)<-c("l","auc","u")
    auc_ci$variable<-c("CMHR","cytotoxic","GZMA","PRF1")
    auc_ci$fenzu<-rep(colnames(yuce_TCGA)[j],4)
    TCGA_AUCCI<-rbind(TCGA_AUCCI,auc_ci)
}

######The figure shows confidence interval
library(ggplot2)
linecol<-c( mypal3,"#666699")
#TCGA
TCGA_CI<-ggplot(data=TCGA_AUCCI,aes(x=heng,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=rep(linecol,4),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",16),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=TCGA_AUCCI$variable)+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
    )+geom_hline(yintercept =0.7,linetype="dashed",size=0.7,color="gray60")

TCGA_CI

#GSE103659
GSE103659_ROC_plot<-list()
GSE103659_AUCCI<-data.frame()
for(j in 1:2){
    roclist <- roc(yuce_GSE103659[,j] ~ CMHR,data= yuce_GSE103659,legacy.axes=T,smooth=T,plot=T,ci=T,auc=T)
    auc1<- round(auc(roclist),2)
    lab<- paste0("AUC_CHMR = ",auc1, "\n")
    p1<-ggroc(roclist,legacy.axes=T,linetype=1,size=1.5)+
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),color="black", linetype="dashed",size=0.5)+
        scale_colour_manual(values = c( mypal3,"#666699"))+
        theme_classic()+
        theme_bw()+
        theme(
            legend.title = element_blank(),
            legend.text = element_text(colour="black",family="Arail",size=14,face="bold"),
            legend.position="none", 
            axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
            axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
            axis.title.y=element_text(family="Arail",size=16,face="bold"), 
            axis.title.x=element_text(family="Arail",size=16,face="bold"),
            plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
            panel.border=element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )+annotate( "text",x = 0.75,y = 0.2,label=lab,size= 3)+ylab(NULL)+xlab(NULL)
    GSE103659_ROC_plot[[j]]<-p1
    GSE103659_AUCCI<-rbind(GSE103659_AUCCI,as.numeric(roclist$ci))
    
}
colnames(GSE103659_AUCCI)<-c("lower","auc","upper")
GSE103659_AUCCI$variable<-c("OS","MGMT")

#GSE60274
GSE60274_ROC_plot<-list()
GSE60274_AUCCI<-data.frame()
for(j in 1:2){
    roclist <- roc(yuce_GSE60274[,j] ~ CMHR,data= yuce_GSE60274,legacy.axes=T,smooth=T,plot=T,ci=T,auc=T)
    auc1<- round(auc(roclist),2)
    lab<- paste0("AUC_CHMR = ",auc1, "\n")
    p1<-ggroc(roclist,legacy.axes=T,linetype=1,size=1.5)+
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),color="black", linetype="dashed",size=0.5)+
        scale_colour_manual(values = c( mypal3,"#666699"))+
        theme_classic()+
        theme_bw()+
        theme(
            legend.title = element_blank(),
            legend.text = element_text(colour="black",family="Arail",size=14,face="bold"),
            legend.position="none", 
            axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
            axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
            axis.title.y=element_text(family="Arail",size=16,face="bold"), 
            axis.title.x=element_text(family="Arail",size=16,face="bold"),
            plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
            panel.border=element_blank(),
            axis.line = element_line(color = "black"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
        )+annotate( "text",x = 0.75,y = 0.2,label=lab,size= 3)+ylab(NULL)+xlab(NULL)
    GSE60274_ROC_plot[[j]]<-p1
    GSE60274_AUCCI<-rbind(GSE60274_AUCCI,as.numeric(roclist$ci))
    
}
colnames(GSE60274_AUCCI)<-c("lower","auc","upper")
GSE60274_AUCCI$variable<-c("OS","MGMT")
GSE_AUC<-rbind(GSE103659_AUCCI,GSE60274_AUCCI)

GSE_AUC$heng<-letters[1:4]
GSE_CI<-ggplot(data=GSE_AUC,aes(x=heng,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=rep(mypal3[1],4),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",4),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=GSE_AUC$variable)+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
    )+geom_hline(yintercept =0.5,linetype="dashed",size=0.7,color="gray60")
all_CI<-plot_grid(TCGA_CI,GSE_CI,ncol=2,rel_widths = c(2,1),align="hv")
all_CI

#######The figure compares the predictive performance of eight CpGct and other CpG sites in previous study
library(ggplot2)
xianse<-c( mypal3,"#666699")
biOS_TCGA<-ggplot(data=TCGAbi[1:11,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[2],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=TCGAbi$zhushi[1:11])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.55,linetype="dashed",size=0.7,color="gray60")
biOS_TCGA

biMGMT_TCGA<-ggplot(data=TCGAbi[12:22,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[3],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=TCGAbi$zhushi[12:22])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.7,linetype="dashed",size=0.7,color="gray60")
biMGMT_TCGA

biOS_GSE10<-ggplot(data=GSE10bi[1:11,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[3],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=GSE10bi$zhushi[1:11])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.55,linetype="dashed",size=0.7,color="gray60")
biOS_GSE10

biMGMT_GSE10<-ggplot(data=GSE10bi[12:22,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[3],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=GSE10bi$zhushi[12:22])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.55,linetype="dashed",size=0.7,color="gray60")
biMGMT_GSE10


biOS_GSE60<-ggplot(data=GSE60bi[1:11,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[3],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=GSE60bi$zhushi[1:11])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.55,linetype="dashed",size=0.7,color="gray60")
biOS_GSE60

biMGMT_GSE60<-ggplot(data=GSE60bi[12:22,],aes(x=site,y=auc))+
    geom_errorbar(aes(ymin=lower,ymax=upper),width=0,color=c(rep(mypal3[1],8),rep(mypal3[3],3)),
                  size=1.25)+
    geom_point(size=3.8, color=rep("black",11),shape=15)+
    theme_classic()+
    scale_x_discrete(labels=GSE60bi$zhushi[12:22])+
    theme_bw()+
    xlab(NULL)+
    ylab(NULL)+
    theme(
        legend.position="none", 
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(angle=90,colour="black",family="Arail",size=8,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=12,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
    )+geom_hline(yintercept =0.55,linetype="dashed",size=0.7,color="gray60")
biMGMT_GSE60


library(cowplot)
plot_grid(biOS_TCGA,biMGMT_TCGA,biOS_GSE10,biMGMT_GSE10,biOS_GSE60,biMGMT_GSE60,
          ncol=2,rel_heights = c(1,1,1),align="hv")

save(all_CI,TCGA_ROC_plot,biOS_GSE10,biOS_GSE60,biOS_TCGA,
     biMGMT_GSE10,biMGMT_GSE60,biMGMT_TCGA,file="CHMR_ROC.Rdata")
