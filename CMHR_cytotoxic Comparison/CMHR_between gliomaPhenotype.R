install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggbreak,ggsignif) 

######The figure show the CHMR score,GZMA expression,PRF1 expression, and cytotoxic score across different glioma phenotypes
his_4 <- ggplot(qif5_1, aes(score,mean, fill = Histoloy)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    geom_errorbar(aes(ymin = mean -0.35*se, ymax = mean +0.35*se), width = 0.17, size = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values =c("#003399","#FF9966","#CC3333","#669999"))+
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
his_4

MGMT_4 <- ggplot(qif5_2, aes(score,mean, fill = MGMT)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    geom_errorbar(aes(ymin = mean -0.35*se, ymax = mean +0.35*se), width = 0.17, size = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values =c(mypal3[1],mypal3[2]))+
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
MGMT_4

IDH_4 <- ggplot(qif5_3, aes(score,mean, fill = IDH)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    geom_errorbar(aes(ymin = mean -0.35*se, ymax = mean +0.35*se), width = 0.17, size = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values =mypal3[c(1,2)])+
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
IDH_4

G60MGMT_plot <- ggplot(MGMT_G60, aes(score,mean, fill = MGMT)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    geom_errorbar(aes(ymin = mean -0.35*se, ymax = mean +0.35*se), width = 0.17, size = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values =mypal3[c(1,2)])+
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
G60MGMT_plot


G10MGMT_plot <- ggplot(MGMT_G10, aes(score,mean, fill = MGMT)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    geom_errorbar(aes(ymin = mean -0.35*se, ymax = mean +0.35*se), width = 0.17, size = 0.2, position = position_dodge(0.7)) +
    scale_fill_manual(values =mypal3[c(1,2)])+
    labs(title = NULL, x = NULL, y = NULL, fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_classic()+theme(legend.position = 0)+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.key.size = unit(0.8,"cm"),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14), 
        axis.text.y=element_text(colour="black",family="Arail",size=14), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
G10MGMT_plot

kruskal.test(as.numeric(yuce_TCGA$PRF1)~yuce_TCGA$Histoloy,data=yuce_TCGA)$p.value
kruskal.test(as.numeric(pheno10$CMHR)~pheno10$MGMT,data=pheno10)$p.value
kruskal.test(as.numeric(pheno60$CMHR)~pheno60$MGMT,data=pheno60)$p.value
library(cowplot)
plot_grid(his_4+theme(legend.position = "none"),IDH_4+theme(legend.position = "none"),MGMT_4+theme(legend.position = "none"),
          MGMTGSE_P+theme(legend.position = "none"),ncol=4)

save(geneCHMRCyt,HriskP,LriskP,gene_boxlist,
     cyt_box,CMHR_cytCor,ScatterCHMR1,ScatterCHMR2,ScatterCHMR2,
     ScatterCHMR3,ScatterCHMR4,his_4,IDH_4,MGMT_4,G10MGMT_plot,G60MGMT_plot,
     file="CHMR_geneExpression.Rdata")
