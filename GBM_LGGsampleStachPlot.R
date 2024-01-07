install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci) 

######The figure shows the sample numbers in different glioma phenotypes
mypal3=c("#003399","#CC3333","#CC6633")

#IDH status
IDH_GBMLGG<-ggplot(IDH_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=mypal3[c(1,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )

#MGMT promoter status
MGMT_GBMLGG<-ggplot(MGMT_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=mypal3[c(1,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )
#histology status
hist_GBMLGG<-ggplot(hist_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=c(mypal3[c(1,3,2)],"#CC99CC"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )
#original phenotype status
orig_GBMLGG<-
    ggplot(orig_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=c(mypal3[c(1,3,2)],"#CC99CC"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )
#TERT promoter status
TERT_GBMLGG<-
    ggplot(TERT_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=mypal3[c(1,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )

#G-CIMP phenotype status
CGIMP_GBMLGG<-
    ggplot(CGIMP_dui, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=mypal3[c(1,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )
#tumor grade
grade_GBMLGG<-ggplot(GL_cancerGrade, aes(x = Var.2, y = value, group = Var.1)) + 
    geom_col(aes(fill=Var.1)) +
    scale_fill_manual(values=mypal3[c(1,3,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )

#sample type
response_GBMLGG<-ggplot(GL_yuanfu, aes(x = Var.2, y = value, group = Var.1)) + 
    geom_col(aes(fill=Var.1)) +
    scale_fill_manual(values=mypal3[c(1,2)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none",
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        axis.line = element_line(color = "black"),
        )

######The figure shows the number of CpG sites or CpGct for their methylaton pattern
CMHC_siteduizhan<-ggplot(CM_celltype, aes(x = cell, y = number, group = type)) + 
    geom_col(aes(fill=type)) +
    scale_fill_manual(values=mypal3[c(1,3)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = number), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="top", 
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        axis.line = element_line(color = "black"),
    )
CMHC_siteduizhan

beta_siteduizhan<-ggplot(beta_celltype, aes(x = cell, y = number, group = type)) + 
    geom_col(aes(fill=type)) +
    scale_fill_manual(values=mypal3[c(1,3)])+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = number), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none", 
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        axis.line = element_line(color = "black"),
    )
beta_siteduizhan

library(cowplot)
plot_grid(beta_siteduizhan,CMHC_siteduizhan,beta_geneposDui,CMHC_geneposDui,
          ncol=2,align="hv",rel_widths = c(1,1))

######The figure shows the genomic location of CpG sites 
CMHC_geneposDui<-ggplot(CMHC_poscount, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=c(mypal3[1],"#FF9933","#0099CC","#FF9999","#336699","#9999CC"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="top", 
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        axis.line = element_line(color = "black"),
    )
CMHC_geneposDui

beta_geneposDui<-ggplot(beta_poscount, aes(x = Var.1, y = value, group = Var.2)) + 
    geom_col(aes(fill=Var.2)) +
    scale_fill_manual(values=c(mypal3[1],"#FF9933","#0099CC","#FF9999","#336699","#9999CC"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="top", 
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        axis.line = element_line(color = "black"),
    )
beta_geneposDui

save(IDH_GBMLGG,MGMT_GBMLGG,hist_GBMLGG,orig_GBMLGG,TERT_GBMLGG,CGIMP_GBMLGG,grade_GBMLGG,response_GBMLGG,file="phenotypeStackPlor.Rdata")
