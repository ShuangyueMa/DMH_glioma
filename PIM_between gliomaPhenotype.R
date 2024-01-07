install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,reshape,ggsignif) 


######The distribution curve shows the PIM score in three NMF subtypes
mypal3=c("#003399","#CC3333","#CC6633")
library(ggsignif)
PIM_distribution<-ggplot(data=PIM_disData, aes(x=PIM, color=NMF_cluster)) +
    geom_density(adjust=1.5, alpha=1,lwd=1.5,linetype = 1) +
    geom_rug()+
    guides(color = guide_legend(title = "NMF Cluster"))+
    theme(
        plot.margin = unit(c(0,0,0,0),"cm"),
        legend.position="none", 
        panel.background = element_blank(),
        legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
        axis.line = element_line(color = "black"),
        axis.text.x=element_text(colour="black",family="Arail",size=10,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=10,face="bold"), 
        axis.title.y=element_text(family="Arail",size=12,face="bold"), 
        axis.title.x=element_text(family="Arail",size=12,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
    )+scale_color_manual(values=mypal3)


######The boxplot shows the association between PIM score and different glioma phenotypes
#PIM differences in tumor types
com_CancerType<-list(c("GBM","LGG"))
GBM_V_LGGpim<-ggplot(gl_clinical,aes(Cancer.Type,PIM.score,col=Cancer.Type))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("cancer type")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_CancerType,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
GBM_V_LGGpim

#PIM differences in different IDH status
com_IDH<-list(c("WT","Mutant"))
IDHMUT_WTpim<-ggplot(IDH_clinical,aes(IDH.status,PIM.score,col=IDH.status))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("IDH status")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_IDH,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
IDHMUT_WTpim

#PIM differences in different MGMT promoter status
com_MGMT<-list(c("Methylated","Unmethylated"))
mgmt_pimbOX<-ggplot(MGMT_clinical,aes(MGMT.promoter.status,PIM.score,col=MGMT.promoter.status))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("MGMT promoter status")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_MGMT,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
mgmt_pimbOX

#PIM differences in different histological phenotypes
com_histo<-list(
    c("astrocytoma","glioblastoma"),
    c("glioblastoma","oligoastrocytoma"),
    c("oligoastrocytoma","oligodendroglioma"),
    c("astrocytoma","oligoastrocytoma"),
    c("glioblastoma","oligodendroglioma"),
    c("astrocytoma","oligodendroglioma")
    
)
histo_pimbOX<-ggplot(histo_clinical,aes(Histology,PIM.score,col=Histology))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=c(mypal3[c(1,3,2)],"#CC99CC"))+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("histology subtype")+
    scale_color_manual(values = c(mypal3[c(1,3,2)],"#CC99CC"))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_histo,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
histo_pimbOX

#PIM differences in different original phenotypes
com_orig<-list(
    c("Classical","IDHmut"),
    c("IDHmut","Mesenchymal"),
    c("Mesenchymal","Proneural"),
    c("Classical","Mesenchymal"),
    c("IDHmut","Proneural"),
    c("Classical","Proneural")
    
)
orig_pimbOX<-ggplot(orig_clinical,aes(Original.Subtype,PIM.score,col=Original.Subtype))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=c(mypal3[c(1,3,2)],"#CC99CC"))+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("origlogy subtype")+
    scale_color_manual(values = c(mypal3[c(1,3,2)],"#CC99CC"))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_orig,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
orig_pimbOX

#PIM differences in different TERT promoter status
com_TERT<-list(c("WT","Mutant"))
TERTMUT_WTpim<-ggplot(TERT_clinical,aes(TERT.promoter.status,PIM.score,col=TERT.promoter.status))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("TERT status")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_TERT,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
TERTMUT_WTpim

#PIM differences in different GCIMP phenotypes
com_GCIMP<-list(c("high","low"))
GCIMPMUT_WTpim<-ggplot(GCIMP_clinical,aes(G.CIMP.phenotype,PIM.score,col=G.CIMP.phenotype))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("GCIMP status")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_GCIMP,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
GCIMPMUT_WTpim

#PIM differences in different tumor grades
tumor_sC<-list(c("G2","G3"),c("G3","G4"),c("G2","G4"))
tumorStg<-ggplot(tumorStage,aes(grade,PIM,col=grade))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3)+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("tumor stage")+
    scale_x_discrete(limits=c("G2","G3","G4"))+
    scale_color_manual(values = mypal3[c(1,3,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=tumor_sC,step_increase = 0.1,map_signif_level = T,y_position=max(pp28$PIM),textsize=5.5,parse=F,color="black",tip_length = 0)
tumorStg

#PIM differences in different sample types
com_sampleType<-list(c("high","low"))
sampleTypepim<-ggplot(sampleTypeData,aes(sampleType,PIM,col=sampleType))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("sampleType")+
    scale_color_manual(values = mypal3[c(1,2)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=com_sampleType,step_increase = 0.1,map_signif_level = T,y_position=0.28,textsize=5.5,parse=F,color="black",tip_length = 0)
sampleTypepim

#PIM differences in different glioma patients treated with Temozolomide
Progre_com<-list(c("Clinical Progressive Disease","Stable Disease"))
TMZ_box<-ggplot(TMZ,aes(measure_of_response,PIM,col=measure_of_response))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=mypal3[c(1,3)])+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("Temozolomide")+
    scale_color_manual(values = mypal3[c(1,3)])+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=Progre_com,step_increase = 0.1,map_signif_level = T,y_position=max(xiang$PIM)+0.06,textsize=5.5,parse=F,color="black",tip_length = 0)
TMZ_box


#PIM differences in different chemotherapeutics
Chemo_Com<-list(c("Avastin","Bevacizumab"),
                c("Bevacizumab","carmustine"),
                c("carmustine","Temozolomide"),
                c("Avastin","carmustine"),
                c("Bevacizumab","Temozolomide"),
                c("Avastin","Temozolomide")
                
)

Chemo_Box<-ggplot(ChemoData,aes(drug_name,PIM,col=drug_name))+
    stat_boxplot(geom="errorbar",width=0.35,size=1.3,position = position_dodge(0.6),color=c(mypal3,"#6E79A7"))+ 
    geom_boxplot(size=1.5,fill="white",outlier.shape = NA,outlier.fill="white",outlier.color="white",alpha=0.65)+
    geom_jitter(stat="identity",width =0.2,shape = 20,size=2.5,alpha=0.65)+
    labs(x="",y="")+
    ggtitle("therapy")+
    scale_color_manual(values = c(mypal3,"#6E79A7"))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
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
    )+geom_signif(comparisons=Chemo_Com,step_increase = 0.1,map_signif_level = T,y_position=max(xiang$PIM)+0.06,textsize=5.5,parse=F,color="black",tip_length = 0)
Chemo_Box

#The figure show the association between PIM score and drug response 
resBar<-
    ggplot(res1, aes(x = Var2, y = Freq, group = Var1)) + 
    geom_col(aes(fill=Var1)) +
    scale_fill_manual(values=c("#CC6633","#CC99CC","pink"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = Freq), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
resBar

nonresBar<-ggplot(non1, aes(x = Var2, y = value, group = Var1)) + 
    geom_col(aes(fill=Var1)) +
    scale_fill_manual(values=c(mypal3[1],mypal3[2]))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label =value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
nonresBar

save(PIM_disData,PIM_distribution,
     gl_clinical,GBM_V_LGGpim,
     IDHMUT_WTpim,IDH_clinical,
     mgmt_pimbOX,MGMT_clinical,
     histo_pimbOX,histo_clinical,
     orig_pimbOX,orig_clinical,
     TERT_clinical,TERTMUT_WTpim,
     GCIMP_clinical,GCIMPMUT_WTpim,
     tumorStage,tumorStg,
     sampleTypeData,sampleTypepim,
     TMZ,TMZ_box,
     ChemoData,Chemo_Box,
     res1,resBar,
     non1,nonresBar,file="PIM_phenotypeBoxplot.Rdata"
)



