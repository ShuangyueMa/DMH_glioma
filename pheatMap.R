#####show colors
install.packages("colourpicker")
library(colourpicker)
mypal3=c("#003399","#CC0000","#CC6633")

######The figure shows the landscape of the immune cell infiltration levels and glioma phenotypes
library(pheatmap)
Figure1_anno_colors<-list(
    Histology=c(astrocytoma="#144E77",glioblastoma="#1D6CA8",oligoastrocytoma="#85AFD0",oligodendroglioma="#C8DDE9"),
    Original.Subtype=c(Classical="#482E86",IDHmut="#5D4592",Mesenchymal="#7982B8",Proneural="#B4B6D6"),
    CancerType=c(GBM="#6B7CBC",LGG="#206CA6"),
    NMF_cluster= c(Cluster1="#003399",Cluster2="#CC0000",Cluster3="#CC6633"),
    IDH.status=c(Mutant="#E1A09E",WT="#0066CC"),
    MGMT.promoter.status=c(Methylated="#206CA6",Unmethylated="#E1A09E"),
    TERT.promoter.status=c(Mutant="#7982B8",WT="#E79F7B"),
    G.CIMP.phenotype=c(high="#83B0D4",low="#B4B6D6"),
    SampleType=c(Primary="#E26F72",Recurrent="#E5A278"),
    tumorgrade=c(G2="#206CA6",G3="#7982B8",G4="skyblue"),
    Response=c(CR="#2597CA",
               PD="#85C6E0",
               PR="#7FAED4",
               SD="#C9DEEA"))
Figure1_pheat<-pheatmap(scale(gbm2),
                        scale="row",
                        annotation_col =chongzhushi,
                        annotation_colors = Figure1_anno_colors,
                        clustering_distance_rows = "euclidean",
                        clustering_method = "complete",
                        cluster_cols = F,
                        cluster_rows = T,
                        fontface="bold",
                        fontfamily= "Arail",
                        fontsize =12.5,
                        fontsize_row =11,
                        width=1500,
                        heigth=1200,
                        show_colnames = F,
                        margin=10,
                        color = c(colorRampPalette(colors = c("#003399","white"))(length(bk)/2),colorRampPalette(colors = c("white","#CC0000"))(length(bk)/2)),
    )

######The figure shows the overlap between NMF subtype and other glioma phenotypes
sankedata
library(ggforce)
sangji<-ggplot(sankedata, aes(x, id = id, split = y, value =intersect)) +
    geom_parallel_sets(aes(fill = NMF_cluster), alpha = 0.85, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.45,fill="white",color="black",size=1) +
    geom_parallel_sets_labels(colour = 'black',angle = 0,size=4) +
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    xlab(NULL)+
    scale_fill_manual(values=mypal3)+
    scale_x_discrete(labels = c("glioblastoma grade","NMF Cluster","Other Subtype")) +
    ggtitle("Overlap of glioblastome NMF Cluster and other Subtype Classification")+
    theme(
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y = element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )

######The figure shows DNA methylation level of Cell-type-related CpG sites
library(scales)
show_col(cellyanse)
anno_cell<-list(cell=c(B=mypal3[1],CD4=mypal3[2],CD8=mypal3[3],Gran="#6495ED",Mono=cellyanse[1],Neut=cellyanse[7]),GSE=c(GSE103541=mypal3[1],GSE166844=mypal3[2],GSE118144=mypal3[3]),SiteType=c(Hyper=mypal3[1],Hypo=mypal3[3]))
bk <- c(seq(-9,-0.1,by=0.01),seq(0,9,by=0.01))
library(pheatmap)
cellsite_ht2<-pheatmap(scale(cellorderEXP),
                       scale="row",
                       annotation_col =group,
                       annotation_row = rowgoup,
                       annotation_colors = anno_cell,
                       clustering_distance_rows = "euclidean",
                       clustering_method = "complete",
                       legend_breaks = c(-10, -5, 0, 5,10),
                       cluster_cols = F,
                       cluster_rows = F,
                       fontface="bold",
                       fontfamily= "Arail",
                       main="cell specifically heterogeneous site",
                       show_colnames = F,
                       show_rownames = F,
                       margin=0.6,
                       color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)))


######The figure shows gene expression of DEG between NMF subtypes
library(gplots)
library(RColorBrewer) 
plot_c1<-c(rep("#297BB8",240),rep("#FAAD60",291),rep("#EF9C92",171))
heatmap.2(as.matrix(DEGrexp),density.info=c("none"),  ColSideColors=plot_c1,  
          dendrogram="col",rowsep = c(172,704),key.xlab = "expression-Zscore",
          colsep = c(240,531),cexRow = 1.3,cexCol = 0.8,main ="DEG between NMF cluster",labCol =FALSE,key.title="expression",
          margins=c(4,17),col=rev(sitecol),colge,scale="row",trace="none",breaks = seq(-0.5,0.5,0.05),Rowv = FALSE,Colv = FALSE)


######The figure shows DNA methylation level of 8 survival-related CpGct for long- and short-survival patients
require(GEOquery)   
require(Biobase)

gy34_data_S<-weidian207[,match(as.character(as.matrix(pheno2074%>%filter(therapy=="34Gy",live=="SS")%>%select(sample))),colnames(weidian207))]
gy34_data_L<-weidian207[,match(as.character(as.matrix(pheno2074%>%filter(therapy=="34Gy",live=="LS")%>%select(sample))),colnames(weidian207))]
gy34_data_S<-gy34_data_S[order(apply(gy34_data_S,1,mean),decreasing = T),]
gy34_data_L<-gy34_data_L[order(apply(gy34_data_L,1,mean),decreasing = T),]

gy34_SP<-Heatmap(as.matrix(gy34_data_S),
                 rect_gp = gpar(col = "white", lwd = 1.2),
                 show_column_names = FALSE,
                 row_order = rownames(gy34_data_S),
                 row_dend_width = unit(2, "cm"),
                 show_heatmap_legend = FALSE,
                 column_dend_reorder = FALSE)
gy34_SP
gy34_LP<-Heatmap(as.matrix(gy34_data_L),
                 rect_gp = gpar(col = "white", lwd = 1.2),
                 row_names_side = "left",
                 show_column_names = FALSE,
                 row_order = rownames(gy34_data_L),
                 row_dend_width = unit(2, "cm"),
                 show_heatmap_legend = FALSE,
                 column_dend_reorder = FALSE)
gy34_LP

gy34_qipao=data.frame(live=c(rep("shortlive",8),rep("longlive",8)),site=c(rownames(gy34_data_S),rownames(gy34_data_L)),mean=c(
    as.numeric(apply(gy34_data_S,1,mean),apply(gy34_data_L,1,mean))
))
se8<-c("#0066CC","#663399","#CC3366","#CC9933","#FF9999","#FF9966","#339933","#0099CC")
gy34_qi<-ggplot(gy34_qipao,aes(live,site))+
    geom_point(data=gy34_qipao,aes(live,site,size=3,color=mean),
               shape=16,show.legend=T)+
    
    scale_color_gradient(low =se8[7] ,high =se8[3])+
    labs(size="",color="average β-value")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        legend.position = "top",
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.text.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
    )

gy34_qi


gy60_data_S<-weidian207[,match(as.character(as.matrix(pheno2074%>%filter(therapy=="60Gy",live=="SS")%>%select(sample))),colnames(weidian207))]
gy60_data_L<-weidian207[,match(as.character(as.matrix(pheno2074%>%filter(therapy=="60Gy",live=="LS")%>%select(sample))),colnames(weidian207))]
gy60_data_S<-gy60_data_S[order(apply(gy60_data_S,1,mean),decreasing = T),]
gy60_data_L<-gy60_data_L[order(apply(gy60_data_L,1,mean),decreasing = T),]

gy60_SP<-Heatmap(as.matrix(gy60_data_S),
                 rect_gp = gpar(col = "white", lwd = 1.2),
                 show_column_names = FALSE,
                 row_order = rownames(gy60_data_S),
                 row_dend_width = unit(2, "cm"),
                 show_heatmap_legend = FALSE,
                 column_dend_reorder = FALSE)
gy60_SP
gy60_LP<-Heatmap(as.matrix(gy60_data_L),
                 rect_gp = gpar(col = "white", lwd = 1.2),
                 row_names_side = "left",
                 show_column_names = FALSE,
                 row_order = rownames(gy60_data_L),
                 row_dend_width = unit(2, "cm"),
                 show_heatmap_legend =FALSE,
                 column_dend_reorder = FALSE)
gy60_LP

gy60_qipao=data.frame(live=c(rep("shortlive",8),rep("longlive",8)),site=c(rownames(gy60_data_S),rownames(gy60_data_L)),mean=c(
    as.numeric(apply(gy60_data_S,1,mean),apply(gy60_data_L,1,mean))
))
se8<-c("#0066CC","#663399","#CC3366","#CC9933","#FF9999","#FF9966","#339933","#0099CC")
gy60_qi<-ggplot(gy60_qipao,aes(live,site))+
    geom_point(data=gy60_qipao,aes(live,site,size=3,color=mean),
               shape=16,show.legend=T)+
    
    scale_color_gradient(low =se8[7] ,high =se8[3])+
    labs(size="",color="average β-value")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        legend.position = "top",
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.text.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
    )

gy60_qi


TCGA_noTMZSS<-yaowu1 %>% filter(treat=="no TMZ") %>% filter(fenzu=="1") #22 10
TCGA_noTMZLS<-yaowu1 %>% filter(treat=="no TMZ") %>% filter(fenzu=="2") #29 10 
noTMZss<-methjisuan[match(weidian8[order(apply(methjisuan[match(weidian8,rownames(methjisuan)),match(TCGA_noTMZSS$sample,colnames(methjisuan))],1,mean),decreasing = T)],rownames(methjisuan))
                    ,match(TCGA_noTMZSS$sample,colnames(methjisuan))]
noTMZls<-methjisuan[match(weidian8[order(apply(methjisuan[match(weidian8,rownames(methjisuan)),match(TCGA_noTMZSS$sample,colnames(methjisuan))],1,mean),decreasing = T)],rownames(methjisuan))
                    ,match(TCGA_noTMZSS$sample,colnames(methjisuan))]


TCGA_noRTSS<-FANG %>% filter(RT=="NO") %>% filter(fenzu2=="1") #55 14
TCGA_noRTLS<-FANG %>% filter(RT=="NO") %>% filter(fenzu2=="2") #53 14
noRTss<-methjisuan[match(weidian8[order(apply(methjisuan[match(weidian8,rownames(methjisuan)),match(TCGA_noRTSS$sample,colnames(methjisuan))],1,mean),decreasing = T)],rownames(methjisuan))
                   ,match(TCGA_noRTSS$sample,colnames(methjisuan))]
noRTls<-methjisuan[match(weidian8[order(apply(methjisuan[match(weidian8,rownames(methjisuan)),match(TCGA_noRTLS$sample,colnames(methjisuan))],1,mean),decreasing = T)],rownames(methjisuan))
                   ,match(TCGA_noRTLS$sample,colnames(methjisuan))]


TCGAnoTMZ_SP<-Heatmap(as.matrix(noTMZss),
                      rect_gp = gpar(col = "white", lwd = 1.2),
                      show_column_names = FALSE,
                      row_order = rownames(noTMZss),
                      row_dend_width = unit(2, "cm"),
                      show_heatmap_legend = FALSE,
                      column_dend_reorder = FALSE)
TCGAnoTMZ_SP
TCGAnoTMZ_LP<-Heatmap(as.matrix(noTMZls),
                      rect_gp = gpar(col = "white", lwd = 1.2),
                      row_names_side = "left",
                      show_column_names = FALSE,
                      row_order = rownames(noTMZls),
                      row_dend_width = unit(2, "cm"),
                      show_heatmap_legend = FALSE,
                      column_dend_reorder = FALSE)
TCGAnoTMZ_LP

TCGAnoRT_SP<-Heatmap(as.matrix(noRTss),
                     rect_gp = gpar(col = "white", lwd = 1.2),
                     show_column_names = FALSE,
                     row_order = rownames(noRTss),
                     row_dend_width = unit(2, "cm"),
                     show_heatmap_legend = FALSE,
                     column_dend_reorder = FALSE)
TCGAnoRT_SP
TCGAnoRT_LP<-Heatmap(as.matrix(noRTls),
                     rect_gp = gpar(col = "white", lwd = 1.2),
                     row_names_side = "left",
                     show_column_names = FALSE,
                     row_order = rownames(noRTls),
                     row_dend_width = unit(2, "cm"),
                     show_heatmap_legend = FALSE,
                     column_dend_reorder = FALSE)
TCGAnoRT_LP


library(ggplot2)
library(reshape)
TCGAnoTMZ_qipao=data.frame(live=c(rep("shortlive",8),rep("longlive",8)),site=c(rownames(noTMZss),rownames(noTMZls)),mean=c(
    as.numeric(apply(noTMZss,1,mean),apply(noTMZls,1,mean))
))
se8<-c("#0066CC","#663399","#CC3366","#CC9933","#FF9999","#FF9966","#339933","#0099CC")
TCGAnoTMZ_qi<-ggplot(TCGAnoTMZ_qipao,aes(live,site))+
    geom_point(data=TCGAnoTMZ_qipao,aes(live,site,size=3,color=mean),
               shape=16,show.legend=T)+
    
    scale_color_gradient(low =se8[7] ,high =se8[3])+
    labs(size="",color="average β-value")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        legend.position = "top",
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.text.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
    )

TCGAnoTMZ_qi


TCGAnoRT_qipao=data.frame(live=c(rep("shortlive",8),rep("longlive",8)),site=c(rownames(noRTss),rownames(noRTls)),mean=c(
    as.numeric(apply(noRTss,1,mean),apply(noRTls,1,mean))
))
se8<-c("#0066CC","#663399","#CC3366","#CC9933","#FF9999","#FF9966","#339933","#0099CC")

TCGAnoRT_qi<-ggplot(TCGAnoRT_qipao,aes(live,site))+
    geom_point(data=TCGAnoRT_qipao,aes(live,site,size=3,color=mean),
               shape=16,show.legend=T)+
    scale_color_gradient(low =se8[7] ,high =se8[3])+
    labs(size="",color="average β-value")+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
        legend.position = "top",
        axis.ticks= element_blank(),
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.text.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
    )

TCGAnoRT_qi





library(cowplot)
library(patchwork)
library(ggplotify)

plot_grid(as.ggplot(TCGAnoTMZ_SP),
          as.ggplot(TCGAnoTMZ_LP),
          as.ggplot(TCGAnoRT_SP),
          as.ggplot(TCGAnoRT_LP),
          as.ggplot(TMZ_SP),
          as.ggplot(TMZ_LP),
          as.ggplot(gy34_SP),
          as.ggplot(gy34_LP),
          as.ggplot(gy60_SP),
          as.ggplot(gy60_LP),
          ncol=2
)

#######The figure shows proportition of CpGct with different methylation changes
TMZdui<-data.frame(zhiliao=c(rep("no",3),rep("yes",3)),leibie=rep(c("down","incre","no"),2),value=c(2,3,3,4,2,2))
TMZduiP<-ggplot(TMZdui, aes(x = zhiliao, y = value, group = leibie)) + 
    geom_col(aes(fill=leibie)) +
    scale_x_discrete(labels=c("no Temozolomide","Temozolomide"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    #coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none", 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
        
    )

RTdui1<-data.frame(zhiliao=c(rep("no",3),rep("yes",3)),leibie=rep(c("down","incre","no"),2),value=c(3,2,3,4,3,1))
RTduiP1<-ggplot(RTdui1, aes(x = zhiliao, y = value, group = leibie)) + 
    geom_col(aes(fill=leibie)) +
    scale_fill_manual(values=c("#CC9933","#CC3366","black"))+
    scale_x_discrete(labels=c("no radiation","radiation 34Gy"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    #coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none", 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
    )


RTdui2<-data.frame(zhiliao=c(rep("no",3),rep("yes",3)),leibie=rep(c("down","incre","no"),2),value=c(3,2,3,3,3,2))
RTduiP2<-ggplot(RTdui2, aes(x = zhiliao, y = value, group = leibie)) + 
    geom_col(aes(fill=leibie)) +
    scale_fill_manual(values=c("#CC9933","#CC3366","black"))+
    scale_x_discrete(labels=c("no radiation","radiation 34Gy"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.25)) +
    geom_text(aes(label = value), position = position_stack(vjust = .5), size = 3)+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    #coord_flip()+
    theme_bw()+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        legend.title = element_blank(),
        legend.position="none", 
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
        
    )


plot_grid(TMZduiP,RTduiP1,RTduiP2,ncol=3,align = "hv",rel_widths = c(1,1,1))
save(
TCGAnoTMZ_SP,TCGAnoTMZ_LP,TCGAnoRT_SP,
TCGAnoRT_LP,TMZ_SP,TMZ_LP,gy34_SP,
gy34_LP,gy60_SP,gy60_LP,TMZduiP,cellsite_ht2,DEGrexp,
RTduiP1,RTduiP2,file="pheatmap.Rdata")
