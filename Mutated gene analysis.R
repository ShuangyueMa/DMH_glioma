library(maftools)
library(ggplot2)
library(ComplexHeatmap)
library(ggsignif)
mypal3=c("#003399","#CC3333","#CC6633")

######The oncoplot shows mutated genes in high-risk and low-risk patients
mut.high<-mut[mut$Tumor_Sample_Barcode %in% shengt[shengt$fenzu=="High",5],]
mut.low<-mut[mut$Tumor_Sample_Barcode %in% shengt[shengt$fenzu=="Low",5],]
dim(mut.high)
dim(mut.low)
high_tubianMAF<-read.maf(mut.high,isTCGA = TRUE)
low_tubianMAF<-read.maf(mut.low,isTCGA = TRUE)
col<-c("#297BB8","#3580BA","#75B789","#F6A6A7","#9E6BAA","skyblue","red","#EF9C92","#FAAD60")
names(col)<-unique(mut$Variant_Classification)

high_tubianMAF@clinical.data$OS<-senlin4Re[match(high_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),1]
high_tubianMAF@clinical.data$tumor.grade<-senlin4Re[match(high_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),24]
high_tubianMAF@clinical.data$MGMT<-senlin4Re[match(high_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),25]
high_tubianMAF@clinical.data$IDH<-senlin4Re[match(high_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),26]
high_tubianMAF@clinical.data$NMF<-senlin4Re[match(high_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),27]
head(high_tubianMAF@clinical.data)
low_tubianMAF@clinical.data$OS<-senlin4Re[match(low_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),1]
low_tubianMAF@clinical.data$tumor.grade<-senlin4Re[match(low_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),24]
low_tubianMAF@clinical.data$MGMT<-senlin4Re[match(low_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),25]
low_tubianMAF@clinical.data$IDH<-senlin4Re[match(low_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),26]
low_tubianMAF@clinical.data$NMF<-senlin4Re[match(low_tubianMAF@clinical.data$Tumor_Sample_Barcode,senlin4Re$bar),27]

OS_col<-c("#297BB8","#75B789")
names(OS_col)<-c("1","0")
grade_col<-c("#666699","#FF9966","#CC6666")
names(grade_col)<-c("G2","G3","G4")
MGMT_col<-c("#3399CC","#CC3333")
names(MGMT_col)<-c("Methylated","Unmethylated")
IDH_col<-c("#9E6BAA","#0066CC")
names(IDH_col)<-c("Mutant","WT")
NMF_col<-mypal3
names(NMF_col)<-c("Cluster1","Cluster2","Cluster3")
fabcolors = list(OS = OS_col,tumor.grade=grade_col,MGMT=MGMT_col,IDH=IDH_col,
                 NMF=NMF_col)


oncoplot(high_tubianMAF,
         top=25,
         logColBar = T,
         colors = col,
         bgCol = "white",
         legendFontSize = 1.2,
         fontSize = 1.15,
         titleFontSize = 1.8,
         clinicalFeatures = c('OS','tumor.grade','MGMT',"IDH","NMF"),
         annotationColor = fabcolors
         
)


oncoplot(low_tubianMAF,
         top=25,
         logColBar = T,
         colors = col,
         bgCol = "white",
         legendFontSize = 1.2,
         fontSize = 1.1,
         titleFontSize = 1.8,
         clinicalFeatures = c('OS','tumor.grade','MGMT',"IDH","NMF"),
         annotationColor = fabcolors
         
)

######The figure show mutated gene expressions for high- and low-risk patients
highcom<-list(c("a1", "a2"),
              c("b1", "b2"),
              c("c1", "c2"),
              c("d1", "d2"),c("e1", "e2"),c("f1", "f2"),c("g1", "g2"),c("h1", "h2"),c("i1", "i2"),c("j1", "j2"),
              c("k1", "k2"),c("l1", "l2"),c("m1", "m2"),c("n1", "n2"),c("p1", "p2"))
HighmutgenePlot<-ggplot(Highdata)+
    geom_boxplot(aes(group, value, fill = fenzu, color = fenzu), 
                 width = 0.5,
                 outlier.shape = NA,
                 size=0,
    )+
    ylab("gene expression (FPKM)")+
    scale_x_discrete(labels = unique(Highdata$gene))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values = mypal3[c(3,1)])+
    scale_color_manual(values = c("black","black"))+
    theme_bw()+
    ylim(0,8.5)+
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          axis.title = element_blank(),
          panel.grid.minor.x  =element_blank(),
          axis.text.x = element_text(angle=45,hjust=0,vjust=0.2),
          axis.line = element_line(color = "black"))+
    geom_signif(aes(group, value),
                comparisons = highcom,
                tip_length = rep(0,30), 
                y_position = rep(7.5,30),
                map_signif_level = T,parse = F)

lowcom<-list(c("a1", "a2"),
             c("b1", "b2"),
             c("c1", "c2"),
             c("d1", "d2"),c("e1", "e2"),c("f1", "f2"),c("g1", "g2"),c("h1", "h2"),c("i1", "i2"),c("j1", "j2"),
             c("k1", "k2"),c("l1", "l2"),c("m1", "m2"),c("n1", "n2"),c("o1","o2"),c("p1", "p2"))
LowmutgenePlot<-ggplot(Lowdata)+
    geom_boxplot(aes(group, value, fill = fenzu, color = fenzu), 
                 width = 0.5,
                 outlier.shape = NA,
                 size=0.3,
    )+
    ylab("gene expression (FPKM)")+
    scale_x_discrete(labels = unique(Lowdata$gene))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values = mypal3[c(3,1)])+
    scale_color_manual(values = c("black","black"))+
    theme_bw()+
    ylim(0,8.5)+
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          panel.grid.minor.x  =element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(angle=45,hjust=0,vjust=0.2),
          axis.line = element_line(color = "black"))+
    geom_signif(aes(group, value),
                comparisons = lowcom,
                tip_length = rep(0,32), 
                y_position = rep(7.5,32),
                map_signif_level = T,parse = F)
gongcom<-list(c("a1", "a2"),
              c("b1", "b2"),
              c("c1", "c2"),
              c("d1", "d2"),c("e1", "e2"),c("f1", "f2"),c("g1", "g2"),c("h1", "h2"),c("i1","i2"))
GongmutgenePlot<-ggplot(gongtudata)+
    geom_boxplot(aes(group, value, fill = fenzu, color = fenzu), 
                 width = 0.5,
                 outlier.shape = NA,
                 size=0.3,
    )+
    ylab("gene expression (FPKM)")+
    scale_x_discrete(labels = unique(gongtudata$gene))+
    scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values = mypal3[c(3,1)])+
    scale_color_manual(values =c("black","black") )+
    theme_bw()+
    ylim(0,7.5)+
    theme(panel.grid = element_blank(), 
          panel.border = element_blank(),
          panel.grid.minor.x  =element_blank(),
          legend.position = "none",
          legend.title = element_blank(),
          #plot.title = element_text(hjust = 0.5),
          axis.title = element_blank(),
          axis.text.x = element_text(angle=45,hjust=0,vjust=0.2),
          axis.line = element_line(color = "black"))+
    geom_signif(aes(group, value),
                comparisons = gongcom,
                tip_length = rep(0,18), 
                y_position = rep(7,18),
                map_signif_level = T,parse = F)


install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 

######The figure shows the CpG site numbers in different immune cell types
mypal3=c("#003399","#CC3333","#CC6633")
cgSiteBar<-ggplot(barr1, aes(x = cell,y =Num,fill=SiteType))+
    geom_bar(stat ="identity",width = 0.65,position = "dodge")+     
    scale_fill_manual(values = rep(c(mypal3[1],mypal3[2]),6))+                
    labs(x = "",y = "Num(x1000)", title = "")+  
    geom_hline(yintercept = c(25,50,75),linetype="dashed",color="gray40")+
    ylim(0,130)+
    scale_y_continuous(expand = expansion(mult=c(0,0)))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.position="top",
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"),
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
cgSiteBar