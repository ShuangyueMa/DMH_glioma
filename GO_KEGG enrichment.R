install.packages("colourpicker")
install.packages("pacman") 
library(pacman)
pacman::p_load(colourpicker,ggplot2,dplyr,viridis,ggpointdensity,ggthemes,ggsci,ggsignif) 
mypal3=c("#FF8C00", "#EE2C2C", "#0000CD")

######The figure shows Go terms enriched with genes related to cell specific CpG sites
data<-godata2%>%select(ONTOLOGY,GOTerm,p.adjust)%>%arrange(ONTOLOGY,-log10(p.adjust))
data$GOTerm<-factor(data$GOTerm,levels = unique(data$GOTerm),ordered = T)
#go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))
GO2<-ggplot(data=Goterm)+
    geom_bar(aes(x=GOTerm,y=-log10(p.adjust),fill=ONTOLOGY),stat = 'identity',position=position_dodge(0.63))+
    scale_fill_manual(values = mypal3)+
    coord_flip()+
    ggtitle("Go Enrichment")+
    xlab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    geom_hline(yintercept =1.30103,linetype="dashed",color="black")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    theme(
        legend.title =element_blank(), 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
GO2<-GO2+annotate("text",label="p.adjust<0.05",family="Arail",
                  parse=F,x=30,y=5,color="black",size=5.5,fontface="italic")

GO2

kegg1<-ggplot(data = kegg,mapping = aes(x=-log10(p.adjust),y=Description))
kegg2<-kegg1+geom_point(mapping = aes(size=Count,colour=-1*log10(p.adjust)))+scale_color_gradient(low ="#0000CD",high = "#EE2C2C")+
    labs(title = "KEGG Pathway Enrichment",x="-log10(p.adjust)",y=" ",colour="-log10(p.adjust)",size="gene number")+
    theme_bw()+
    theme(
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"),
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),        
        panel.grid.minor = element_blank()
    )
kegg2


######The figure shows sample numbers in three datasets
sample_bar<- ggplot(dataSample, aes(x = data,y =number,fill=dat))+
    geom_bar(stat ="identity",width = 0.55,position = "dodge")+     
    scale_fill_manual(values =mypal3[c(3,1)])+                 
    labs(x = "",y = "patients", title = "")+  
    scale_y_continuous(expand = expansion(mult=c(0,0)))+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    theme(
        plot.margin = unit(c(0.1,0.1,0.1,0.1),"cm"),
        legend.key.size = unit(0.8,"cm"),
        legend.title = element_blank(),
        legend.position="top", 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
sample_bar

######The figure shows EWAS Go terms of cell-type-specific CpG sites
GOcmhcPlot<-ggplot(data=GOcmhc)+
    geom_bar(aes(x=heng,y=-log10(p.Value),fill=Ontology.Type),stat = 'identity',position=position_dodge(0.63))+
    scale_fill_manual(values = mypal3)+
    coord_flip()+
    xlab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    geom_hline(yintercept =1.30103,linetype="dashed",color="black")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    scale_x_discrete(labels=rev(GOcmhc$Term[c(14:1,18:15,24:19)]))+
    theme(
        legend.position="top", 
        legend.title =element_blank(), 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
GOcmhcPlot

######The figure shows EWAS Go terms of CpG sites related to prognosis
GObetaPlot<-ggplot(data=GObeta)+
    geom_bar(aes(x=heng,y=-log10(p.Value),fill=Ontology.Type),stat = 'identity',position=position_dodge(0.63))+
    scale_fill_manual(values = mypal3)+
    coord_flip()+
    xlab("")+
    theme_classic()+theme(legend.position = 0)+
    theme_bw()+
    geom_hline(yintercept =1.30103,linetype="dashed",color="black")+
    scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
    scale_x_discrete(labels=rev(GObeta$Term[c(14:1,18:15,24:19)]))+
    theme(
        legend.position="top", 
        legend.title =element_blank(), 
        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
    )
GObetaPlot
