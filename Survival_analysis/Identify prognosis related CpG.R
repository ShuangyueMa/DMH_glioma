library(survival)
library(survminer)
library(forestplot)
mypal3=c("#003399","#CC0000","#CC6633")

######The KM survival curves for prognosis-realted CpG sites (β-value) 
TCGA_cgbetaKM<-list()
for(i in prognosis_Cpg){
cgValue<-betaTCGA[,match(cg,colnames(betaTCGA))]
betaTCGA<-beta_TCGA[order(cgValue),]
betaTCGA$group<-c(rep("Low-β",272),rep("high-β",272))
fit_TCGA1 <- survfit(Surv(OS.time,OS) ~group , 
                     data = betaTCGA)
TCGAsur_beta<-ggsurvplot(fit_TCGA1, 
                         data =betaTCGA,  
                         conf.int = F, 
                         surv.median.line = "hv",  
                         pval.coord=c(4500,0.60),
                         ggtheme = theme_bw(base_family = "Arail")+theme(
                             legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                             legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                             legend.position="none", 
                             axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                             axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                            axis.title.x = element_blank(),
                             axis.title.y = element_blank(),
                             plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                             panel.border =element_blank(),                      
                             axis.line = element_line(color = "black"),
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank()),
                         pval.method.size=TRUE,
                         censor=T,
                         risk.table = F, 
                         xlab = "Overall survival (days)", 
                         legend = c(0.65,0.83), 
                         palette = c(mypal3[1],mypal3[3]),
                         legend.title="β",
                         legend.labs=c("Low","High"),
                         fontsize=4,
                         log.rank.weights="1",
                         risk.table.fontsize=5,
                         font.main = c(16, "bold", "black"),
                         break.x.by = 1000,) 
TCGA_cgBetaKM[[i]]<-"TCGAsur_beta"
print(i)
}

GSE103659_cgbetaKM<-list()
for(i in prognosis_Cpg){
    cgValue<-betaGSE103659[,match(cg,colnames(betaGSE103659))]
    betaGSE103659<-beta_GSE103659[order(cgValue),]
    betaGSE103659$group<-c(rep("Low-β",90),rep("high-β",90))
    fit_GSE1036591 <- survfit(Surv(OS.time,OS) ~group , 
                              data = betaGSE103659)
    GSE103659sur_beta<-ggsurvplot(fit_GSE1036591, 
                                  data =betaGSE103659,  
                                  conf.int = F, 
                                  surv.median.line = "hv",  
                                  pval.coord=c(4500,0.60),
                                  ggtheme = theme_bw(base_family = "Arail")+theme(
                                      legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                                      legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                                      legend.position="none", 
                                      axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                      axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                                      panel.border =element_blank(),                      
                                      axis.line = element_line(color = "black"),
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank()),
                                  pval.method.size=TRUE,
                                  censor=T,
                                  risk.table = F, 
                                  xlab = "Overall survival (days)", 
                                  legend = c(0.65,0.83), 
                                  palette = c(mypal3[1],mypal3[3]),
                                  legend.title="β",
                                  legend.labs=c("Low","High"),
                                  fontsize=4,
                                  log.rank.weights="1",
                                  risk.table.fontsize=5,
                                  font.main = c(16, "bold", "black"),
                                  break.x.by = 1000,) 
    GSE103659_cgBetaKM[[i]]<-"GSE103659sur_beta"
    print(i)
}

GSE60274_cgbetaKM<-list()
for(i in prognosis_Cpg){
    cgValue<-betaGSE60274[,match(cg,colnames(betaGSE60274))]
    betaGSE60274<-beta_GSE60274[order(cgValue),]
    betaGSE60274$group<-c(rep("Low-β",31),rep("high-β",31))
    fit_GSE602741 <- survfit(Surv(OS.time,OS) ~group , 
                             data = betaGSE60274)
    GSE60274sur_beta<-ggsurvplot(fit_GSE602741, 
                                 data =betaGSE60274,  
                                 conf.int = F, 
                                 surv.median.line = "hv",  
                                 pval.coord=c(4500,0.60),
                                 ggtheme = theme_bw(base_family = "Arail")+theme(
                                     legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                                     legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                                     legend.position="none", 
                                     axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                     axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                                     panel.border =element_blank(),                      
                                     axis.line = element_line(color = "black"),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank()),
                                 pval.method.size=TRUE,
                                 censor=T,
                                 risk.table = F, 
                                 xlab = "Overall survival (days)", 
                                 legend = c(0.65,0.83), 
                                 palette = c(mypal3[1],mypal3[3]),
                                 legend.title="β",
                                 legend.labs=c("Low","High"),
                                 fontsize=4,
                                 log.rank.weights="1",
                                 risk.table.fontsize=5,
                                 font.main = c(16, "bold", "black"),
                                 break.x.by = 1000,) 
    GSE60274_cgBetaKM[[i]]<-"GSE60274sur_beta"
    print(i)
}


######The KM survival curves for prognosis-realted CpGct (CHMR score)
TCGA_cgCHMRKM<-list()
for(i in prognosis_Cpg){
    cgValue<-CHMRTCGA[,match(cg,colnames(CHMRTCGA))]
    CHMRTCGA<-CHMR_TCGA[order(cgValue),]
    CHMRTCGA$group<-c(rep("Low-CHMR",272),rep("high-CHMR",272))
    fit_TCGA1 <- survfit(Surv(OS.time,OS) ~group , 
                         data = CHMRTCGA)
    TCGAsur_CHMR<-ggsurvplot(fit_TCGA1, 
                             data =CHMRTCGA,  
                             conf.int = F, 
                             surv.median.line = "hv",  
                             pval.coord=c(4500,0.60),
                             ggtheme = theme_bw(base_family = "Arail")+theme(
                                 legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                                 legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                                 legend.position="none", 
                                 axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                 axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                 axis.title.x = element_blank(),
                                 axis.title.y = element_blank(),
                                 plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                                 panel.border =element_blank(),                      
                                 axis.line = element_line(color = "black"),
                                 panel.grid.major = element_blank(), 
                                 panel.grid.minor = element_blank()),
                             pval.method.size=TRUE,
                             censor=T,
                             risk.table = F, 
                             xlab = "Overall survival (days)", 
                             legend = c(0.65,0.83), 
                             palette = c(mypal3[1],mypal3[3]),
                             legend.title="CHMR",
                             legend.labs=c("Low","High"),
                             fontsize=4,
                             log.rank.weights="1",
                             risk.table.fontsize=5,
                             font.main = c(16, "bold", "black"),
                             break.x.by = 1000,) 
    TCGA_cgCHMRKM[[i]]<-"TCGAsur_CHMR"
    print(i)
}

GSE103659_cgCHMRKM<-list()
for(i in prognosis_Cpg){
    cgValue<-CHMRGSE103659[,match(cg,colnames(CHMRGSE103659))]
    CHMRGSE103659<-CHMR_GSE103659[order(cgValue),]
    CHMRGSE103659$group<-c(rep("Low-CHMR",90),rep("high-CHMR",90))
    fit_GSE1036591 <- survfit(Surv(OS.time,OS) ~group , 
                              data = CHMRGSE103659)
    GSE103659sur_CHMR<-ggsurvplot(fit_GSE1036591, 
                                  data =CHMRGSE103659,  
                                  conf.int = F, 
                                  surv.median.line = "hv",  
                                  pval.coord=c(4500,0.60),
                                  ggtheme = theme_bw(base_family = "Arail")+theme(
                                      legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                                      legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                                      legend.position="none", 
                                      axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                      axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank(),
                                      plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                                      panel.border =element_blank(),                      
                                      axis.line = element_line(color = "black"),
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank()),
                                  pval.method.size=TRUE,
                                  censor=T,
                                  risk.table = F, 
                                  xlab = "Overall survival (days)", 
                                  legend = c(0.65,0.83), 
                                  palette = c(mypal3[1],mypal3[3]),
                                  legend.title="CHMR",
                                  legend.labs=c("Low","High"),
                                  fontsize=4,
                                  log.rank.weights="1",
                                  risk.table.fontsize=5,
                                  font.main = c(16, "bold", "black"),
                                  break.x.by = 1000,) 
    GSE103659_cgCHMRKM[[i]]<-"GSE103659sur_CHMR"
    print(i)
}

GSE60274_cgCHMRKM<-list()
for(i in prognosis_Cpg){
    cgValue<-CHMRGSE60274[,match(cg,colnames(CHMRGSE60274))]
    CHMRGSE60274<-CHMR_GSE60274[order(cgValue),]
    CHMRGSE60274$group<-c(rep("Low-CHMR",31),rep("high-CHMR",31))
    fit_GSE602741 <- survfit(Surv(OS.time,OS) ~group , 
                             data = CHMRGSE60274)
    GSE60274sur_CHMR<-ggsurvplot(fit_GSE602741, 
                                 data =CHMRGSE60274,  
                                 conf.int = F, 
                                 surv.median.line = "hv",  
                                 pval.coord=c(4500,0.60),
                                 ggtheme = theme_bw(base_family = "Arail")+theme(
                                     legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                                     legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                                     legend.position="none", 
                                     axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                     axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                                     panel.border =element_blank(),                      
                                     axis.line = element_line(color = "black"),
                                     panel.grid.major = element_blank(), 
                                     panel.grid.minor = element_blank()),
                                 pval.method.size=TRUE,
                                 censor=T,
                                 risk.table = F, 
                                 xlab = "Overall survival (days)", 
                                 legend = c(0.65,0.83), 
                                 palette = c(mypal3[1],mypal3[3]),
                                 legend.title="CHMR",
                                 legend.labs=c("Low","High"),
                                 fontsize=4,
                                 log.rank.weights="1",
                                 risk.table.fontsize=5,
                                 font.main = c(16, "bold", "black"),
                                 break.x.by = 1000,) 
    GSE60274_cgCHMRKM[[i]]<-"GSE60274sur_CHMR"
    print(i)
}



