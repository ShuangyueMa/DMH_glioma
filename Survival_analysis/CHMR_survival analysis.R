library(survival)
library(survminer)
library(forestplot)
mypal3=c("#003399","#CC0000","#CC6633")

######The KM survival curves for high- and low-risk patients in TCGA, GSE103659, and GSE60274
#KM analysis of CHMR score for TCGA
modelTCGA<-as.numeric(summary(coxph(formula = Surv(OS.time, OS) ~ cg04848682+cg23757461+cg04359840+cg16269755+cg03180426+cg25578728+cg24035035+cg24770798, data = SurvivalTCGA, 
                                    method = "efron"))$coef[,1])
SurvivalTCGA$CMHR<-model1[1]*SurvivalTCGA[,grep("cg04848682",colnames(SurvivalTCGA))]+
    model1[2]*SurvivalTCGA[,grep("cg23757461",colnames(SurvivalTCGA))]+
    model1[3]*SurvivalTCGA[,grep("cg04359840",colnames(SurvivalTCGA))]+
    model1[4]*SurvivalTCGA[,grep("cg16269755",colnames(SurvivalTCGA))]+
    model1[5]*SurvivalTCGA[,grep("cg03180426",colnames(SurvivalTCGA))]+
    model1[6]*SurvivalTCGA[,grep("cg25578728",colnames(SurvivalTCGA))]+
    model1[7]*SurvivalTCGA[,grep("cg24035035",colnames(SurvivalTCGA))]+
    model1[8]*SurvivalTCGA[,grep("cg24770798",colnames(SurvivalTCGA))]
SurvivalTCGA$CMHR[order(SurvivalTCGA$CMHR)][272]
SurvivalTCGA$group<-ifelse(SurvivalTCGA$CMHR<=SurvivalTCGA$CMHR[order(SurvivalTCGA$CMHR)][272],"Low","High")
fit_riskTCGA <- survfit(Surv(OS.time,OS) ~group ,  
                        data = SurvivalTCGA)
surv <- Surv(time = SurvivalTCGA$OS.time, event = SurvivalTCGA$OS);data.survdiff=survdiff(surv~group,data=SurvivalTCGA)
prisk_TCGA = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
prisk_TCGA
KM_TCGArisk<-ggsurvplot(fit_riskTCGA, 
                        data =SurvivalTCGA, 
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
                        legend.title="CMHR",
                        legend.labs=c("Low","High"),
                        fontsize=4,
                        log.rank.weights="1",
                        risk.table.fontsize=5,
                        font.main = c(16, "bold", "black"),
                        break.x.by = 1000) 
KM_TCGArisk

#KM analysis of CHMR score for GSE103659
modelGSE103659<-as.numeric(summary(coxph(formula = Surv(OS.time, OS) ~ cg04848682+cg23757461+cg04359840+cg16269755+cg03180426+cg25578728+cg24035035+cg24770798, data = SurvivalGSE103659, 
                                         method = "efron"))$coef[,1])
SurvivalGSE103659$CMHR<-model1[1]*SurvivalGSE103659[,grep("cg04848682",colnames(SurvivalGSE103659))]+
    model1[2]*SurvivalGSE103659[,grep("cg23757461",colnames(SurvivalGSE103659))]+
    model1[3]*SurvivalGSE103659[,grep("cg04359840",colnames(SurvivalGSE103659))]+
    model1[4]*SurvivalGSE103659[,grep("cg16269755",colnames(SurvivalGSE103659))]+
    model1[5]*SurvivalGSE103659[,grep("cg03180426",colnames(SurvivalGSE103659))]+
    model1[6]*SurvivalGSE103659[,grep("cg25578728",colnames(SurvivalGSE103659))]+
    model1[7]*SurvivalGSE103659[,grep("cg24035035",colnames(SurvivalGSE103659))]+
    model1[8]*SurvivalGSE103659[,grep("cg24770798",colnames(SurvivalGSE103659))]
SurvivalGSE103659$CMHR[order(SurvivalGSE103659$CMHR)][90]
SurvivalGSE103659$group<-ifelse(SurvivalGSE103659$CMHR<=SurvivalGSE103659$CMHR[order(SurvivalGSE103659$CMHR)][90],"Low","High")
fit_riskGSE103659 <- survfit(Surv(OS.time,OS) ~group ,  
                             data = SurvivalGSE103659)
surv <- Surv(time = SurvivalGSE103659$OS.time, event = SurvivalGSE103659$OS);data.survdiff=survdiff(surv~group,data=SurvivalGSE103659)
prisk_GSE103659 = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
prisk_GSE103659
KM_GSE103659risk<-ggsurvplot(fit_riskGSE103659, 
                             data =SurvivalGSE103659, 
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
                             legend.title="CMHR",
                             legend.labs=c("Low","High"),
                             fontsize=4,
                             log.rank.weights="1",
                             risk.table.fontsize=5,
                             font.main = c(16, "bold", "black"),
                             break.x.by = 1000) 
KM_GSE103659risk

#KM analysis of CHMR score for GSE60274
modelGSE60274<-as.numeric(summary(coxph(formula = Surv(OS.time, OS) ~ cg04848682+cg23757461+cg04359840+cg16269755+cg03180426+cg25578728+cg24035035+cg24770798, data = SurvivalGSE60274, 
                                        method = "efron"))$coef[,1])
SurvivalGSE60274$CMHR<-model1[1]*SurvivalGSE60274[,grep("cg04848682",colnames(SurvivalGSE60274))]+
    model1[2]*SurvivalGSE60274[,grep("cg23757461",colnames(SurvivalGSE60274))]+
    model1[3]*SurvivalGSE60274[,grep("cg04359840",colnames(SurvivalGSE60274))]+
    model1[4]*SurvivalGSE60274[,grep("cg16269755",colnames(SurvivalGSE60274))]+
    model1[5]*SurvivalGSE60274[,grep("cg03180426",colnames(SurvivalGSE60274))]+
    model1[6]*SurvivalGSE60274[,grep("cg25578728",colnames(SurvivalGSE60274))]+
    model1[7]*SurvivalGSE60274[,grep("cg24035035",colnames(SurvivalGSE60274))]+
    model1[8]*SurvivalGSE60274[,grep("cg24770798",colnames(SurvivalGSE60274))]
SurvivalGSE60274$CMHR[order(SurvivalGSE60274$CMHR)][31]
SurvivalGSE60274$group<-ifelse(SurvivalGSE60274$CMHR<=SurvivalGSE60274$CMHR[order(SurvivalGSE60274$CMHR)][31],"Low","High")
fit_riskGSE60274 <- survfit(Surv(OS.time,OS) ~group ,  
                            data = SurvivalGSE60274)
surv <- Surv(time = SurvivalGSE60274$OS.time, event = SurvivalGSE60274$OS);data.survdiff=survdiff(surv~group,data=SurvivalGSE60274)
prisk_GSE60274 = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
prisk_GSE60274
KM_GSE60274risk<-ggsurvplot(fit_riskGSE60274, 
                            data =SurvivalGSE60274, 
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
                            legend.title="CMHR",
                            legend.labs=c("Low","High"),
                            fontsize=4,
                            log.rank.weights="1",
                            risk.table.fontsize=5,
                            font.main = c(16, "bold", "black"),
                            break.x.by = 1000) 
KM_GSE60274risk


######The forestplot for CHMR score in three datasets: TCGA GSE103659 GSE60274
#TCGA cox model created by CMHR score, age, tumor grade, IDH status, and MGMT promoter status
coxCHMR_TCGA<- coxph(Surv(OS.time, OS)~CMHR+Age+TumorGrade+IDH_status+MGMT_status,data=TCGACoxData)
forest_TCGA<-ggforest(coxCHMR_TCGA.cox,main="hazard ratio",cpositions=c(0.02,0.22,0.4),fontsize=1.4,refLabel="reference",noDigits=2)
Forest4_TCGA<-forestplot(TCGA_Forest[,c(1,2,6,7)], 
                         mean=TCGA_Forest[,3],   
                         lower=TCGA_Forest[,4],  
                         upper=TCGA_Forest[,5], 
                         zero=1,            
                         boxsize=0.5,      
                         graph.pos= "right" ,
                         hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                         "2" = gpar(lty=2),
                                         "7"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                         graphwidth = unit(.25,"npc"),
                         xlab="AIC: 1023.48 C-index:0.65 logRank:2.80e-06",
                         xticks=c(0.4,1,3,5,7,10) ,
                         is.summary=c(T,F,F,F,F,F),
                         txt_gp=fpTxtGp(
                             label = list(
                                 title=gpar(fontfamily = "Arail", cex=1.2),
                                 gpar(fontfamily = "Times", fontface="italic"),
                                 gpar(fontfamily = "",
                                      col = "black")),
                             ticks = gpar(fontfamily = "Arail", cex=1),
                             xlab  = gpar(fontfamily = "Arail", fontface="bold",cex=1.0)),
                         lwd.zero=1,
                         lwd.xaxis=1, 
                         lty.ci=1.3,
                         ci.vertices =T,
                         ci.vertices.height=0.2, 
                         clip=c(0.1,8),
                         ineheight=unit(8, 'mm'), 
                         line.margin=unit(8, 'mm'),
                         colgap=unit(6, 'mm'),
                         fn.ci_norm="fpDrawDiamondCI", 
                         col=fpColors(box =mypal3[1], 
                                      lines =mypal3[1], 
                                      zero = "black"),
                         
                         
)
Forest4_TCGA


#GSE103659 cox model created by CMHR score, age, and MGMT promoter status 
coxCHMR_GSE10<- coxph(Surv(OS.time, OS)~CHMR+Age+MGMT_status,data=GSE103659CoxData)
forest_GSE10<-ggforest(coxCHMR_GSE10.cox,main="hazard ratio",cpositions=c(0.02,0.22,0.4),fontsize=1.4,refLabel="reference",noDigits=2)
Forest_GSE10<-forestplot(GSE103659_Forest[,c(1,2,6,7)], 
                            mean=GSE103659_Forest[,3],   
                            lower=GSE103659_Forest[,4],  
                            upper=GSE103659_Forest[,5], 
                            zero=1,            
                            boxsize=0.5,      
                            graph.pos= "right" ,
                            hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                            "2" = gpar(lty=2),
                                            "7"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                            graphwidth = unit(.25,"npc"),
                            xlab="AIC: 1023.48 C-index:0.65 logRank:2.80e-06",
                            xticks=c(0.4,1,3,5,7,10) ,
                            is.summary=c(T,F,F,F,F,F),
                            txt_gp=fpTxtGp(
                                label = list(
                                    title=gpar(fontfamily = "Arail", cex=1.2),
                                    gpar(fontfamily = "Times", fontface="italic"),
                                    gpar(fontfamily = "",
                                         col = "black")),
                                ticks = gpar(fontfamily = "Arail", cex=1),
                                xlab  = gpar(fontfamily = "Arail", fontface="bold",cex=1.0)),
                            lwd.zero=1,
                            lwd.xaxis=1, 
                            lty.ci=1.3,
                            ci.vertices =T,
                            ci.vertices.height=0.2, 
                            clip=c(0.1,8),
                            ineheight=unit(8, 'mm'), 
                            line.margin=unit(8, 'mm'),
                            colgap=unit(6, 'mm'),
                            fn.ci_norm="fpDrawDiamondCI", 
                            col=fpColors(box =mypal3[1], 
                                         lines =mypal3[1], 
                                         zero = "black"),
                            
                            
)
Forest4_GSE10

#GSE60274 cox model created by CMHR score, age, and MGMT promoter status 
coxCHMR_GSE60<- coxph(Surv(OS.time, OS)~CMHR+Age+TumorGrade+IDH_status+MGMT_status,data=GSE60CoxData)
forest_GSE60<-ggforest(coxCHMR_GSE60.cox,main="hazard ratio",cpositions=c(0.02,0.22,0.4),fontsize=1.4,refLabel="reference",noDigits=2)
Forest4_GSE60<-forestplot(GSE60_Forest[,c(1,2,6,7)], 
                          mean=GSE60_Forest[,3],   
                          lower=GSE60_Forest[,4],  
                          upper=GSE60_Forest[,5], 
                          zero=1,            
                          boxsize=0.5,      
                          graph.pos= "right" ,
                          hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                          "2" = gpar(lty=2),
                                          "7"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                          graphwidth = unit(.25,"npc"),
                          xlab="AIC: 1023.48 C-index:0.65 logRank:2.80e-06",
                          xticks=c(0.4,1,3,5,7,10) ,
                          is.summary=c(T,F,F,F,F,F),
                          txt_gp=fpTxtGp(
                              label = list(
                                  title=gpar(fontfamily = "Arail", cex=1.2),
                                  gpar(fontfamily = "Times", fontface="italic"),
                                  gpar(fontfamily = "",
                                       col = "black")),
                              ticks = gpar(fontfamily = "Arail", cex=1),
                              xlab  = gpar(fontfamily = "Arail", fontface="bold",cex=1.0)),
                          lwd.zero=1,
                          lwd.xaxis=1, 
                          lty.ci=1.3,
                          ci.vertices =T,
                          ci.vertices.height=0.2, 
                          clip=c(0.1,8),
                          ineheight=unit(8, 'mm'), 
                          line.margin=unit(8, 'mm'),
                          colgap=unit(6, 'mm'),
                          fn.ci_norm="fpDrawDiamondCI", 
                          col=fpColors(box =mypal3[1], 
                                       lines =mypal3[1], 
                                       zero = "black"),
                          
                          
)
Forest4_GSE60



