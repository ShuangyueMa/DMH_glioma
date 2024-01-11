library(survival)
library(survminer)
library(forestplot)
mypal3=c("#003399","#CC0000","#CC6633")

######The KM survival curves for NMF subtypes
fit <- survfit(Surv(OS.time,OS) ~NMF_cluster ,  
               data = NMF_survivalData)
nmf_OS<-ggsurvplot(fit,
                   data =NMF_survivalData,  
                   conf.int = F, 
                   pval = TRUE, 
                   surv.median.line = "hv",  
                   pval.coord=c(2000,0.75),
                   ggtheme = theme_bw(base_family = "Arail")+theme(
                       legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                       legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                       legend.position="none", 
                       axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                       axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                       axis.title.y=element_text(family="Arail",size=16,face="bold"), 
                       axis.title.x=element_text(family="Arail",size=16,face="bold"), 
                       plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                       panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),                       
                       panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank()),
                   pval.method.size=TRUE,
                   censor=T,
                   risk.table = F, 
                   xlab = "Overall Survival(days)", 
                   legend = c(0.8,0.75), 
                   palette = mypal3,
                   legend.labs=c("NMF_Cluster1","NMF_Cluster2","NMF_Cluster3"),
                   fontsize=4,
                   log.rank.weights="1",
                   risk.table.fontsize=5,
                   legend.title = "", 
                   break.x.by = 1000,)  

fit1 <- survfit(Surv(PFI.time,PFI) ~NMF_cluster ,  
                data = NMF_survivalData)

nmf_pfi<-ggsurvplot(fit1,
                    data =NMF_survivalData,  
                    conf.int = F, 
                    pval = TRUE, 
                    surv.median.line = "hv",  
                    pval.coord=c(2000,0.75),
                    ggtheme = theme_bw(base_family = "Arail")+theme(
                        legend.title =  element_text(colour="black",family="Arail",size=10,face="bold"),
                        legend.text = element_text(colour="black",family="Arail",size=10,face="bold"),
                        legend.position="none", 
                        axis.text.x=element_text(colour="black",family="Arail",size=14,face="bold"), 
                        axis.text.y=element_text(colour="black",family="Arail",size=14,face="bold"), 
                        axis.title.y=element_text(family="Arail",size=16,face="bold"), 
                        axis.title.x=element_text(family="Arail",size=16,face="bold"), 
                        plot.title = element_text(family="Arail",size=16,face="bold",hjust = 0.5,lineheight=0.2), 
                        panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"),                       #axis.line = element_line(color = "black", size = 1.5),
                        panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank()),
                    pval.method.size=TRUE,
                    censor=T,
                    risk.table = F, 
                    xlab = "Progrssion Free Interval (days)", 
                    legend = c(0.8,0.75), 
                    palette = mypal3,
                    legend.labs=c("NMF_Cluster1","NMF_Cluster2","NMF_Cluster3"),
                    fontsize=4,
                    log.rank.weights="1",
                    risk.table.fontsize=5,
                    legend.title = "", 
                    break.x.by = 1000,)  
nmf_pfi

######The forestplot for survival analysis of PIM score
forest2<-ggforest(PIM.cox,main="hazard ratio",cpositions=c(0.02,0.22,0.4),fontsize=1.4,refLabel="reference",noDigits=2)
Forest_PIM<-forestplot(Forest4Result[,c(1,2,6,7)], 
                       mean=Forest4Result[,3],   
                       lower=Forest4Result[,4],  
                       upper=Forest4Result[,5], 
                       zero=1,            
                       boxsize=0.5,      
                       graph.pos= "right" ,
                       hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                       "2" = gpar(lty=2),
                                       "17"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                       graphwidth = unit(.25,"npc"),
                       xlab="AIC: 938.99 C-index:0.84 logRank:8.10e-21",
                       xticks=c(0.4,1,3,5,7,10) ,
                       is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F),
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

