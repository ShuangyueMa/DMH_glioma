library(survival)
library(survminer)
library(forestplot)
mypal3=c("#003399","#CC0000","#CC6633")

######KM survival analysis for three CpG sites realted to radiation therapy in GSE195684 dataset
cg=c("cg03180426","cg23757461","cg24770798")
RZ_KMPlot<-list()
for(i in 1:length(cg))
RZ_sample<-filter(phenoRad,treatment=="RT"){
RZ_sample$value<-as.numeric(BetaRad[match(cg[i],rownames(BetaRad)),match(RZ_sample$sample,colnames(BetaRad))])
RZ_sample1<-RZ_sample[order(RZ_sample$value),]
RZ_sample1$group<-c(rep("Low",36),rep("High",37))
fit_GSE195684_RT <- survfit(Surv(OS.time,OS) ~group ,  
                        data = RZ_sample1)
surv <- Surv(time = RZ_sample1$OS.time, event = RZ_sample1$OS);data.survdiff=survdiff(surv~group,data=RZ_sample1)
p1rt = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
p1rt
p1_rtPlot<-ggsurvplot(fit_GSE195684_RT,
                      data =RZ_sample1,  
                      conf.int = F, 
                      surv.median.line = "hv",  
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
                      legend.title=cg[i],
                      legend.labs=c("Low","High"),
                      fontsize=4,
                      log.rank.weights="1",
                      risk.table.fontsize=5,
                      font.main = c(16, "bold", "black"),
                      break.x.by = 200) 
p1_rtPlot
RZ_KMPlot[[i]]<-p1_rtPlot
print(i)
}


######The forestplot for survival analysis of eight CpGct in GSE195684 dataset
cox_19data<-as.data.frame(cbind(phenoRad,guodu1ge))
cgcox_result<- summary(coxph(Surv(OS.time, OS)~age+gender+cg04848682+cg04359840+
                        cg03180426+cg23757461+cg16269755+cg25578728+cg24770798+cg24035035, data=cox_19data))
CpG_Forest_rad<-forestplot(cgcox_result[,c(1,2,6,7)], 
                          mean=cgcox_result[,3],   
                          lower=cgcox_result[,4],  
                          upper=cgcox_result[,5], 
                          zero=1,            
                          boxsize=0.5,      
                          graph.pos= "right" ,
                          hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                          "2" = gpar(lty=2),
                                          "14"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                          graphwidth = unit(.25,"npc"),
                          xlab="AIC: 849.88 C-index: 0.63 log-rank p: 0.002",
                          xticks=c(0.01,1,2,4,6,8,10,12) ,
                          is.summary=c(T,F,F,F,F,F,F,F,F,F,F,F,F),
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
CpG_Forest_rad
