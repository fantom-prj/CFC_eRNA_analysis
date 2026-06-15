library(dplyr) 
library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(ggalt)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)
library(tidyr)


###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig2=paste0(primary_folder,"fig2/out/")
path_fig2_data=paste0(primary_folder,"fig2/data/")

#=======
theme1=theme(  panel.background = element_rect(fill = "white", color = "white", linewidth = 0.25, linetype = "solid"),
               panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', color = "grey90"), 
               axis.text.x = element_text(color ="black"),
               axis.text.y = element_text(color ="black"),
               plot.title = element_text(hjust = 0.5,  color ="black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
               plot.subtitle = element_text(color="black", hjust=0.5),
               text = element_text(size=6),
               strip.text = element_text(size=6),
               legend.background = element_rect(fill="white", color="black", linewidth=0.25),
               legend.title = element_text(hjust=0.5),
               legend.text = element_text(lineheight = 0.6, margin = margin(l = 1, unit = "pt")),
               legend.key.size = unit(0.2, 'cm'),
               legend.margin = margin(0.02,0.02,0.02,0.02, "cm"),
               axis.line = element_line(linewidth = 0.25, colour = "black"),
               axis.ticks = element_line(linewidth = 0.25,colour = "black"),
               strip.background =element_rect(fill="grey90"),
               strip.text.x = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               strip.text.y = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               legend.key = element_blank(),
               legend.position = "right",
               legend.box.margin = margin(t=-10, r=0, b=0, l=-10),
               plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

#===============================================================================
#Needed
#f2b
RLE.CRE=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.RLE.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

RLE.CRE=log10(RLE.CRE+0.01)
RLE.CRE=data.frame(RLE.CRE)
RLE.CRE[RLE.CRE>3]=3
RLE.CREp=RLE.CRE[which(rownames(RLE.CRE)%in% CREanno$CREID[which(CREanno$promoter_type=="promoter-like")]),]
RLE.CREe=RLE.CRE[which(rownames(RLE.CRE)%in% CREanno$CREID[which(CREanno$promoter_type=="enhancer-like")]),]
RLE.CREu=RLE.CRE[which(rownames(RLE.CRE)%in% CREanno$CREID[which(CREanno$promoter_type=="unclassed")]),]
# sample 10% for visualization
RLE.CREp1=RLE.CREp[sample(nrow(RLE.CREp),nrow(RLE.CREp)/10),]
RLE.CREe1=RLE.CREe[sample(nrow(RLE.CREe),nrow(RLE.CREe)/10),]
RLE.CREu1=RLE.CREu[sample(nrow(RLE.CREu),nrow(RLE.CREu)/10),]

library(ComplexHeatmap)
library(circlize)
library("RColorBrewer")

col_fun <- colorRamp2(
  seq(min(RLE.CREp1, na.rm = TRUE), max(RLE.CREp, na.rm = TRUE), length.out = 25),
  colorRampPalette(c("white", "navy"))(25))

out1 <- Heatmap(
  RLE.CREp1,
  name = "Legand",
  col = col_fun,
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
  row_dend_width = unit(8, "mm"),
  height = unit(1.53, "cm"),
  use_raster = TRUE, row_dend_gp = gpar(lwd = 0.25))
  
col_fun2 <- colorRamp2(
  seq(min(RLE.CREe1, na.rm = TRUE), max(RLE.CREp, na.rm = TRUE), length.out = 25),
  colorRampPalette(c("white", "darkred"))(25))

out2 <- Heatmap(
  RLE.CREe1,
  name = "Legand",
  col = col_fun2,
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
  row_dend_width = unit(8, "mm"),
  height = unit(1.62, "cm"),
  use_raster = TRUE, row_dend_gp = gpar(lwd = 0.25))

col_fun3 <- colorRamp2(
  seq(min(RLE.CREu1, na.rm = TRUE), max(RLE.CREp, na.rm = TRUE), length.out = 25),
  colorRampPalette(c("white", "black"))(25))

out3 <- Heatmap(
  RLE.CREu1,
  name = "Legand",
  col = col_fun3,
  clustering_method_rows = "ward.D2",
  cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
  row_dend_width = unit(8, "mm"),
  height = unit(0.5, "cm"),
  use_raster = TRUE, row_dend_gp = gpar(lwd = 0.25))

ht_list <- out1 %v% out2 %v% out3

pdf(paste0(path_fig2,"f2b.heatmap.pdf"), width = 2.1, height = 3)
draw(ht_list)
dev.off()

#=====================
#Needed
#f2c
cpg=read.delim(paste0(path_fig2_data,"tCRE_CpG_island_nooverlap.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpg$anno_region=factor(cpg$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
cpg=cpg[which(cpg$orientation == "Others"),]
cpg$signalID=gsub("n","min_",cpg$signalID)
cpg$signalID=factor(cpg$signalID, levels=c("min_20","min_40","min_60","min_80"))

tata=read.delim(paste0(path_fig2_data,"tCRE_TATA_nooverlap.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
tata$anno_region=factor(tata$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
tata=tata[which(tata$orientation == "Others"),]
tata$signalID=gsub("n","min_",tata$signalID)
tata$signalID=factor(tata$signalID, levels=c("min_1","min_2","min_3","min_4"))

#=========
C11=ggplot()+
  scale_color_npg()+
  geom_line(data = cpg[which(cpg$anno_region %in% c("promoter-like","enhancer-like")),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.9,0.8), legend.direction = "vertical")
C12=ggplot()+
  scale_color_npg()+
  geom_line(data = tata[which(tata$anno_region %in% c("promoter-like","enhancer-like")),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="TATA box")+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.02,0.04,0.06,0.08))+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.8,0.8), legend.direction = "vertical")

pdf(paste0(path_fig2,"f2c.tCRE_nonoverlap_CpG_TATA.pdf"), width = 1.8, height = 2)
grid.arrange(C11,C12, ncol=2, nrow = 1, widths = c(1.1,1))
dev.off()

#=====================
#Needed
#f2d
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE=CREanno[which(CREanno$representative == "Yes"),]

CRE$express="THP-1"
CRE$express[which(CRE$NSC>0 | CRE$Neuron>0)]="NSC or Neuron"
CRE$express[which(CRE$iPS>0)]="iPSC"

CRE$group="Un-marked"
CRE$group[which(CRE$K27ME3_iPS=="Yes")]="Repressed"
CRE$group[which(CRE$K27Ac_iPS=="Yes")]="Active"
CRE$group[which(CRE$K27ME3_iPS=="Yes" & CRE$K27Ac_iPS=="Yes")]="Co-marked"

CREab2=CRE[which(CRE$promoter_type == "enhancer-like" & CRE$CGI == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREab3=CRE[which(CRE$promoter_type == "enhancer-like" & CRE$TATA == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREab4=CRE[which(CRE$promoter_type == "enhancer-like" & CRE$CGInap == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREab5=CRE[which(CRE$promoter_type == "enhancer-like" & CRE$CGIap == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

CREab2$group2="CGI"
CREab3$group2="TATA"
CREab4$group2="CGInap"
CREab5$group2="CGIap"

CREab1=rbind(CREab2,CREab3,CREab4,CREab5)
CREab1$express=factor(CREab1$express, levels=c("THP-1","NSC or Neuron","iPSC"))
CREab2=CREab1%>%group_by(group2, group)%>%dplyr::summarise(tCRE=sum(count))
CREab2$label=paste0("n=",CREab2$tCRE)
CREab1$group=factor(CREab1$group,levels=c("Un-marked","Repressed","Active","Co-marked"))

f2d=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enhancer transcription with\niPSC chromatin state", fill=NULL) +
  scale_fill_manual(values=c("white","grey","black"), guide = guide_legend(reverse = TRUE) )+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.35), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(rows=vars(group2))+
  geom_bar(data=CREab1, mapping=aes(x=percent, y=group, fill=express), linewidth=0.2, color = "black", stat="identity", alpha=1) + 
  geom_text(data=CREab2, mapping=aes(x=1.02, y=group, label=label), size=1.8, hjust=0)+
  theme1 + theme(legend.position = "bottom" )
pdf(paste0(path_fig2,"f2d.tCRE_Expression_repressed_enhacer_dCGI.pdf"), width = 1.7, height = 2.1)
print(f2d)
dev.off()

#=====================
#Needed
#f2e
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE=CREanno[which(CREanno$representative == "Yes"),]

CRE$group="Un-marked"
CRE$group[which(CRE$K27ME3_iPS=="Yes")]="Repressed"
CRE$group[which(CRE$K27Ac_iPS=="Yes")]="Active"
CRE$group[which(CRE$K27ME3_iPS=="Yes" & CRE$K27Ac_iPS=="Yes")]="Co-marked"

CREa=reshape2::melt(CRE[which(CRE$promoter_type == "enhancer-like"),c("group","CGI","CGIap","CGInap","TATA")], id=1)
CREab=CREa%>%group_by(group, variable, value)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREab$value=factor(CREab$value, levels=c("Yes","No"))
CREab$label=paste0(CREab$group,"_",CREab$value)

CRE$'Un-marked'="No"
CRE$'Un-marked'[which(CRE$group == "Un-marked")]="Yes"
CRE$Active="No"
CRE$Active[which(CRE$group == "Active")]="Yes"
CRE$'Co-marked'="No"
CRE$'Co-marked'[which(CRE$group == "Co-marked")]="Yes"
CRE$Repressed="No"
CRE$Repressed[which(CRE$group == "Repressed")]="Yes"

CREac=reshape2::melt(CRE[which(CRE$promoter_type == "enhancer-like"),c("Repressed","Active","Co-marked","Un-marked","CGI","CGIap","CGInap","TATA")], id=c(1:4))
colnames(CREac)[c(5,6)]=c("variable2","value2")
CREac=reshape2::melt(CREac, id=c(5:6))
CREac$result=paste0(CREac$value,"_",CREac$value2)
CREac1=CREac%>%group_by(variable, variable2, result)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

CREad=spread(CREac1[,c(1:4)],key=3, value=4)
colnames(CREad)=c("group","variable","no_no","no_yes","yes_no","yes_yes")
for(i in 1:nrow(CREad)){
  GSEATasting <- matrix(c(CREad$no_no[i], CREad$no_yes[i], CREad$yes_no[i], CREad$yes_yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  CREad$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  CREad$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
CREad$label="***"
CREad$label[which(CREad$p.val>=0.001)]="**"
CREad$label[which(CREad$p.val>=0.01)]="*"
CREad$label[which(CREad$p.val>=0.05)]="n.s."
CREad$label1=paste0("OR=",signif(CREad$OR,2),CREad$label)
CREad$group=factor(CREad$group,levels=c("Un-marked","Repressed","Active","Co-marked"))

f2e=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enhancer structure with\niPSC chromatin state", fill="CGI/ CGIap/\nCGInap/ TATA") +
  scale_fill_manual(values=c("black","white"))+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.35), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(rows=vars(variable))+
  geom_bar(data=CREab, mapping=aes(x=percent, y=group, fill=value), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=CREad, mapping=aes(y=group, x=1.01, label=label1), size=1.8, color="black", hjust=0, vjust=0.5)+
  theme1+theme(legend.position=c(0.25,0.75))
pdf(paste0(path_fig2,"f2e.tCRE_repressed_enhancer_dCGI.pdf"), width = 1.8, height = 2.1)
print(f2e)
dev.off()

#=====================
#Needed
#f2f
CRE=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

dat0=CRE[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like"),]%>%group_by(promoter_type,SE_all)%>%dplyr::summarise(count=n())
dat1=CRE[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like"),]%>%group_by(promoter_type,SE_iPSC)%>%dplyr::summarise(count=n())
dat2=CRE[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like"),]%>%group_by(promoter_type,SE_NSC)%>%dplyr::summarise(count=n())
dat3=CRE[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like"),]%>%group_by(promoter_type,SE_Neuron)%>%dplyr::summarise(count=n())
dat4=CRE[which(CRE$promoter_type == "promoter-like" | CRE$promoter_type == "enhancer-like"),]%>%group_by(promoter_type,orientation)%>%dplyr::summarise(count=n())
dat0$group="Neuron-series"
dat1$group="iPSC-specific"
dat2$group="NSC-specific"
dat3$group="Neuron-specific"
dat4$group="Directionality"
colnames(dat0)[2]="variable"
colnames(dat1)[2]="variable"
colnames(dat2)[2]="variable"
colnames(dat3)[2]="variable"
colnames(dat4)[2]="variable"
dat=rbind(dat0,dat1,dat2,dat3,dat4)
dat$group2="Super enhancer"
dat$group2[which(dat$group=="Directionality")]="Directionality"
dat$group[which(dat$group=="Directionality")]=dat$promoter_type[which(dat$group=="Directionality")]
dat$group=factor(dat$group,levels=c("promoter-like","enhancer-like","iPSC-specific","NSC-specific","Neuron-specific","Neuron-series"))
dat$variable=factor(dat$variable, levels=c("Others","1D","2D","SE","TE"))

f2f=ggplot(dat, aes(x=group, y=count, fill=variable, group=variable)) + 
  scale_fill_manual(values=c("white","#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF"))+
  labs(x=NULL , y ="Number of tCRE" , title ="Directionality & super enhancer",  fill=NULL)+
  facet_grid(cols=vars(group2),  scales="free_x", space="free_x")+
  geom_bar(alpha=1, linewidth=0.25, stat="identity", color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2f.tCRE.number.directionality.SE.pdf"), width = 1.8, height = 1.4)
print(f2f)
dev.off() 

#=====================
#Needed
#f2g
data2=read.delim(paste0(path_fig2_data,"filtered_tCRE_for_conservation.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2$group4=factor(data2$group4,levels=c("enhancer-like","promoter-like","iPSC","NSC","Neuron"))

data4=data2%>%group_by(group2, group4, feature)%>%dplyr::summarise(median=median(mean_1001,na.rm=T),count=n())
data1wilcox3=data2%>%group_by(group2,group4)%>%dplyr::summarise(p=wilcox.test(mean_1001 ~feature, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f2g=ggplot() + 
  labs(x=NULL, y = "Mean PhastCon score",title= "Conservation of tCREs") +
  scale_fill_npg(guide=NULL)+
  facet_wrap(vars(group4), ncol=5, scale="free_x")+
  coord_cartesian(ylim=c(0,0.9))+
  ggdist::stat_halfeye(data=data2, mapping=aes(y=mean_1001, x=feature, fill=group2),adjust = .75, width = .6, .width = 0, justification = -.4, point_colour = NA, alpha=0.4)+
  geom_boxplot(data=data2, mapping=aes(y=mean_1001, x=feature, fill=group2), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.3, outlier.shape = NA) + 
  geom_text(data=data1wilcox3, mapping=aes(y=0.82, x=1.5, label=label),size=2.4)+
  geom_segment(data=data1wilcox3, mapping=aes(x=1,xend=2,y=0.8,yend=0.8), linewidth=0.25)+
  theme1+theme(legend.position="bottom",legend.box.margin = margin(t=-10, r=0, b=0, l=-10))
pdf(paste0(path_fig2,"f2g.enhancer.phastcon.pdf"), width = 3.3, height = 1.4)
print(f2g)
dev.off() 



#=====================
#Needed
#f2h
summaryRE=read.delim(paste0(path_fig2_data,"promoter_RE_FEresult.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
summaryRE$FE_logOR=log(summaryRE$OR)
summaryRE$sig_level="ns"
summaryRE$sig_level[which(summaryRE$p.val<0.05)]="*"
summaryRE$sig_level[which(summaryRE$p.val<0.01)]="**"
summaryRE$sig_level[which(summaryRE$p.val<0.001)]="***"

summaryRE$promoter_type=factor(summaryRE$promoter_type, levels=c("enhancer_2D","enhancer_1D","promoter_2D","promoter_1D", "CTCF-alone", "unclassed","enhancer-like","promoter-like"))
summaryRE$RE=factor(summaryRE$RE, levels=c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat"))
summaryRE$sig_level=factor(summaryRE$sig_level, levels=c("ns","*","**","***"))

f2h=ggplot() + 
  labs(x=NULL, y = "Promoter type",title= "Enrichment of repetitive element in tCREs") +
  geom_point(data=summaryRE, mapping=aes(y=promoter_type, x=RE, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.2,0.6,1.2,1.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2h.FE_promoter_type_tCRE_RE.pdf"), width = 2.5, height = 1.5)
print(f2h)
dev.off() 

#=====================
#Needed
#f2i
summarySTR=read.delim(paste0(path_fig2_data,"promoter_STR_FEresult.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
summarySTR$FE_logOR=log(summarySTR$OR)
summarySTR$sig_level="ns"
summarySTR$sig_level[which(summarySTR$p.val<0.05)]="*"
summarySTR$sig_level[which(summarySTR$p.val<0.01)]="**"
summarySTR$sig_level[which(summarySTR$p.val<0.001)]="***"

summarySTR$promoter_type=factor(summarySTR$promoter_type, levels=c("CTCF-alone", "unclassed","enhancer-like","promoter-like"))
summarySTR$STR=factor(summarySTR$STR, levels=c("CGG","CCG","CT","AG","GT","AC","CCT","AGG"))
summarySTR$sig_level=factor(summarySTR$sig_level, levels=c("ns","*","**","***"))

f2i=ggplot() + 
  labs(x=NULL, y = "Promoter type",title= "Enrichment of STR in tCRE summits") +
  geom_point(data=summarySTR, mapping=aes(y=promoter_type, x=STR, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.2,0.6,1.2,1.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"f2i.FE_promoter_type_tCRE_STR.pdf"), width = 2.5, height = 1.5)
print(f2i)
dev.off() 


