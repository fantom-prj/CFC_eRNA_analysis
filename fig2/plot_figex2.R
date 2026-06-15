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
               legend.background = element_rect(fill="white", color="black", size=0.25),
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

#============================================================================================
#Needed
#ex2b
atac3=read.delim(paste0(path_fig2_data,"Neuron_tCRE_ATAC.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

atac3$ATAC_supported=factor(atac3$ATAC_supported, levels=c("No_ATAC","ATAC"))
atac3$promoter_type=factor(atac3$promoter_type, levels=c("promoter-like","enhancer-like","unclassed"))
atac4=atac3[which(atac3$variable == "Neuron-series"),]

ex2b=ggplot() + 
  scale_fill_manual(values=c("white","grey"))+
  labs(x=NULL, y ="count", title="tCRE from Neuron series (n = 65,646)", fill=NULL)+
  geom_bar(data=atac4, mapping=aes(x=promoter_type, y=count, fill=ATAC_supported), stat="identity", size=0.25, color="black")+
  geom_text(data=atac4[which(atac4$ATAC_supported == "ATAC"),], mapping=aes(x=promoter_type, y=count-200, label=label), vjust=1, size=2, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2b.Neuron_tCRE_ATAC.pdf"), width = 2.2, height = 1.5)
print(ex2b)
dev.off() 

#============================================================================================
#Needed
#ex2cd
CRE1=read.delim(paste0(path_fig2_data,"tCRE_unclass_Neuron.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

ss1=CRE1%>%group_by(ATAC,typeStr,cell_type)%>%dplyr::summarise(count=n())
ss1$cell_type=factor(ss1$cell_type, levels=c("GENCODE","iPSC-specific","NSC-specific","Neuron-specific","non-specific"))
ss1$ATAC=factor(ss1$ATAC, levels=c("withATAC","noATAC"))
ss1$typeStr[which(ss1$typeStr == "gene_tss")]="Annotated"
ss1$typeStr[which(ss1$typeStr == "unanno_tss")]="Un-annotated"
ss1=ss1%>%group_by(ATAC,typeStr)%>%dplyr::mutate(percent=count/sum(count))
ss1$label=paste0(signif(ss1$percent*100,2),"%")
ss1$label[which(ss1$typeStr == "Annotated")]=NA
ss1$label[which(ss1$ATAC == "noATAC")]=NA
ex2c=ggplot() + 
  labs(x=NULL, y = "Number of tCRE",title= "Source of unclassed tCRE", fill=NULL) +
  scale_fill_npg()+
  facet_grid(cols=vars(ATAC))+
  geom_bar(data=ss1, mapping=aes(y=count, x=typeStr, fill=cell_type), size=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=ss1, mapping=aes(y=count, x=typeStr, group=cell_type, label=label), size=1.8, color="black", position = position_stack(vjust = 0.5, reverse=F))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2c.unclassed_tCRE_source.pdf"), width = 2.3, height = 1.5)
print(ex2c)
dev.off()


ss2=CRE1[which(CRE1$typeStr == "unanno_tss"),]%>%group_by(ATAC,typeStr,class2,ChromHMM)%>%dplyr::summarise(count=n())
ss2$ATAC=factor(ss2$ATAC, levels=c("withATAC","noATAC"))
ss2$typeStr[which(ss2$typeStr == "unanno_tss")]="Un-annotated"
ss2=ss2%>%group_by(ATAC,ChromHMM,)%>%dplyr::mutate(percent=count/sum(count))
ss2$label=paste0(signif(ss2$percent*100,2),"%")
ss2$label[which(ss2$ChromHMM == "ChroHMM")]=NA
ss2$label[which(ss2$ATAC == "withATAC")]=NA
ss3=ss2%>%group_by(ATAC,ChromHMM,)%>%dplyr::summarise(count=sum(count))
ss3$label=paste0("n=", ss3$count)

ex2d=ggplot() + 
  labs(x=NULL, y = "% of tCRE",title= "TSS summit position", fill=NULL) +
  scale_fill_npg()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1.1), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(ATAC))+
  geom_bar(data=ss2, mapping=aes(y=percent, x=ChromHMM, fill=class2), size=0.25, color = "black", stat="identity") + 
  geom_text(data=ss2, mapping=aes(y=percent, x=ChromHMM, group=class2, label=label), size=1.8, color="white", position = position_stack(vjust = 0.5, reverse=F))+
  geom_text(data=ss3, mapping=aes(y=1.01, x=ChromHMM, label=label), size=1.8, color="black", vjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2d.unclassed_tCRE_source2.pdf"), width = 2.3, height = 1.5)
print(ex2d)
dev.off()

#=====================
#Needed
#ex2e

gc2=read.delim(paste0(path_fig2_data,"tCRE_both_GCcontent_nooverlap_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gc2$promoter_type=factor(gc2$promoter_type, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
gc2$group[which(gc2$group == "random")]="background"
gc2$group[which(gc2$group == "anno_region")]="tested\nregion"
gc2$group=factor(gc2$group, levels=c("tested\nregion", "background"))

phastcon=read.delim(paste0(path_fig2_data,"tCRE_nooverlap_phastcon_all_way_plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
phastcon$variable=factor(phastcon$variable, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
phastcon=phastcon[which(phastcon$signalID %in% c("4way","17way","30way","100way")),]
phastcon$signalID=factor(phastcon$signalID, levels=c("4way","17way","30way","100way"))

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

inr=read.delim(paste0(path_fig2_data,"tCRE_INR_nooverlap_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
inr$anno_region=factor(inr$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
inr=inr[which(inr$orientation == "Others"),]
inr$signalID=gsub("n","min_",inr$signalID)
inr$signalID=factor(inr$signalID, levels=c("min_1","min_2","min_3","min_4"))
#=========
#labs(x="Distance relative to cluster summit (nt)", y="% of tCRE overlap signal peaks")
C9=ggplot(gc2[which(gc2$promoter_type != "CTCF-alone"),], aes(x=position, y=GCprecent/100, color=group)) + 
  scale_color_npg()+
  scale_y_continuous(labels = scales::percent)+
  labs(x="tCRE summit" , y ="% of GC" , title ="GC content",  color=NULL)+
  facet_grid(rows=vars(promoter_type) ,  scales="free_y")+
  geom_line(alpha=0.75, size=0.25)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C10=ggplot(phastcon[which(phastcon$variable != "CTCF-alone"),], aes(x=positionV3, y=score, color=signalID, group=signalID)) + 
  scale_color_npg()+
  coord_cartesian(ylim = c(0.1, NA))+
  labs(x=NULL , y ="Mean phastCon score" , title ="Conservation",  color=NULL)+
  facet_grid(rows=vars(variable) ,  scales="free_y")+
  geom_line(alpha=0.75, size=0.25)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C11=ggplot()+
  scale_color_npg()+
  geom_line(data = cpg[which(cpg$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C12=ggplot()+
  scale_color_npg()+
  geom_line(data = tata[which(tata$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="TATA box")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C13=ggplot()+
  scale_color_npg()+
  geom_line(data = inr[which(inr$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="Initiator")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")

pdf(paste0(path_fig2,"ex2e.tCRE_nonoverlap_CpG_TATA_inia_phase.pdf"), width = 4.6, height = 3)
grid.arrange(C9,C10, C11,C12, C13, ncol=5, nrow = 1, widths = c(1.1,1.1,1,1,1))
dev.off()

#=====================
#Needed
#ex2f
cpg=read.delim(paste0(path_fig2_data,"tCRE_CpG_island_nooverlap.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpg=cpg[which(cpg$anno_region %in% c("enhancer-AP","enhancer-NAP")),]
cpg$label="Major"
cpgm=read.delim(paste0(path_fig2_data,"tCRE_CpG_island_minorstrand.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpgm=cpgm[which(cpgm$anno_region %in% c("enhancer-AP","enhancer-NAP")),]
cpgm$label="Minor"

cpg=rbind(cpg,cpgm)
cpg$signalID=gsub("n","min_",cpg$signalID)
cpg$signalID=factor(cpg$signalID, levels=c("min_20","min_40","min_60","min_80"))

C11=ggplot()+
  scale_color_npg()+
  geom_line(data = cpg, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_wrap( ~label+anno_region, scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
pdf(paste0(path_fig2,"ex2f1.tCRE_minorstrand_CpG_TATA.pdf"), width = 1.8, height = 2.2)
print(C11)
dev.off()
#==
cpgm=read.delim(paste0(path_fig2_data,"tCRE_CpG_island_minorstrand.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpgm=cpgm[which(cpgm$anno_region == "enhancer-like" & cpgm$orientation == "Others"),]
cpgm$signalID=gsub("n","min_",cpgm$signalID)
cpgm$signalID=factor(cpgm$signalID, levels=c("min_20","min_40","min_60","min_80"))

tatam=read.delim(paste0(path_fig2_data,"tCRE_TATA_minorstrand.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
tatam=tatam[which(tatam$anno_region == "enhancer-like" & tatam$orientation == "Others"),]
tatam$signalID=gsub("n","min_",tatam$signalID)
tatam$signalID=factor(tatam$signalID, levels=c("min_1","min_2","min_3","min_4"))
C12=ggplot()+
  scale_color_npg()+
  geom_line(data = cpgm, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C13=ggplot()+
  scale_color_npg()+
  geom_line(data = tatam, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, size=0.25)+
  labs(color=NULL, x=NULL, y=NULL, title="TATA box")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")

pdf(paste0(path_fig2,"ex2f2.tCRE_minorstrand_CpG_TATA.pdf"), width = 1.8, height = 1.2)
grid.arrange(C12,C13, ncol=2, nrow = 1, widths = c(1,1))
dev.off()

#=====================
#Needed
#ex2g
CRE1b=read.delim(paste0(path_fig2_data,"tCRE_minor_strand_summary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

k1=CRE1b%>%group_by(antisense_tCRE,promoter_type, CpGTATA_main,CpGTATA_minor)%>%dplyr::summarise(count1=n())%>%dplyr::mutate(percent=count1/sum(count1))
k2=k1[grep("CGI",k1$CpGTATA_main),]%>%group_by(antisense_tCRE,promoter_type,CpGTATA_minor)%>%dplyr::mutate(group2="CGI")
k3=k1[grep("Null",k1$CpGTATA_main),]%>%group_by(antisense_tCRE,promoter_type,CpGTATA_minor)%>%dplyr::mutate(group2="Null")
k4=k1[grep("TATA",k1$CpGTATA_main),]%>%group_by(antisense_tCRE,promoter_type,CpGTATA_minor)%>%dplyr::mutate(group2="TATA")
k5=rbind(k2,k3,k4)

k5$CpGTATA_minor=factor(k5$CpGTATA_minor, levels=c("CGI","Null","Mix","TATA"))
k5$antisense_tCRE=factor(k5$antisense_tCRE, levels=c("No","Yes"))
k5$group2=factor(k5$group2, levels=c("CGI","Null","TATA"))

ex2g=ggplot()+
  geom_bar(data=k5[which(k5$promoter_type=="enhancer-like"),], mapping=aes(x=antisense_tCRE, y=count1, fill=CpGTATA_minor), stat="identity", color="black", linewidth=0.2)+
  labs(x="Presence of minor strand", y="count", fill="Class of\nminor strand", title="Enhancer-like tCRE")+
  facet_grid(cols=vars(group2), scales = "free_y")+
  scale_fill_npg()+
  theme1
pdf(paste0(path_fig2,"ex2g.tCRE_CpGTATA_antisense_number.pdf"), width = 2, height = 1.2)
print(ex2g)
dev.off() 

#=====================
#Needed
#ex2i
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE=CREanno[which(CREanno$representative == "Yes" & CREanno$promoter_type == "enhancer-like"),]

CRE$Null="Yes"
CRE$Null[which(CRE$CGI == "Yes" | CRE$TATA == "Yes")]="No"

CREa=reshape2::melt(CRE[,c("Distance_PLScCRE","CGI","CGIap","CGInap","TATA","Null")], id=1)
CREa=CREa[which(CREa$value == "Yes"),]

CREa$variable=factor(CREa$variable, levels=c("CGI","CGIap","CGInap","Null","TATA"))
CREab=CREa%>%group_by(variable)%>%dplyr::summarise(count=n(),median=median(Distance_PLScCRE),mean=mean(Distance_PLScCRE))

ex2i=ggplot() + 
  labs(x=NULL, y = "Minimum distance (nt)",title= "Distance of enhancer\nfrom promoter cCRE") +
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"),guide=NULL)+
  ggdist::stat_halfeye(data=CREa, mapping=aes(y=Distance_PLScCRE+1, x=variable, fill=variable),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=CREa, mapping=aes(y=Distance_PLScCRE+1, x=variable, fill=variable), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  #geom_point(data=CREab, mapping=aes(y=mean, x=variable, color=variable), size=0.3, shape=19, color="red")+
  geom_text(data=CREab, mapping=aes(y=median+1, x=variable, label=paste0(signif(median,3),"nt (n=",count,")")), size=2, angle=90, vjust=-0.7, hjust=0.5, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex2i.tCRE_CpGTATA_PLS_distance.pdf"), width = 1.5, height = 1.8)
print(ex2i)
dev.off() 

#=====================
#Needed
#ex2j
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE=CREanno[which(CREanno$representative == "Yes"),]

need=c("iPS","NSC","Neuron")
CRE1=CRE[,c(grep(need[1],colnames(CRE)),which(colnames(CRE)%in% c("promoter_type","CGI","CGIap","CGInap","TATA")))]
CRE2=CRE[,c(grep(need[2],colnames(CRE)),which(colnames(CRE)%in% c("promoter_type","CGI","CGIap","CGInap","TATA")))]
CRE3=CRE[,c(grep(need[3],colnames(CRE)),which(colnames(CRE)%in% c("promoter_type","CGI","CGIap","CGInap","TATA")))]
CRE3=CRE3[,c(2:10)]
CRE1$cell="iPSC"
CRE2$cell="NSC"
CRE3$cell="Neuron"
colnames(CRE1)=c("tCRE_count","SE","K27me3","K27ac","promoter_type","CGI","TATA","CGInap","CGIap","cell")
colnames(CRE2)=colnames(CRE1)
colnames(CRE3)=colnames(CRE1)
CRE1=rbind(CRE1,CRE2,CRE3)
CRE1=CRE1[which(CRE1$promoter_type == "enhancer-like" & CRE1$tCRE_count>0),]
CRE1$Repressed="No"
CRE1$Repressed[which(CRE1$K27me3 == "Yes" & CRE1$K27ac == "No")]="Yes"
CRE1$Active="No"
CRE1$Active[which(CRE1$K27me3 == "No" & CRE1$K27ac == "Yes")]="Yes"
CRE1$'Co-marked'="No"
CRE1$'Co-marked'[which(CRE1$K27me3 == "Yes" & CRE1$K27ac == "Yes")]="Yes"
CRE1$'Un-marked'="No"
CRE1$'Un-marked'[which(CRE1$K27me3 == "No" & CRE1$K27ac == "No")]="Yes"
CRE1$group="Un-marked"
CRE1$group[which(CRE1$K27me3=="Yes")]="Repressed"
CRE1$group[which(CRE1$K27ac=="Yes")]="Active"
CRE1$group[which(CRE1$K27me3=="Yes" & CRE1$K27ac=="Yes")]="Co-marked"

#for percentage
CREa1=reshape2::melt(CRE1[,c("cell","group","CGI","CGIap","CGInap","TATA")], id=c(1,2))
CREab1=CREa1%>%group_by(cell, group, variable, value)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREab1$value=factor(CREab1$value, levels=c("Yes","No"))
CREab1$label=paste0(CREab1$group,"_",CREab1$value)

CREak=reshape2::melt(CRE1[,c("cell","Repressed","Active","Co-marked","Un-marked","CGI","CGIap","CGInap","TATA")], id=c(1:5))
colnames(CREak)[c(6,7)]=c("variable2","value2")
CREak=reshape2::melt(CREak, id=c(1,6:7))
CREak$result=paste0(CREak$value,"_",CREak$value2)
CREak1=CREak%>%group_by(cell, variable, variable2, result)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

CREak2=spread(CREak1[,c(1:5)],key=4, value=5)
colnames(CREak2)=c("cell","group","variable","no_no","no_yes","yes_no","yes_yes")
for(i in 1:36){
  GSEATasting <- matrix(c(CREak2$no_no[i], CREak2$no_yes[i], CREak2$yes_no[i], CREak2$yes_yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  CREak2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  CREak2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
CREak2$label="***"
CREak2$label[which(CREak2$p.val>=0.001)]="**"
CREak2$label[which(CREak2$p.val>=0.01)]="*"
CREak2$label[which(CREak2$p.val>=0.05)]=""
CREak2$label1=paste0("OR=",signif(CREak2$OR,2),CREak2$label)
CREab1$group=factor(CREab1$group,levels=c("Un-marked","Repressed","Active","Co-marked"))
CREab1$cell=factor(CREab1$cell, levels=c("iPSC", "NSC", "Neuron"))
CREak2$cell=factor(CREak2$cell, levels=c("iPSC", "NSC", "Neuron"))

ex2j=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enrichment with chromatin state", fill="CGI/ CGIap/\nCGInap/ TATA") +
  scale_fill_manual(values=c("black","white"))+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.42), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(cell),rows=vars(variable))+
  geom_bar(data=CREab1, mapping=aes(x=percent, y=group, fill=value), linewidth=0.2, color = "black", stat="identity", alpha=1) + 
  geom_text(data=CREak2, mapping=aes(y=group, x=1.01, label=label1), size=1.8, color="black", hjust=0, vjust=0.5)+
  theme1+theme(legend.position=c(0.1,0.5))
pdf(paste0(path_fig2,"ex2j.tCRE_cell_type3_expressing_poised_enhace_dCGI.pdf"), width = 4.3, height = 1.9)
print(ex2j)
dev.off()





