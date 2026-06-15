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
library(ggrastr)

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


#============================================================================================
#Needed
#ex3bc
combine=read.delim(paste0(path_fig2_data,"rose_SE_info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

combine$group=factor(combine$group, levels=c("iPSC","NSC","Neuron"))
cutoff=data.frame(cbind(group=c("iPSC","NSC","Neuron"),cutoff=c(16537,8480,13252), SE_region=c(1932,1278,1790)))
cutoff$label=paste0("cutoff = ",cutoff$cutoff,"\nSE regions = ",cutoff$SE_region)
cutoff$group=factor(cutoff$group, levels=c("iPSC","NSC","Neuron"))

ex3b=ggplot() + 
  labs(x="Ranked region", y = "Count",title= "Identification of super enhancer regions") +
  scale_color_npg(guide="none")+
  coord_cartesian(ylim=c(0,200000))+
  geom_line(data=combine, mapping=aes(y=total_count, x=rank, color=factor(group)), linewidth=0.25)+
  geom_hline(data=combine, mapping=aes(yintercept=cutoff), linetype = "dashed", linewidth=0.25, color="black")+
  geom_text(data=cutoff, mapping=aes(y=160000, x=50, label=label), hjust=0, vjust=0, size=2)+
  facet_grid(cols=vars(group), scale="free_x")+
  theme1
pdf(paste0(path_fig2,"ex3b.super.enhancer.27ac.pdf"), width = 3.5, height = 1.5)
print(ex3b)
dev.off() 

combine$size=combine$V4-combine$V3
combine1=combine[which(combine$total_count > combine$cutoff),]
combine2=combine1%>%group_by(group)%>%dplyr::summarise(median=median(size))

ex3c=ggplot() + 
  labs(y="size (bp)",title= "Size of SE regions", x=NULL) +
  scale_fill_npg(guide="none")+
  coord_cartesian(ylim=c(0,200000))+
  geom_boxplot(data=combine1, mapping=aes(y=size, x=group, fill=group), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.8, outlier.shape = NA) + 
  geom_text(data=combine2, mapping=aes(y=median+3000, x=group, label=signif(median,3)), size=2, vjust=0, hjust=0.5, color="black")+
  theme1
pdf(paste0(path_fig2,"ex3c.super.enhancer.27ac3.size.pdf"), width = 1.5, height = 1.3)
print(ex3c)
dev.off() 

#=====================
#Needed
#ex3d
phastcon=read.delim(paste0(path_fig2_data,"tCRE_nooverlap_phastcon_directionality_all_way_plot.tsv.gz"),header=T, stringsAsFactors=F, check.names = F)
phastcon$promoter_type=factor(phastcon$promoter_type, levels=c("promoter-like","enhancer-like","CTCF-alone"))
phastcon=phastcon[which(phastcon$signalID %in% c("4way","17way","30way","100way")),]
phastcon$signalID=factor(phastcon$signalID, levels=c("4way","17way","30way","100way"))

ex3d=ggplot(phastcon[which(phastcon$promoter_type %in% c("promoter-like","enhancer-like")),], aes(x=positionV3, y=score, color=variable, group=variable)) + 
  scale_color_npg()+
  coord_cartesian(ylim = c(0.1, NA), xlim=c(-1000,1000))+
  labs(x="Stranded distance relative to tCRE summit" , y ="Mean phastCon score" , title ="Conservation of bidirectional tCRE",  color=NULL)+
  facet_grid(rows=vars(promoter_type), cols=vars(signalID),  scales="free_y")+
  geom_vline(xintercept=c(-500, 500), linetype="dashed", linewidth=0.15, color="black")+
  geom_line(alpha=0.75, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.85,0.87))
pdf(paste0(path_fig2,"ex3d.directionality.phastcon.pdf"), width = 3, height = 2)
print(ex3d)
dev.off() 

#=====================
#Needed
#ex3e
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

CREanno1=CREanno[which(CREanno$promoter_type %in% c("promoter-like","enhancer-like") & CREanno$region1001_rep == "Yes"),]
CREanno1=CREanno1[which(CREanno1$orientation %in% c("2D","1D")),]
CREanno1$group=paste0(substr(CREanno1$promoter_type,1,1),"_",CREanno1$orientation)

CREanno1$group2="Null"
CREanno1$group2[which(CREanno1$any_CpG_island == "Yes")]="CGI"
CREanno1$group2[which(CREanno1$TATA_box == "Yes")]="TATA"
CREanno1$group2[which(CREanno1$any_CpG_island == "Yes" & CREanno1$TATA_box == "Yes")]="Mix"
CREanno1$group2=factor(CREanno1$group2, levels=c("Null","TATA","Mix","CGI"))

CREanno2=CREanno1%>%group_by(group, group2)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
CREanno2$label=paste0(signif(CREanno2$percent,2)*100,"%")
CREanno2$label[which(CREanno2$group2 %in% c("Null","Mix"))]=NA

ex3e=ggplot(CREanno2, aes(x=as.factor(group), y=percent, fill=group2)) + 
  scale_fill_manual(values=c("white","#4DBBD5FF","darkorchid3","#E64B35FF"))+
  labs(fill=NULL, x=NULL, y ="% of tCREs", title="Regulatory elements")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85,alpha=0.7)+
  geom_text(aes(label=label), size=1.8, position=position_stack(vjust=0.5), color="black", angle=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.32,0.8))
pdf(paste0(path_fig2,"ex3e.CPG_TATA_1D_2D.pdf"), width = 1.25, height = 2)
print(ex3e)
dev.off() 

#=====================
#Needed
#ex3f
data2=read.delim(paste0(path_fig2_data,"filtered_tCRE_for_conservation.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2=data2[which(data2$group2=="Directionality"),] #keep only 1D 2D
data2$group=paste0(substr(data2$group4,1,1),"_",data2$feature)

data2a=reshape2::melt(data2[,c(6,8,3,10)], id=c(3,4))
data2a$variable=factor(data2a$variable, levels=c("TATA","CGI"))
data2b=data2a%>%group_by(group,variable, value)%>%dplyr::summarise(median=median(phastCon17_mean_1001, na.rm=T),count=n())
data2c=data2a%>%group_by(group,variable)%>%dplyr::summarise(p=wilcox.test(phastCon17_mean_1001 ~value, alternative = "two.sided")$p.value)
data2c$label="***"
data2c$label[which(data2c$p>=0.001)]="**"
data2c$label[which(data2c$p>=0.01)]="*"
data2c$label[which(data2c$p>=0.05)]="n.s."

ex3f=ggplot() + 
  labs(x=NULL, y = "Mean PhastCon score",title= "Conservation from regulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF","#7E6148FF"),guide=NULL)+
  facet_grid(rows=vars(variable), cols=vars(group), scale="free_x")+
  coord_cartesian(ylim=c(0,0.95))+
  ggdist::stat_halfeye(data=data2a, mapping=aes(y=phastCon17_mean_1001, x=value, fill=variable),adjust = .75, width = .5, .width = 0, justification = -.4, point_colour = NA, alpha=0.3)+
  geom_boxplot(data=data2a, mapping=aes(y=phastCon17_mean_1001, x=value, fill=variable), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA, alpha=0.7) + 
  geom_text(data=data2b, mapping=aes(y=median, x=value, label=signif(median,3)), angle=90, size=2, vjust=(-0.8), hjust=0.5, color="black")+
  geom_text(data=data2c, mapping=aes(y=0.9, x=1.5, label=label),size=2.4)+
  geom_segment(data=data2c, mapping=aes(x=1,xend=2,y=0.85,yend=0.85), linewidth=0.25)+
  theme1
pdf(paste0(path_fig2,"ex3f.PhastCon_CpG_TATA_1D_2D.pdf"), width = 2.9, height = 2)
print(ex3f)
dev.off() 

#=====================
#Needed
#ex3g
data0=read.delim(paste0(path_fig2_data,"all.CAGE.summary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data0$total_read=data0$no_raw_read-data0$raw_read_in_cluster
data0$non_q=data0$no_raw_read-data0$passed_read_num
data0$q=data0$passed_read_num-data0$CTSS_unG
data0$total_read2=data0$passed_read_num-data0$Qualified_read_in_cluster
data0$non_q_cluster=data0$raw_read_in_cluster-data0$Qualified_read_in_cluster
data0$non_q_outside_cluster=data0$non_q-data0$non_q_cluster
colnames(data0)[c(5,7,8,10,11,12,13,14,15)]=c("Qualified Read w/ unencodedG","Inside TSScluster","q_Inside TSScluster", "Outside TSScluster","Non-qualified Read","Qualified Read","q_Outside TSScluster","nonq_Inside TSScluster","nonq_Outside TSScluster")

data0b=reshape2::melt(data0[,c(9,1,5,11,12)], id=c(1,2))
data0b=data0b%>%group_by(seq,lib)%>%dplyr::mutate(percent=value/sum(value))
data0b$lib=factor(data0b$lib, levels=c("iPSC","NSC","Neuron","THP-1","dTHP-1"))
data0b$variable=factor(data0b$variable, levels=c("Non-qualified Read","Qualified Read","Qualified Read w/ unencodedG"))
data0b$seq=factor(data0b$seq, levels=c("CFC-seq","CAGE"))

ex3g=ggplot(data0b, aes(x=lib, y=percent, fill=variable)) + 
  scale_fill_manual(values=c("white","grey","black"))+
  labs(fill=NULL, x=NULL, y ="% of reads", title="Read quality & un-encoded G")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85)+
  facet_grid(cols=vars(seq), scale="free", space="free")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = "bottom", legend.direction = "vertical")
pdf(paste0(path_fig2,"ex3g.SCAFE.stepwise.read.pdf"), width = 2, height = 1.8)
print(ex3g)
dev.off() 

#=====================
#Needed
#ex3h
summary0=read.delim(paste0(path_fig2_data,"CTSS_1nt_CAGE_CFC_compare.summary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

summary0$group=factor(summary0$group, levels= c("iPSC","NSC","Neuron","Neuron-series"))
summary0$group2=factor(summary0$group2, levels=c("CFC-\nspecific","both","CAGE-\nspecific"))
summary0$class[which(summary0$class == "CTSS")]="QualifiedRead" 
summary0$class[which(summary0$class == "unG")]="Un-encoded G"  
summary0$label=paste0(signif(summary0$percent,3)*100,"%")
summary0$label[which(summary0$group2!="both")]=NA

ex3h=ggplot(summary0, aes(x=as.factor(group), y=percent, fill=group2)) + 
  scale_fill_manual(values=c("#E64B35FF","darkorchid3","#4DBBD5FF"))+
  labs(fill=NULL, x=NULL, y ="% of read", title="Concordant TSS position")+
  facet_grid(cols=vars(class))+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85, alpha=0.5)+
  geom_text(aes(label=label), size=1.8, position=position_stack(vjust=0.5), color="black", angle=90)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex3h.CTSS.compare.pdf"), width = 2, height = 1.8)
print(ex3h)
dev.off() 

#=====================
#Needed
#ex3i
count_m=read.delim(paste0(path_fig2_data,"tCRE.identification.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

count_m1=count_m%>%group_by(CFC1,CAGE1)%>%dplyr::summarise(count=n(), group=1)
count_m2=count_m%>%group_by(CFC2,CAGE2)%>%dplyr::summarise(count=n(), group=2)
count_m3=count_m%>%group_by(CFC3,CAGE3)%>%dplyr::summarise(count=n(), group=3)
count_m1$tech=c("CAGE-\nspecific","CFC-\nspecific","Both")
count_m2$tech=c("CAGE-\nspecific","CFC-\nspecific","Both")
count_m3$tech=c("None","CAGE-\nspecific","CFC-\nspecific","Both")
count_m=rbind(count_m1[,c(4,5,3)],count_m2[,c(4,5,3)],count_m3[,c(4,5,3)])
count_m$tech=factor(count_m$tech, levels=c("CFC-\nspecific","Both","CAGE-\nspecific","None"))
count_m=count_m%>%group_by(group)%>%dplyr::mutate(percent=count/sum(count))
count_m$label=paste0(signif(count_m$percent,3)*100,"%")
count_m$label[which(count_m$tech!="Both")]=NA


ex3i=ggplot(count_m, aes(x=as.factor(group), y=count, fill=tech)) + 
  scale_fill_manual(values=c("#E64B35FF","darkorchid3","#4DBBD5FF","white"))+
  labs(fill=NULL, x="Read per tCRE cutoff", y ="Number of tCREs", title="Joint SCAFE\ntCREs identification")+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85, alpha=0.5)+
  geom_text(aes(label=label), size=1.8, position=position_stack(vjust=0.5), color="black", angle=90)+
  theme1
pdf(paste0(path_fig2,"ex3i.tCRE.identification.pdf"), width = 1.4, height = 1.8)
print(ex3i)
dev.off() 

#=====================
#Needed
#ex3j
summary1=read.delim(paste0(path_fig2_data,"tCRE_quantification_CFC_CAGE_RLE_CPM.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
summary1$cell=factor(summary1$cell, levels=c("iPSC","NSC","Neuron"))
summary2=summary1%>%group_by(cell)%>%dplyr::summarise(pearson=cor(CFC, CAGE, method="pearson"), count=n())
summary2$label=paste0("r = ", signif(summary2$pearson,3))
ex3j=ggplot(summary1, aes(x=CAGE, y=CFC)) + 
  labs(x="CAGE" , y ="CFC-seq",  title="Joint SCAFE\ntCRE quantification (log10CPM)",color=NULL)+
  geom_text(data=summary2, aes(x=-1, y=3.5, label=label), vjust=0, hjust=0, size=2, color="black")+
  facet_wrap(~cell , ncol=2 , scales="fixed")+
  geom_point_rast(shape=21, size=0.1, color="black", alpha=0.05)+
  theme1
pdf(paste0(path_fig2,"ex3j.tCRE_quantification_CAGE_CFC.pdf"), width = 1.6, height = 1.9)
print(ex3j)
dev.off() 


#=====================
#Needed
#ex3k
asummary1=read.delim(paste0(path_fig2_data,"RE_MAPQ_compare.NA_removed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# MAPQ -> 1 2 3  4
# CFC  -> 0 5 10 20
# CAGE -> 0 2 3  225
colnames(asummary1)[c(8:10)]=c("CFC-seq","CAGE_PE","CAGE_SE")
asummary2=reshape2::melt(asummary1[,c(1,7,8,9,10)], id=c(1,2))
asummary2=asummary2[which(asummary2$repClass %in% c("DNA","LTR","LINE","SINE")),]
asummary2$repClass=factor(asummary2$repClass, levels=c("DNA","LTR","LINE","SINE"))
ex3k=ggplot(asummary2, aes(x=max_MQ, y=value, color=variable, group=variable)) + 
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("0","2|5","3|10","225|20"))+
  labs(x="MAPQ threshold (CAGE|CFC)" , y ="% of TSS cluster",  title="TSS clusters at repeat elements",color=NULL)+
  facet_wrap(~repClass , ncol=7 , scales="fixed")+
  geom_line(linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position=c(0.9,0.2))
pdf(paste0(path_fig2,"ex3k.RE_TSS_cluster.pdf"), width = 2.4, height = 1.7)
print(ex3k)
dev.off() 

#=====================
#Needed
#ex3l
allcluster1=read.delim(paste0(path_fig2_data,"ont_ss.Neuronalone.cluster.final.rm_na.result.tsv.gz"),header=T, stringsAsFactors = T, check.names = F)
allcluster0=allcluster1[which( !is.na(allcluster1$repFamily)),]
allcluster0=tidyr::separate_rows(allcluster0, repName, repClass, repFamily, repSize, milliDiv, sep=";")
clusterrep=allcluster0%>%group_by(V4,repClass)%>%dplyr::slice_max(milliDiv) #each cluster only has one entry in each category
clusterrep$group5="others"
clusterrep$group5[which(clusterrep$ontCAGE_max_MAPQ>=20 & clusterrep$SEssCAGE_max_MAPQ <225)]="CFC"
clusterrep=clusterrep[which(clusterrep$repClass %in% c("DNA","LTR","LINE","SINE")),]
clusterrep$repClass=factor(clusterrep$repClass, levels=c("DNA","LTR","LINE","SINE"))

repcluster4=clusterrep%>%group_by(repClass)%>%dplyr::summarise(p=wilcox.test(as.numeric(milliDiv) ~group5, alternative = "two.sided")$p.value)
repcluster4$label="***"
repcluster4$label[which(repcluster4$p>=0.001)]="**"
repcluster4$label[which(repcluster4$p>=0.01)]="*"
repcluster4$label[which(repcluster4$p>=0.05)]="n.s."
clusterrep$group5=factor(clusterrep$group5, levels=c("CFC","others"))

ex3l=ggplot() + 
  geom_boxplot(data=clusterrep, mapping=aes(y=as.numeric(milliDiv)/1000, x=group5, fill=group5), linewidth=0.25, color = "black", notch = FALSE, width = 0.6, outlier.shape = NA, alpha=0.8) + 
  labs(x="CFC-seq versus others", y = "Divergence",title= "Evolutionary age of repeat elements", fill=NULL) +
  scale_fill_manual(values=c("#E64B35FF","#00A087FF"),guide=NULL)+
  facet_grid(cols=vars(repClass))+
  geom_text(data=repcluster4, mapping=aes(y=0.5, x=1.5, label=label), size=2.2, color="black")+
  annotate("segment", x=1,xend=2,y=0.45,yend=0.45, linewidth=0.25)+
  coord_cartesian(ylim=c(0,0.55))+
  theme1
pdf(paste0(path_fig2,"ex3l.RE_MAPQ_compare.NA_removed.evolutionary_age.pdf"), width = 2.5, height = 1.7)
print(ex3l)
dev.off() 

#=====================
#Needed
#ex3m
CREanno2=read.delim(paste0(path_fig2_data,"tCRE_RE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

CREanno2$RE=factor(CREanno2$RE, levels=c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat", "Null", "Multiple"))
CREanno2$promoter_type=factor(CREanno2$promoter_type, levels=c("promoter-like","enhancer-like","CTCF-alone","unclassed"))
CREanno3=CREanno2%>%group_by(promoter_type, RE)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count))
CREanno4=CREanno2%>%group_by(RE)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count))

ex3m=ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15), breaks=c(0,0.25,0.5,0.75,1))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of tCRE" , title ="Presence of repeat elements in tCRE",  fill=NULL)+
  geom_bar(data=CREanno3, mapping=aes(x=promoter_type, y=percent, fill=RE), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=CREanno3[which(CREanno3$RE=="Null"),], mapping=aes(x=promoter_type, y=0.3, label=paste0(signif(percent*100,3),"%")), size=1.8, color="black")+
  geom_text(data=CREanno3[which(CREanno3$RE=="Null"),], mapping=aes(x=promoter_type, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"ex3m.Presence_of_RE_inside_tCRE.pdf"), width = 2, height = 1.7)
print(ex3m)
dev.off() 




