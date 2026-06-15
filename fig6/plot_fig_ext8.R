library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)
library(tidyr)
library(dplyr) 
library(ggrastr)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#============
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig6=paste0(primary_folder,"fig6/out/")
path_fig6_data=paste0(primary_folder,"fig6/data/")
setwd(path_fig6_data)

#============
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
#Figext8a

data2=read.delim(paste0(path_fig6_data,"end5_cluster_SNP_number_commononly.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data5=reshape2::melt(data2[,c(1:4,11,5:10)], id=c(1:5))
data5$group2="Ex5_cluster"
data5$group2[grep("exon",data5$variable)]="Exon"
data5$group3="eQTL"
data5$group3[grep("GWAS",data5$variable)]="GWAS"
data5$group3[grep("SNP",data5$variable)]="all_SNP"
data5=data5[which(data5$CpGTATA != "Others"),]
data6=data5[which(data5$ex5cluster_class%in%c("p_ncRNA","e_ncRNA","mRNA") & data5$group2=="Exon"),]%>%group_by(ex5cluster_class, group2, group3, CpGTATA,variable)%>%dplyr::summarise(mean=mean(value),sd=sd(value),per1000=mean(value/exon_coverage*1000, na.rm=T))
data7=data5[which(data5$ex5cluster_class%in%c("p_ncRNA","e_ncRNA","mRNA") & data5$group2=="Ex5_cluster"),]%>%group_by(ex5cluster_class, group2, group3, CpGTATA,variable)%>%dplyr::summarise(mean=mean(value),sd=sd(value),per1000=mean(value/ex5_length*1000))
data8=rbind(data6,data7)
data8$CpGTATA=factor(data8$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))

ex8a=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  scale_alpha_manual(values=c(0.5,1,1), guide=NULL)+
  labs(x=NULL , y ="Mean number of SNP per 1kb" , title ="SNP count from enhancer structure",  fill=NULL)+
  facet_grid(cols=vars(group2), rows=vars(group3), scales="free_y")+
  geom_bar(data=data8[which(data8$ex5cluster_class == "e_ncRNA"),], mapping=aes(x=CpGTATA, y=per1000, fill=CpGTATA, alpha=group3), linewidth=0.25, stat="identity")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8a.SNP_count_from_enhancer_structure.pdf"), width = 1.8, height = 1.7)
print(ex8a)
dev.off() 

#============
#separate celltype
#Figext8b&c

data3=read.delim(paste0(path_fig6_data,"HiC5kb_q001_ex5cluster_ATAC_GC.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
#this file contain a 5kb-bin in each row, bins contain redundant features (TATA, CGI, etc) were removed
n5_interact2=read.delim(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_summit_cell_collapsed.tsv.gz"), header=T, stringsAsFactors = F, check.names = T)
#this file contain all Hi-C intersect with ex5_cluster

i1=n5_interact2[grep("iPSC",n5_interact2$cell),]%>%group_by(source_group, CpGTATA_source, source_ex5)%>%dplyr::summarise(codingGene_connectivity=n())
s1=n5_interact2[grep("NSC",n5_interact2$cell),]%>%group_by(source_group, CpGTATA_source, source_ex5)%>%dplyr::summarise(codingGene_connectivity=n())
n1=n5_interact2[grep("NRN",n5_interact2$cell),]%>%group_by(source_group, CpGTATA_source, source_ex5)%>%dplyr::summarise(codingGene_connectivity=n())
datai1=left_join(data3[which(data3$iPSC > 0),c(1:3,5:6)], i1[,c(3:4)], by=c("n5_string"="source_ex5"), copy=F)
datas1=left_join(data3[which(data3$NSC > 0),c(1:3,5,7)], s1[,c(3:4)], by=c("n5_string"="source_ex5"), copy=F)
datan1=left_join(data3[which(data3$Neuron > 0),c(1:3,5,8)], n1[,c(3:4)], by=c("n5_string"="source_ex5"), copy=F)
datai1$cell="iPSC"
datas1$cell="NSC"
datan1$cell="Neuron"
colnames(datai1)[5]="ATAC_CPM"
colnames(datas1)[5]="ATAC_CPM"
colnames(datan1)[5]="ATAC_CPM"
data2=rbind(datai1,datas1,datan1)
data2$codingGene_connectivity[which(is.na(data2$codingGene_connectivity))]=0
data2=data2[which(data2$ex5cluster_class == "e_ncRNA"),]
data2$cell=factor(data2$cell, levels=c("iPSC","NSC","Neuron"))
data2$CpGTATA=factor(data2$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))

q2=data2%>%group_by(cell, CpGTATA)%>%dplyr::summarise(median=median(codingGene_connectivity), count=n(), ATACmedian=median(ATAC_CPM))

#p-value
data4a=data2[which(data2$CpGTATA %in% c("CGIap","Null")),c("cell","CpGTATA","codingGene_connectivity","ATAC_CPM")]
data4b=data2[which(data2$CpGTATA %in% c("CGInap","Null")),c("cell","CpGTATA","codingGene_connectivity","ATAC_CPM")]
data4c=data2[which(data2$CpGTATA %in% c("TATA","Null")),c("cell","CpGTATA","codingGene_connectivity","ATAC_CPM")]
data4d=data2[which(data2$CpGTATA %in% c("CGInap","TATA")),c("cell","CpGTATA","codingGene_connectivity","ATAC_CPM")]

data4a$group2="CGIap"
data4b$group2="CGInap"
data4c$group2="TATA"
data4d$group2="CGInap_TATA"
data4=rbind(data4a,data4b,data4c,data4d)

data5=data4%>%group_by(cell,group2)%>%dplyr::summarise(p=wilcox.test(codingGene_connectivity ~CpGTATA, alternative = "two.sided")$p.value,
                                                                        p2=wilcox.test(ATAC_CPM ~CpGTATA, alternative = "two.sided")$p.value)
data5$label="***"
data5$label[which(data5$p>=0.001)]="**"
data5$label[which(data5$p>=0.01)]="*"
data5$label[which(data5$p>=0.05)]="n.s."
data5$label2="***"
data5$label2[which(data5$p2>=0.001)]="**"
data5$label2[which(data5$p2>=0.01)]="*"
data5$label2[which(data5$p2>=0.05)]="n.s."

ex8b=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="Hi-C connection per ex5_cluster" , title ="Coding gene promoter connectivity from enhancer",  fill=NULL)+
  coord_cartesian(ylim=c(0,30))+
  facet_grid(cols=vars(cell), scales="free_y")+
  geom_boxplot(data=data2, mapping=aes(x=CpGTATA, y=codingGene_connectivity, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=q2, mapping=aes(y=median, x=CpGTATA, label=median), size=1.8, vjust=-0.2, hjust=0.5, color="black")+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=26.25, x=2, label=label),size=2.4)+
  annotate("segment", x=1,xend=3,y=26,yend=26, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=24.75, x=2.5, label=label),size=2.4)+
  annotate("segment",x=2,xend=3,y=24.5,yend=24.5, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=27.75, x=3.5, label=label),size=2.4)+
  annotate("segment",x=3,xend=4,y=27.5,yend=27.5, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=29.25, x=3, label=label),size=2.4)+
  annotate("segment",x=2,xend=4,y=29,yend=29, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8b.HiC5kb_q001_ex5cluster_CGI_TATA_3_cells.pdf"), width = 2.5, height = 1.7)
print(ex8b)
dev.off() 

#===correlation between ATAC CPM and connectivity
data2k=data2[which(!is.na(data2$codingGene_connectivity)),]
data3k=data2k%>%group_by(cell)%>%dplyr::summarise(count=n(),spearman=cor.test(ATAC_CPM,codingGene_connectivity, method="spearman")$estimate, p=cor.test(ATAC_CPM,codingGene_connectivity, method="spearman")$p.value)

ex8c=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="ATAC CPM of 5kb Hi-C bin" , title ="Chromatin accessibility",  fill=NULL)+
  coord_cartesian(ylim=c(0,240))+
  facet_grid(cols=vars(cell), scales="free_y")+
  geom_boxplot(data=data2, mapping=aes(x=CpGTATA, y=ATAC_CPM, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.3, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=q2, mapping=aes(y=ATACmedian, x=CpGTATA, label=signif(ATACmedian,3)), angle=90, size=1.8, vjust=-0.5, hjust=0.5, color="black")+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=175, x=2, label=label2),size=2.4)+
  annotate("segment", x=1,xend=3,y=170,yend=170, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=155, x=2.5, label=label2),size=2.4)+
  annotate("segment", x=2,xend=3,y=150,yend=150, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=195, x=3.5, label=label2),size=2.4)+
  annotate("segment", x=3,xend=4,y=190,yend=190, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=215, x=2.5, label=label2),size=2.4)+
  annotate("segment", x=1,xend=4,y=210,yend=210, linewidth=0.25)+
  geom_text(data=data3k, mapping=aes(y=235, x=2.5, label=paste0("rho=",signif(spearman, 3),"***")), size=1.8)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8c.HiC5kb_ATACCPM_ex5cluster_CGI_TATA_exclude_zero_HiC.pdf"), width = 2.5, height = 1.7)
print(ex8c)
dev.off() 

#============
#ext8d
data3=read.delim(paste0(path_fig6_data,"HiC5kb_q001_ex5cluster_ATAC_GC.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
#this file contain a 5kb-bin in each row, bins contain redundant features (TATA, CGI, etc) were removed

data3$ex5cluster_class=factor(data3$ex5cluster_class, levels=c("p_ncRNA","e_ncRNA"))
data3$CpGTATA=factor(data3$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))

#p-value
data4a=data3[which(data3$ex5cluster_class == "e_ncRNA" & data3$CpGTATA %in% c("CGIap","Null")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]
data4b=data3[which(data3$ex5cluster_class == "e_ncRNA" & data3$CpGTATA %in% c("CGInap","Null")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]
data4c=data3[which(data3$CpGTATA %in% c("TATA","Null")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]
data4d=data3[which(data3$ex5cluster_class == "e_ncRNA" & data3$CpGTATA %in% c("CGInap","TATA")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]
data4e=data3[which(data3$ex5cluster_class != "e_ncRNA" & data3$CpGTATA %in% c("CGI","Null")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]
data4f=data3[which(data3$ex5cluster_class != "e_ncRNA" & data3$CpGTATA %in% c("CGI","TATA")),c("ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent")]

data4a$group2="CGIap"
data4b$group2="CGInap"
data4c$group2="TATA"
data4d$group2="CGInap_TATA"
data4e$group2="CGI"
data4f$group2="CGI_TATA"
data4=rbind(data4a,data4b,data4c,data4d,data4e,data4f)

data5=data4%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(codingGene_connectivity ~CpGTATA, alternative = "two.sided")$p.value,
                                                                   p2=wilcox.test(bin_5kb_GCcontent ~CpGTATA, alternative = "two.sided")$p.value)

data5$label="***"
data5$label[which(data5$p>=0.001)]="**"
data5$label[which(data5$p>=0.01)]="*"
data5$label[which(data5$p>=0.05)]="n.s."
data5$label2="***"
data5$label2[which(data5$p2>=0.001)]="**"
data5$label2[which(data5$p2>=0.01)]="*"
data5$label2[which(data5$p2>=0.05)]="n.s."

#===correlation between GC content and connectivity
data3h=data3[which(!is.na(data3$codingGene_connectivity)),]
data4h=data3h%>%group_by(ex5cluster_class)%>%dplyr::summarise(count=n(),spearman=cor.test(bin_5kb_GCcontent,codingGene_connectivity, method="spearman")$estimate, p=cor.test(bin_5kb_GCcontent,codingGene_connectivity, method="spearman")$p.value)

ex8d=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="GC content of 5kb Hi-C bin" , title ="GC content",  fill=NULL)+
  coord_cartesian(ylim=c(0.25,1))+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(cols=vars(ex5cluster_class), scales="free_x", space="free_x")+
  geom_boxplot(data=data3, mapping=aes(x=CpGTATA, y=bin_5kb_GCcontent, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=0.76, x=2, label=label2),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGI"),], mapping=aes(y=0.76, x=2, label=label2),size=2.4)+
  annotate("segment", x=1,xend=3,y=0.75,yend=0.75, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=0.71, x=2.5, label=label2),size=2.4)+
  annotate("segment", x=2,xend=3,y=0.7,yend=0.7, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=0.81, x=3.5, label=label2),size=2.4)+
  annotate("segment", x=3,xend=4,y=0.8,yend=0.8, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=0.86, x=3, label=label2),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGI_TATA"),], mapping=aes(y=0.86, x=3, label=label2),size=2.4)+
  annotate("segment", x=2,xend=4,y=0.85,yend=0.85, linewidth=0.25)+
  geom_text(data=data4h, mapping=aes(y=0.95, x=2.5, label=paste0("rho=",signif(spearman,3),"***")), size=1.8)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8d.HiC5kb_GCcontent_ex5cluster_CGI_TATA_exclude_zero_HiC.pdf"), width = 2.4, height = 1.7)
print(ex8d)
dev.off() 

#===========
#ex8e
#SE and 2D
data3=read.delim(paste0(path_fig6_data,"HiC5kb_q001_ex5cluster_ATAC_GC.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
colnames(data3)[c(12:13)]=c("2D","SE")

data2b=reshape::melt(data3[,c(2,4, 12, 13)], id=c(1,2))
data2b=data2b[which(!is.na(data2b$value)),]
data2b=data2b[which(data2b$value!="Others"),]
data2b$value[which(data2b$value %in% c("2D","SE"))]="Yes"
data2b$value[which(data2b$value %in% c("1D","TE"))]="No"
data2b$value=factor(data2b$value, levels=c("Yes","No"))

data2b$ex5cluster_class=factor(data2b$ex5cluster_class, levels=c("p_ncRNA","e_ncRNA","other_ncRNA"))
q2=data2b%>%group_by(ex5cluster_class, variable, value)%>%dplyr::summarise(median=median(codingGene_connectivity), count=n(),atleast1=length(which(codingGene_connectivity>0)))

#p-value
data5b=data2b%>%group_by(ex5cluster_class,variable)%>%dplyr::summarise(p=wilcox.test(codingGene_connectivity ~value, alternative = "two.sided")$p.value)
data5b$label="***"
data5b$label[which(data5b$p>=0.001)]="**"
data5b$label[which(data5b$p>=0.01)]="*"
data5b$label[which(data5b$p>=0.05)]="n.s."

ex8e=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF"), guide=NULL)+
  labs(x=NULL , y ="Hi-C connection per ex5_cluster" , title ="Coding gene promoter connectivity",  fill=NULL)+
  coord_cartesian(ylim=c(0,25))+
  facet_grid(cols=vars(variable), scales="free")+
  geom_boxplot(data=data2b[which(data2b$ex5cluster_class == "e_ncRNA"),], mapping=aes(x=value, y=codingGene_connectivity, fill=variable), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=q2[which(q2$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=median, x=value, label=median), size=1.8, vjust=-0.2, hjust=0.5, color="black")+
  geom_text(data=data5b[which(data5b$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=24, x=1.5, label=label),size=2.4)+
  annotate("segment", x=1, xend=2, y=22, yend=22, linewidth=0.25)+
  theme1
pdf(paste0(path_fig6,"ex8e.HiC5kb_q001_ex5cluster_CGI_TATA.SE_2D.pdf"), width = 1.6, height = 1.7)
print(ex8e)
dev.off() 

#============
#ex8f
k22=read.delim(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_FE_structure.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#this data contain only e_ncRNA

k22$source_CpGTATA=sapply(strsplit(k22$group,"_"),"[",1)
k22$target_CpGTATA=sapply(strsplit(k22$group,"_"),"[",2)
k22$target_CpGTATA=paste0("mRNA:",k22$target_CpGTATA)
k22$source_CpGTATA=factor(k22$source_CpGTATA,levels=c("CGIap","CGInap","Null","TATA"))
k22$label2="***"
k22$label2[which(k22$p.val>=0.001)]="**"
k22$label2[which(k22$p.val>=0.01)]="*"
k22$label2[which(k22$p.val>=0.05)]="n.s."

ex8f=ggplot(k22, aes(x=source_CpGTATA, y=OR, fill=source_CpGTATA)) + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x="Regulatory elements from source bin" , y ="Odds ratio" , title ="Enrichment with Hi-C connection", fill=NULL)+
  scale_y_continuous(limits=c(0,1.5), breaks=c(0,0.5,1,1.5))+
  facet_grid(cols=vars(target_CpGTATA))+
  geom_bar(alpha=0.7, linewidth=0.25, stat="identity", color="black")+
  geom_text(mapping=aes(y=OR+0.02, x=source_CpGTATA, label=label2),size=2.2, hjust=0.5, vjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8f.HiC.FE_structure.e_ncRNA.pdf"), width = 2, height = 1.7)
print(ex8f)
dev.off() 

#============
#ex8g
#ABC model
ABC=read.delim(paste0(path_fig6_data,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone_summary.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
ABC=ABC[which(ABC$CpGTATA != "Others"),]

ABC$cell=factor(ABC$cell, levels=c("iPSC","NSC","Neuron"))
ABC$CpGTATA=factor(ABC$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))
ABC00=ABC%>%group_by(cell,CpGTATA)%>%dplyr::summarise(count=median(count), ABC.Score=median(ABC.Score))

#p-value
data4a=ABC[which(ABC$CpGTATA %in% c("CGIap","Null")),c("cell","CpGTATA","count","ABC.Score")]
data4b=ABC[which(ABC$CpGTATA %in% c("CGInap","Null")),c("cell","CpGTATA","count","ABC.Score")]
data4c=ABC[which(ABC$CpGTATA %in% c("TATA","Null")),c("cell","CpGTATA","count","ABC.Score")]
data4d=ABC[which(ABC$CpGTATA %in% c("CGInap","TATA")),c("cell","CpGTATA","count","ABC.Score")]

data4a$group2="CGIap"
data4b$group2="CGInap"
data4c$group2="TATA"
data4d$group2="CGInap_TATA"
data4=rbind(data4a,data4b,data4c,data4d)

data5=data4%>%group_by(cell,group2)%>%dplyr::summarise(p=wilcox.test(count ~CpGTATA, alternative = "two.sided")$p.value,
                                                       p2=wilcox.test(ABC.Score ~CpGTATA, alternative = "two.sided")$p.value)
data5$label="***"
data5$label[which(data5$p>=0.001)]="**"
data5$label[which(data5$p>=0.01)]="*"
data5$label[which(data5$p>=0.05)]="n.s."
data5$label2="***"
data5$label2[which(data5$p2>=0.001)]="**"
data5$label2[which(data5$p2>=0.01)]="*"
data5$label2[which(data5$p2>=0.05)]="n.s."

ex8g=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="# of connection by ABC model" , title ="Gene connectivity from enhancers",  fill=NULL)+
  coord_cartesian(ylim=c(0,27))+
  facet_grid(cols=vars(cell), scales="free_y")+
  geom_boxplot(data=ABC, mapping=aes(x=CpGTATA, y=count, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=ABC00, mapping=aes(y=count, x=CpGTATA, label=count), size=1.8, vjust=-0.2, hjust=0.5, color="black")+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=22.25, x=2, label=label),size=2.4)+
  annotate("segment", x=1,xend=3,y=22,yend=22, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=20.25, x=2.5, label=label),size=2.4)+
  annotate("segment",x=2,xend=3,y=20,yend=20, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=24.25, x=3.5, label=label),size=2.4)+
  annotate("segment",x=3,xend=4,y=24,yend=24, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=26.25, x=3, label=label),size=2.4)+
  annotate("segment",x=2,xend=4,y=26,yend=26, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8g.ABC_002_ex5cluster_CGI_TATA_3_cells.pdf"), width = 2.5, height = 1.7)
print(ex8g)
dev.off() 

#===============
#ex8h
#same data processing from ex8g

ex8h=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="ABC score per enhancer" , title ="Overall enhancer activity",  fill=NULL)+
  coord_cartesian(ylim=c(0,0.2))+
  facet_grid(cols=vars(cell), scales="free_y")+
  geom_boxplot(data=ABC, mapping=aes(x=CpGTATA, y=ABC.Score, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.35, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=ABC00, mapping=aes(y=ABC.Score, x=CpGTATA, label=signif(ABC.Score,2)), size=1.8, angle=90, vjust=-0.55, hjust=0.5, color="black")+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=0.145, x=2, label=label2),size=2.4)+
  annotate("segment", x=1,xend=3,y=0.14,yend=0.14, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=0.125, x=2.5, label=label2),size=2.4)+
  annotate("segment",x=2,xend=3,y=0.12,yend=0.12, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=0.165, x=3.5, label=label2),size=2.4)+
  annotate("segment",x=3,xend=4,y=0.16,yend=0.16, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=0.185, x=3, label=label2),size=2.4)+
  annotate("segment",x=2,xend=4,y=0.18,yend=0.18, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex8h.ABC_score_ex5cluster_CGI_TATA_3_cells.pdf"), width = 2.7, height = 1.7)
print(ex8h)
dev.off() 

