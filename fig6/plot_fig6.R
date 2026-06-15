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
#Needed
#Fig.6a&b
data=read.delim(paste0(path_fig6_data,"eQTL_GWAS_table5_exon_acceptor_donor.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)

data1=data%>%group_by(region, group, transcript_novelty)%>%dplyr::summarise(count=n())
data1$transcript_novelty=factor(data1$transcript_novelty, levels=c("Transcript from novel gene","Novel isoform","ENST"))
data1$region=factor(data1$region, levels=c("Exon","Donor","Acceptor"))

f6b=ggplot(data1[which(data1$group %in% c("eQTL_all","GWAS")),], aes(x=group, y=(count/100), fill=transcript_novelty)) + 
  scale_fill_manual(values=c("#E64B35FF","black","white"))+
  scale_x_discrete(labels=c("eQTL","GWAS"))+
  labs(x=NULL , y ="Number of SNP (hundreds)" , title ="Annotated SNPs in SALA transcriptome",  fill=NULL)+
  facet_wrap(vars(region), ncol=4,  scales="free")+
  geom_bar(alpha=1, linewidth=0.25, stat="identity", color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6b.eQTL_GWAS_table5_exon_acceptor_donor.pdf"), width = 3.1, height = 1.7)
print(f6b)
dev.off() 


#reduce group and calculate the %
data1$group2="novel"
data1$group2[which(data1$transcript_novelty=="ENST")]="ENST"
data2=data1%>%group_by(region, group, group2)%>%dplyr::summarise(count=sum(count))%>%dplyr::mutate(base=max(count))
data2$percent=data2$count/data2$base
data2$label=paste0(signif(data2$percent*100,2),"%")
data2$label[which(data2$label=="100%")]=NA
data2$group2=factor(data2$group2, levels=c("novel","ENST"))
data2$region=factor(data2$region, levels=c("Exon","Donor","Acceptor"))

f6a=ggplot(data2[which(data2$group %in% c("eQTL_all","GWAS")),], aes(x=group, y=percent, fill=group2)) + 
  scale_fill_manual(values=c("#E64B35FF","white"))+
  scale_x_discrete(labels=c("eQTL","GWAS"))+
  labs(x=NULL , y ="Number of SNP" , title ="Annotated SNPs in SALA transcriptome",  fill=NULL)+
  facet_wrap(vars(region), ncol=4,  scales="fixed")+
  geom_bar(alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data2[which(data2$group %in% c("eQTL_all","GWAS")),], mapping=aes(x=group, y=1, group=group2, label=label), size=1.8, color="black", vjust=0)+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6a.eQTL_GWAS_table5_exon_acceptor_donor_percent.pdf"), width = 2.2, height = 1.7)
print(f6a)
dev.off() 

#============
#fig6c
data2=read.delim(paste0(path_fig6_data,"end5_cluster_SNP_number_commononly.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data3=reshape2::melt(data2[,c(1:2,11:15)], id=c(1:3))
data3$group2="Ex5_cluster"
data3$group2[grep("exon",data3$variable)]="Exon"
data3$group3="eQTL"
data3$group3[grep("GWAS",data3$variable)]="GWAS"
data3=data3[which(data3$CpGTATA != "Others"),]
data4=data3[which(data3$ex5cluster_class%in%c("e_ncRNA")),]%>%group_by(ex5cluster_class, group2, group3, CpGTATA,variable)%>%dplyr::summarise(mean=mean(value, na.rm =T),sd=sd(value, na.rm =T))

data2a=data3[which(data3$CpGTATA %in% c("CGIap","Null")),]
data2b=data3[which(data3$CpGTATA %in% c("CGInap","Null")),]
data2c=data3[which(data3$CpGTATA %in% c("TATA","Null")),]
data2d=data3[which(data3$CpGTATA %in% c("CGInap","TATA")),]
data2a$group4="CGIap"
data2b$group4="CGInap"
data2c$group4="TATA"
data2d$group4="CGInap_TATA"
data2e=rbind(data2a,data2b,data2c,data2d)
data2e=data2e[which(data2e$ex5cluster_class %in% c("e_ncRNA")),]

data2f=data2e%>%group_by(ex5cluster_class,group2, group3,group4)%>%dplyr::summarise(p=wilcox.test(value ~CpGTATA, alternative = "two.sided")$p.value)
data2f$label="***"
data2f$label[which(data2f$p>=0.001)]="**"
data2f$label[which(data2f$p>=0.01)]="*"
data2f$label[which(data2f$p>=0.05)]="n.s."

data3$CpGTATA=factor(data3$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))

f6c=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="% of SNP from background" , title ="Annotated SNPs from enhancer structure",  fill=NULL)+
  scale_y_continuous(labels = scales::percent)+
  facet_grid(cols=vars(group2), rows=vars(group3), scales="free_y")+
  coord_cartesian(ylim=c(0,0.08))+
  geom_violin(data=data3[which(data3$ex5cluster_class == "e_ncRNA"),], mapping=aes(x=CpGTATA, y=value, fill=CpGTATA), linewidth=0.25, alpha = 0.8)+
  geom_point(data=data4[which(data4$ex5cluster_class == "e_ncRNA"),], mapping=aes(x=CpGTATA, y=mean), size=0.2, color="red")+
  geom_text(data=data2f[which(data2f$ex5cluster_class == "e_ncRNA" & data2f$group4 == "CGIap"),], mapping=aes(y=0.06, x=2, label=label),size=2.4)+
  geom_text(data=data2f[which(data2f$ex5cluster_class == "e_ncRNA" & data2f$group4 == "CGInap"),], mapping=aes(y=0.05, x=2.5, label=label),size=2.4)+
  geom_text(data=data2f[which(data2f$ex5cluster_class == "e_ncRNA" & data2f$group4 == "TATA"),], mapping=aes(y=0.07, x=3.5, label=label),size=2.4)+
  geom_text(data=data2f[which(data2f$ex5cluster_class == "e_ncRNA" & data2f$group4 == "CGInap_TATA"),], mapping=aes(y=0.08, x=3, label=label),size=2.4)+
  annotate("segment", x=1.1, xend=2.9, y=0.058, yend=0.058, linewidth=0.25)+
  annotate("segment", x=2.1, xend=2.9, y=0.048, yend=0.048, linewidth=0.25)+
  annotate("segment", x=3.1, xend=3.9, y=0.068, yend=0.068, linewidth=0.25)+
  annotate("segment", x=2.1, xend=3.9, y=0.078, yend=0.078, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6c.rate_from_commonSNPs.eQTL_GWAS_CGI_TATA.e_ncRNA.pdf"), width = 2, height = 1.7)
print(f6c)
dev.off() 

#============
#Fig6d
#Needed
data3=read.delim(paste0(path_fig6_data,"HiC5kb_q001_ex5cluster_ATAC_GC.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
#this file contain a 5kb-bin in each row, bins contain redundant features (TATA, CGI, etc) were removed

data3$ex5cluster_class=factor(data3$ex5cluster_class, levels=c("p_ncRNA","e_ncRNA"))
data3$CpGTATA=factor(data3$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))
q2=data3%>%group_by(ex5cluster_class, CpGTATA)%>%dplyr::summarise(median=median(codingGene_connectivity,na.rm=T), count=n(),atleast1=length(which(!is.na(codingGene_connectivity))))

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

f6d=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="Hi-C connection per ex5_cluster" , title ="Coding gene promoter connectivity",  fill=NULL)+
  coord_cartesian(ylim=c(0,40))+
  facet_grid(cols=vars(ex5cluster_class), scales="free_x", space="free_x")+
  geom_boxplot(data=data3, mapping=aes(x=CpGTATA, y=codingGene_connectivity, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=q2, mapping=aes(y=median, x=CpGTATA, label=median), size=1.8, vjust=-0.2, hjust=0.5, color="black")+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=33.5, x=2, label=label),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGI"),], mapping=aes(y=33.5, x=2, label=label),size=2.4)+
  annotate("segment", x=1,xend=3,y=33,yend=33, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=30.5, x=2.5, label=label),size=2.4)+
  annotate("segment", x=2,xend=3,y=30,yend=30, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=36.5, x=3.5, label=label),size=2.4)+
  annotate("segment", x=3,xend=4,y=36,yend=36, linewidth=0.25)+
  geom_text(data=data5[which(data5$group2 == "CGInap_TATA"),], mapping=aes(y=39.5, x=3, label=label),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGI_TATA"),], mapping=aes(y=39.5, x=3, label=label),size=2.4)+
  annotate("segment", x=2,xend=4,y=39,yend=39, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6d.HiC5kb_q001_ex5cluster_CGI_TATA_exclude_zero_HiC.pdf"), width = 1.7, height = 1.7)
print(f6d)
dev.off() 

#============
#Needed
#f6e
data8=read.delim(paste0(path_fig6_data,"CpGTATA_repeat_element_FE.n5_exon.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)

data8$feature1=factor(data8$feature1, levels=c("TATA","Null","CGInap","CGIap"))
data8$feature2=factor(data8$feature2, levels=c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat"))
data8$sig_level=factor(data8$sig_level, levels=c("ns","*","**","***"))

f6e=ggplot() + 
  labs(x=NULL, y = "Enhancer regulatory group",title= "Enrichment with repetitive elements") +
  facet_grid(rows=vars(group))+
  geom_point(data=data8[which(data8$gene_group %in% c("e_ncRNA")),], mapping=aes(y=feature1, x=feature2, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradientn(colors = c("#3C5488FF", "white", "#DC0000FF"), values = scales::rescale(c(-3.95, 0, 1.7)), limits = c(-3.95, 1.7),)+ 
  scale_size_manual(values=c(0.2,0.6,1.2,1.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6e.e_ncRNA_FE_CpGTATA_repeat_element_ex5cluster_exon.pdf"), width = 2.3, height = 1.7)
print(f6e)
dev.off() 


#==================
#Needed
#f6f
data=read.delim(paste0(path_fig6_data,"ex5_cluster_LTR_distribution_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data$anno_region=factor(data$anno_region, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
data$CpGTATA=factor(data$CpGTATA, levels=c("CGI","TATA","Null"))
data1=data%>%group_by(anno_region,CpGTATA)%>%dplyr::summarise(size=unique(size))
data1$label=paste0(data1$CpGTATA,": ",data1$size)
f6f <- ggplot()+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF","grey"))+
  coord_cartesian(xlim=c(-2000,2000), ylim=c(0,0.52))+
  geom_line(data = data[which(data$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=V4, y=V5, color=CpGTATA, group=CpGTATA), linewidth=0.25)+
  geom_text(data=data1[which(data1$CpGTATA == "CGI" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2000, y=0.4, label=label), vjust=0, hjust=0, size=1.8, color="#00A087FF")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2000, y=0.5, label=label), vjust=0, hjust=0, size=1.8, color="#3C5488FF")+
  geom_text(data=data1[which(data1$CpGTATA == "Null" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2000, y=0.45, label=label), vjust=0, hjust=0, size=1.8, color="grey")+
  facet_grid(cols=vars(anno_region))+
  labs(color=NULL, x="Distance from the stranded ex5_cluster summit", y="% of ex5_cluster", title="LTR distribution from different regulatory elements")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6f.LTR_distribution.pdf"), width = 3.3, height = 1.7)
print(f6f)
dev.off() 

#============
#f6g
library(ComplexHeatmap)
library(circlize)
library(grid)
RLE.ex5=read.delim(paste0(path_fig6_data,"ex5_cluster_Neuron_THP1.RLE.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
data2a=read.delim(paste0(path_fig6_data,"TE_group_CpGTATA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2a$CpGTATA[grep("CGI",data2a$CpGTATA)]="CGI"
data2a=data2a[which(data2a$ex5cluster_class %in% c("e_ncRNA") & data2a$ CpGTATA %in% c("CGI","Null","TATA")),]

RLE.ex5=log10(RLE.ex5+0.01)
RLE.ex5=data.frame(RLE.ex5)
RLE.ex5[RLE.ex5>3]=3

RLE.ex5t=RLE.ex5[which(rownames(RLE.ex5)%in% data2a$n5_string[which(data2a$CpGTATA=="TATA")]),]
RLE.ex5tltr=RLE.ex5t[which(rownames(RLE.ex5t)%in% data2a$n5_string[which(data2a$group_ex5cluster == "LTR")]),]
RLE.ex5txltr=RLE.ex5t[-which(rownames(RLE.ex5t)%in% rownames(RLE.ex5tltr)),]

col_fun5 <- colorRamp2(c(min(RLE.ex5tltr), max(RLE.ex5tltr)), c("white", "#DC0000FF"))
col_fun6 <- colorRamp2(c(min(RLE.ex5txltr), max(RLE.ex5txltr)), c("white", "#7E6148FF"))

out5 <- Heatmap(as.matrix(RLE.ex5tltr),
                name = "value",
                col = col_fun5,
                cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE, row_dend_gp = gpar(lwd = 0.3),
                clustering_method_rows = "ward.D2",
                row_dend_width = unit(10, "mm"))
out6 <- Heatmap(as.matrix(RLE.ex5txltr),
                name = "value",
                col = col_fun6,
                cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE, row_dend_gp = gpar(lwd = 0.3),
                clustering_method_rows = "ward.D2",
                row_dend_width = unit(10, "mm"))

out5@matrix_param$height <- unit(nrow(RLE.ex5tltr), "null")
out6@matrix_param$height <- unit(nrow(RLE.ex5txltr), "null")
pdf(paste0(path_fig6,"f6g.CpGTATA_ex5cluster_heatmap_LTR.pdf"), width = 1.9, height = 1.6)
draw(out5 %v% out6)
dev.off()

#============
#f6h
data7=read.delim(paste0(path_fig6_data,"LTR_NFY_TEAD_positive_number.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data7$label[which(data7$label == "None")]="Others"
data7$variable=factor(data7$variable, levels=c("non_iPSC","all_iPSC","Co_marked","Active","Repressed","Un_marked"))
data7$label=factor(data7$label, levels=c("NF-Y","Both","TEAD","Others"))
data7$percent[which(data7$label %in% c("Both"))]=NA

f6h=ggplot(data=data7[which(data7$CpGTATA == "TATA"),], mapping=aes(y=count, x=variable, fill=label)) + 
  labs(y="# of ex5_cluster", x = NULL, title= "TF binding in LTR TATA-box ex5_cluster", fill="TF\nbinding") +
  facet_wrap(vars(ex5cluster_class),ncol=2, scale="free_y")+
  geom_bar(mapping=aes (alpha=variable), linewidth=0.2, stat="identity", color="black")+
  geom_text(mapping=aes(label=percent), position = position_stack(vjust = 0.5), size=1.8, alpha=1, color="black")+
  scale_alpha_manual(values=c(0.25,0.8,0.8,0.8,0.8,0.8,0.8), guide=NULL)+
  scale_fill_manual(values=c("Others"="grey", "NF-Y"="orange", "Both"="purple", "TEAD"="red"))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6h.NFY.TEAD.LTR_TATA.number.pdf"), width = 3.1, height = 2.1)
print(f6h)
dev.off() 

#=======
data5=read.delim(paste0(path_fig6_data,"LTR_ex5_cluster_ATAC_CPM.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data5=data5[which(data5$CpGTATA == "TATA" & data5$TFbind != "Both"),]
data5$group="Non-iPSC"
data5$group[which(data5$iPSC_t == "Yes")]="All-iPSC"
data5$TFbind[which(data5$group == "Non-iPSC")]="Non-iPSC"
data5$TFbind=factor(data5$TFbind, levels=c("Non-iPSC","Others","NF-Y","TEAD"))
data5a=data5%>%group_by(ex5cluster_class,TFbind)%>%dplyr::summarise(median=median(ATACcpm_iPSC),count=n())

f6i=ggplot() + 
  labs(x=NULL , y ="iPSC ATAC CPM" , title ="Chromatin accessibility of\nLTR-TATA ex5_clusters",  fill=NULL)+
  coord_cartesian(ylim=c(0,80))+
  facet_grid(rows=vars(ex5cluster_class), scales="free_y")+
  geom_boxplot(data=data5, mapping=aes(x=TFbind, y=ATACcpm_iPSC, fill=TFbind, alpha=TFbind), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=data5a, mapping=aes(y=median, x=TFbind, label=signif(median,2)), angle=90, size=1.8, vjust=-1.1, hjust=0.5, color="black")+
  geom_text(data=data5a, mapping=aes(y=75, x=TFbind, label=paste0("n=",count)), size=1.8, color="black")+
  scale_alpha_manual(values=c(0.25,0.8,0.8,0.8), guide=NULL)+
  scale_fill_manual(values=c("Non-iPSC"="grey","Others"="grey", "NF-Y"="orange", "TEAD"="red"), guide=NULL)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"f6i.NFY.TEAD.LTR_TATA.ATACcpm.pdf"), width =1.4, height = 2.1)
print(f6i)
dev.off() 


