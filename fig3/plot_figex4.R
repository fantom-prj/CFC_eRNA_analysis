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
path_fig3=paste0(primary_folder,"fig3/out/")
path_fig3_data=paste0(primary_folder,"fig3/data/")

#=====================
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

#======================
#Needed
#ex4bc
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
table5$transcript_novelty[which(table5$transcript_novelty == "Transcript from novel gene")]="Transcript from novel TU"
table5$polyA[which(table5$polyA == "Yes")]="p(A)"
table5$polyA[which(table5$polyA == "No")]="non-p(A)"
table5$polyA=factor(table5$polyA, levels=c("p(A)","non-p(A)"))

n5support=table5%>%group_by(transcript_novelty,n5_support)%>%dplyr::summarise(count=n())
n5support$group="All"

ex4b=ggplot(n5support, aes(y=transcript_novelty, x=count, fill=n5_support)) + 
  scale_fill_manual(values=c("SCAFE"="#4DBBD5FF","SCAFE & GENCODE"="#E64B35FF"))+
  labs(fill=NULL, x=NULL, y =NULL, title="Support of Ex5_clusters from SALA Final")+
  facet_grid(rows=vars(group))+
  scale_x_continuous(labels = label_comma())+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85)+
  theme1+theme(legend.box.margin = margin(t=-10, r=0, b=0, l=-10),legend.position = "right", legend.direction = "vertical",
               plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"))

n3support=table5%>%group_by(transcript_novelty,n3_support,polyA)%>%dplyr::summarise(count=n())
n3support$n3_support=factor(n3support$n3_support, levels=c("no_support","n3_cluster","GENCODE & n3_cluster"))

ex4c=ggplot(n3support, aes(y=transcript_novelty, x=count, fill=n3_support)) + 
  scale_x_continuous(breaks=c(0,30000,60000), labels = label_comma())+
  scale_fill_manual(values=c("no_support"="white", "n3_cluster"="#4DBBD5FF","GENCODE & n3_cluster"="#E64B35FF"))+
  labs(fill=NULL, x=NULL, y =NULL, title="Support of Ex3_clusters from SALA Final")+
  facet_grid(rows=vars(polyA))+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85)+
  theme1+theme(legend.box.margin = margin(t=-10, r=0, b=0, l=-10),legend.position = "right", legend.direction = "vertical")

pdf(paste0(path_fig3,"ex4bc.end_support_transcript_base.pdf"), width = 2.8, height = 1.6)
grid.arrange(arrangeGrob(ex4b, ex4c, ncol=1, nrow = 2, heights = c(1.2,1.8)))
dev.off() 

#======================
#Needed
#ex4e
library(VennDiagram)
all=read.delim(paste0(path_fig3_data,"venn_ENST_SALA_TALON_Iso.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)

n12=length(intersect(all$geneID[which(all$analysis == "SALA")], all$geneID[which(all$analysis == "TALON")]))
n23=length(intersect(all$geneID[which(all$analysis == "TALON")], all$geneID[which(all$analysis == "IsoQuant")]))
n13=length(intersect(all$geneID[which(all$analysis == "SALA")], all$geneID[which(all$analysis == "IsoQuant")]))
n123=length(intersect(all$geneID[which(all$analysis == "TALON")], intersect(all$geneID[which(all$analysis == "SALA")], all$geneID[which(all$analysis == "IsoQuant")])))

grid.newpage()
ex4e=draw.triple.venn(length(which(all$analysis == "SALA")), length(which(all$analysis == "TALON")), length(which(all$analysis == "IsoQuant")), n12, n23, n13, n123, category=c("SALA","TALON","Isoquant"),
                    cex = 0.7, cat.cex = 0.75, cat.col = "black", lwd=c(0.25,0.25,0.25))
pdf(paste0(path_fig3,"ex4e.SALA_TALON_iso_ENST.venn.pdf"), width = 1.1, height = 1.1)
grid.draw(ex4e)
dev.off()

#==
missENST=read.delim(paste0(path_fig3_data,"venn_ENST_SALA_missed.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
missENST$structural_category=gsub("full-splice_match","FSM",missENST$structural_category)
missENST$structural_category=gsub("incomplete-splice_match","ISM",missENST$structural_category)
missENST$structural_category[which(is.na(missENST$structural_category))]="Null"
missENST1=missENST%>%group_by(group, structural_category)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
missENST1$group=factor(missENST1$group,levels=c("IsoQuant alone","IsoQuant & TALON","TALON alone"))

ex4e2=ggplot(missENST1, aes(x=group, y=percent, fill=structural_category)) + 
  #facet_grid(cols=vars(SALA_partial))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF","white"))+
  labs(fill=NULL, x=NULL, y ="% of transcript", title= "Associated model\nof absent ENST")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex4e2.SQANTI3.pdf"), width = 1, height = 1.3)
print(ex4e2)
dev.off() 

ex4e3=ggplot() + 
  labs(x = "Difference of TSS + TES  (nt)",title= NULL, y=NULL,color=NULL, fill=NULL) +
  geom_density(data=missENST[which(missENST$structural_category == "FSM"),], mapping=aes(x=AbsSumEnd, fill=structural_category), color="#E64B35FF", linewidth=0.25)+
  scale_fill_manual(values=c("#E64B35FF"))+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme1+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major = element_blank(),legend.position=c(0.8,0.5))
pdf(paste0(path_fig3,"ex4e3.distance.SQANTI3.pdf"), width = 1, height = 0.5)
print(ex4e3)
dev.off()

#======================
#Needed
#ex4f
SQANTI3=read.delim(paste0(path_fig3_data,"SQANTI3.txt.gz"), header=T, check.names = F, stringsAsFactors = F)
SQANTI3$structural_category[which(SQANTI3$structural_category=="full-splice_match")]="FSM"
SQANTI3$structural_category[which(SQANTI3$structural_category=="incomplete-splice_match")]="ISM"
SQANTI3$structural_category[which(SQANTI3$structural_category=="novel_in_catalog")]="NIC"
SQANTI3$structural_category[which(SQANTI3$structural_category=="novel_not_in_catalog")]="NNC"
colnames(SQANTI3)[c(5,6,7,9,11)]=c("IsoQuant_Sensitive", "IsoQuant_Default", "SALA_Default", "SALA_Final","TALON_Read_filtered")

SQANTI3a=reshape2::melt(SQANTI3[,c(1,9,7,11,5,6)], id=1)
SQANTI3a=SQANTI3a%>%group_by(variable)%>%dplyr::mutate(percent=value/sum(value))
SQANTI3a$structural_category=factor(SQANTI3a$structural_category,levels=c("FSM","ISM","NIC","NNC","genic","genic_intron","antisense","intergenic","fusion"))

ex4f=ggplot(SQANTI3a, aes(x=variable, y=percent, fill=structural_category)) + 
  scale_fill_npg()+
  labs(fill="Category", x=NULL, y ="% of transcript", title= "Category of novel transcripts")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex4f.SQANTI3.pdf"), width = 1.5, height = 1.7)
print(ex4f)
dev.off() 

#======================
#Needed
#ex4g
datab1=read.delim(paste0(path_fig3_data,"transcript.group.count.tsv.gz"), header=T, stringsAsFactors=F, check.names=F)
datab1=datab1[which(datab1$name %in% c("SALA_Final","SALA_Default","TALON_Read_filtered","IsoQuant_Sensitive","IsoQuant_Default")),]

datab1$analysis=sapply(strsplit(datab1$name,"_"),"[",1)
datab1$analysis=factor(datab1$analysis, levels=c("SALA","TALON","IsoQuant"))
datab1$name=factor(datab1$name, levels=c("SALA_Final","SALA_Default","TALON_Read_filtered","IsoQuant_Sensitive","IsoQuant_Default"))
datab1$V7[which(datab1$V7 == "Transcript_from_novel_gene")]="Transcript\nfrom novel TU"
datab1$V7[which(datab1$V7 == "ENST")]="GENCODE\ntranscript"
datab1$V7[which(datab1$V7 == "Novel_isoform")]="Novel isoform\nfrom ENSG"
ex4g=ggplot() + 
  labs(x=NULL, y ="Number of transcripts", title="Transcript discovered by different assemblers", fill="Assember")+
  geom_bar(data=datab1, mapping=aes(x=name, fill=analysis, y=count), color="black", stat = "identity", linewidth=0.25, width=0.8)+
  facet_wrap(vars(V7),ncol=3, scale="free")+
  scale_y_continuous(labels = label_number(suffix = "K", scale = 0.001))+
  scale_fill_npg()+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1))

pdf(paste0(path_fig3,"ex4g.transcript_discovery.pdf"), width = 3.6, height = 1.7)
print(ex4g)
dev.off() 

#======================
#Needed
#ex4h
CRE.hit.rate=read.delim(paste0(path_fig3_data,"transcript.hit.rate.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE.hit.rate$Others=CRE.hit.rate$novel.Tx.total-CRE.hit.rate$novel.Tx.hit
colnames(CRE.hit.rate)[6]="Supported"
CRE.hit.rate$group2="Raw"
CRE.hit.rate$group2[-grep("Raw",CRE.hit.rate$name)]="Filtered"
CRE.hit.rate=CRE.hit.rate[which(CRE.hit.rate$name %in% c("SALA_Final","SALA_Default","TALON_Read_filtered","IsoQuant_Sensitive","IsoQuant_Default")),]

CRE.hit.rate1=reshape2::melt(CRE.hit.rate[,c("group2","group","platform","name","Supported","Others")], id=c(1:4))
CRE.hit.rate1=CRE.hit.rate1%>%group_by(group2,group,platform,name)%>%dplyr::mutate(percent=value/sum(value))
CRE.hit.rate1$name=factor(CRE.hit.rate1$name, levels=c("SALA_Raw","SALA_Read_filtered","SALA_Final","SALA_Default","TALON_Raw","TALON_Read_filtered","IsoQuant_Sensitive","IsoQuant_Default"))
CRE.hit.rate1$variable=factor(CRE.hit.rate1$variable, levels=c("Others","Supported"))
CRE.hit.rate1$platform=factor(CRE.hit.rate1$platform, levels=c("SALA","TALON","IsoQuant"))
CRE.hit.rate1$group2=factor(CRE.hit.rate1$group2, levels=c("Raw","Filtered"))

ex4h=ggplot() + 
  labs(x=NULL, y ="% of transcripts", title="Transcript from novel transcriptional unit", fill="Analysis", alpha=NULL)+
  geom_bar(data=CRE.hit.rate1[which(CRE.hit.rate1$group %in% c("cCRE","ATAC","TSS_cluster")),], mapping=aes(x=name, fill=platform, alpha=variable,  y=percent), color="black", stat = "identity", linewidth=0.25, width=0.8)+
  facet_grid(cols = vars(group), scale="free", space="fixed")+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  scale_fill_npg(guide = NULL)+
  scale_alpha_manual(values=c(0,1))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1), legend.position = "bottom")

pdf(paste0(path_fig3,"ex4h.compare_SALA_n5.pdf"), width = 2.8, height = 2)
print(ex4h)
dev.off() 

#======================
#Needed
#for figure 1
CRE.hit.rate=read.delim(paste0(path_fig3_data,"transcript.hit.rate.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE.hit.rate=CRE.hit.rate[which(CRE.hit.rate$group == "ATAC"),]
CRE.hit.rate=CRE.hit.rate[which(CRE.hit.rate$name %in% c("SALA_Final","TALON_Read_filtered","IsoQuant_Default")),]
CRE.hit.rate$platform=factor(CRE.hit.rate$platform, levels=c("SALA","TALON","IsoQuant"))
path_fig1="/osc-fs_home/yip/CFC_seq_paper_fig_data/fig1/"
fig1_SALA=ggplot() + 
  labs(x="ATAC support", y ="# Novel transcript", title=NULL, color=NULL)+
  geom_point(data=CRE.hit.rate, mapping=aes(x=novel.rate, color=platform, y=novel.Tx.total),size=0.5)+
  scale_x_continuous(limits=c(0,1), breaks=c(0, 0.5,1))+
  scale_y_log10(limits=c(1, 1000000), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  scale_color_npg()+
  theme1+theme(panel.grid.major = element_blank(),legend.position = c(0.2,0.4))

pdf(paste0(path_fig1,"fig1_SALA.compare_SALA_n5.pdf"), width = 0.8, height = 1)
print(fig1_SALA)
dev.off() 

#======================
#Needed
#ex4i
sum1=read.delim(paste0(path_fig3_data,"RNA_class_summary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#generated from SALA_Final.R
sum1$count[which(sum1$class %in% c("GENCODE:others","GENCODE:protein_coding","Novel:potential_coding"))]=NA

ex4i=ggplot() + 
  scale_fill_npg()+
  labs(x=NULL, y =NULL, title="RNA class from SALA transcriptome", fill=NULL)+
  geom_bar(data=sum1, mapping=aes(x=group, fill=class, y=percent), linewidth=0.2, color="black", stat = "identity", alpha=0.7)+
  geom_text(data=sum1, mapping=aes(x=group, fill=class, y=percent, label=count), position = position_stack(vjust = 0.5), size=1.8)+
  scale_y_continuous(labels = scales::percent, limits=c(0,1.1), breaks=c(0,0.25,0.5,0.75,1))+  
  theme1+ theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex4i.RNA.class.pdf"), width = 1.8, height = 2)
print(ex4i)
dev.off() 

#======================
#Needed
#ex4j
sum3a=read.delim(paste0(path_fig3_data,"ncRNA.class.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

sum3a$T4_Novel_geneClass=gsub("GENCODE:lncRNA","GENCODE\ngene:lncRNA",sum3a$T4_Novel_geneClass)
sum3a$T4_Novel_geneClass=gsub("Novel:lncRNA","Novel gene:\nlncRNA",sum3a$T4_Novel_geneClass)
sum3a$T4_Novel_geneClass=gsub("Novel:short_ncRNA","Novel gene:\nshort_ncRNA",sum3a$T4_Novel_geneClass)

sum3a$label=sum3a$T_count
sum3a$label[which(sum3a$include == "No")]=NA

ex4j=ggplot() + 
  scale_fill_manual(values = c("Yes"="grey", "No"="white"))+
  labs(x=NULL, y ="Number of transcripts", title="Source of non-coding RNA transcript", fill="include")+
  geom_bar(data=sum3a, mapping=aes(x=Novel_transcriptClass, fill=include, y=T_count), linewidth=0.25, color="black", stat = "identity", alpha=1)+
  geom_text(data=sum3a, mapping=aes(x=Novel_transcriptClass, y=T_count, label=label), angle=45, hjust=0.2, vjust=-0.2, size=1.8)+
  scale_y_log10(limits=c(1, 1000000), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  facet_grid(cols=vars(T4_Novel_geneClass), scale="free", space="free")+
  theme1+
  theme(axis.text.x = element_text(color ="black", angle=45, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex4j.ncRNA.class.pdf"), width = 2.6, height = 2)
print(ex4j)
dev.off() 


