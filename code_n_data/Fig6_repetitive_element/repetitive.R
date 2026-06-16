library(dplyr) 
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

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
output_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/Tx_n_ex5_cluster/")
path_fig6_data=paste0(primary_folder,"fig6/data/")
MEME_input2=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/MEME_input2/")
MEME_output2=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/MEME_output2/")
NFY_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/NFY_ChIP/")
TEAD_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/TEAD_ChIP/")


#===============================================================================
#repeat element with n5_cluster and transcript model
options(scipen=999)

#=====================
#bash
#bedtools intersect with repeat, intersecting all table5, later filtered to our ex5_cluster alone from paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")
#repeatmasker downloaded from http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
setwd(output_folder)
system("bedtools intersect -wa -wb -a ex5_cluster_table5.bed.gz -b UCSC.hg38.repeat.bed.gz | gzip > ex5_cluster.rmsk.bed.gz")
system("bedtools intersect -wa -wb -a table5.final.detected.alone.bed6.bed.gz -b UCSC.hg38.repeat.bed.gz | gzip > exon.rmsk.bed.gz")

#=====================
#parse and collapse the intersect result
n5_repeat=read.delim(paste0(output_folder,"ex5_cluster.rmsk.bed.gz"), header=F, stringsAsFactors = F)
n5_repeat$n5_length=n5_repeat$V3-n5_repeat$V2
n5_repeat$repeat_length=n5_repeat$V9-n5_repeat$V8
n5_repeat$overlap_start=n5_repeat$V8
n5_repeat$overlap_start[which(n5_repeat$V8<n5_repeat$V2)]=n5_repeat$V2[which(n5_repeat$V8<n5_repeat$V2)]
n5_repeat$overlap_end=n5_repeat$V9
n5_repeat$overlap_end[which(n5_repeat$V9>n5_repeat$V3)]=n5_repeat$V3[which(n5_repeat$V9>n5_repeat$V3)]
n5_repeat$overlap_length=n5_repeat$overlap_end-n5_repeat$overlap_start
n5_repeat$overlap_percent=n5_repeat$overlap_length/n5_repeat$n5_length
k=n5_repeat%>%group_by(V14)%>%dplyr::summarise(count=n())
n5_repeat$group="others"
n5_repeat$group[which(n5_repeat$V14 %in% k$V14[which(k$count > 150)])]=n5_repeat$V14[which(n5_repeat$V14 %in% k$V14[which(k$count > 150)])]
write.table(n5_repeat, paste0(output_folder,"ex5_cluster.rmsk.bed.gz"), col.names=T, row.names=F, sep="\t", quote=F)

exon_repeat=read.delim(paste0(output_folder,"exon.rmsk.bed.gz"), header=F, stringsAsFactors = F)
exon_repeat$n5_length=exon_repeat$V3-exon_repeat$V2
exon_repeat$repeat_length=exon_repeat$V9-exon_repeat$V8
exon_repeat$overlap_start=exon_repeat$V8
exon_repeat$overlap_start[which(exon_repeat$V8<exon_repeat$V2)]=exon_repeat$V2[which(exon_repeat$V8<exon_repeat$V2)]
exon_repeat$overlap_end=exon_repeat$V9
exon_repeat$overlap_end[which(exon_repeat$V9>exon_repeat$V3)]=exon_repeat$V3[which(exon_repeat$V9>exon_repeat$V3)]
exon_repeat$overlap_length=exon_repeat$overlap_end-exon_repeat$overlap_start
exon_repeat$overlap_percent=exon_repeat$overlap_length/exon_repeat$n5_length
k=exon_repeat%>%group_by(V14)%>%dplyr::summarise(count=n())
exon_repeat$group="others"
exon_repeat$group[which(exon_repeat$V14 %in% c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat"))]=exon_repeat$V14[which(exon_repeat$V14 %in% c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat"))]
write.table(exon_repeat, paste0(output_folder,"exon.rmsk.bed.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#collapse the ex5_cluster RE into per ex5_cluster base, with overlap length filter
n5_repeat=n5_repeat[which(n5_repeat$overlap_length >=6),]
n5_repeat=n5_repeat%>%group_by(V4)%>%dplyr::mutate(RE_group_6nt=paste(unique(group), collapse=";"),RE_subgroup_6nt=paste(unique(V15), collapse=";"))
n5_repeat1=n5_repeat%>%group_by(V4)%>%dplyr::slice_max(overlap_percent,with_ties = FALSE)
n5_repeat1=n5_repeat1[,c(4,13:15,22,16,20,21,23,24)]
write.table(n5_repeat1, gzfile(paste0(output_folder,"ex5_cluster.rmsk.collapse_n5cluster.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#collapse the exon RE into per ex5_cluster base, with overlap length filter
exon_repeat1=exon_repeat%>%group_by(V4,V13,V14,V15,group)%>%dplyr::summarise(hit_exon=n(), overlap_length=sum(overlap_length))
exon_repeat1=exon_repeat1[which(exon_repeat1$overlap_length >= 200),]
exon_repeat1=right_join(table5[,c(1,62,57)], exon_repeat1, by=c("model_ID"="V4"),copy=F)
exon_repeat1=exon_repeat1%>%group_by(n5_string, group, V15)%>%dplyr::slice_max(overlap_length,with_ties = FALSE)
exon_repeat1=exon_repeat1%>%group_by(n5_string)%>%dplyr::mutate(RE_group_200nt=paste(unique(group), collapse=";"),RE_subgroup_200nt=paste(unique(V15), collapse=";"))
exon_repeat1=exon_repeat1%>%group_by(n5_string)%>%dplyr::slice_max(overlap_length,with_ties = FALSE)
write.table(exon_repeat1, gzfile(paste0(output_folder,"exon.rmsk.bed.collapse_n5cluster.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#all files located in [primary_folder]/code_n_data/Fig6_repetitive_element/Tx_n_ex5_cluster
#===============================================================================


#===============================================================================
#only for exon, find coverage
#% of exon coverage, number of exon covered.
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
exon_repeat1=read.delim(paste0(output_folder,"exon.rmsk.bed.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
exon_repeat1=inner_join(exon_repeat1, data1[,c(1,9,85)], by=c("n5_string"), copy=F)
exon_repeat1$overlap_percent=exon_repeat1$overlap_length/exon_repeat1$transcript_length
exon_repeat13=exon_repeat1[which(!is.na(exon_repeat1$ex5cluster_class)),]

write.table(exon_repeat13,gzfile(paste0(path_fig6_data,"RE_ex5_cluster.exon_coverage_nexon.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#file located in [primary_folder]/fig6/data
#for -> fig ex8i

#===============================================================================
####prepare fisher exact
n5_repeat1=read.delim(paste0(output_folder,"ex5_cluster.rmsk.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5_repeat1$value="Yes"
n5_repeat2=spread(n5_repeat1[,c(1,5,11)], key=2, value=3)
n5_repeat2[is.na(n5_repeat2)]="No"

path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$value="Yes"
data2=spread(data1[,c("ex5cluster_class","n5_string","CpGTATA","value")], key=3, value=4)
data2=left_join(data2, n5_repeat2, by=c("n5_string"="V4"), copy=F)
data2[is.na(data2)]="No"
data2=data2[,-c(7,13)]
data5=data.frame(matrix(nrow=0, ncol=5))
colnames(data5)=c("gene_group","feature1","feature2","count","label")

for (i in 1:5){
  for (j in 1:7){
    data3=data2[,c(1,2+i,7+j)]
    colnames(data3)=c("gene_group","feature1", "feature2")
    data4=data3%>%group_by(gene_group,feature1,feature2)%>%dplyr::summarise(count=n())
    data4$label=paste0(data4$feature1,"_",data4$feature2)
    data4$feature1=colnames(data2[2+i])
    data4$feature2=colnames(data2[7+j])
    data5=rbind(data5,data4)}}
data6=spread(data5, key=5, value=4)
data6[is.na(data6)]=0

data6=data6[-which(data6$gene_group != "e_ncRNA" & data6$feature1 == "CGIap"),]
data6=data6[-which(data6$gene_group != "e_ncRNA" & data6$feature1 == "CGInap"),]
data6=data6[-which(data6$gene_group == "e_ncRNA" & data6$feature1 == "CGI"),]

for(i in 1:nrow(data6)){
  GSEATasting <- matrix(c(data6$Yes_Yes[i], data6$Yes_No[i], data6$No_Yes[i], data6$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  data6$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  data6$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
data6$pv=paste0("p = ",signif(data6$p.val,3))
data6$sig_level="ns"
data6$sig_level[which(data6$p.val<0.05)]="*"
data6$sig_level[which(data6$p.val<0.01)]="**"
data6$sig_level[which(data6$p.val<0.001)]="***"
data6$FE_logOR=log(data6$OR)

#
exon_repeat1=read.delim(paste0(output_folder,"exon.rmsk.bed.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
exon_repeat1$value="Yes"
exon_repeat2=spread(exon_repeat1[,c(2,7,12)], key=2, value=3)
exon_repeat2[is.na(exon_repeat2)]="No"

path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$value="Yes"
data2=spread(data1[,c("ex5cluster_class","n5_string","CpGTATA","value")], key=3, value=4)
data2=left_join(data2, exon_repeat2, by="n5_string", copy=F)
data2[is.na(data2)]="No"
data2=data2[,-c(7,13)]
data5=data.frame(matrix(nrow=0, ncol=5))
colnames(data5)=c("gene_group","feature1","feature2","count","label")
for (i in 1:5){
  for (j in 1:7){
    data3=data2[,c(1,2+i,7+j)]
    colnames(data3)=c("gene_group","feature1", "feature2")
    data4=data3%>%group_by(gene_group,feature1,feature2)%>%dplyr::summarise(count=n())
    data4$label=paste0(data4$feature1,"_",data4$feature2)
    data4$feature1=colnames(data2[2+i])
    data4$feature2=colnames(data2[7+j])
    data5=rbind(data5,data4)}}
data7=spread(data5, key=5, value=4)
data7[is.na(data7)]=0

data7=data7[-which(data7$gene_group != "e_ncRNA" & data7$feature1 == "CGIap"),]
data7=data7[-which(data7$gene_group != "e_ncRNA" & data7$feature1 == "CGInap"),]
data7=data7[-which(data7$gene_group == "e_ncRNA" & data7$feature1 == "CGI"),]

for(i in 1:nrow(data7)){
  GSEATasting <- matrix(c(data7$Yes_Yes[i], data7$Yes_No[i], data7$No_Yes[i], data7$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  data7$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  data7$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
data7$pv=paste0("p = ",signif(data7$p.val,3))
data7$sig_level="ns"
data7$sig_level[which(data7$p.val<0.05)]="*"
data7$sig_level[which(data7$p.val<0.01)]="**"
data7$sig_level[which(data7$p.val<0.001)]="***"
data7$FE_logOR=log(data7$OR)

data6$group="Ex5_cluster"
data7$group="Exon"
data8=rbind(data6, data7)
write.table(data8, gzfile(paste0(path_fig6_data,"CpGTATA_repeat_element_FE.n5_exon.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#file located in [primary_folder]/fig6/data
#for -> fig 6e
#===============================================================================


#===============================================================================
#actual occurence
n5_repeat1=read.delim(paste0(output_folder,"ex5_cluster.rmsk.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
exon_repeat1=read.delim(paste0(output_folder,"exon.rmsk.bed.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2a=left_join(data1[,c("ex5cluster_class","n5_string","CpGTATA")], n5_repeat1[,c("V4","group")], by=c("n5_string"="V4"), copy=F)
data2a=left_join(data2a, exon_repeat1[,c("n5_string","group")], by="n5_string", copy=F, suffix=c("_ex5cluster","_exon"))
data2a[is.na(data2a)]="Null"
write.table(data2a, gzfile(paste0(path_fig6_data,"TE_group_CpGTATA.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data
#for -> fig ex9a
#fig ex9a shows distribution in e_ncRNA, additional plots below show p_ncRNA, other_ncRNA, mRNA & CTCF_ncRNA

data3a=data2a%>%group_by(ex5cluster_class,CpGTATA,group_ex5cluster)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="ex5_cluster")
data4a=data2a%>%group_by(ex5cluster_class,CpGTATA,group_exon)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="exon")
colnames(data3a)[3]="REgroup"
colnames(data4a)[3]="REgroup"
data5a=rbind(data3a, data4a)
data5a$CpGTATA=factor(data5a$CpGTATA, levels=c("CGI","Null","TATA","Others"))
data5a$REgroup=factor(data5a$REgroup, levels=c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat", "Null", "others"))
data6a=unique(data5a[,c(1,2,6,7)])

ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Presence of repeat elements inside ex5_cluster \nor transcript model of p_ncRNA",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("p_ncRNA") & data5a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=percent, fill=REgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("p_ncRNA") & data6a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Presence of repeat elements inside ex5_cluster \nor transcript model of other_ncRNA",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("other_ncRNA") & data5a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=percent, fill=REgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("other_ncRNA") & data6a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Presence of repeat elements inside ex5_cluster \nor transcript model of mRNA",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("mRNA") & data5a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=percent, fill=REgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("mRNA")  & data6a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Presence of repeat elements inside ex5_cluster \nor transcript model of mRNA",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("CTCF_ncRNA") & data5a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=percent, fill=REgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("CTCF_ncRNA")  & data6a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))


#===============================================================================
#actual occurrence for family
n5_repeat1=read.delim(paste0(output_folder,"ex5_cluster.rmsk.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
exon_repeat1=read.delim(paste0(output_folder,"exon.rmsk.bed.collapse_n5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data2a=left_join(data1[,c("ex5cluster_class","n5_string","CpGTATA")], n5_repeat1[,c(1,4,5)], by=c("n5_string"="V4"), copy=F)
data2a=left_join(data2a, exon_repeat1[,c(2,6,7)], by="n5_string", copy=F, suffix=c("_ex5cluster","_exon"))
data2a[is.na(data2a)]="Null"
write.table(data2a,gzfile(paste0(path_fig6_data,"family_of_LTR.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data
#for -> fig ex9b

#===============================================================================
#heatmap on ex5_cluster TPM to show cell type specificity of enhancer , divided into dCGI, CGI, Null and TATA
#count matrix from SCAFE counting CTSS inside genuine TSS cluster per ex5_cluster
#counting refer to [primary_folder]/code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count_n5cluster/00_count.sh
#this contain all ex5, filter later for visualization
count_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count_n5cluster/output/count_matrix/")
count_ex5=read.delim(paste0(count_path,"ontCAGE.Neuron_THP1.count.txt"), header=T, stringsAsFactors = F, check.names = F)
count_ex5$iPSC=rowSums(count_ex5[,c(22:23)])
count_ex5$NSC=rowSums(count_ex5[,c(2:3)])
count_ex5$Neuron=rowSums(count_ex5[,c(4:5)])
count_ex5$THP1=rowSums(count_ex5[,c(6:13)])
count_ex5$dTHP1=rowSums(count_ex5[,c(14:21)])
rownames(count_ex5)=count_ex5$CREID
count_ex51=count_ex5[,c(24:28)]

d <- DGEList(counts=count_ex51)
RLE <- calcNormFactors(d, method="RLE")
RLE.ex5=cpm(RLE, normalized.lib.sizes=TRUE)
write.table(RLE.ex5,paste0(path_fig6_data,"ex5_cluster_Neuron_THP1.RLE.tsv.gz"), col.names=T, row.names=T, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data
#for -> fig 6g & ex9e
#===============================================================================



#===============================================================================
#prepare data for MEME motif enrichment
path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2a=read.delim(paste0(path_fig6_data,"TE_group_CpGTATA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

# add repeat element to data1
data1=left_join(data1, data2a[,c("n5_string","group_ex5cluster","group_exon")], by="n5_string", copy=F)
colnames(data1)[c(129,130)]=c("repeat_ex5cluster","repeat_exon")
write.table(data1,gzfile(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

data1$TFstart=data1$summit_start-400
data1$TFstart[which(data1$strand == "-")]=data1$summit_start[which(data1$strand == "-")]-100
data1$TFend=data1$summit_end+100
data1$TFend[which(data1$strand == "-")]=data1$summit_end[which(data1$strand == "-")]+400

library(GenomicRanges)
gr <- GRanges(seqnames = data1$chr, ranges = IRanges(start = data1$TFstart+1, end = data1$TFend))
reduced_gr <- GenomicRanges::reduce(gr, min.gapwidth=1)
overlap_hits <- findOverlaps(gr, reduced_gr)
data1$TF501_ID = subjectHits(overlap_hits)

data1$CpGTATA[which(data1$CpGTATA == "Null")]=NA
data2a=data1%>%group_by(TF501_ID)%>%dplyr::summarise(CpGTATA=paste(unique(na.omit(CpGTATA)),collapse=";"),
                                                     promoter_type=paste(unique(promoter_type), collapse=";"))

data3=data2a[-grep(";",data2a$promoter_type),]
data3=data3[-grep(";",data3$CpGTATA),]
data3=data3[which(data3$CpGTATA != "Others"),]

data1a=data1[which(data1$TF501_ID %in% unique(data3$TF501_ID)),]
data1b=data1a%>%group_by(TF501_ID)%>%dplyr::slice_max(count) #take highest count
data1b=data1b%>%group_by(TF501_ID)%>%dplyr::slice_min(ex5_length) # take the shorter length of ex5
data1c=data1b%>%group_by(TF501_ID)%>%dplyr::slice_sample(n=1) #take random one if more than one per TF501_ID

data1b=data1b[which(data1b$ex5cluster_class %in% c("p_ncRNA","e_ncRNA","other_ncRNA","mRNA")),]
data1b$CpGTATA[which(is.na(data1b$CpGTATA))]="Null"
write.table(data1b[,c("n5_string","chr","TFstart", "TFend", "TF501_ID")], gzfile(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

ak=data1b%>%group_by(ex5cluster_class, CpGTATA)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

f1=c("p_ncRNA","e_ncRNA","other_ncRNA","mRNA")
f2=c("dCGI","uCGI","Null","TATA")
for (i in 1:nrow(ak)){
    data1c=data1b[which(data1b$ex5cluster_class == ak$ex5cluster_class[i] & data1b$CpGTATA == ak$CpGTATA[i]),c("chr","TFstart","TFend","n5_string","count","strand")]
    write.table(data1c[order(data1c$chr,data1c$TFstart),],gzfile(paste0(MEME_input2,ak$ex5cluster_class[i],"_",ak$CpGTATA[i],".bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)}

data1c=data1b[which(data1b$CpGTATA %in% c("CGIap","CGInap") & data1b$ex5cluster_class == "e_ncRNA"),c("chr","TFstart","TFend","n5_string","count","strand")]
write.table(data1c[order(data1c$chr,data1c$TFstart),],gzfile(paste0(MEME_input2,"e_ncRNA_CGI.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

eTATA=read.delim(paste0(MEME_input2,"e_ncRNA_TATA.bed.gz"), header=F)
eTATA1=eTATA[which(eTATA$V4 %in% data1$n5_string[which(data1$repeat_ex5cluster == "LTR")]),]
write.table(eTATA1,gzfile(paste0(MEME_input2,"e_ncRNA_TATALTR.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

uTATA=read.delim(paste0(MEME_input2,"other_ncRNA_TATA.bed.gz"), header=F)
uTATA1=uTATA[which(uTATA$V4 %in% data1$n5_string[which(data1$repeat_ex5cluster == "LTR")]),]
write.table(uTATA1,gzfile(paste0(MEME_input2,"other_ncRNA_TATALTR.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#bed files located in [primary_folder]/code_n_data/Fig6_repetitive_element/MEME_input2


#===============================================================================
#bash
#MEME
setwd(MEME_input2)
system("for file in *.bed.gz; do bedtools getfasta -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed \"$file\" -fo \"${file%.bed.gz}.fa\"; done")

code1 <- paste(
"for file in *.fa; do sea",
"--p \"$file\"",
"--m JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt",
"--thresh 10000000",
"--text",
"--noseqs >  \"../MEME_output2/${file%.fa}.MEME.tsv\"; done")
system(code1)

#MEME results located in [primary_folder]/code_n_data/Fig6_repetitive_element/MEME_output2

#===============================================================================
#parse MEME result
files=list.files(path=MEME_output2, pattern=".MEME.tsv")
files.names=gsub(".MEME.tsv","",files)
meme0=data.frame()
for (i in 1:length(files)){
  meme=read.delim(paste0(MEME_output2,files[i]), header=T, stringsAsFactors = F, check.names = F, nrow=879)
  meme$logQval=-meme$LOG_QVALUE
  meme$Set=files.names[i]
  meme0=rbind(meme0,meme)
}
write.table(meme0[,c("ID","ALT_ID","ENR_RATIO","logQval","Set")],gzfile(paste0(MEME_output2,"CFC.SEAall.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)
meme=read.delim(paste0(MEME_output2,"CFC.SEAall.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
meme1=meme[-grep("mRNA",meme$Set),]
need=unique(meme1$ID[which(meme1$logQval>2 & meme1$ENR_RATIO > 12)])
meme1a=meme1[which(meme1$ID %in% need),]
meme1b=spread(meme1a[,c(1,2,5,3)], key=3, value=4)
meme1b$geneName=toupper(meme1b$ALT_ID)

#====
#filter for expressed TF gene, use CAGE data
cage250=read.delim(paste0(MEME_input2,"count_matrix_no_rRNA.join_tech_rep.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cage250$iPSC=rowSums(cage250[,c(25,26)])
cage250$NSC=rowSums(cage250[,c(29,30)])
cage250$Neuron=rowSums(cage250[,c(27,28)])
cage250$THP1_0=rowSums(cage250[,c(31,32)])
cage250$THP1_24=rowSums(cage250[,c(33,34)])
cage250$THP1_96=rowSums(cage250[,c(35,36)])
d <- DGEList(counts=cage250[c(37:42)])
RLE <- calcNormFactors(d, method="RLE")
RLE.cage250=data.frame(cpm(RLE, normalized.lib.sizes=TRUE))
RLE.cage250$geneID=rownames(RLE.cage250)
RLE.cage250$MaxCPM=apply(RLE.cage250[,c(1:6)], 1, max)
#==

gene_list=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript_to_gene.tsv"), header=F, stringsAsFactors = F)
gene_list=unique(gene_list[,c(2,5)])
meme1b=left_join(meme1b, gene_list, by=c("geneName"="V5"), copy=F)
meme1b=left_join(meme1b, RLE.cage250[,c(7,8)], by=c("V2"="geneID"), copy=F)
meme1b=meme1b[which(meme1b$MaxCPM>1),] #at least one cell type with CPM>1
meme1b=meme1b[which(meme1b$ID != "MA1537.2"),] #remove redundant name with different motif ID
rownames(meme1b)=meme1b$ALT_ID
meme1b=meme1b[,c(2:15)]
write.table(meme1b, gzfile(paste0(path_fig6_data,"meme_result.tsv.gz")), row.names=T, col.names=T, sep="\t", quote=F)

#==
#count the number of ex5_cluster included
content=data.frame(matrix(nrow=16, ncol=2))
colnames(content)=c("set","count")

files=list.files(path=MEME_input2, pattern=".bed.gz")
files.names=gsub(".bed.gz","",files)
for (i in 1:length(files)){
  bed=read.delim(paste0(MEME_input2,files[i]), header=F)
  content$set[i]=files.names[i]
  content$count[i]=nrow(bed)}
write.table(content,paste0(path_fig6_data,"meme_result.count.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
#files located in [primary_folder]/fig6/data

#===============================================================================
# parse the NFYB binding results
# ENCSR146UIC is a ChIP-seq of NFYB on WTC11 iPSC line
# the conservative IDR thresholded peaks were used
# overlapping performed in ./NFY_ex5.R
NFY_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/NFY_ChIP/")
NFY=read.delim(paste0(NFY_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.NFY.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

# -400 to 100 around summit is considered binding positive
NFY1=NFY[-which(NFY$locE<4500 | NFY$locS>5100),]
data1$NFY_bind_iPSC="No"
data1$NFY_bind_iPSC[which(data1$n5_string %in% unique(NFY1$n5_string))]="Yes"

#===============================================================================
# parse the TEAD4 binding results
# FANTOM6 ChIP-seq of TEAD4 on WTC11 iPSC line 
# the merged peaks were used
# overlapping performed in ./TEAD_ex5.R

TEAD=read.delim(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD_merge.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

# -400 to 100 around summit is considered binding positive
TEAD1=TEAD[-which(TEAD$locE<4500 | TEAD$locS>5100),]
data1$TEAD_bind_iPSC="No"
data1$TEAD_bind_iPSC[which(data1$n5_string %in% unique(TEAD1$n5_string))]="Yes"

write.table(data1,paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# subset for non-overlap analysis
subset1=read.delim(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=data1[which(data1$n5_string %in% subset1$n5_string),]
data1$CpGTATA[grep("CGI",data1$CpGTATA)]="CGI"

data1$LTR="No"
data1$LTR[which(data1$repeat_ex5cluster == "LTR")]="Yes"

data3=data1[which(data1$iPSC>0),]%>%group_by(ex5cluster_class,CpGTATA,LTR,NFY_bind_iPSC)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
data3$label=paste0(data3$LTR,"_",data3$NFY_bind_iPSC)
data3a=spread(data3[,c(1,2,7,5)],key=3,value=4)
data3a[is.na(data3a)]=0

data3t=data1[which(data1$iPSC>0),]%>%group_by(ex5cluster_class,CpGTATA,LTR,TEAD_bind_iPSC)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
data3t$label=paste0(data3t$LTR,"_",data3t$TEAD_bind_iPSC)
data3b=spread(data3t[,c(1,2,7,5)],key=3,value=4)
data3b[is.na(data3b)]=0

data3a$group="NF-Y"
data3b$group="TEAD"
data3a=rbind(data3a,data3b)

for(i in 1:nrow(data3a)){
  GSEATasting <- matrix(c(data3a$Yes_Yes[i], data3a$Yes_No[i], data3a$No_Yes[i], data3a$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  data3a$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  data3a$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
data3a$pv=paste0("p = ",signif(data3a$p.val,3))
data3a$sig_level="ns"
data3a$sig_level[which(data3a$p.val<0.05)]="*"
data3a$sig_level[which(data3a$p.val<0.01)]="**"
data3a$sig_level[which(data3a$p.val<0.001)]="***"
data3a$FE_logOR=log(data3a$OR)
data3a$comparison="LTR_TFbind"
write.table(data3a,gzfile(paste0(path_fig6_data,"NFY_TEAD_LTR.n5_FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#file located in [primary_folder]/fig6/data
#for -> fig ex9h

#===================================
# enrichment between TF binding and H3K27 chromatin state
data1$Un_marked="No"
data1$Un_marked[which(data1$K27Ac_iPS == "No" & data1$K27ME3_iPS == "No")]="Yes"
data1$Active="No"
data1$Active[which(data1$K27Ac_iPS == "Yes" & data1$K27ME3_iPS == "No")]="Yes"
data1$Co_marked="No"
data1$Co_marked[which(data1$K27Ac_iPS == "Yes" & data1$K27ME3_iPS == "Yes")]="Yes"
data1$Repressed="No"
data1$Repressed[which(data1$K27Ac_iPS == "No" & data1$K27ME3_iPS == "Yes")]="Yes"

#NFY
data4a=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Un_marked,NFY_bind_iPSC)%>%dplyr::summarise(count=n(), group="Un-marked")
data4b=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Active,NFY_bind_iPSC)%>%dplyr::summarise(count=n(), group="Active")
data4c=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Repressed,NFY_bind_iPSC)%>%dplyr::summarise(count=n(), group="Repressed")
data4d=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Co_marked,NFY_bind_iPSC)%>%dplyr::summarise(count=n(), group="Co-marked")

colnames(data4a)[3]="K27"
colnames(data4b)=colnames(data4a)
colnames(data4c)=colnames(data4a)
colnames(data4d)=colnames(data4a)

data4=rbind(data4a,data4b,data4c,data4d)
data4$label=paste0(data4$K27,"_",data4$NFY_bind_iPSC)
data4a=spread(data4[,c(1,2,6,7,5)],key=4,value=5)
data4a[is.na(data4a)]=0

#TEAD
data5a=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Un_marked,TEAD_bind_iPSC)%>%dplyr::summarise(count=n(), group="Un-marked")
data5b=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Active,TEAD_bind_iPSC)%>%dplyr::summarise(count=n(), group="Active")
data5c=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Repressed,TEAD_bind_iPSC)%>%dplyr::summarise(count=n(), group="Repressed")
data5d=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$iPSC>0 & data1$repeat_ex5cluster == "LTR"),]%>%group_by(ex5cluster_class,CpGTATA,Co_marked,TEAD_bind_iPSC)%>%dplyr::summarise(count=n(), group="Co-marked")

colnames(data5a)[3]="K27"
colnames(data5b)=colnames(data5a)
colnames(data5c)=colnames(data5a)
colnames(data5d)=colnames(data5a)

data5=rbind(data5a,data5b,data5c,data5d)
data5$label=paste0(data5$K27,"_",data5$TEAD_bind_iPSC)
data5a=spread(data5[,c(1,2,6,7,5)],key=4,value=5)
data5a[is.na(data5a)]=0

data4a$group2="NF-Y"
data5a$group2="TEAD"
data4a=rbind(data4a,data5a)

for(i in 1:nrow(data4a)){
  GSEATasting <- matrix(c(data4a$Yes_Yes[i], data4a$Yes_No[i], data4a$No_Yes[i], data4a$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  data4a$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  data4a$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
data4a$pv=paste0("p = ",signif(data4a$p.val,3))
data4a$sig_level="ns"
data4a$sig_level[which(data4a$p.val<0.05)]="*"
data4a$sig_level[which(data4a$p.val<0.01)]="**"
data4a$sig_level[which(data4a$p.val<0.001)]="***"
data4a$FE_logOR=log(data4a$OR)
data4a$comparison="K27_TFbind"

write.table(data4a,gzfile(paste0(path_fig6_data,"NFY_TEAD_K27_FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#file located in [primary_folder]/fig6/data
#for -> fig ex9i

#===============================================================================
RLE.ex5=read.delim(paste0(path_fig6_data,"ex5_cluster_Neuron_THP1.RLE.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
RLE.ex5$n5_string=rownames(RLE.ex5)
data5=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$CpGTATA != "CGI" & data1$repeat_ex5cluster == "LTR"),]
data5=left_join(data5, RLE.ex5[,c(1,6)], by="n5_string", suffix=c("_count","_CPM"), copy=F)
data6=reshape2::melt(data5[,c(1,9,85,138,131,132,134:137)], id=1:6)
data6=data6[which(data6$value == "Yes"),]
data6a=data6[which(data6$iPSC > 0),]%>%group_by(ex5cluster_class,CpGTATA,variable,NFY_bind_iPSC, TEAD_bind_iPSC)%>%dplyr::summarise(CPM=median(iPSC_CPM), count=n())
data6b=data6[which(data6$iPSC ==0),]%>%group_by(ex5cluster_class,CpGTATA,NFY_bind_iPSC, TEAD_bind_iPSC)%>%dplyr::summarise(CPM=median(iPSC_CPM), count=n(), variable="non_iPSC")
data6c=data6[which(data6$iPSC >0),]%>%group_by(ex5cluster_class,CpGTATA,NFY_bind_iPSC, TEAD_bind_iPSC)%>%dplyr::summarise(CPM=median(iPSC_CPM), count=n(), variable="all_iPSC")
data7=bind_rows(data6a,data6b,data6c)
data7$label="None"
data7$label[which(data7$NFY_bind_iPSC=="Yes")]="NF-Y"
data7$label[which(data7$TEAD_bind_iPSC=="Yes")]="TEAD"
data7$label[which(data7$TEAD_bind_iPSC=="Yes" & data7$NFY_bind_iPSC=="Yes")]="Both"
data7=data7%>%group_by(ex5cluster_class,CpGTATA,variable)%>%mutate(percent=paste0(signif(count/sum(count)*100,2),"%"))

write.table(data7, gzfile(paste0(path_fig6_data,"LTR_NFY_TEAD_positive_number.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#file located in [primary_folder]/fig6/data
#for -> fig 6h

#===============================================================================
#ATAC count for non-redundant ex5_cluster 501 bp
#take the scATAC count as bulk ATAC count
#source scATAC data can be found from DDBJ: DRA019608 (DRR618513- DRR618515)
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)
library(chromVAR)
subset1=read.delim(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
write.table(subset1[order(subset1$chr,subset1$TFstart),c(2,3,4,1)], gzfile(paste0(NFY_folder,"ex5_cluster_used_motif.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

peak_ex5cluster501 <- getPeaks(paste0(NFY_folder,"ex5_cluster_used_motif.bed.gz"), sort_peaks = TRUE)
IDList <- c("iPS","NSC","NRN")
md.List <- list()
for(i in IDList){
  tmp <- read.table(
    file = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/singlecell.csv"), 
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1)[-1, ]
  md.List[[i]] <- tmp}

for(i in IDList){md.List[[i]] <- md.List[[i]][md.List[[i]]$passed_filters > 500, ]}

## create fragment objects
frags.List <- list()
for(i in IDList){
  tmp <- CreateFragmentObject(
    path = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/fragments.tsv.gz"),
    cells = rownames(md.List[[i]]))
  frags.List[[i]] <- tmp}
plan("multicore", workers = 6)
options(future.globals.maxSize = 50 * 1024^3) 
counts.List <- list()
for(i in IDList){
  tmp <-  FeatureMatrix(
    fragments = frags.List[[i]],
    features = peak_ex5cluster501,
    cells = rownames(md.List[[i]]))
  counts.List[[i]] <- tmp}

c1=as.data.frame(rowSums(as.matrix(counts.List[["iPS"]])))
c2=as.data.frame(rowSums(as.matrix(counts.List[["NSC"]])))
c3=as.data.frame(rowSums(as.matrix(counts.List[["NRN"]])))
c1$peakID=rownames(c1)
c2$peakID=rownames(c2)
c3$peakID=rownames(c3)
c0=full_join(c1,c2, by="peakID",copy=F)
c0=full_join(c0,c3, by="peakID",copy=F)
colnames(c0)[c(1,3,4)]=c("ATACcount_iPSC","ATACcount_NSC","ATACcount_Neuron")
c0$chr=sapply(strsplit(c0$peakID,"-"),"[",1)
c0$start=as.numeric(sapply(strsplit(c0$peakID,"-"),"[",2))-1
c0$end=as.numeric(sapply(strsplit(c0$peakID,"-"),"[",3))
c0=left_join(c0, subset1[,c(1:4)],by=c("chr"="chr","start"="TFstart","end"="TFend"),copy=F)

write.table(c0[,c(5:8,1,3,4)],gzfile(paste0(NFY_folder,"ex5_cluster501_ATAC_celltype.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#get CPM
TF_ATAC=read.delim(paste0(NFY_folder,"ex5_cluster501_ATAC_celltype.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

rownames(TF_ATAC)=TF_ATAC$n5_string
TF_ATAC1=TF_ATAC[,-c(1:4)]
cpm_TF_ATAC = cpm(TF_ATAC1)
write.table(cpm_TF_ATAC,gzfile(paste0(NFY_folder,"ex5_cluster501_ATAC_celltype.tsv.gz")), row.names=T, col.names=T, sep="\t", quote=F)

#===============================================================================
#extract ATAC signal into needed ex5_cluster
cpm_TF_ATAC=read.delim(paste0(NFY_folder,"ex5_cluster501_ATAC_celltype.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpm_TF_ATAC$n5_string=rownames(cpm_TF_ATAC)

subset1=read.delim(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=data1[which(data1$n5_string %in% subset1$n5_string),]
data1$CpGTATA[grep("CGI",data1$CpGTATA)]="CGI"

data1$LTR="No"
data1$LTR[which(data1$repeat_ex5cluster == "LTR")]="Yes"
data1$iPSC_t="No"
data1$iPSC_t[which(data1$iPSC>0)]="Yes"
data1$TFbind="Others"
data1$TFbind[which(data1$NFY_bind_iPSC == "Yes")]="NF-Y"
data1$TFbind[which(data1$TEAD_bind_iPSC == "Yes")]="TEAD"
data1$TFbind[which(data1$NFY_bind_iPSC == "Yes"& data1$TEAD_bind_iPSC == "Yes")]="Both"
data5=data1[which(data1$ex5cluster_class %in% c("e_ncRNA","other_ncRNA") & data1$CpGTATA != "CGI" & data1$repeat_ex5cluster == "LTR"),]
data5=left_join(data5, cpm_TF_ATAC, by="n5_string", copy=F)

data5=data5[,c("n5_string","ex5cluster_class","CpGTATA","iPSC_t","TFbind","ATACcount_iPSC")]
colnames(data5)[6]="ATACcpm_iPSC"
write.table(data5,gzfile(paste0(path_fig6_data,"LTR_ex5_cluster_ATAC_CPM.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data
#for -> fig 6i

#===============================================================================
# prepare supp Table S19
data1=read.delim(paste0(primary_folder,"fig4/data/features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5ee=table5[which(table5$ex5cluster_class == "e_ncRNA"),]
table5ee=table5ee[which(table5ee$n5_string %in% data1$n5_string),]
table5ee1=unique(table5ee[,c("n5_string","T4_gene_ID","T4_gene_name")])
length(unique(table5ee1$n5_string))#28354
length(unique(table5ee1$T4_gene_ID))#23854
table5ee1=left_join(table5ee1,data1[,c("n5_string","SE_source","CpGTATA","repeat_ex5cluster","repeat_exon","NFY_bind_iPSC","TEAD_bind_iPSC")], by="n5_string", copy=F)

gene4=read.delim(paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/all.match.genetogene_short_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(gene4)=c("CFCseq_geneID","FANTOM_CAT_geneID","LncBook_geneID","Refseq_geneID","GENCODEv47_geneID")
ref=separate_rows(gene4[which(!is.na(gene4$Refseq_geneID)),c(1,4)],Refseq_geneID, sep=";")
ref=ref[which(ref$CFCseq_geneID %in% table5ee1$T4_gene_ID),]

v47=separate_rows(gene4[which(!is.na(gene4$GENCODEv47_geneID)),c(1,5)],GENCODEv47_geneID, sep=";")
gencode.info=fread(paste0("/analysisdata/fantom6/Interactome/resources/gencode_v47/transcript_to_gene.tsv"), header=F)
gencode.info=unique(gencode.info[,c(2,5)])
v47=left_join(v47, gencode.info, by=c("GENCODEv47_geneID"="V2"),copy=F)
v47=v47[which(v47$CFCseq_geneID %in% table5ee1$T4_gene_ID),]

EVL=read.delim(paste0(primary_folder,"code_n_data/Fig6_repetitive_element/EVLncRNAs3/EVLncRNAs3_function_human.txt"),header=T, stringsAsFactors = F, check.names=F)
EVL1=EVL[which(EVL$`LncRNA name` %in% table5ee1$T4_gene_name),]
EVL1=left_join(EVL1, unique(table5ee1[,c(2,3)]), by=c("LncRNA name"="T4_gene_name"), copy=F)
EVL1a=EVL1[,c(26,2,6,8,23)]%>%group_by(T4_gene_ID)%>%dplyr::summarise(Link="Direct",`LncRNA name`=paste(unique(`LncRNA name`),collapse=";"),Function=paste(unique(`Molecular fucntions`),collapse=";"), Disease=paste(unique(Disease),collapse=";"),PMID=paste(unique(PMID),collapse=";"))

EVL2=EVL[which(EVL$`LncRNA name` %in% ref$Refseq_geneID),]
EVL2=left_join(EVL2, ref, by=c("LncRNA name"="Refseq_geneID"), copy=F)
EVL2a=EVL2[,c(26,2,6,8,23)]%>%group_by(CFCseq_geneID)%>%dplyr::summarise(Link="Thu_Refseq_gene",`LncRNA name`=paste(unique(`LncRNA name`),collapse=";"),Function=paste(unique(`Molecular fucntions`),collapse=";"), Disease=paste(unique(Disease),collapse=";"),PMID=paste(unique(PMID),collapse=";"))

EVL3=EVL[which(EVL$`LncRNA name` %in% v47$V5),]
EVL3=EVL3[-which(EVL3$`LncRNA name` %in% unique(EVL1$`LncRNA name`)),]
EVL3=left_join(EVL3, v47[,c(1,3)], by=c("LncRNA name"="V5"), copy=F)
EVL3=EVL3[-which(EVL3$CFCseq_geneID %in% unique(EVL2$CFCseq_geneID)),]
EVL3a=EVL3[,c(26,2,6,8,23)]%>%group_by(CFCseq_geneID)%>%dplyr::summarise(Link="Thu_GCv47_gene",`LncRNA name`=paste(unique(`LncRNA name`),collapse=";"),Function=paste(unique(`Molecular fucntions`),collapse=";"), Disease=paste(unique(Disease),collapse=";"),PMID=paste(unique(PMID),collapse=";"))

colnames(EVL2a)=colnames(EVL1a)
colnames(EVL3a)=colnames(EVL1a)
EVL1a=rbind(EVL1a,EVL2a,EVL3a)
table5ee1=left_join(table5ee1, EVL1a, by="T4_gene_ID",copy=F)

write.table(table5ee1,paste0(primary_folder,"supplementary_table/TableS18_eRNA_features_EVLncRNA.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

