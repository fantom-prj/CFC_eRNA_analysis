#output stored in [primary_folder]/fig2/data

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
path_fig2_data=paste0(primary_folder,"fig2/data/")
TATA_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/TATAbox/")
TATAex5_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/TATAbox/")

#===============================================================================
#bash
#Genome wide TATA box

setwd(TATA_path)
system("scanMotifGenomeWide.pl MA0108.3.motif hg38 -bed -keepAll -p 10 >hg38.TATAbox.bed")

#==============================
#R

options(scipen=999)
tata=fread("hg38.TATAbox.bed", header=F, stringsAsFactors = F)
tata=tata[which(nchar(tata$V1)<=5),] #main chromaosome only
tata$V2=tata$V2-1
tata$V4=paste0("TBP",1:nrow(tata))
write.table(tata, gzfile("hg38.TATAbox.main.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

cluster=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.bed.bgz"), header=F, stringsAsFactors = F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
cluster=cluster[which(cluster$V4 %in% table5$n5_string),]
write.table(cluster[order(cluster$V1,cluster$V7),c(1,7,8,4,5,6)], paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), row.names=F, col.names=F, sep="\t", quote=F)
cluster_summit=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), header=F, stringsAsFactors = F)
cluster_summit1=cluster_summit
cluster_summit1$V2=cluster_summit1$V2-50
cluster_summit1$V3=cluster_summit1$V3+50
cluster_summit1$V2[which(cluster_summit1$V2 <0)]=0
write.table(cluster_summit1[order(cluster_summit1$V1,cluster_summit1$V2),], gzfile(paste0(TATAex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#=================================
#bedtools
#bash
setwd(TATA_path)
system("rm hg38.TATAbox.bed")
setwd(TATAex5_path)
system(paste0("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.50bp_extend.bed.gz -b ",TATA_path,"hg38.TATAbox.main.bed.gz -s | gzip > Neuron_THP1.S3.end5.summit.table5.50bp_extend.TBP.bed.gz"))
system("rm hg38.TATAbox.main.bed.gz")

#=================================
#collapse intersect into metagene plot
#R

cluster_summit1=read.delim(paste0(TATAex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.bed.gz"), header=F, stringsAsFactors = F)

TATA=read.delim(paste0(TATAex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.TBP.bed.gz"), header=F, stringsAsFactors = F)
length(unique(TATA$V4))#25943
TATA=left_join(TATA, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
TATA$locS=TATA$V2-TATA$V2_ori-50
TATA$locE=TATA$V3-TATA$V2_ori-50
TATA1=TATA[which(TATA$V6 == "+" & TATA$V2==TATA$V8),c(1,14,15,4,10,11)]
TATA2=TATA[which(TATA$V6 == "-" & TATA$V3==TATA$V9),c(1,15,14,4,10,11)]
TATA2$locE=TATA2$locE * (-1)
TATA2$locS=TATA2$locS * (-1)
TATA2$V1="chr1"
TATA1$V1="chr1"
colnames(TATA2)=colnames(TATA1)
TATA1$locS=TATA1$locS+50
TATA1$locE=TATA1$locE+50
TATA2$locS=TATA2$locS+51
TATA2$locE=TATA2$locE+51
TATA=rbind(TATA1,TATA2)
colnames(TATA)[c(4,5,6)]=c("n5_string","name","motif_score")
TATA=left_join(TATA, unique(table5[,c("n5_string","CREID")]),by="n5_string",copy=F)
write.table(TATA[order(TATA$locS),],paste0(TATAex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.TBP.fakebed.tsv"), row.names=F, col.names=T, sep="\t", quote=F)

TATA$locE=TATA$locS+4

CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
np=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")])]
ne=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")])]
nu=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")])]
nc=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")])]
np2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")])]
np1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")])]
ne2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")])]
ne1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")])]

TATAp=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]),]
TATAe=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]),]
TATAu=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")]),]
TATAc=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]),]
TATAp2D=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")]),]
TATAp1D=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")]),]
TATAe2D=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")]),]
TATAe1D=TATA[which(TATA$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")]),]
need=c(1,2,3,4)
for (i in 1:length(need)){
  TATAc1=TATAc[which(TATAc$motif_score >= need[i]),]
  TATAp1=TATAp[which(TATAp$motif_score >= need[i]),]
  TATAp1D1=TATAp1D[which(TATAp1D$motif_score >= need[i]),]
  TATAp2D1=TATAp2D[which(TATAp2D$motif_score >= need[i]),]
  TATAe1=TATAe[which(TATAe$motif_score >= need[i]),]
  TATAe1D1=TATAe1D[which(TATAe1D$motif_score >= need[i]),]
  TATAe2D1=TATAe2D[which(TATAe2D$motif_score >= need[i]),]
  TATAu1=TATAu[which(TATAu$motif_score >= need[i]),]
  #
  TATAc1=TATAc1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAp1=TATAp1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAp1D1=TATAp1D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAp2D1=TATAp2D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAe1=TATAe1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAe1D1=TATAe1D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAe2D1=TATAe2D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  TATAu1=TATAu1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  path2=paste0(TATAex5_path,"n",need[i],"/")
  write.table(TATAc1[order(TATAc1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1[order(TATAp1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1[order(TATAe1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAu1[order(TATAu1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp2D1[order(TATAp2D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAp2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1D1[order(TATAp1D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAp1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe2D1[order(TATAe2D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAe2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1D1[order(TATAe1D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.TATAe1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:100),c(1:101),c(-50:50))
write.table(fake2,gzfile(paste0(TATAex5_path,"TATA.location101.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#================================
#bash
#bedtools count
setwd(paste0(TATAex5_path,"n1/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",TATAex5_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(TATAex5_path,"n2/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",TATAex5_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(TATAex5_path,"n3/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",TATAex5_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(TATAex5_path,"n4/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",TATAex5_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

#===============================
#R
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=TATAex5_path, pattern=".result.bed.gz", recursive=T)
files=files[grep("end5",files)]
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(TATAex5_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(path1, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("TATAe",data$group)]="enhancer-like"
data$anno_region[grep("TATAu",data$group)]="unclassed"
data$anno_region[grep("TATAc",data$group)]="CTCF-alone"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="TATA_box"
write.table(data, gzfile(paste0(path_fig2_data,"end5_cluster_TATA_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)




