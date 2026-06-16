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
INR_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/INR/")
INRex5_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/INR/")

#===========================
#Genome wide initiator
#bash
setwd(INR_path)
system("scanMotifGenomeWide.pl INR_oldJaspar_10100.motif hg38 -bed -keepAll -p 10 >hg38.INR.bed")
# the motif file need to be generated, all separators are tab 

#===============================================
cluster=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.bed.bgz"), header=F, stringsAsFactors = F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
cluster=cluster[which(cluster$V4 %in% table5$n5_string),]
write.table(cluster[order(cluster$V1,cluster$V7),c(1,7,8,4,5,6)], paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), row.names=F, col.names=F, sep="\t", quote=F)
cluster_summit=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), header=F, stringsAsFactors = F)
cluster_summit1=cluster_summit
cluster_summit1$V2=cluster_summit1$V2-50
cluster_summit1$V3=cluster_summit1$V3+50
cluster_summit1$V2[which(cluster_summit1$V2 <0)]=0
write.table(cluster_summit1[order(cluster_summit1$V1,cluster_summit1$V2),], gzfile(paste0(INRex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================
#bedtools
#bash
setwd(INRex5_path)
system("sh bed_intersect2.sh")  
system("rm hg38.INR.bed")
system("rm part_*")

#==========================
#R
options(scipen=999)
cluster_summit1=read.delim(paste0(INRex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.bed.gz"), header=F, stringsAsFactors = F)

files=list.files(path=INRex5_path, pattern="bed.gz")
files=files[1:10]

INR=read.delim(paste0(INRex5_path,files[1]), header=F, stringsAsFactors = F)
INR=left_join(INR, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
INR=INR[which(INR$V11>=48),]
for (i in 2:length(files)){
  INR1=read.delim(paste0(INRex5_path,files[i]), header=F, stringsAsFactors = F)
  INR1=left_join(INR1, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
  INR1=INR1[which(INR1$V11>=48),]
  INR=rbind(INR,INR1)}

length(unique(INR$V4))#71012
INR$locS=INR$V2-INR$V2_ori-50
INR$locE=INR$V3-INR$V2_ori-50
INR1=INR[which(INR$V6 == "+" & INR$V2==INR$V8),c(1,14,15,4,10,11)]
INR2=INR[which(INR$V6 == "-" & INR$V3==INR$V9),c(1,15,14,4,10,11)]
INR2$locE=INR2$locE * (-1)
INR2$locS=INR2$locS * (-1)
INR2$V1="chr1"
INR1$V1="chr1"
colnames(INR2)=colnames(INR1)
INR1$locS=INR1$locS+50
INR1$locE=INR1$locE+50
INR2$locS=INR2$locS+51
INR2$locE=INR2$locE+51
INR=rbind(INR1,INR2)

colnames(INR)[c(4,5,6)]=c("n5_string","name","motif_score")
INR=left_join(INR, unique(table5[,c("n5_string","CREID")]),by="n5_string",copy=F)
write.table(INR[order(INR$locS),],gzfile(paste0(INRex5_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.INR.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
INR$locE=INR$locS+2

CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
np=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")])]
ne=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")])]
nu=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")])]
nc=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")])]
np2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")])]
np1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")])]
ne2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")])]
ne1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")])]

INRp=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]),]
INRe=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]),]
INRu=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")]),]
INRc=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]),]
INRp2D=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")]),]
INRp1D=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")]),]
INRe2D=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")]),]
INRe1D=INR[which(INR$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")]),]
need=c(48.5,49,49.5,50)
for (i in 1:length(need)){
  INRc1=INRc[which(INRc$motif_score >= need[i]),]
  INRp1=INRp[which(INRp$motif_score >= need[i]),]
  INRp1D1=INRp1D[which(INRp1D$motif_score >= need[i]),]
  INRp2D1=INRp2D[which(INRp2D$motif_score >= need[i]),]
  INRe1=INRe[which(INRe$motif_score >= need[i]),]
  INRe1D1=INRe1D[which(INRe1D$motif_score >= need[i]),]
  INRe2D1=INRe2D[which(INRe2D$motif_score >= need[i]),]
  INRu1=INRu[which(INRu$motif_score >= need[i]),]
  #
  INRc1=INRc1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRp1=INRp1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRp1D1=INRp1D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRp2D1=INRp2D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRe1=INRe1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRe1D1=INRe1D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRe2D1=INRe2D1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  INRu1=INRu1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
  path2=paste0(INRex5_path,"n",i,"/")
  write.table(INRc1[order(INRc1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRp1[order(INRp1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRe1[order(INRe1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRu1[order(INRu1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRp2D1[order(INRp2D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRp2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRp1D1[order(INRp1D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRp1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRe2D1[order(INRe2D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRe2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(INRe1D1[order(INRe1D1$locS),c(1:6)],gzfile(paste0(path2,"end5.summit_50bp_extend.INRe1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:100),c(1:101),c(-50:50))
write.table(fake2,gzfile(paste0(INRex5_path,"INR.location101.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#================================
#bedtools count
#bash
setwd(paste0(INRex5_path,"n1/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",INRex5_path,"INR.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(INRex5_path,"n2/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",INRex5_path,"INR.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(INRex5_path,"n3/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",INRex5_path,"INR.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(INRex5_path,"n4/"))
system(paste0("for file in end5.summit_50bp_extend*; do bedtools intersect -wa -a ",INRex5_path,"INR.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

#===============================
#R
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=INRex5_path, pattern=".result.bed.gz", recursive=T)
files=files[grep("end5",files)]
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(INRex5_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(INRex5_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("INRe",data$group)]="enhancer-like"
data$anno_region[grep("INRu",data$group)]="unclassed"
data$anno_region[grep("INRc",data$group)]="CTCF-alone"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="INR"
write.table(data, gzfile(paste0(path_fig2_data,"end5_cluster_INR_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

