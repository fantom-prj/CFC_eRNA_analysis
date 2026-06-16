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
CGI_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/CGI/")
CGIex5_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/CGI/")

#===================================================
#bash
setwd(CGIex5_path)
cmd <- paste0(
  "wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz ",
  "| gunzip -c ",
  "| awk 'BEGIN{ OFS=\"\\t\" }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }' ",
  "| sort-bed - ",
  "> cpgIslandExt.hg38.bed")
system(cmd)

#===================================================
#R
setwd(CGIex5_path)
cpg=read.delim("cpgIslandExt.hg38.bed", header=F, stringsAsFactors = F)
cpg=cpg[which(nchar(cpg$V1)<=5),]
write.table(cpg, gzfile("cpgIslandExt.hg38.main_chr.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===================================================
#n5 cluster base (from sala)
setwd(CGIex5_path)
cluster_summit1=read.delim("Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz", header=F, stringsAsFactors = F)

#===================================================
#bedtools
#bash
setwd(CGI_path)
system("rm cpgIslandExt.hg38.bed")
system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed .gz-b cpgIslandExt.hg38.main_chr.bed.gz | gzip > Neuron_THP1.S3.end5.summit.table5.5000bp_extend.CpG.bed.gz")

#==========================
#R
setwd(CGI_path)
options(scipen=999)
CPG=read.delim("Neuron_THP1.S3.end5.summit.table5.5000bp_extend.CpG.bed.gz", header=F, stringsAsFactors = F)
length(unique(CPG$V4))#46992
CPG=left_join(CPG, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
CPG$locS=CPG$V2-CPG$V2_ori-5000
CPG$locE=CPG$V3-CPG$V2_ori-5000
CPG1=CPG[which(CPG$V6 == "+"),c(1,18,19,4,10,12,16)]
CPG2=CPG[which(CPG$V6 == "-"),c(1,19,18,4,10,12,16)]
CPG2$locE=CPG2$locE * (-1)
CPG2$locS=CPG2$locS * (-1)
CPG2$V1="chr1"
CPG1$V1="chr1"
colnames(CPG2)=colnames(CPG1)
CPG1$locS=CPG1$locS+5000
CPG1$locE=CPG1$locE+5000
CPG2$locS=CPG2$locS+5001
CPG2$locE=CPG2$locE+5001
CPG=rbind(CPG1,CPG2)
colnames(CPG)[c(4,5,6,7)]=c("n5_string","name","cpgNum","obsExp")
CPG=left_join(CPG, unique(table5[,c("n5_string","CREID")]),by="n5_string",copy=F)
write.table(CPG[order(CPG$locS),],gzfile("Neuron_THP1.S3.end5.summit.table5.5000bp_extend.CpG.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
np=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")])]
ne=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")])]
nu=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")])]
nc=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")])]
np2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")])]
np1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")])]
ne2D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")])]
ne1D=cluster_summit1$V4[which(cluster_summit1$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")])]

CPGp=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]),]
CPGe=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]),]
CPGu=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")]),]
CPGc=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]),]
CPGp2D=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")]),]
CPGp1D=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")]),]
CPGe2D=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")]),]
CPGe1D=CPG[which(CPG$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")]),]

need=c(20,40,60,80)
for (i in 1:length(need)){
  CPGc1=CPGc[which(CPGc$cpgNum >= need[i]),]
  CPGp1=CPGp[which(CPGp$cpgNum >= need[i]),]
  CPGp1D1=CPGp1D[which(CPGp1D$cpgNum >= need[i]),]
  CPGp2D1=CPGp2D[which(CPGp2D$cpgNum >= need[i]),]
  CPGe1=CPGe[which(CPGe$cpgNum >= need[i]),]
  CPGe1D1=CPGe1D[which(CPGe1D$cpgNum >= need[i]),]
  CPGe2D1=CPGe2D[which(CPGe2D$cpgNum >= need[i]),]
  CPGu1=CPGu[which(CPGu$cpgNum >= need[i]),]
  path2=paste0(CGIex5_path,"n",need[i],"/")
  write.table(CPGc1[order(CPGc1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGc.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1[order(CPGp1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGp.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1[order(CPGe1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGe.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGu1[order(CPGu1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGu.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp2D1[order(CPGp2D1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGp2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1D1[order(CPGp1D1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGp1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe2D1[order(CPGe2D1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGe2D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1D1[order(CPGe1D1$locS),c(1:4)],gzfile(paste0(path2,"end5.summit_5kb_extend.CpGe1D.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile("CpG.location10001.fakebed.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#================================
#bedtools count
#bash
setwd(paste0(CGIex5_path,"n20/"))
system(paste0("for file in *end5.summit_5kb_extend; do bedtools intersect -wa -a ",CGIex5_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(CGIex5_path,"n40/"))
system(paste0("for file in *end5.summit_5kb_extend; do bedtools intersect -wa -a ",CGIex5_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(CGIex5_path,"n60/"))
system(paste0("for file in *end5.summit_5kb_extend; do bedtools intersect -wa -a ",CGIex5_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))

setwd(paste0(CGIex5_path,"n80/"))
system(paste0("for file in *end5.summit_5kb_extend; do bedtools intersect -wa -a ",CGIex5_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed}.result.bed.gz\"; done"))


#===============================
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=CGIex5_path, pattern=".result.bed", recursive=T)
files=files[grep("end5",files)]
files=files[-grep("_ud",files)]
files=files[-grep("_d",files)]
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed","",files.names)
data=read.delim(paste0(CGIex5_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(CGIex5_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("CpGe",data$group)]="enhancer-like"
data$anno_region[grep("CpGu",data$group)]="unclassed"
data$anno_region[grep("CpGc",data$group)]="CTCF-alone"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="CpG_island"
write.table(data, gzfile(paste0(path_fig2_data,"end5_cluster_CpG_island_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

