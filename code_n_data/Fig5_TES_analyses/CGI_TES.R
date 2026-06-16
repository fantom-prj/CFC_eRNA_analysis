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
library(readxl)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#####################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig5_data=paste0(primary_folder,"fig5/data/")
CGItes_path=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/")

#===================================================
#bash
setwd(CGItes_path)
cmd <- paste0(
  "wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz ",
  "| gunzip -c ",
  "| awk 'BEGIN{ OFS=\"\\t\" }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }' ",
  "| sort-bed - ",
  "> cpgIslandExt.hg38.bed")
system(cmd)

#===================================================
#R
setwd(CGItes_path)
cpg=read.delim("cpgIslandExt.hg38.bed", header=F, stringsAsFactors = F)
cpg=cpg[which(nchar(cpg$V1)<=5),]
write.table(cpg, gzfile("cpgIslandExt.hg38.main_chr.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===================================================
# use 3' end as reference -> position 0
setwd(CGItes_path)
options(scipen=999)
transcriptTES1=read.delim(paste0(CGItes_path,"Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz"),header=F, stringsAsFactors = F)

#==========================
#bash
#bedtools
setwd(CGItes_path)
system("bedtools intersect -wb -a Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz -b cpgIslandExt.hg38.main_chr.bed.gz | gzip > Neuron_THP1.S3.TES.table5.5000bp_extend.CpG.bed.gz")

#==========================
setwd(CGItes_path)
CPG=read.delim("Neuron_THP1.S3.TES.table5.5000bp_extend.CpG.bed.gz", header=F, stringsAsFactors = F)
length(unique(CPG$V4))#50066
CPG=left_join(CPG, transcriptTES1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
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
colnames(CPG)[c(4,5,6,7)]=c("TESID","name","cpgNum","obsExp")

table5b=read.delim("TESID_restricted_to_n3_string.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
table5a2=table5b[which(table5b$TEScount >= 3),]

CPG=left_join(CPG, table5a2[,c(1:4,6)],by="TESID",copy=F)
CPG=CPG[which(!is.na(CPG$ex5cluster_class)),]
write.table(CPG[order(CPG$locS),],gzfile("Neuron_THP1.S3.TES.table5.5000bp_extend.CpG.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

CPGmA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "mRNA")]
CPGmN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "mRNA")]
CPGeA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "e_ncRNA")]
CPGeN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA")]
CPGeN1=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "CGIap")]
CPGeN2=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "CGInap")]
CPGeN3=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "Null")]
CPGeN4=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "TATA")]
CPGpA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "p_ncRNA")]
CPGpN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "p_ncRNA")]
CPGpN1=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "p_ncRNA" & table5a2$CpGTATA == "CGI")]

need=c(20,40,60,80)
for (i in 1:length(need)){
  CPGmAa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGmA),]
  CPGmNa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGmN),]
  CPGeAa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeA),]
  CPGeNa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeN),]
  CPGeN1a=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeN1),]
  CPGeN2a=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeN2),]
  CPGeN3a=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeN3),]
  CPGeN4a=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGeN4),]
  CPGpAa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGpA),]
  CPGpNa=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGpN),]
  CPGpN1a=CPG[which(CPG$cpgNum >= need[i] & CPG$TESID %in% CPGpN1),]
  
  path2=paste0(CGItes_path,"n",need[i],"/")
  write.table(CPGmAa[order(CPGmAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGmA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGmNa[order(CPGmNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGmN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeAa[order(CPGeAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeNa[order(CPGeNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeN1a[order(CPGeN1a$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeN1.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeN2a[order(CPGeN2a$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeN2.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeN3a[order(CPGeN3a$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeN3.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeN4a[order(CPGeN4a$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGeN4.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGpAa[order(CPGpAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGpA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGpNa[order(CPGpNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGpN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGpN1a[order(CPGpN1a$locS),c(1:4)],gzfile(paste0(path2,"TES_5kb_extend.CPGpN1.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}
  

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile("CpG.location10001.fakebed.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)
#================================
#bedtools count
#bash
setwd(paste0(CGItes_path,"n20/"))
system(paste0("for file in TES_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ",CGItes_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGItes_path,"n40/"))
system(paste0("for file in TES_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ",CGItes_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGItes_path,"n60/"))
system(paste0("for file in TES_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ",CGItes_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGItes_path,"n80/"))
system(paste0("for file in TES_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ",CGItes_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

#===============================
#R
size1=c(length(CPGeA),length(CPGeN), length(CPGeN1), length(CPGeN2), length(CPGeN3), length(CPGeN4), length(CPGmA), length(CPGmN), length(CPGpA), length(CPGpN), length(CPGpN1))
size1=rep(size1,4)
files=list.files(path=CGItes_path, pattern=".result.bed.gz", recursive=T)
files=files[grep("TES_5kb_extend",files)]

files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(CGItes_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size1[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(CGItes_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size1[i]
  data=rbind(data,data1)}

data$anno_region="mRNA_pA"
data$anno_region[grep("CPGmN",data$group)]="mRNA_pN"
data$anno_region[grep("CPGeA",data$group)]="eRNA_pA"
data$anno_region[grep("CPGeN",data$group)]="eRNA_pN"
data$anno_region[grep("CPGeN1",data$group)]="eRNA_pN_CGIap"
data$anno_region[grep("CPGeN2",data$group)]="eRNA_pN_CGInap"
data$anno_region[grep("CPGeN3",data$group)]="eRNA_pN_Null"
data$anno_region[grep("CPGeN4",data$group)]="eRNA_pN_TATA"
data$anno_region[grep("CPGpA",data$group)]="p_ncRNA_pA"
data$anno_region[grep("CPGpN",data$group)]="p_ncRNA_pN"
data$anno_region[grep("CPGpN1",data$group)]="p_ncRNA_pN_CGI"

data$group2="CpG_island"

data$polyA="polyA"
data$polyA[grep("pN",data$anno_region)]="non-polyA"
data$ex5_cluster=sapply(strsplit(data$anno_region,"_p"),"[",1)
data$ex5_cluster[which(data$ex5_cluster=="eRNA")]="e_ncRNA"
write.table(data, gzfile(paste0(path_fig5_data,"TES_CpG_island_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# stored in [primary_folder]/fig5/data

