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

#===============================================================================
#bash
#Genome wide TATA box

setwd(TATA_path)
system("scanMotifGenomeWide.pl MA0108.3.motif hg38 -bed -keepAll -p 10 >hg38.TATAbox.bed")

#==============================
#R

setwd(TATA_path)
options(scipen=999)
tata=fread("hg38.TATAbox.bed", header=F, stringsAsFactors = F)
tata=tata[which(nchar(tata$V1)<=5),] #main chromaosome only
tata$V2=tata$V2-1
tata$V4=paste0("TBP",1:nrow(tata))
write.table(tata, gzfile("hg38.TATAbox.main.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

CREsummit=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz"), header=F, stringsAsFactors = F)
CREsummit$V2=CREsummit$V2-50
CREsummit$V3=CREsummit$V3+50
CREsummit$V2[which(CREsummit$V2 <0)]=0
write.table(CREsummit, gzfile("ontCAGE.Neuron_THP1.CRE.stranded_summit_50bp_extend.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#==============================
#bedtools
#bash
setwd(TATA_path)
system("rm hg38.TATAbox.bed")
system("bedtools intersect -wb -a ontCAGE.Neuron_THP1.CRE.stranded_summit_50bp_extend.bed.gz -b hg38.TATAbox.main.bed.gz -s | gzip > ontCAGE.Neuron_THP1.CRE.stranded_summit_50bb_extend.TBP.bed.gz")
system("rm hg38.TATAbox.main.bed.gz")

#==========================
#collapse intersect into metagene plot
#R
setwd(TATA_path)
TATA=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_50bb_extend.TBP.bed.gz", header=F, stringsAsFactors = F)
length(unique(TATA$V4))#26770
TATA=left_join(TATA, CREsummit[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
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
colnames(TATA)[c(4,5,6)]=c("CREID","name","motif_score")
write.table(TATA[order(TATA$locS),],gzfile("ontCAGE.Neuron_THP1.CRE.stranded_summit_50bb_extend.TBP.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)
TATA$locE=TATA$locS+4

TATA=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_50bb_extend.TBP.fakebed.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno=CREanno[which(CREanno$representative == "Yes"),] #take major strand alone

np=CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]
ne=CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]
nu=CREanno$CREID[which(CREanno$promoter_type == "unclassed")]
nc=CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]
np2D=CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="2D")]
np1D=CREanno$CREID[which(CREanno$promoter_type == "promoter-like" & CREanno$orientation =="1D")]
ne2D=CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="2D")]
ne1D=CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$orientation =="1D")]
neNAP=CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$Distance_PLScCRE >= 2000)]
neAP=CREanno$CREID[which(CREanno$promoter_type == "enhancer-like" & CREanno$Distance_PLScCRE < 2000)]

TATAp=TATA[which(TATA$CREID %in% np),]
TATAe=TATA[which(TATA$CREID %in% ne),]
TATAu=TATA[which(TATA$CREID %in% nu),]
TATAc=TATA[which(TATA$CREID %in% nc),]
TATAp2D=TATA[which(TATA$CREID %in% np2D),]
TATAp1D=TATA[which(TATA$CREID %in% np1D),]
TATAe2D=TATA[which(TATA$CREID %in% ne2D),]
TATAe1D=TATA[which(TATA$CREID %in% ne1D),]
TATAeNAP=TATA[which(TATA$CREID %in% neNAP),]
TATAeAP=TATA[which(TATA$CREID %in% neAP),]

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
  TATAeNAP1=TATAeNAP[which(TATAeNAP$motif_score >= need[i]),]
  TATAeAP1=TATAeAP[which(TATAeAP$motif_score >= need[i]),]
  #
  TATAc1=TATAc1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp1=TATAp1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp1D1=TATAp1D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp2D1=TATAp2D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe1=TATAe1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe1D1=TATAe1D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe2D1=TATAe2D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAu1=TATAu1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAeNAP1=TATAeNAP1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAeAP1=TATAeAP1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  
  path2=paste0(TATA_path,"n",need[i],"/")
  write.table(TATAc1[order(TATAc1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAc.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1[order(TATAp1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1[order(TATAe1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAu1[order(TATAu1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAu.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp2D1[order(TATAp2D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp2D.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1D1[order(TATAp1D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp1D.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe2D1[order(TATAe2D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe2D.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1D1[order(TATAe1D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe1D.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAeNAP1[order(TATAeNAP1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAeNAP.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAeAP1[order(TATAeAP1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAeAP.nooverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:100),c(1:101),c(-50:50))
write.table(fake2,gzfile("TATA.location101.fakebed.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#================================
#bash
#bedtools count
setwd(paste0(TATA_path,"n1/"))
system(paste0("for file in *nooverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",TATA_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(TATA_path,"n2/"))
system(paste0("for file in *nooverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",TATA_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(TATA_path,"n3/"))
system(paste0("for file in *nooverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",TATA_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(TATA_path,"n4/"))
system(paste0("for file in *nooverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",TATA_path,"TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

#===============================
#R
setwd(TATA_path)
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(neAP),length(neNAP),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=TATA_path, pattern=".nooverlap.result.bed.gz", recursive=T)
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(TATA_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(TATA_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("TATAe",data$group)]="enhancer-like"
data$anno_region[grep("TATAu",data$group)]="unclassed"
data$anno_region[grep("TATAc",data$group)]="CTCF-alone"
data$anno_region[grep("TATAeAP",data$group)]="enhancer-AP"
data$anno_region[grep("TATAeNAP",data$group)]="enhancer-NAP"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="TATA_box"
write.table(data, gzfile(paste0(path_fig2_data,"tCRE_TATA_nooverlap.result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#minorstrand
setwd(TATA_path)
system("bedtools intersect -wb -a ontCAGE.Neuron_THP1.CRE.minorstrand.complete.summit50.bed.gz -b hg38.TATAbox.main.bed.gz -s | gzip > ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_50bb_extend.TBP.bed.gz")
system("rm hg38.TATAbox.main.bed.gz")

#==========================
CREsummit=read.delim("ontCAGE.Neuron_THP1.CRE.minorstrand.complete.summit50.bed.gz", header=F, stringsAsFactors = F)
TATA=read.delim("ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_50bb_extend.TBP.bed.gz", header=F, stringsAsFactors = F)
length(unique(TATA$V4))#4753
TATA=left_join(TATA, CREsummit[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
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
colnames(TATA)[c(4,5,6)]=c("CREID","name","motif_score")
write.table(TATA[order(TATA$locS),],gzfile("ontCAGE.Neuron_THP1.CRE.minorstrand_summit_50bb_extend.TBP.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

CREanno=read.delim(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno1=CREanno[which(CREanno$CREID %in% CREsummit$V4),] #take minor strand alone

np=CREanno1$CREID[which(CREanno1$promoter_type == "promoter-like")]
ne=CREanno1$CREID[which(CREanno1$promoter_type == "enhancer-like")]
nu=CREanno1$CREID[which(CREanno1$promoter_type == "unclassed")]
nc=CREanno1$CREID[which(CREanno1$promoter_type == "CTCF-alone")]
np2D=CREanno1$CREID[which(CREanno1$promoter_type == "promoter-like" & CREanno1$orientation =="2D")]
np1D=CREanno1$CREID[which(CREanno1$promoter_type == "promoter-like" & CREanno1$orientation =="1D")]
ne2D=CREanno1$CREID[which(CREanno1$promoter_type == "enhancer-like" & CREanno1$orientation =="2D")]
ne1D=CREanno1$CREID[which(CREanno1$promoter_type == "enhancer-like" & CREanno1$orientation =="1D")]
neNAP=CREanno1$CREID[which(CREanno1$promoter_type == "enhancer-like" & CREanno1$Distance_PLScCRE >= 2000)]
neAP=CREanno1$CREID[which(CREanno1$promoter_type == "enhancer-like" & CREanno1$Distance_PLScCRE < 2000)]

TATAp=TATA[which(TATA$CREID %in% np),]
TATAe=TATA[which(TATA$CREID %in% ne),]
TATAu=TATA[which(TATA$CREID %in% nu),]
TATAc=TATA[which(TATA$CREID %in% nc),]
TATAp2D=TATA[which(TATA$CREID %in% np2D),]
TATAp1D=TATA[which(TATA$CREID %in% np1D),]
TATAe2D=TATA[which(TATA$CREID %in% ne2D),]
TATAe1D=TATA[which(TATA$CREID %in% ne1D),]
TATAeNAP=TATA[which(TATA$CREID %in% neNAP),]
TATAeAP=TATA[which(TATA$CREID %in% neAP),]

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
  TATAeNAP1=TATAeNAP[which(TATAeNAP$motif_score >= need[i]),]
  TATAeAP1=TATAeAP[which(TATAeAP$motif_score >= need[i]),]
  #
  TATAc1=TATAc1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp1=TATAp1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp1D1=TATAp1D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAp2D1=TATAp2D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe1=TATAe1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe1D1=TATAe1D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAe2D1=TATAe2D1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAu1=TATAu1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAeNAP1=TATAeNAP1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  TATAeAP1=TATAeAP1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
  
  path2=paste0(TATA_path,"n",need[i],"/")
  write.table(TATAc1[order(TATAc1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAc.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1[order(TATAp1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1[order(TATAe1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAu1[order(TATAu1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAu.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp2D1[order(TATAp2D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp2D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAp1D1[order(TATAp1D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAp1D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe2D1[order(TATAe2D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe2D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAe1D1[order(TATAe1D1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAe1D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAeNAP1[order(TATAeNAP1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAeNAP.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(TATAeAP1[order(TATAeAP1$locS),c(1:6)],gzfile(paste0(path2,"summit_50bp_extend.TATAeAP.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

#===============================================================================
#bash
#bedtools count
setwd(paste0(TATA_path,"n1/"))
system(paste0("for file in *.minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ", TATA_path, "TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))
setwd(paste0(TATA_path,"n2/"))
system(paste0("for file in *.minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ", TATA_path, "TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))
setwd(paste0(TATA_path,"n3/"))
system(paste0("for file in *.minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ", TATA_path, "TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))
setwd(paste0(TATA_path,"n4/"))
system(paste0("for file in *.minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ", TATA_path, "TATA.location101.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

#===============================================================================
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(neAP),length(neNAP),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=TATA_path, pattern=".minorstrand.result.bed.gz", recursive=T)

files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".minorstrand.result.bed.gz","",files.names)
data=read.delim(paste0(TATA_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(TATA_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("TATAe",data$group)]="enhancer-like"
data$anno_region[grep("TATAu",data$group)]="unclassed"
data$anno_region[grep("TATAc",data$group)]="CTCF-alone"
data$anno_region[grep("TATAeAP",data$group)]="enhancer-AP"
data$anno_region[grep("TATAeNAP",data$group)]="enhancer-NAP"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="TATA_box"
write.table(data, gzfile(paste0(path_fig2_data,"tCRE_TATA_minorstrand.result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)



