#output stored in [primary_folder]/fig2/data

#53-196: major strand summit
#197-318: minor strand summit

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

#===================================================
#bash
setwd(CGI_path)
cmd <- paste0(
  "wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz ",
  "| gunzip -c ",
  "| awk 'BEGIN{ OFS=\"\\t\" }{ print $2, $3, $4, $5$6, $7, $8, $9, $10, $11, $12 }' ",
  "| sort-bed - ",
  "> cpgIslandExt.hg38.bed")
system(cmd)

#===================================================
#R
setwd(CGI_path)
cpg=read.delim("cpgIslandExt.hg38.bed.gz", header=F, stringsAsFactors = F)
cpg=cpg[which(nchar(cpg$V1)<=5),]
write.table(cpg, gzfile("cpgIslandExt.hg38.main_chr.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===================================================
setwd(CGI_path)
CREsummit=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz"), header=F, stringsAsFactors = F)
CREsummit$V2=CREsummit$V2-5000
CREsummit$V3=CREsummit$V3+5000
CREsummit$V2[which(CREsummit$V2 <0)]=0
write.table(CREsummit, gzfile("ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===================================================
#bedtools
#bash
setwd(CGI_path)
system("rm cpgIslandExt.hg38.bed")
system("bedtools intersect -wb -a ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.bed.gz -b cpgIslandExt.hg38.main_chr.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.CpG.bed.gz")

#===================================================
#major strand summit
#R
setwd(CGI_path)
options(scipen=999)
CREsummit10001=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.bed.gz", header=F, stringsAsFactors = F)
CPG=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.CpG.bed.gz", header=F, stringsAsFactors = F)
length(unique(CPG$V4))#41856
CPG=left_join(CPG, CREsummit10001[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
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
colnames(CPG)[c(4,5,6,7)]=c("CREID","name","cpgNum","obsExp")
write.table(CPG[order(CPG$locS),],gzfile("ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.CpG.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

CPG=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.CpG.fakebed.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

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

CPGp=CPG[which(CPG$CREID %in% np),]
CPGe=CPG[which(CPG$CREID %in% ne),]
CPGu=CPG[which(CPG$CREID %in% nu),]
CPGc=CPG[which(CPG$CREID %in% nc),]
CPGp2D=CPG[which(CPG$CREID %in% np2D),]
CPGp1D=CPG[which(CPG$CREID %in% np1D),]
CPGe2D=CPG[which(CPG$CREID %in% ne2D),]
CPGe1D=CPG[which(CPG$CREID %in% ne1D),]
CPGeNAP=CPG[which(CPG$CREID %in% neNAP),]
CPGeAP=CPG[which(CPG$CREID %in% neAP),]

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
  CPGeNAP1=CPGeNAP[which(CPGeNAP$cpgNum >= need[i]),]
  CPGeAP1=CPGeAP[which(CPGeAP$cpgNum >= need[i]),]
  
  path2=paste0(CGI_path,"n",need[i],"/")
  write.table(CPGc1[order(CPGc1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGc.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1[order(CPGp1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1[order(CPGe1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGu1[order(CPGu1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGu.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp2D1[order(CPGp2D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp2D.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1D1[order(CPGp1D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp1D.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe2D1[order(CPGe2D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe2D.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1D1[order(CPGe1D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe1D.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeNAP1[order(CPGeNAP1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGeNAP.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeAP1[order(CPGeAP1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGeAP.noverlap.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile("CpG.location10001.fakebed.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)
#================================
#bedtools count
#bash
setwd(paste0(CGI_path,"n20/"))
system(paste0("for file in *noverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n40/"))
system(paste0("for file in *noverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n60/"))
system(paste0("for file in *noverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n80/"))
system(paste0("for file in *noverlap.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

#===============================
#R
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(neAP),length(neNAP),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=CGI_path, pattern="noverlap.result.bed.gz", recursive=T)
files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(CGI_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(CGI_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("CpGe",data$group)]="enhancer-like"
data$anno_region[grep("CpGu",data$group)]="unclassed"
data$anno_region[grep("CpGc",data$group)]="CTCF-alone"
data$anno_region[grep("CpGeAP",data$group)]="enhancer-AP"
data$anno_region[grep("CpGeNAP",data$group)]="enhancer-NAP"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="CpG_island"
write.table(data, gzfile(paste0(path_fig2_data,"tCRE_CpG_island_nooverlap.result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#minor strand
setwd(CGI_path)
system("bedtools intersect -wb -a ontCAGE.Neuron_THP1.CRE.minorstrand.complete.summit5000.bed.gz -b cpgIslandExt.hg38.main_chr.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_5kb_extend.CpG.bed.gz")
CREsummit=read.delim("ontCAGE.Neuron_THP1.CRE.minorstrand.complete.summit5000.bed.gz", header=F, stringsAsFactors = F)

CPG=read.delim("ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_5kb_extend.CpG.bed.gz", header=F, stringsAsFactors = F)
length(unique(CPG$V4))#11627
CPG=left_join(CPG, CREsummit[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
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
colnames(CPG)[c(4,5,6,7)]=c("CREID","name","cpgNum","obsExp")
write.table(CPG[order(CPG$locS),],gzfile("ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_5kb_extend.CpG.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

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

CPGp=CPG[which(CPG$CREID %in% np),]
CPGe=CPG[which(CPG$CREID %in% ne),]
CPGu=CPG[which(CPG$CREID %in% nu),]
CPGc=CPG[which(CPG$CREID %in% nc),]
CPGp2D=CPG[which(CPG$CREID %in% np2D),]
CPGp1D=CPG[which(CPG$CREID %in% np1D),]
CPGe2D=CPG[which(CPG$CREID %in% ne2D),]
CPGe1D=CPG[which(CPG$CREID %in% ne1D),]
CPGeNAP=CPG[which(CPG$CREID %in% neNAP),]
CPGeAP=CPG[which(CPG$CREID %in% neAP),]

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
  CPGeNAP1=CPGeNAP[which(CPGeNAP$cpgNum >= need[i]),]
  CPGeAP1=CPGeAP[which(CPGeAP$cpgNum >= need[i]),]
  
  path2=paste0(CGI_path,"n",need[i],"/")
  
  write.table(CPGc1[order(CPGc1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGc.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1[order(CPGp1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1[order(CPGe1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGu1[order(CPGu1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGu.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp2D1[order(CPGp2D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp2D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGp1D1[order(CPGp1D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGp1D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe2D1[order(CPGe2D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe2D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGe1D1[order(CPGe1D1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGe1D.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeNAP1[order(CPGeNAP1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGeNAP.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(CPGeAP1[order(CPGeAP1$locS),c(1:4)],gzfile(paste0(path2,"summit_5kb_extend.CpGeAP.minorstrand.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

#================================
#bedtools count
#bash
setwd(paste0(CGI_path,"n20/"))
system(paste0("for file in *minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n40/"))
system(paste0("for file in *minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n60/"))
system(paste0("for file in *minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

setwd(paste0(CGI_path,"n80/"))
system(paste0("for file in *minorstrand.fakebed.bed.gz; do bedtools intersect -wa -a ",CGI_path,"CpG.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done"))

#===============================
size=c(length(nc),length(ne),length(ne1D),length(ne2D),length(neAP),length(neNAP),length(np),length(np1D),length(np2D),length(nu))
size=rep(size,4)
files=list.files(path=CGI_path, pattern="minorstrand.result.bed.gz", recursive=T)

files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".minorstrand.result.bed.gz","",files.names)
data=read.delim(paste0(CGI_path, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(CGI_path, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size[i]
  data=rbind(data,data1)}

data$anno_region="promoter-like"
data$anno_region[grep("CpGe",data$group)]="enhancer-like"
data$anno_region[grep("CpGu",data$group)]="unclassed"
data$anno_region[grep("CpGc",data$group)]="CTCF-alone"
data$anno_region[grep("CpGeAP",data$group)]="enhancer-AP"
data$anno_region[grep("CpGeNAP",data$group)]="enhancer-NAP"
data$orientation="Others"
data$orientation[grep("2D",data$group)]="2D"
data$orientation[grep("1D",data$group)]="1D"
data$group="CpG_island"
write.table(data, gzfile(paste0(path_fig2_data,"CpG_island_minorstrand.result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================

C11=ggplot()+
  scale_color_npg()+
  geom_line(data = data[which(data$orientation=="Others"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")



