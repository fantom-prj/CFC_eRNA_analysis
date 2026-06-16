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
path_fig6_data=paste0(primary_folder,"fig6/data/")
LTR_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/LTR_result/")

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)


#===============================================================================
#n5 cluster base (from sala) and restricted to ex5_cluster here: paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz"), header=F, stringsAsFactors = F)
write.table(cluster_summit1,gzfile(paste0(LTR_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
cluster_summit1=left_join(cluster_summit1, unique(data1[,c(1,9,27,85)]),by=c("V4"="n5_string"),copy=F)

#===============================================================================
#bash
#bedtools

setwd(LTR_folder)
system("bedtools merge -i UCSC.hg38.LTR.bed.gz -s -c 4,5,6 -o distinct,distinct,distinct | gzip > UCSC.hg38.LTR.merged.bed.gz")
system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz -b UCSC.hg38.LTR.merged.bed.gz -s | gzip > Neuron_THP1.S3.end5.summit.table5.5000bp_extend.LTR.bed.gz")
# -> strand specific overlap, no need marjor strand

#==========================
LTR=read.delim(paste0(LTR_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.LTR.bed.gz"), header=F, stringsAsFactors = F)
length(unique(LTR$V4))#30568
LTR=left_join(LTR, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
LTR$locS=LTR$V2-LTR$V2_ori-5000
LTR$locE=LTR$V3-LTR$V2_ori-5000
LTR1=LTR[which(LTR$V6 == "+"),c(1,14,15,4,10,11,12)]
LTR2=LTR[which(LTR$V6 == "-"),c(1,15,14,4,10,11,12)]
LTR2$locE=LTR2$locE * (-1)
LTR2$locS=LTR2$locS * (-1)
LTR2$V1="chr1"
LTR1$V1="chr1"
colnames(LTR2)=colnames(LTR1)
LTR1$locS=LTR1$locS+5000
LTR1$locE=LTR1$locE+5000
LTR2$locS=LTR2$locS+5001
LTR2$locE=LTR2$locE+5001
LTR=rbind(LTR1,LTR2)
colnames(LTR)[c(4,5,6,7)]=c("n5_string","subfamily","name","strand")
LTR=left_join(LTR, unique(data1[,c(1,9,27,83)]),by="n5_string",copy=F)
write.table(LTR[order(LTR$locS),],paste0(LTR_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.LTR.fakebed.tsv"), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#grouping according to promoter type and regulatory element
LTR=read.delim(paste0(LTR_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.LTR.fakebed.tsv"), header=T, stringsAsFactors = F, check.names = F)
ntata=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]

nCGI=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA %in% c("CGI","CGIap","CGInap"))]
nCGI_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA %in% c("CGI"))]
nCGI_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA %in% c("CGIap","CGInap"))]
nCGI_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA %in% c("CGI"))]
nCGI_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA %in% c("CGI"))]

nNull=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA %in% c("Null"))]
nNull_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA %in% c("Null"))]
nNull_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA %in% c("Null"))]
nNull_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA %in% c("Null"))]
nNull_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA %in% c("Null"))]

LTRtata=LTR[which(LTR$n5_string %in% ntata),]
LTRtata_p=LTR[which(LTR$n5_string %in% ntata_p),]
LTRtata_e=LTR[which(LTR$n5_string %in% ntata_e),]
LTRtata_m=LTR[which(LTR$n5_string %in% ntata_m),]
LTRtata_o=LTR[which(LTR$n5_string %in% ntata_o),]
LTRCGI=LTR[which(LTR$n5_string %in% nCGI),]
LTRCGI_p=LTR[which(LTR$n5_string %in% nCGI_p),]
LTRCGI_e=LTR[which(LTR$n5_string %in% nCGI_e),]
LTRCGI_m=LTR[which(LTR$n5_string %in% nCGI_m),]
LTRCGI_o=LTR[which(LTR$n5_string %in% nCGI_o),]
LTRNull=LTR[which(LTR$n5_string %in% nNull),]
LTRNull_p=LTR[which(LTR$n5_string %in% nNull_p),]
LTRNull_e=LTR[which(LTR$n5_string %in% nNull_e),]
LTRNull_m=LTR[which(LTR$n5_string %in% nNull_m),]
LTRNull_o=LTR[which(LTR$n5_string %in% nNull_o),]

path1=paste0(LTR_folder,"fakebed/")
write.table(LTRtata[order(LTRtata$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRtata.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRtata_p[order(LTRtata_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRtata_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRtata_e[order(LTRtata_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRtata_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRtata_m[order(LTRtata_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRtata_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRtata_o[order(LTRtata_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRtata_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRCGI[order(LTRCGI$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRCGI.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRCGI_p[order(LTRCGI_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRCGI_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRCGI_e[order(LTRCGI_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRCGI_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRCGI_m[order(LTRCGI_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRCGI_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRCGI_o[order(LTRCGI_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRCGI_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRNull[order(LTRNull$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRNull.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRNull_p[order(LTRNull_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRNull_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRNull_e[order(LTRNull_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRNull_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRNull_m[order(LTRNull_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRNull_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(LTRNull_o[order(LTRNull_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.LTRNull_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile(paste0(path1,"LTR.location10001.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools count

setwd(path1)
system("for file in end5.summit_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a LTR.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================
#parse the counting result
path1=paste0(LTR_folder,"fakebed/")
size=c(length(nCGI_e),length(nCGI_m),length(nCGI_o),length(nCGI_p),length(nCGI),length(nNull_e),length(nNull_m),length(nNull_o),length(nNull_p),length(nNull),length(ntata_e),length(ntata_m),length(ntata_o),length(ntata_p),length(ntata))

files=list.files(path=path1, pattern=".result.bed.gz", recursive=T)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(path1, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$V5=data$V5/size[1]
data$size=size[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(path1, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$V5=data1$V5/size[i]
  data1$size=size[i]
  data=rbind(data,data1)}

data$anno_region="All"
data$anno_region[grep("_e",data$group)]="e_ncRNA"
data$anno_region[grep("_p",data$group)]="p_ncRNA"
data$anno_region[grep("_m",data$group)]="mRNA"
data$anno_region[grep("_o",data$group)]="other_ncRNA"
data$CpGTATA="CGI"
data$CpGTATA[grep("tata",data$group)]="TATA"
data$CpGTATA[grep("Null",data$group)]="Null"

data$group2="LTR"
write.table(data,gzfile(paste0(path_fig6_data,"ex5_cluster_LTR_distribution_result.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
