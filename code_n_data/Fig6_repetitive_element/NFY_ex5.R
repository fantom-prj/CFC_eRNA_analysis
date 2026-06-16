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
NFY_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/NFY_ChIP/")

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
data1=read.delim(paste0(primary_folder,"fig4/data/features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#data1_iPS=data1[which(data1$iPSC>0),]

#===============================================================================
#use 500bp non-overlap region region, defined from ./repetitive.R
subset1=read.delim(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=data1[which(data1$n5_string %in% subset1$n5_string),]

#===============================================================================
#n5 cluster base (from sala), intersect with all, filter to subset1 later
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz"), header=F, stringsAsFactors = F)
write.table(cluster_summit1,gzfile(paste0(NFY_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools

setwd(NFY_folder)
# ENCSR146UIC is a ChIP-seq of NFYB on WTC11 iPSC line
# the conservative IDR thresholded peaks were used

system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz -b ENCFF765MOP.bed.gz | gzip > Neuron_THP1.S3.end5.summit.table5.5000bp_extend.NFY.bed.gz")

#==========================
NFY=read.delim(paste0(NFY_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.NFY.bed.gz"), header=F, stringsAsFactors = F)
length(unique(NFY$V4))#18396
NFY=left_join(NFY, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
NFY$locS=NFY$V2-NFY$V2_ori-5000
NFY$locE=NFY$V3-NFY$V2_ori-5000
NFY1=NFY[which(NFY$V6 == "+"),c(1,18,19,4,11,12)]
NFY2=NFY[which(NFY$V6 == "-"),c(1,19,18,4,11,12)]
NFY2$locE=NFY2$locE * (-1)
NFY2$locS=NFY2$locS * (-1)
NFY2$V1="chr1"
NFY1$V1="chr1"
colnames(NFY2)=colnames(NFY1)
NFY1$locS=NFY1$locS+5000
NFY1$locE=NFY1$locE+5000
NFY2$locS=NFY2$locS+5001
NFY2$locE=NFY2$locE+5001
NFY=rbind(NFY1,NFY2)
colnames(NFY)[c(4,5,6)]=c("n5_string","score","strand")
NFY=left_join(NFY, unique(data1[,c(1,9,27,85)]),by="n5_string",copy=F)
NFY=NFY[which(!is.na(NFY$ex5cluster_class)),]
write.table(NFY[order(NFY$locS),],gzfile(paste0(NFY_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.NFY.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#grouping according to promoter type and regulatory element
NFY=read.delim(paste0(NFY_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.NFY.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz"), header=F, stringsAsFactors = F)
cluster_summit1=left_join(cluster_summit1, unique(data1b[,c(1,9,27,77,80,85,111,117,129)]),by=c("V4"="n5_string"),copy=F)
cluster_summit1=cluster_summit1[which(!is.na(cluster_summit1$ex5cluster_class)),]
cluster_summit1=cluster_summit1[which(cluster_summit1$iPSC > 0),] #take iPSC expressing alone

ntata=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]

ntata_pl=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_el=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_ol=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]

ntata_pla=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_ela=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_ola=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_plno=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_elno=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_olno=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]


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

NFYtata=NFY[which(NFY$n5_string %in% ntata),]
NFYtata_p=NFY[which(NFY$n5_string %in% ntata_p),]
NFYtata_e=NFY[which(NFY$n5_string %in% ntata_e),]
NFYtata_m=NFY[which(NFY$n5_string %in% ntata_m),]
NFYtata_o=NFY[which(NFY$n5_string %in% ntata_o),]
NFYtata_pl=NFY[which(NFY$n5_string %in% ntata_pl),]
NFYtata_el=NFY[which(NFY$n5_string %in% ntata_el),]
NFYtata_ol=NFY[which(NFY$n5_string %in% ntata_ol),]

NFYtata_pla=NFY[which(NFY$n5_string %in% ntata_pla),]
NFYtata_ela=NFY[which(NFY$n5_string %in% ntata_ela),]
NFYtata_ola=NFY[which(NFY$n5_string %in% ntata_ola),]
NFYtata_plno=NFY[which(NFY$n5_string %in% ntata_plno),]
NFYtata_elno=NFY[which(NFY$n5_string %in% ntata_elno),]
NFYtata_olno=NFY[which(NFY$n5_string %in% ntata_olno),]

NFYCGI=NFY[which(NFY$n5_string %in% nCGI),]
NFYCGI_p=NFY[which(NFY$n5_string %in% nCGI_p),]
NFYCGI_e=NFY[which(NFY$n5_string %in% nCGI_e),]
NFYCGI_m=NFY[which(NFY$n5_string %in% nCGI_m),]
NFYCGI_o=NFY[which(NFY$n5_string %in% nCGI_o),]
NFYNull=NFY[which(NFY$n5_string %in% nNull),]
NFYNull_p=NFY[which(NFY$n5_string %in% nNull_p),]
NFYNull_e=NFY[which(NFY$n5_string %in% nNull_e),]
NFYNull_m=NFY[which(NFY$n5_string %in% nNull_m),]
NFYNull_o=NFY[which(NFY$n5_string %in% nNull_o),]

path1=paste0(NFY_folder,"fakebed/")
write.table(NFYtata[order(NFYtata$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_p[order(NFYtata_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_e[order(NFYtata_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_m[order(NFYtata_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_o[order(NFYtata_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_pl[order(NFYtata_pl$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_pl.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_el[order(NFYtata_el$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_el.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_ol[order(NFYtata_ol$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_ol.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_pla[order(NFYtata_pla$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_pla.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_ela[order(NFYtata_ela$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_ela.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_ola[order(NFYtata_ola$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_ola.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_plno[order(NFYtata_plno$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_plno.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_elno[order(NFYtata_elno$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_elno.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYtata_olno[order(NFYtata_olno$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYtata_olno.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

write.table(NFYCGI[order(NFYCGI$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYCGI.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYCGI_p[order(NFYCGI_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYCGI_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYCGI_e[order(NFYCGI_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYCGI_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYCGI_m[order(NFYCGI_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYCGI_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYCGI_o[order(NFYCGI_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYCGI_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYNull[order(NFYNull$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYNull.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYNull_p[order(NFYNull_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYNull_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYNull_e[order(NFYNull_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYNull_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYNull_m[order(NFYNull_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYNull_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(NFYNull_o[order(NFYNull_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.NFYNull_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile(paste0(path1,"NFY.location10001.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools count

setwd(path1)
system("for file in end5.summit_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a NFY.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================
#parse the counting result
path1=paste0(NFY_folder,"fakebed/")
size=c(length(nCGI_e),length(nCGI_m),length(nCGI_o),length(nCGI_p),length(nCGI),length(nNull_e),length(nNull_m),length(nNull_o),length(nNull_p),length(nNull),
       length(ntata_e),length(ntata_el),length(ntata_ela),length(ntata_elno),
       length(ntata_m),length(ntata_o),length(ntata_ol),length(ntata_ola),length(ntata_olno),
       length(ntata_p),length(ntata_pl),length(ntata_pla),length(ntata_plno),length(ntata))

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
data$CpGTATA[grep("pl",data$group)]="TATA_LTR"
data$CpGTATA[grep("el",data$group)]="TATA_LTR"
data$CpGTATA[grep("ol",data$group)]="TATA_LTR"

data$CpGTATA[grep("pla",data$group)]="TATA_LTR_active"
data$CpGTATA[grep("ela",data$group)]="TATA_LTR_active"
data$CpGTATA[grep("ola",data$group)]="TATA_LTR_active"
data$CpGTATA[grep("plno",data$group)]="TATA_LTR_noMark"
data$CpGTATA[grep("elno",data$group)]="TATA_LTR_noMark"
data$CpGTATA[grep("olno",data$group)]="TATA_LTR_noMark"

data$group2="NFY"
write.table(data,gzfile(paste0(path_fig6_data,"ex5_cluster_NFY_distribution_result.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================


