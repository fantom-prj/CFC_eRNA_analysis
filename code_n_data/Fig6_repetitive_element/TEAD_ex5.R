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
TEAD_folder=paste0(primary_folder,"code_n_data/Fig6_repetitive_element/TEAD_ChIP/")

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
data1=read.delim(paste0(primary_folder,"fig4/data/features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

#===============================================================================
#use 500bp non-overlap region region, defined from ./repetitive.R
subset1=read.delim(paste0(path_fig6_data,"ex5_cluster_used_motif.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=data1[which(data1$n5_string %in% subset1$n5_string),]

#===============================================================================
#n5 cluster base (from sala)
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed"), header=F, stringsAsFactors = F)
write.table(cluster_summit1,gzfile(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools

setwd(TEAD_folder)
# FANTOM6 ChIP-seq of TEAD4 on WTC11 iPSC line 
# the IDR thresholded peaks were used

system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz -b TEAD4_iPSC_IDR_p005.bed.gz | gzip > Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD.bed.gz")
system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed.gz -b TEAD4_2c_output_merge.bed.gz | gzip > Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD_merge.bed.gz")

#==========================
TEAD=read.delim(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD.bed.gz"), header=F, stringsAsFactors = F)
length(unique(TEAD$V4))#5533
TEAD=left_join(TEAD, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
TEAD$locS=TEAD$V2-TEAD$V2_ori-5000
TEAD$locE=TEAD$V3-TEAD$V2_ori-5000
TEAD1=TEAD[which(TEAD$V6 == "+"),c(1,28,29,4,11,12)]
TEAD2=TEAD[which(TEAD$V6 == "-"),c(1,29,28,4,11,12)]
TEAD2$locE=TEAD2$locE * (-1)
TEAD2$locS=TEAD2$locS * (-1)
TEAD2$V1="chr1"
TEAD1$V1="chr1"
colnames(TEAD2)=colnames(TEAD1)
TEAD1$locS=TEAD1$locS+5000
TEAD1$locE=TEAD1$locE+5000
TEAD2$locS=TEAD2$locS+5001
TEAD2$locE=TEAD2$locE+5001
TEAD=rbind(TEAD1,TEAD2)
colnames(TEAD)[c(4,5,6)]=c("n5_string","score","strand")
TEAD=left_join(TEAD, unique(data1[,c(1,9,27,85)]),by="n5_string",copy=F)
TEAD=TEAD[which(!is.na(TEAD$ex5cluster_class)),]
write.table(TEAD[order(TEAD$locS),],gzfile(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#grouping according to promoter type and regulatory element
TEAD=read.delim(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed"), header=F, stringsAsFactors = F)
cluster_summit1=left_join(cluster_summit1, unique(data1b[,c(1,9,27,77,80,85,111,117,129)]),by=c("V4"="n5_string"),copy=F)
cluster_summit1=cluster_summit1[which(!is.na(cluster_summit1$ex5cluster_class)),]

ntata=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_pl=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_el=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_ol=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_pli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_eli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_oli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_plni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]
ntata_elni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]
ntata_olni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]

ntata_plia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_elia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_olia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_plino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_elino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_olino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]


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

TEADtata=TEAD[which(TEAD$n5_string %in% ntata),]
TEADtata_p=TEAD[which(TEAD$n5_string %in% ntata_p),]
TEADtata_e=TEAD[which(TEAD$n5_string %in% ntata_e),]
TEADtata_m=TEAD[which(TEAD$n5_string %in% ntata_m),]
TEADtata_o=TEAD[which(TEAD$n5_string %in% ntata_o),]
TEADtata_pl=TEAD[which(TEAD$n5_string %in% ntata_pl),]
TEADtata_el=TEAD[which(TEAD$n5_string %in% ntata_el),]
TEADtata_ol=TEAD[which(TEAD$n5_string %in% ntata_ol),]
TEADtata_pli=TEAD[which(TEAD$n5_string %in% ntata_pli),]
TEADtata_eli=TEAD[which(TEAD$n5_string %in% ntata_eli),]
TEADtata_oli=TEAD[which(TEAD$n5_string %in% ntata_oli),]
TEADtata_plni=TEAD[which(TEAD$n5_string %in% ntata_plni),]
TEADtata_elni=TEAD[which(TEAD$n5_string %in% ntata_elni),]
TEADtata_olni=TEAD[which(TEAD$n5_string %in% ntata_olni),]
TEADtata_plia=TEAD[which(TEAD$n5_string %in% ntata_plia),]
TEADtata_elia=TEAD[which(TEAD$n5_string %in% ntata_elia),]
TEADtata_olia=TEAD[which(TEAD$n5_string %in% ntata_olia),]
TEADtata_plino=TEAD[which(TEAD$n5_string %in% ntata_plino),]
TEADtata_elino=TEAD[which(TEAD$n5_string %in% ntata_elino),]
TEADtata_olino=TEAD[which(TEAD$n5_string %in% ntata_olino),]

TEADCGI=TEAD[which(TEAD$n5_string %in% nCGI),]
TEADCGI_p=TEAD[which(TEAD$n5_string %in% nCGI_p),]
TEADCGI_e=TEAD[which(TEAD$n5_string %in% nCGI_e),]
TEADCGI_m=TEAD[which(TEAD$n5_string %in% nCGI_m),]
TEADCGI_o=TEAD[which(TEAD$n5_string %in% nCGI_o),]
TEADNull=TEAD[which(TEAD$n5_string %in% nNull),]
TEADNull_p=TEAD[which(TEAD$n5_string %in% nNull_p),]
TEADNull_e=TEAD[which(TEAD$n5_string %in% nNull_e),]
TEADNull_m=TEAD[which(TEAD$n5_string %in% nNull_m),]
TEADNull_o=TEAD[which(TEAD$n5_string %in% nNull_o),]

path1=paste0(TEAD_folder,"fakebed/")
write.table(TEADtata[order(TEADtata$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_p[order(TEADtata_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_e[order(TEADtata_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_m[order(TEADtata_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_o[order(TEADtata_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_pl[order(TEADtata_pl$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_pl.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_el[order(TEADtata_el$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_el.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_ol[order(TEADtata_ol$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_ol.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_pli[order(TEADtata_pli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_pli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_eli[order(TEADtata_eli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_eli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_oli[order(TEADtata_oli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_oli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plni[order(TEADtata_plni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elni[order(TEADtata_elni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olni[order(TEADtata_olni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plia[order(TEADtata_plia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elia[order(TEADtata_elia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olia[order(TEADtata_olia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plino[order(TEADtata_plino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elino[order(TEADtata_elino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olino[order(TEADtata_olino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

write.table(TEADCGI[order(TEADCGI$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_p[order(TEADCGI_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_e[order(TEADCGI_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_m[order(TEADCGI_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_o[order(TEADCGI_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull[order(TEADNull$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_p[order(TEADNull_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_e[order(TEADNull_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_m[order(TEADNull_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_o[order(TEADNull_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile(paste0(path1,"TEAD.location10001.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools count

setwd(path1)
system("for file in end5.summit_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a TEAD.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================
#parse the counting result
path1=paste0(TEAD_folder,"fakebed/")
size=c(length(nCGI_e),length(nCGI_m),length(nCGI_o),length(nCGI_p),length(nCGI),length(nNull_e),length(nNull_m),length(nNull_o),length(nNull_p),length(nNull),
       length(ntata_e),length(ntata_el),length(ntata_eli),length(ntata_elia),length(ntata_elino),length(ntata_elni),
       length(ntata_m),length(ntata_o),length(ntata_ol),length(ntata_oli),length(ntata_olia),length(ntata_olino),length(ntata_olni),
       length(ntata_p),length(ntata_pl),length(ntata_pli),length(ntata_plia),length(ntata_plino),length(ntata_plni),length(ntata))

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
data$CpGTATA[grep("pli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("eli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("oli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("plni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("elni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("olni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("plia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("elia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("olia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("plino",data$group)]="TATA_LTR_iPS_noMark"
data$CpGTATA[grep("elino",data$group)]="TATA_LTR_iPS_noMark"
data$CpGTATA[grep("olino",data$group)]="TATA_LTR_iPS_noMark"

data$group2="TEAD"
write.table(data,gzfile(paste0(path_fig6_data,"ex5_cluster_TEAD_distribution_result.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================

data=read.delim(paste0(path_fig6_data,"ex5_cluster_TEAD_distribution_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data$anno_region=factor(data$anno_region, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
data=data[which(data$CpGTATA %in% c("TATA_LTR","TATA_LTR_iPS","TATA_LTR_iPS_active","TATA_LTR_iPS_noMark","Null")),]
data$CpGTATA=factor(data$CpGTATA, levels=c("TATA_LTR","TATA_LTR_iPS","TATA_LTR_iPS_active","TATA_LTR_iPS_noMark","Null"))
data1=data%>%group_by(anno_region,CpGTATA)%>%dplyr::summarise(size=unique(size))
data1$label=paste0(data1$CpGTATA,": ",data1$size)
ggplot()+
  scale_color_manual(values=c("black", "#3C5488FF","green","purple","red","grey"))+
  coord_cartesian(xlim=c(-2500,2500), ylim=c(0,0.32))+
  geom_line(data = data[which(data$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=V4, y=V5, color=CpGTATA, group=CpGTATA), linewidth=0.25)+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.225, label=label), vjust=0, hjust=0, size=2.8, color="black")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.25, label=label), vjust=0, hjust=0, size=2.8, color="#3C5488FF")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS_active" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.275, label=label), vjust=0, hjust=0, size=2.8, color="green")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS_noMark" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.3, label=label), vjust=0, hjust=0, size=2.8, color="purple")+
  #geom_text(data=data1[which(data1$CpGTATA == "Null" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.5, label=label), vjust=0, hjust=0, size=2.8, color="red")+
  facet_grid(cols=vars(anno_region))+
  labs(color=NULL, x="Distance from the stranded ex5_cluster summit", y="% of ex5_cluster", title="LTR distribution from different regulatory elements")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))


#==========================
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed"), header=F, stringsAsFactors = F)

TEAD=read.delim(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD_merge.bed.gz"), header=F, stringsAsFactors = F)
length(unique(TEAD$V4))#10776
TEAD=left_join(TEAD, cluster_summit1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
TEAD$locS=TEAD$V2-TEAD$V2_ori-5000
TEAD$locE=TEAD$V3-TEAD$V2_ori-5000
TEAD1=TEAD[which(TEAD$V6 == "+"),c(1,14,15,4)]
TEAD2=TEAD[which(TEAD$V6 == "-"),c(1,15,14,4)]
TEAD2$locE=TEAD2$locE * (-1)
TEAD2$locS=TEAD2$locS * (-1)
TEAD2$V1="chr1"
TEAD1$V1="chr1"
colnames(TEAD2)=colnames(TEAD1)
TEAD1$locS=TEAD1$locS+5000
TEAD1$locE=TEAD1$locE+5000
TEAD2$locS=TEAD2$locS+5001
TEAD2$locE=TEAD2$locE+5001
TEAD=rbind(TEAD1,TEAD2)
colnames(TEAD)[c(4)]=c("n5_string")
TEAD=left_join(TEAD, unique(data1[,c(1,9,27,85)]),by="n5_string",copy=F)
TEAD=TEAD[which(!is.na(TEAD$ex5cluster_class)),]
write.table(TEAD[order(TEAD$locS),],gzfile(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD_merge.fakebed.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#grouping according to promoter type and regulatory element
TEAD=read.delim(paste0(TEAD_folder,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.TEAD_merge.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cluster_summit1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.5000bp_extend.bed"), header=F, stringsAsFactors = F)
cluster_summit1=left_join(cluster_summit1, unique(data1b[,c(1,9,27,77,80,85,111,117,129)]),by=c("V4"="n5_string"),copy=F)
cluster_summit1=cluster_summit1[which(!is.na(cluster_summit1$ex5cluster_class)),]

ntata=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_p=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_e=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_m=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("mRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_o=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA")]
ntata_pl=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_el=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_ol=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR")]
ntata_pli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_eli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_oli=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0)]
ntata_plni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]
ntata_elni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]
ntata_olni=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC ==0)]

ntata_plia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_elia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_olia=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "Yes" & cluster_summit1$K27ME3_iPS == "No")]
ntata_plino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("p_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_elino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("e_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]
ntata_olino=cluster_summit1$V4[which(cluster_summit1$ex5cluster_class %in% c("other_ncRNA") & cluster_summit1$CpGTATA == "TATA" & cluster_summit1$repeat_ex5cluster == "LTR" & cluster_summit1$iPSC >0 & cluster_summit1$K27Ac_iPS == "No" & cluster_summit1$K27ME3_iPS == "No")]


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

TEADtata=TEAD[which(TEAD$n5_string %in% ntata),]
TEADtata_p=TEAD[which(TEAD$n5_string %in% ntata_p),]
TEADtata_e=TEAD[which(TEAD$n5_string %in% ntata_e),]
TEADtata_m=TEAD[which(TEAD$n5_string %in% ntata_m),]
TEADtata_o=TEAD[which(TEAD$n5_string %in% ntata_o),]
TEADtata_pl=TEAD[which(TEAD$n5_string %in% ntata_pl),]
TEADtata_el=TEAD[which(TEAD$n5_string %in% ntata_el),]
TEADtata_ol=TEAD[which(TEAD$n5_string %in% ntata_ol),]
TEADtata_pli=TEAD[which(TEAD$n5_string %in% ntata_pli),]
TEADtata_eli=TEAD[which(TEAD$n5_string %in% ntata_eli),]
TEADtata_oli=TEAD[which(TEAD$n5_string %in% ntata_oli),]
TEADtata_plni=TEAD[which(TEAD$n5_string %in% ntata_plni),]
TEADtata_elni=TEAD[which(TEAD$n5_string %in% ntata_elni),]
TEADtata_olni=TEAD[which(TEAD$n5_string %in% ntata_olni),]
TEADtata_plia=TEAD[which(TEAD$n5_string %in% ntata_plia),]
TEADtata_elia=TEAD[which(TEAD$n5_string %in% ntata_elia),]
TEADtata_olia=TEAD[which(TEAD$n5_string %in% ntata_olia),]
TEADtata_plino=TEAD[which(TEAD$n5_string %in% ntata_plino),]
TEADtata_elino=TEAD[which(TEAD$n5_string %in% ntata_elino),]
TEADtata_olino=TEAD[which(TEAD$n5_string %in% ntata_olino),]

TEADCGI=TEAD[which(TEAD$n5_string %in% nCGI),]
TEADCGI_p=TEAD[which(TEAD$n5_string %in% nCGI_p),]
TEADCGI_e=TEAD[which(TEAD$n5_string %in% nCGI_e),]
TEADCGI_m=TEAD[which(TEAD$n5_string %in% nCGI_m),]
TEADCGI_o=TEAD[which(TEAD$n5_string %in% nCGI_o),]
TEADNull=TEAD[which(TEAD$n5_string %in% nNull),]
TEADNull_p=TEAD[which(TEAD$n5_string %in% nNull_p),]
TEADNull_e=TEAD[which(TEAD$n5_string %in% nNull_e),]
TEADNull_m=TEAD[which(TEAD$n5_string %in% nNull_m),]
TEADNull_o=TEAD[which(TEAD$n5_string %in% nNull_o),]

path1=paste0(TEAD_folder,"fakebed_merge/")
write.table(TEADtata[order(TEADtata$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_p[order(TEADtata_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_e[order(TEADtata_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_m[order(TEADtata_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_o[order(TEADtata_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_pl[order(TEADtata_pl$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_pl.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_el[order(TEADtata_el$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_el.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_ol[order(TEADtata_ol$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_ol.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_pli[order(TEADtata_pli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_pli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_eli[order(TEADtata_eli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_eli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_oli[order(TEADtata_oli$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_oli.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plni[order(TEADtata_plni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elni[order(TEADtata_elni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olni[order(TEADtata_olni$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olni.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plia[order(TEADtata_plia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elia[order(TEADtata_elia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olia[order(TEADtata_olia$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olia.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_plino[order(TEADtata_plino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_plino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_elino[order(TEADtata_elino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_elino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADtata_olino[order(TEADtata_olino$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADtata_olino.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

write.table(TEADCGI[order(TEADCGI$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_p[order(TEADCGI_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_e[order(TEADCGI_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_m[order(TEADCGI_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADCGI_o[order(TEADCGI_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADCGI_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull[order(TEADNull$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_p[order(TEADNull_p$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_p.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_e[order(TEADNull_e$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_e.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_m[order(TEADNull_m$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_m.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(TEADNull_o[order(TEADNull_o$locS),c(1:4)],gzfile(paste0(path1,"end5.summit_5kb_extend.TEADNull_o.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile(paste0(path1,"TEAD.location10001.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools count

setwd(path1)
system("for file in end5.summit_5kb_extend*.fakebed.bed.gz; do bedtools intersect -wa -a TEAD.location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#===============================================================================
#parse the counting result
path1=paste0(TEAD_folder,"fakebed_merge/")
size=c(length(nCGI_e),length(nCGI_m),length(nCGI_o),length(nCGI_p),length(nCGI),length(nNull_e),length(nNull_m),length(nNull_o),length(nNull_p),length(nNull),
       length(ntata_e),length(ntata_el),length(ntata_eli),length(ntata_elia),length(ntata_elino),length(ntata_elni),
       length(ntata_m),length(ntata_o),length(ntata_ol),length(ntata_oli),length(ntata_olia),length(ntata_olino),length(ntata_olni),
       length(ntata_p),length(ntata_pl),length(ntata_pli),length(ntata_plia),length(ntata_plino),length(ntata_plni),length(ntata))

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
data$CpGTATA[grep("pli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("eli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("oli",data$group)]="TATA_LTR_iPS"
data$CpGTATA[grep("plni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("elni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("olni",data$group)]="TATA_LTR_notiPS"
data$CpGTATA[grep("plia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("elia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("olia",data$group)]="TATA_LTR_iPS_active"
data$CpGTATA[grep("plino",data$group)]="TATA_LTR_iPS_noMark"
data$CpGTATA[grep("elino",data$group)]="TATA_LTR_iPS_noMark"
data$CpGTATA[grep("olino",data$group)]="TATA_LTR_iPS_noMark"

data$group2="TEAD"
write.table(data,gzfile(paste0(path_fig6_data,"ex5_cluster_TEAD_merge_distribution_result.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================

data=read.delim(paste0(path_fig6_data,"ex5_cluster_TEAD_merge_distribution_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data$anno_region=factor(data$anno_region, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
data=data[which(data$CpGTATA %in% c("TATA_LTR","TATA_LTR_iPS","TATA_LTR_iPS_active","TATA_LTR_iPS_noMark","Null")),]
data$CpGTATA=factor(data$CpGTATA, levels=c("TATA_LTR","TATA_LTR_iPS","TATA_LTR_iPS_active","TATA_LTR_iPS_noMark","Null"))
data1=data%>%group_by(anno_region,CpGTATA)%>%dplyr::summarise(size=unique(size))
data1$label=paste0(data1$CpGTATA,": ",data1$size)
ggplot()+
  scale_color_manual(values=c("black", "#3C5488FF","green","purple","red","grey"))+
  coord_cartesian(xlim=c(-2500,2500), ylim=c(0,0.32))+
  geom_line(data = data[which(data$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=V4, y=V5, color=CpGTATA, group=CpGTATA), linewidth=0.25)+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.225, label=label), vjust=0, hjust=0, size=2.8, color="black")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.25, label=label), vjust=0, hjust=0, size=2.8, color="#3C5488FF")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS_active" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.275, label=label), vjust=0, hjust=0, size=2.8, color="green")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA_LTR_iPS_noMark" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.3, label=label), vjust=0, hjust=0, size=2.8, color="purple")+
  #geom_text(data=data1[which(data1$CpGTATA == "Null" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-2500, y=0.5, label=label), vjust=0, hjust=0, size=2.8, color="red")+
  facet_grid(cols=vars(anno_region))+
  labs(color=NULL, x="Distance from the stranded ex5_cluster summit", y="% of ex5_cluster", title="LTR distribution from different regulatory elements")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))





