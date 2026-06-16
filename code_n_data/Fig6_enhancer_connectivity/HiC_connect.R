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
output_folder=paste0(primary_folder,"code_n_data/Fig6_enhancer_connectivity/5kb/")
path_fig6_data=paste0(primary_folder,"fig6/data/")

#use ex5_cluster to intersect with hi-C region
#for 5kb resolution
#=============================================================================================================================
#Hi-C raw data -> DRA019572 (DRR614942- DRR614953)
#parse Hi-C GOTHIC result 
#additional filter: FDR < 0.01; read count >=5
path11="/analysisdata/fantom6/Interactome/HiC_KI/xufeng/GOTHIC/"
options(scipen=999)
iPS.hic=read.delim(paste0(path11,"IPS_sig_0.05_cis_HiC_5KB.bedpe"), header=F, stringsAsFactors = F)
NSC.hic=read.delim(paste0(path11,"NSC_sig_0.05_cis_HiC_5KB.bedpe"), header=F, stringsAsFactors = F)
NRN.hic=read.delim(paste0(path11,"NEU_sig_0.05_cis_HiC_5KB.bedpe"), header=F, stringsAsFactors = F)
iPS.hica=unique(iPS.hic[,c(1:3)])
iPS.hicb=unique(iPS.hic[,c(4:6)])
NSC.hica=unique(NSC.hic[,c(1:3)])
NSC.hicb=unique(NSC.hic[,c(4:6)])
NRN.hica=unique(NRN.hic[,c(1:3)])
NRN.hicb=unique(NRN.hic[,c(4:6)])
iPS.hica$key=paste0("n",1:nrow(iPS.hica))
iPS.hicb$key=paste0("n",1:nrow(iPS.hicb))
NSC.hica$key=paste0("n",1:nrow(NSC.hica))
NSC.hicb$key=paste0("n",1:nrow(NSC.hicb))
NRN.hica$key=paste0("n",1:nrow(NRN.hica))
NRN.hicb$key=paste0("n",1:nrow(NRN.hicb))
iPS.hic=left_join(iPS.hic,iPS.hica, by=c("V1","V2","V3"), copy=F)
iPS.hic=left_join(iPS.hic,iPS.hicb, by=c("V4","V5","V6"), copy=F, suffix=c("_left","_right"))
NSC.hic=left_join(NSC.hic,NSC.hica, by=c("V1","V2","V3"), copy=F)
NSC.hic=left_join(NSC.hic,NSC.hicb, by=c("V4","V5","V6"), copy=F, suffix=c("_left","_right"))
NRN.hic=left_join(NRN.hic,NRN.hica, by=c("V1","V2","V3"), copy=F)
NRN.hic=left_join(NRN.hic,NRN.hicb, by=c("V4","V5","V6"), copy=F, suffix=c("_left","_right"))

iPS.fdr=fread("GOTHIC_5KB_0.05_sig_full_list.tsv", header=T, select=c(1,2,3))
colnames(iPS.fdr)[c(2,3)]=c("IPS_q_value","count")
iPS.fdr=iPS.fdr[which(iPS.fdr$IPS_q_value < 0.01 & iPS.fdr$count >=5),]
iPS.hic$Contacts=paste0(iPS.hic$V1,"_",iPS.hic$V2+2500,"_",iPS.hic$V5+2500)
iPS.hic=right_join(iPS.hic,iPS.fdr,by="Contacts",copy=F)

NSC.fdr=fread("GOTHIC_5KB_0.05_sig_full_list.tsv", header=T, select=c(1,4,5))
colnames(NSC.fdr)[c(2,3)]=c("NSC_q_value","count")
NSC.fdr=NSC.fdr[which(NSC.fdr$NSC_q_value < 0.01 & NSC.fdr$count >=5),]
NSC.hic$Contacts=paste0(NSC.hic$V1,"_",NSC.hic$V2+2500,"_",NSC.hic$V5+2500)
NSC.hic=right_join(NSC.hic,NSC.fdr,by="Contacts",copy=F)

NRN.fdr=fread("GOTHIC_5KB_0.05_sig_full_list.tsv", header=T, select=c(1,6,7))
colnames(NRN.fdr)[c(2,3)]=c("NEU_q_value","count")
NRN.fdr=NRN.fdr[which(NRN.fdr$NEU_q_value < 0.01 & NRN.fdr$count >=5),]
NRN.hic$Contacts=paste0(NRN.hic$V1,"_",NRN.hic$V2+2500,"_",NRN.hic$V5+2500)
NRN.hic=right_join(NRN.hic,NRN.fdr,by="Contacts",copy=F)

iPS.hic$IPS_q_value= format(iPS.hic$IPS_q_value, scientific = T)
NSC.hic$NSC_q_value= format(NSC.hic$NSC_q_value, scientific = T)
NRN.hic$NEU_q_value= format(NRN.hic$NEU_q_value, scientific = T)

write.table(iPS.hic,gzfile(paste0(output_folder,"iPSC_GOTHIC_5KB_intermediate.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(NSC.hic,gzfile(paste0(output_folder,"NSC_GOTHIC_5KB_intermediate.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(NRN.hic,gzfile(paste0(output_folder,"NRN_GOTHIC_5KB_intermediate.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

write.table(iPS.hica,gzfile(paste0(output_folder,"iPSCa_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(iPS.hicb,gzfile(paste0(output_folder,"iPSCb_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(NSC.hica,gzfile(paste0(output_folder,"NSCa_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(NSC.hicb,gzfile(paste0(output_folder,"NSCb_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(NRN.hica,gzfile(paste0(output_folder,"NRNa_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(NRN.hicb,gzfile(paste0(output_folder,"NRNb_001.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#intermediate Hi-C result located in [primary_folder]/code_n_data/Fig6_enhancer_connectivity/5kb
NRN.hic=fread(paste0(output_folder,"NRN_GOTHIC_5KB_intermediate.tsv.gz"), header=T)

#===============================================================================
#bash
#intersect ex5_cluster summit with the 5kb bin region

setwd(output_folder)
system("for file in *_001.bed.gz; do bedtools intersect -wa -wb -a Neuron_THP1.S3.end5.summit.table5.bed.gz -b \"$file\" | gzip > \"${file%.bed.gz}.n5_string_summit.bed.gz\"; done")

#intermediate Hi-C result located in [primary_folder]/code_n_data/Fig6_enhancer_connectivity/5kb

#===============================================================================
#summarize the intersected files from each chromosome
setwd(output_folder)
data1=fread(paste0(output_folder,"iPSCa_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
data2=fread(paste0(output_folder,"iPSCb_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
iPS.hic1=iPS.hic[which(iPS.hic$key_left %in% unique(data1$V10) & iPS.hic$key_right %in% unique(data2$V10)),]
need=unique(iPS.hic1$V1)

for (i in 1: length(need)){
  data1a=data1[grep(need[i],data1$V1),]
  data2a=data2[grep(need[i],data2$V1),]
  iPS.hica=iPS.hic1[which(iPS.hic1$V1 == need[i]),c(7:10)]
  iPS.hica=left_join(iPS.hica, data1a[,c(2:4)], by=c("key_left"="V10"),copy=F, relationship = "many-to-many")
  iPS.hica=left_join(iPS.hica, data2a[,c(2:4)], by=c("key_right"="V10"),copy=F, relationship = "many-to-many")
  colnames(iPS.hica)[c(5,6,7,8)]=c("start_left","ex5_cluster_left","start_right","ex5_cluster_right")
  iPS.hica$together=paste0(iPS.hica$ex5_cluster_left,"|",iPS.hica$ex5_cluster_right)
  dup=unique(iPS.hica$together[which(duplicated(iPS.hica$together))])
  
  if (length(dup)==0){ iPS.hica2=iPS.hica
  iPS.hica1=iPS.hica[-which(iPS.hica$together %in% dup),]} else {
    iPS.hica1=iPS.hica[-which(iPS.hica$together %in% dup),]
    iPS.hica2=iPS.hica[which(iPS.hica$together %in% dup),]
    iPS.hica2=iPS.hica2%>%group_by(together)%>%slice_min(IPS_q_value)}
  iPS.hica=rbind(iPS.hica1,iPS.hica2)
  write.table(iPS.hica[order(iPS.hica$start_left,iPS.hica$start_right),c(3:8)],gzfile(paste0(need[i],"_iPSC_001_n5_string_summit_HiC_contact5k.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)
}
setwd(output_folder)
data1=fread(paste0(output_folder,"NSCa_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
data2=fread(paste0(output_folder,"NSCb_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
NSC.hic1=NSC.hic[which(NSC.hic$key_left %in% unique(data1$V10) & NSC.hic$key_right %in% unique(data2$V10)),]
for (i in 1: length(need)){
  data1a=data1[grep(need[i],data1$V1),]
  data2a=data2[grep(need[i],data2$V1),]
  NSC.hica=NSC.hic1[which(NSC.hic1$V1 == need[i]),c(7:10)]
  NSC.hica=left_join(NSC.hica, data1a[,c(2:4)], by=c("key_left"="V10"),copy=F, relationship = "many-to-many")
  NSC.hica=left_join(NSC.hica, data2a[,c(2:4)], by=c("key_right"="V10"),copy=F, relationship = "many-to-many")
  colnames(NSC.hica)[c(5,6,7,8)]=c("start_left","ex5_cluster_left","start_right","ex5_cluster_right")
  NSC.hica$together=paste0(NSC.hica$ex5_cluster_left,"|",NSC.hica$ex5_cluster_right)
  dup=unique(NSC.hica$together[which(duplicated(NSC.hica$together))])
  if (length(dup)==0){ NSC.hica2=NSC.hica
  NSC.hica1=NSC.hica[-which(NSC.hica$together %in% dup),]} else {
    NSC.hica1=NSC.hica[-which(NSC.hica$together %in% dup),]
    NSC.hica2=NSC.hica[which(NSC.hica$together %in% dup),]
    NSC.hica2=NSC.hica2%>%group_by(together)%>%slice_min(IPS_q_value)}
  NSC.hica=rbind(NSC.hica1,NSC.hica2)
  write.table(NSC.hica[order(NSC.hica$start_left,NSC.hica$start_right),c(3:8)],gzfile(paste0(need[i],"_NSC_001_n5_string_summit_HiC_contact5k.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)
}
setwd(output_folder)
data1=fread(paste0(output_folder,"NRNa_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
data2=fread(paste0(output_folder,"NRNb_001.n5_string_summit.bed.gz"),header=F, select=c(1,2,4,10))
NRN.hic1=NRN.hic[which(NRN.hic$key_left %in% unique(data1$V10) & NRN.hic$key_right %in% unique(data2$V10)),]
for (i in 1: length(need)){
  data1a=data1[grep(need[i],data1$V1),]
  data2a=data2[grep(need[i],data2$V1),]
  NRN.hica=NRN.hic1[which(NRN.hic1$V1 == need[i]),c(7:10)]
  NRN.hica=left_join(NRN.hica, data1a[,c(2:4)], by=c("key_left"="V10"),copy=F, relationship = "many-to-many")
  NRN.hica=left_join(NRN.hica, data2a[,c(2:4)], by=c("key_right"="V10"),copy=F, relationship = "many-to-many")
  colnames(NRN.hica)[c(5,6,7,8)]=c("start_left","ex5_cluster_left","start_right","ex5_cluster_right")
  NRN.hica$together=paste0(NRN.hica$ex5_cluster_left,"|",NRN.hica$ex5_cluster_right)
  dup=unique(NRN.hica$together[which(duplicated(NRN.hica$together))])
  if (length(dup)==0){ NRN.hica2=NRN.hica
  NRN.hica1=NRN.hica[-which(NRN.hica$together %in% dup),]} else {
    NRN.hica1=NRN.hica[-which(NRN.hica$together %in% dup),]
    NRN.hica2=NRN.hica[which(NRN.hica$together %in% dup),]
    NRN.hica2=NRN.hica2%>%group_by(together)%>%slice_min(IPS_q_value)}
  NRN.hica=rbind(NRN.hica1,NRN.hica2)
  write.table(NRN.hica[order(NRN.hica$start_left,NRN.hica$start_right),c(3:8)],gzfile(paste0(need[i],"_NRN_001_n5_string_summit_HiC_contact5k.tsv.gz")),row.names=F, col.names=T, sep="\t", quote=F)
}

rm(iPS.hic,NSC.hic,NRN.hic)

#===============================================================================

n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
coding_n5_string=unique(table5$n5_string[which(table5$Gencode_transcriptClass2 == "protein_coding")])

i=1
files=list.files(path=output_folder, pattern="001_n5_string_summit_HiC_contact5k")
files.names1=sapply(strsplit(files,"_"),"[",1)
files.names2=sapply(strsplit(files,"_"),"[",2)
n5_interact=read.delim(paste0(output_folder, files[i]),header=T, stringsAsFactors = F, check.names = F)
n5_interact$distance=n5_interact$start_right-n5_interact$start_left
n5_interact=left_join(n5_interact, data1[,c(1,9)], by=c("ex5_cluster_left"="n5_string"), copy=F)
n5_interact=left_join(n5_interact, data1[,c(1,9)], by=c("ex5_cluster_right"="n5_string"), copy=F, suffix=c("_left","_right"))
n5_interact$ex5cluster_class_left[which(n5_interact$ex5_cluster_left %in% coding_n5_string)]="mRNA"
n5_interact$ex5cluster_class_right[which(n5_interact$ex5_cluster_right %in% coding_n5_string)]="mRNA"
n5_interact=n5_interact[which(!is.na(n5_interact$ex5cluster_class_left) & !is.na(n5_interact$ex5cluster_class_right)),]
n5_interact=n5_interact[-which(n5_interact$ex5cluster_class_left == "mRNA" & n5_interact$ex5cluster_class_right == "mRNA"),]
n5_interact=n5_interact[-which(n5_interact$ex5cluster_class_left != "mRNA" & n5_interact$ex5cluster_class_right != "mRNA"),]
n5_interact$left_bin=paste0(sapply(strsplit(n5_interact$Contacts,"_"),"[",1),"_",sapply(strsplit(n5_interact$Contacts,"_"),"[",2))
n5_interact$right_bin=paste0(sapply(strsplit(n5_interact$Contacts,"_"),"[",1),"_",sapply(strsplit(n5_interact$Contacts,"_"),"[",3))
colnames(n5_interact)[2]="q_value"
n5_interact$q_value= format(n5_interact$q_value, scientific = T)
n5_interact$source_group=n5_interact$ex5cluster_class_left
n5_interact$source_group[which(n5_interact$ex5cluster_class_left=="mRNA")]=n5_interact$ex5cluster_class_right[which(n5_interact$ex5cluster_class_left=="mRNA")]
n5_interact$source_ex5=n5_interact$ex5_cluster_left
n5_interact$source_ex5[which(n5_interact$ex5cluster_class_left=="mRNA")]=n5_interact$ex5_cluster_right[which(n5_interact$ex5cluster_class_left=="mRNA")]
n5_interact$source_bin=n5_interact$left_bin
n5_interact$source_bin[which(n5_interact$ex5cluster_class_left=="mRNA")]=n5_interact$right_bin[which(n5_interact$ex5cluster_class_left=="mRNA")]
n5_interact$target_ex5=n5_interact$ex5_cluster_left
n5_interact$target_ex5[which(n5_interact$ex5cluster_class_left!="mRNA")]=n5_interact$ex5_cluster_right[which(n5_interact$ex5cluster_class_left!="mRNA")]
n5_interact$target_bin=n5_interact$left_bin
n5_interact$target_bin[which(n5_interact$ex5cluster_class_left!="mRNA")]=n5_interact$right_bin[which(n5_interact$ex5cluster_class_left!="mRNA")]

n5_interact=left_join(n5_interact, n5cluster[,c(4,80)], by=c("source_ex5"="n5_string"),copy=F)
n5_interact=left_join(n5_interact, n5cluster[,c(4,80)], by=c("target_ex5"="n5_string"),copy=F, suffix=c("_source","_target"))
n5_interact$cell=files.names2[i]
for (i in 2: length(files)){
  n5_interact1=read.delim(paste0(output_folder, files[i]),header=T, stringsAsFactors = F, check.names = F)
  n5_interact1$distance=n5_interact1$start_right-n5_interact1$start_left
  n5_interact1=left_join(n5_interact1, data1[,c(1,9)], by=c("ex5_cluster_left"="n5_string"), copy=F)
  n5_interact1=left_join(n5_interact1, data1[,c(1,9)], by=c("ex5_cluster_right"="n5_string"), copy=F, suffix=c("_left","_right"))
  n5_interact1$ex5cluster_class_left[which(n5_interact1$ex5_cluster_left %in% coding_n5_string)]="mRNA"
  n5_interact1$ex5cluster_class_right[which(n5_interact1$ex5_cluster_right %in% coding_n5_string)]="mRNA"
  n5_interact1=n5_interact1[which(!is.na(n5_interact1$ex5cluster_class_left) & !is.na(n5_interact1$ex5cluster_class_right)),]
  n5_interact1=n5_interact1[-which(n5_interact1$ex5cluster_class_left == "mRNA" & n5_interact1$ex5cluster_class_right == "mRNA"),]
  n5_interact1=n5_interact1[-which(n5_interact1$ex5cluster_class_left != "mRNA" & n5_interact1$ex5cluster_class_right != "mRNA"),]
  n5_interact1$left_bin=paste0(sapply(strsplit(n5_interact1$Contacts,"_"),"[",1),"_",sapply(strsplit(n5_interact1$Contacts,"_"),"[",2))
  n5_interact1$right_bin=paste0(sapply(strsplit(n5_interact1$Contacts,"_"),"[",1),"_",sapply(strsplit(n5_interact1$Contacts,"_"),"[",3))
  colnames(n5_interact1)[2]="q_value"
  n5_interact1$source_group=n5_interact1$ex5cluster_class_left
  n5_interact1$source_group[which(n5_interact1$ex5cluster_class_left=="mRNA")]=n5_interact1$ex5cluster_class_right[which(n5_interact1$ex5cluster_class_left=="mRNA")]
  n5_interact1$source_ex5=n5_interact1$ex5_cluster_left
  n5_interact1$source_ex5[which(n5_interact1$ex5cluster_class_left=="mRNA")]=n5_interact1$ex5_cluster_right[which(n5_interact1$ex5cluster_class_left=="mRNA")]
  n5_interact1$source_bin=n5_interact1$left_bin
  n5_interact1$source_bin[which(n5_interact1$ex5cluster_class_left=="mRNA")]=n5_interact1$right_bin[which(n5_interact1$ex5cluster_class_left=="mRNA")]
  n5_interact1$target_ex5=n5_interact1$ex5_cluster_left
  n5_interact1$target_ex5[which(n5_interact1$ex5cluster_class_left!="mRNA")]=n5_interact1$ex5_cluster_right[which(n5_interact1$ex5cluster_class_left!="mRNA")]
  n5_interact1$target_bin=n5_interact1$left_bin
  n5_interact1$target_bin[which(n5_interact1$ex5cluster_class_left!="mRNA")]=n5_interact1$right_bin[which(n5_interact1$ex5cluster_class_left!="mRNA")]
  
  n5_interact1=left_join(n5_interact1, n5cluster[,c(4,80)], by=c("source_ex5"="n5_string"),copy=F)
  n5_interact1=left_join(n5_interact1, n5cluster[,c(4,80)], by=c("target_ex5"="n5_string"),copy=F, suffix=c("_source","_target"))
  n5_interact1$cell=files.names2[i]
  n5_interact=rbind(n5_interact,n5_interact1)}

length(unique(n5_interact$Contacts)) # 138040
n5_interact$ex5_contact_ID=paste0(n5_interact$source_ex5,"_",n5_interact$target_ex5)
length(unique(n5_interact$ex5_contact_ID)) # 340179
n5_interact$CpGTATA_target[grep("CGI",n5_interact$CpGTATA_target)]="CGI"
write.table(n5_interact,gzfile(paste0(output_folder,"hiC_result_allcell_ex5_cluster_summit.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#collapse cell type
n5_interact=read.delim(paste0(output_folder,"hiC_result_allcell_ex5_cluster_summit.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5_interact1=n5_interact[,c(19,20,2)]%>%group_by(ex5_contact_ID)%>%dplyr::summarise(cell=paste(cell,collapse=";"),q_value=paste(q_value, collapse=";"))
n5_interact2=unique(n5_interact[,c(1,7,12:18,20)])
n5_interact2=left_join(n5_interact2,n5_interact1, by="ex5_contact_ID",copy=F)
table5k=table5%>%group_by(n5_string)%>%dplyr::summarise(geneID=paste(unique(T4_gene_ID),collapse=";"),geneName=paste(unique(T4_gene_name),collapse=";"))
n5_interact2=left_join(n5_interact2, table5k, by=c("source_ex5"="n5_string"),copy=F)
n5_interact2=left_join(n5_interact2, table5k, by=c("target_ex5"="n5_string"),copy=F, suffix=c("_source","_target"))
n5_interact2$target_group="coding_gene"
n5_interact2=n5_interact2[which(n5_interact2$source_ex5 %in% data1$n5_string),]
write.table(n5_interact2,gzfile(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_summit_cell_collapsed.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# Hi-C connectivity file located in [primary_folder]/fig6/data

#===============================================================================
#find the source number, GC content and ATAC hit
setwd(output_folder)
options(scipen=999)
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$bin_5kb=paste0(data1$chr,"_",ceiling(data1$end/5000)*5000-2500)
data1$bin_5kb_start=ceiling(data1$end/5000)*5000-5000
data1$bin_5kb_end=ceiling(data1$end/5000)*5000
write.table(unique(data1[order(data1$chr,data1$bin_5kb_start),c("chr","bin_5kb_start","bin_5kb_end","bin_5kb")]), gzfile("hiC_hit_source_bin.bed.gz"),col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#getfasta and GC content
setwd(output_folder)
system("bedtools getfasta -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed hiC_hit_source_bin.bed.gz | gzip > hiC_hit_source_bin.tsv.gz")
cmd = "zcat hiC_hit_source_bin.tsv.gz | awk '{
    seq = $5;
    g = gsub(/G/, \"\", seq);
    c = gsub(/C/, \"\", seq);
    gc_content = (g + c) / 5000;
    OFS=\"\\t\"; print $0, gc_content;
}' | gzip > hiC_hit_source_bin_with_gc.tsv.gz"
system(cmd)
#===============================================================================
#link GC content of 5kb bin region to the ex5_cluster 
CG=fread(paste0(output_folder,"hiC_hit_source_bin_with_gc.tsv.gz"), header=F, select=c(4,6))
colnames(CG)[2]="bin_5kb_GCcontent"
data1=left_join(data1, CG, by=c("bin_5kb"="V4"),copy=F)

#===============================================================================
#bash
#find the total count from HiC of the source 5kb bin
#source data of Hi-C not provided, please download from DDBJ
setwd("/analysisdata/fantom6/Interactome/HiC_KI")
system("zcat iPSC_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $2, $3-1, $3}' | sort -k1,1 -k2,2n > iPSC_contacts1.bed")
system("zcat iPSC_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $4, $5-1, $5}' | sort -k1,1 -k2,2n > iPSC_contacts2.bed")
system("zcat NSC_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $2, $3-1, $3}' | sort -k1,1 -k2,2n > NSC_contacts1.bed")
system("zcat NSC_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $4, $5-1, $5}' | sort -k1,1 -k2,2n > NSC_contacts2.bed")
system("zcat NEU_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $2, $3-1, $3}' | sort -k1,1 -k2,2n > NRN_contacts1.bed")
system("zcat NEU_contacts.pairs.gz | awk 'BEGIN{OFS=\"\t\"} {print $4, $5-1, $5}' | sort -k1,1 -k2,2n > NRN_contacts2.bed")

#split
setwd("/analysisdata/fantom6/Interactome/HiC_KI/hiC_count2")
system("for file in *_contacts*; do split -l 100000000 \"$file\" \"${file%.bed}_\"; done")

#bedtools count
setwd("/analysisdata/fantom6/Interactome/HiC_KI/hiC_count2")
system(paste0("for file in *_contacts*; do bedtools intersect -wa -a ",output_folder,"hiC_hit_source_bin.bed.gz -b \"$file\" -c > \"hiC_hit_source_bin.${file}.bed\"; done"))

#===============================================================================
path12="/analysisdata/fantom6/Interactome/HiC_KI/hiC_count2/"
files=list.files(path=path12, pattern="bed")
need=c("iPSC","NSC","NRN")
for (j in 1:3){
  files1=files[grep(need[j],files)]
  count1=read.delim(paste0(path12,files1[1]), header=F, stringsAsFactors = F)
  for (i in 2:length(files1)){
    count2=read.delim(paste0(path12,files1[i]), header=F, stringsAsFactors = F)
    count1=rbind(count1, count2)
    count1=count1%>%group_by(V1, V2, V3, V4)%>%dplyr::summarise(V5=sum(V5))}
  write.table(count1,paste0(path12,"hiC_count_",need[j],".tsv"),col.names=F, row.names=F, sep="\t", quote=F)}

files=list.files(path=path12, pattern="tsv")
count1=read.delim(paste0(path12, files[1]), header=F, stringsAsFactors = F)
for (i in 2:3){
  count2=read.delim(paste0(path12, files[i]), header=F, stringsAsFactors = F)
  count1=left_join(count1, count2[,c(4,5)], by="V4", copy=F)}
colnames(count1)[c(5:7)]=c("HiCrawCount_iPSC","HiCrawCount_NSC","HiCrawCount_Neuron")

data1=left_join(data1, count1[,c(4:7)], by=c("bin_5kb"="V4"),copy=F)
write.table(count1,gzfile(paste0(output_folder,"hiC_hit_source_bin_count.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/code_n_data/Fig6_enhancer_connectivity/5kb

#===============================================================================
#take the scATAC count as bulk ATAC count
#source scATAC data can be found from DDBJ
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)
library(chromVAR)
peak_hiC <- getPeaks(paste0(output_folder,"hiC_hit_source_bin.bed.gz"), sort_peaks = TRUE)
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
    features = peak_hiC,
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
c0$start=sapply(strsplit(c0$peakID,"-"),"[",2)
c0$bin_5kb=paste0(c0$chr,"_",as.numeric(c0$start)+2499)
write.table(c0[,c(7,1,3,4)],gzfile(paste0(output_folder,"hiC_hit_source_bin_ATAC_celltype.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#file located in [primary_folder]/code_n_data/Fig6_enhancer_connectivity/5kb

#===============================================================================
#add ATAC CPM for Hi-C source bin
ATACcount=read.delim(paste0(output_folder,"hiC_hit_source_bin_ATAC_celltype.tsv.gz"), row.names=1, header=T, stringsAsFactors = F, check.names = F)
d <- DGEList(counts=ATACcount)
TMM <- calcNormFactors(d, method="TMM")
TMM.ATAC=cpm(TMM, normalized.lib.sizes=TRUE)
write.table(TMM.ATAC,paste0(output_folder,"hiC_hit_source_bin_ATAC_TMM_celltype.tsv.gz"), col.names=T, row.names=T, sep="\t", quote=F)

ATACCPM=read.delim(paste0(output_folder,"hiC_hit_source_bin_ATAC_TMM_celltype.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ATACCPM$bin_5kb=rownames(ATACCPM)
colnames(ATACCPM)[1:3]=c("ATACCPM_iPSC","ATACCPM_NSC","ATACCPM_Neuron")
data1=left_join(data1, ATACCPM, by="bin_5kb",copy=F)
write.table(data1,paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), col.names=T, row.names=F, quote=F, sep="\t")

#===============================================================================
#add connectivity
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5_interact2=read.delim(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_summit_cell_collapsed.tsv.gz"), header=T, stringsAsFactors = F, check.names = T)

q2=unique(n5_interact2)%>%group_by(source_ex5)%>%dplyr::summarise(codingGene_connectivity=length(unique(target_bin)))
##one Hi-C contact only count for once, adjust for target bin alone
##adjustment from source was done next section, remove if a bin contain both TATA and CpG

data1=left_join(data1, q2, by=c("n5_string"="source_ex5"), copy=F)
write.table(data1,paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), col.names=T, row.names=F, quote=F, sep="\t")

#===============================================================================
#extract features from ex5_cluster annotation
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

exclude0=unique(data1$bin_5kb[which(data1$CpGTATA == "Others")])
data2=data1[-which(data1$bin_5kb %in% exclude0),] #Others refer to multiple regulatory element, bins contain others were removed.
data2$CpGTATA[which(data2$CpGTATA == "Null")]=NA
data2a=data2%>%group_by(bin_5kb)%>%dplyr::summarise(CpGTATA=paste(unique(na.omit(CpGTATA)),collapse=";"),count=n())
exclude1=data2a$bin_5kb[grep(";",data2a$CpGTATA)]
data2b=data2[-which(data2$bin_5kb %in% exclude1),] #remove the bins that contain multiple regulatory elements
data2b$CpGTATA[which(is.na(data2b$CpGTATA))]="Null"
data2b$rank=1
data2b$rank[which(data2b$CpGTATA == "Null")]=0
data2b=data2b%>%group_by(bin_5kb)%>%slice_max(rank)
data3=data2b[which(data2b$ex5cluster_class %in% c("p_ncRNA","e_ncRNA")),]
data3=data3%>%group_by(bin_5kb,CpGTATA)%>%slice_min(ex5cluster_class)
data3=data3%>%group_by(bin_5kb,CpGTATA)%>%slice_sample(n=1)

data3=data3[which(!is.na(data3$codingGene_connectivity)),] #only consider those ex5_clusters showing at least one connection
write.table(data3[,c("n5_string","ex5cluster_class","CpGTATA","codingGene_connectivity","bin_5kb_GCcontent","ATACCPM_iPSC","ATACCPM_NSC","ATACCPM_Neuron","iPSC","NSC","Neuron","orientation","SE_all")],gzfile(paste0(path_fig6_data,"HiC5kb_q001_ex5cluster_ATAC_GC.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data

#===============================================================================
#are the source and target ex5_cluster enriched with same regulatory element (CGI/TATA)?

#define CGI TATA from the mRNA side
#data1$all_CGI="No"
#data1$all_CGI[which(data1$upstream_CpG_island == "Yes" | data1$downstream_CpG_island == "Yes")]="Yes"

#this contain ex5_cluster of all mRNA

n5_interact2=read.delim(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_summit_cell_collapsed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5_interact2=n5_interact2[which(n5_interact2$source_group == "e_ncRNA"),]
#only allow source bin with one regulatory element
n5_interact2s=unique(n5_interact2[,c("source_bin","CpGTATA_source")])
n5_interact2s$CpGTATA_source[which(n5_interact2s$CpGTATA_source == "Null")]=NA
n5_interact2s=n5_interact2s%>%group_by(source_bin)%>%dplyr::summarise(CpGTATA_source=paste(unique(na.omit(CpGTATA_source)),collapse=";"))
exclude1=n5_interact2s$source_bin[grep(";",n5_interact2s$CpGTATA_source)]
exclude0=n5_interact2s$source_bin[grep("Others", n5_interact2s$CpGTATA_source)] #Others refer to multiple regulatory element, bins contain others were removed.
n5_interact2s=n5_interact2s[-which(n5_interact2s$source_bin %in% union(exclude1, exclude0)),]
n5_interact2s$CpGTATA_source[which(n5_interact2s$CpGTATA_source == "")]="Null"
n5_interact2s$value="Yes"
n5_interact3s=spread(n5_interact2s,key=2,value=3)
n5_interact3s[is.na(n5_interact3s)]="No"

#only allow target bin with one regulatory element
n5_interact2t=unique(n5_interact2[,c("target_bin","CpGTATA_target")])
n5_interact2t$CpGTATA_target[which(n5_interact2t$CpGTATA_target == "Null")]=NA
n5_interact2t=n5_interact2t%>%group_by(target_bin)%>%dplyr::summarise(CpGTATA_target=paste(unique(na.omit(CpGTATA_target)),collapse=";"))
exclude1=n5_interact2t$target_bin[grep(";",n5_interact2t$CpGTATA_target)]
exclude0=n5_interact2t$target_bin[grep("Others", n5_interact2t$CpGTATA_target)] #Others refer to multiple regulatory element, bins contain others were removed.
n5_interact2t=n5_interact2t[-which(n5_interact2t$target_bin %in% union(exclude1, exclude0)),]
n5_interact2t$CpGTATA_target[which(n5_interact2t$CpGTATA_target == "")]="Null"
n5_interact2t$value="Yes"
n5_interact3t=spread(n5_interact2t,key=2,value=3)
n5_interact3t[is.na(n5_interact3t)]="No"

final=left_join(unique(n5_interact2[,c("Contacts","source_bin","target_bin")]),n5_interact3s, by="source_bin",copy=F)
final=left_join(final,n5_interact3t, by="target_bin",copy=F, suffix=c("_source","_target"))
colnames(final)[c(4,5)]=paste0(colnames(final)[c(4,5)],"_source")
colnames(final)[8]=paste0(colnames(final)[8],"_target")
final=final[which(!is.na(final$CGIap_source) & !is.na(final$CGI_target)),]

k13=final%>%group_by(CGIap_source,CGI_target)%>%dplyr::summarise(count=n(),group="CGIap_CGI")
k13$label=paste0(k13$CGIap_source,"_",k13$CGI_target)
k14=final%>%group_by(CGIap_source,TATA_target)%>%dplyr::summarise(count=n(),group="CGIap_TATA")
k14$label=paste0(k14$CGIap_source,"_",k14$TATA_target)
k15=final%>%group_by(CGInap_source,CGI_target)%>%dplyr::summarise(count=n(),group="CGInap_CGI")
k15$label=paste0(k15$CGInap_source,"_",k15$CGI_target)
k16=final%>%group_by(CGInap_source,TATA_target)%>%dplyr::summarise(count=n(),group="CGInap_TATA")
k16$label=paste0(k16$CGInap_source,"_",k16$TATA_target)
k17=final%>%group_by(Null_source,CGI_target)%>%dplyr::summarise(count=n(),group="Null_CGI")
k17$label=paste0(k17$Null_source,"_",k17$CGI_target)
k18=final%>%group_by(Null_source,TATA_target)%>%dplyr::summarise(count=n(),group="Null_TATA")
k18$label=paste0(k18$Null_source,"_",k18$TATA_target)
k19=final%>%group_by(TATA_source,CGI_target)%>%dplyr::summarise(count=n(),group="TATA_CGI")
k19$label=paste0(k19$TATA_source,"_",k19$CGI_target)
k20=final%>%group_by(TATA_source,TATA_target)%>%dplyr::summarise(count=n(),group="TATA_TATA")
k20$label=paste0(k20$TATA_source,"_",k20$TATA_target)

k21=rbind(k13,k14,k15,k16,k17,k18,k19,k20)
k22=spread(k21[3:5],key=3,value=1)

for(i in 1:nrow(k22)){
  GSEATasting <- matrix(c(k22$Yes_Yes[i], k22$Yes_No[i], k22$No_Yes[i], k22$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  k22$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  k22$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
k22$pv=paste0("p = ",signif(k22$p.val,3))
write.table(k22, gzfile(paste0(path_fig6_data,"hiC_result_allcell_ex5_cluster_FE_structure.tsv.gz")), col.names=T, row.names = F, quote = F, sep="\t")

#result located in [primary_folder]/fig6/data
# -> for fig. ex8f



