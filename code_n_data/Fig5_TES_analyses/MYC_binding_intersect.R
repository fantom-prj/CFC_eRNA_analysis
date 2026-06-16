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
path_fig5_data=paste0(primary_folder,"fig5/data/")
path_myc=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/chipseq_MYC_humanES/")
path_TES=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/")
CGItes_path=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/")

# intersect with myc chip-seq from human ES
#=====================================================
#version one, using the defined peaks from GSM1505809
options(scipen=999)

transcriptTES1=read.delim(paste0(path_TES,"CGI/Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz"),header=F, stringsAsFactors = F)
transcriptTES1$V2=as.numeric(sapply(strsplit(transcriptTES1$V4,"_"),"[",2))
transcriptTES1$V3=as.numeric(sapply(strsplit(transcriptTES1$V4,"_"),"[",3))

write.table(transcriptTES1[order(transcriptTES1$V1, transcriptTES1$V2),],gzfile(paste0(path_myc,"Neuron_THP1.S3.TES.table5.1bp.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

myc=read.delim(paste0(path_myc,"GSM1505809_cMyc_110812_h64.bed.peak.txt"), header=F, stringsAsFactors = F)
myc$V1=paste0("chr",myc$V1)
myc$V6="."
write.table(myc[order(myc$V1,myc$V2),],paste0(path_myc,"GSM1505809_cMyc_110812_h64.bed.peak.bed"), col.names=F, row.names=F, sep="\t", quote=F)

#==================
#liftover
#bash
setwd(path_myc)
system("liftOver GSM1505809_cMyc_110812_h64.bed.peak.bed hg19ToHg38.over.chain.gz GSM1505809_cMyc_110812_h64.bed.peak.hg38.bed unmapped.bed")
#==================
myc=read.delim("GSM1505809_cMyc_110812_h64.bed.peak.hg38.bed", header=F, stringsAsFactors = F)
write.table(myc[order(myc$V1,myc$V2),],"GSM1505809_cMyc_110812_h64.bed.peak.hg38.bed", col.names=F, row.names=F, sep="\t", quote=F)

#==================
#bedtools
#bash
setwd(path_myc)
system("bedtools closest -a Neuron_THP1.S3.TES.table5.1bp.bed.gz -b GSM1505809_cMyc_110812_h64.bed.peak.hg38.bed -D a | gzip > myc_chip_hES_Neuron_THP1.S3.TES.table5.1bp.bed.gz")

#==================
table5b=read.delim(paste0(path_TES,"CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
myc_result=read.delim(paste0(path_myc,"myc_chip_hES_Neuron_THP1.S3.TES.table5.1bp.bed.gz"), header=F, stringsAsFactors = F)
myc_result$MYC_TES0="No"
myc_result$MYC_TES0[which(myc_result$V13==0)]="Yes"
length(which(myc_result$MYC_TES0=="Yes"))
myc_result=myc_result%>%group_by(V4)%>%dplyr::slice_max(V11)
myc_result=myc_result%>%group_by(V4)%>%dplyr::slice_max(V10)
myc_result=left_join(table5b, myc_result,  by=c("TESID"="V4"),copy=F)
length(unique(myc_result$TESID)) #84955
write.table(myc_result,gzfile(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#TES position with human ES MYC chip-seq data
#=============================================
#bedtools intersect
#bash
setwd(paste0(path_myc,"position/"))
system("bedtools intersect -wb -a Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz -b ../GSM1505809_cMyc_110812_h64.bed.peak.hg38.bed | gzip > Neuron_THP1.S3.TES.table5.5000bp_extend.MYC_chip.bed.gz")
#=============================================
#R
options(scipen=999)

setwd(paste0(path_myc,"position/"))
transcriptTES1=read.delim(paste0(path_TES,"CGI/Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz"),header=F, stringsAsFactors = F)
myc=read.delim("Neuron_THP1.S3.TES.table5.5000bp_extend.MYC_chip.bed.gz", header=F, stringsAsFactors = F)
length(unique(myc$V4))#43755
myc=left_join(myc, transcriptTES1[,c(2,4)], by="V4", copy=F, suffix=c("","_ori"))
myc$locS=myc$V2-myc$V2_ori-5000
myc$locE=myc$V3-myc$V2_ori-5000
myc1=myc[which(myc$V6 == "+" & myc$V2==myc$V8),c(1,14,15,4,10,11)]
myc2=myc[which(myc$V6 == "-" & myc$V3==myc$V9),c(1,15,14,4,10,11)]
myc2$locE=myc2$locE * (-1)
myc2$locS=myc2$locS * (-1)
myc2$V1="chr1"
myc1$V1="chr1"
colnames(myc2)=colnames(myc1)
myc1$locS=myc1$locS+5000
myc1$locE=myc1$locE+5000
myc2$locS=myc2$locS+5001
myc2$locE=myc2$locE+5001
myc=rbind(myc1,myc2)
colnames(myc)[c(4,5,6)]=c("TESID","name","score")

table5b=read.delim(paste0(path_TES,"CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5a2=table5b[which(table5b$TEScount >= 3),]
table5a2[which(table5a2$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA")),]%>%group_by(ex5cluster_class, polyA)%>%dplyr::summarise(count=n())
table5a2[which(table5a2$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA")),]%>%group_by(ex5cluster_class, polyA, downstream_CpG_island)%>%dplyr::summarise(count=n())
#this file determine the total number
myc=left_join(myc, table5a2[,c(1:4,6)],by="TESID",copy=F)
myc=myc[which(!is.na(myc$ex5cluster_class)),]
write.table(myc[order(myc$locS),],gzfile("Neuron_THP1.S3.TES.table5.5000bp_extend.MYC_chip.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

myc=read.delim("Neuron_THP1.S3.TES.table5.5000bp_extend.MYC_chip.fakebed.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
mycmA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "mRNA")]
mycmN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "mRNA")]
myceA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "e_ncRNA")]
myceN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA")]
myceN1=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "CGIap")]
myceN2=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "CGInap")]
myceN3=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "Null")]
myceN4=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "e_ncRNA" & table5a2$CpGTATA == "TATA")]
mycpA=table5a2$TESID[which(table5a2$polyA == "Yes" & table5a2$ex5cluster_class == "p_ncRNA")]
mycpN=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "p_ncRNA")]
mycpN1=table5a2$TESID[which(table5a2$polyA == "No" & table5a2$ex5cluster_class == "p_ncRNA" & table5a2$CpGTATA == "CGI")]

need=c(50,100,150)            
path1=paste0(path_myc,"position/")
for (i in 1:length(need)){
  mycmAa=myc[which(myc$score >= need[i] & myc$TESID %in% mycmA),]
  mycmNa=myc[which(myc$score >= need[i] & myc$TESID %in% mycmN),]
  myceAa=myc[which(myc$score >= need[i] & myc$TESID %in% myceA),]
  myceNa=myc[which(myc$score >= need[i] & myc$TESID %in% myceN),]
  myceN1a=myc[which(myc$score >= need[i] & myc$TESID %in% myceN1),]
  myceN2a=myc[which(myc$score >= need[i] & myc$TESID %in% myceN2),]
  myceN3a=myc[which(myc$score >= need[i] & myc$TESID %in% myceN3),]
  myceN4a=myc[which(myc$score >= need[i] & myc$TESID %in% myceN4),]
  mycpAa=myc[which(myc$score >= need[i] & myc$TESID %in% mycpA),]
  mycpNa=myc[which(myc$score >= need[i] & myc$TESID %in% mycpN),]
  mycpN1a=myc[which(myc$score >= need[i] & myc$TESID %in% mycpN1),]

  path2=paste0(path1,"n",need[i],"/")
  write.table(mycmAa[order(mycmAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.mycmA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(mycmNa[order(mycmNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.mycmN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceAa[order(myceAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceNa[order(myceNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceN1a[order(myceN1a$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceN1.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceN2a[order(myceN2a$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceN2.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceN3a[order(myceN3a$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceN3.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(myceN4a[order(myceN4a$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.myceN4.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(mycpAa[order(mycpAa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.mycpA.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(mycpNa[order(mycpNa$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.mycpN.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(mycpN1a[order(mycpN1a$locS),c(1:4)],gzfile(paste0(path2,"TES_5000b_extend.mycpN1.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)}

fake2=cbind("chr1",c(0:10000),c(1:10001),c(-5000:5000))
write.table(fake2,gzfile(paste0(path1,"location10001.fakebed.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#================================
#bedtools count
setwd(paste0(path_myc,"position/n50/"))
system("for file in TES_5000b_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ../location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

setwd(paste0(path_myc,"position/n100/"))
system("for file in TES_5000b_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ../location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

setwd(paste0(path_myc,"position/n150/"))
system("for file in TES_5000b_extend*.fakebed.bed.gz; do bedtools intersect -wa -a ../location10001.fakebed.bed.gz -b \"$file\" -c | gzip > \"${file%.fakebed.bed.gz}.result.bed.gz\"; done")

#================================
#R
path1=paste0(path_myc,"position/")
size1=c(length(myceA),length(myceN), length(myceN1), length(myceN2), length(myceN3), length(myceN4),length(mycmA), length(mycmN), length(mycpA), length(mycpN), length(mycpN1))
size1=rep(size1,3)
files=list.files(path=path1, pattern=".result.bed.gz", recursive=T)

files.index=sapply(strsplit(files, "\\/"),"[",1)
files.names=sapply(strsplit(files, "extend."),"[",2)
files.names=gsub(".result.bed.gz","",files.names)
data=read.delim(paste0(path1, files[1]), header=F, stringsAsFactors = F)
data$group=files.names[1]
data$signalID=files.index[1]
data$V5=data$V5/size1[1]
for (i in 2:length(files)){
  data1=read.delim(paste0(path1, files[i]), header=F, stringsAsFactors = F)
  data1$group=files.names[i]
  data1$signalID=files.index[i]
  data1$V5=data1$V5/size1[i]
  data=rbind(data,data1)}

data$anno_region="mRNA_pA"
data$anno_region[grep("mycmN",data$group)]="mRNA_pN"
data$anno_region[grep("myceA",data$group)]="eRNA_pA"
data$anno_region[grep("myceN",data$group)]="eRNA_pN"
data$anno_region[grep("myceN1",data$group)]="eRNA_pN_CGIap"
data$anno_region[grep("myceN2",data$group)]="eRNA_pN_CGInap"
data$anno_region[grep("myceN3",data$group)]="eRNA_pN_Null"
data$anno_region[grep("myceN4",data$group)]="eRNA_pN_TATA"
data$anno_region[grep("mycpA",data$group)]="p_ncRNA_pA"
data$anno_region[grep("mycpN",data$group)]="p_ncRNA_pN"
data$anno_region[grep("mycpN1",data$group)]="p_ncRNA_pN_CGI"

data$group2="MYC_ChIP"
data$polyA="polyA"
data$polyA[grep("pN",data$anno_region)]="non-polyA"
data$ex5_cluster=sapply(strsplit(data$anno_region,"_p"),"[",1)
data$ex5_cluster[which(data$ex5_cluster=="eRNA")]="e_ncRNA"
write.table(data, gzfile(paste0(path_fig5_data,"TES_ChIP_MYC_position.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#result for ex7k
#===============================================================================















# NOT used in the manuscript
#===============================================================================
#use the union from 4 human ES cells from ChIP atlas, using fdr=1e-5
#merge
#bash
setwd(paste0(path_myc,"other_myc_chipseq/"))
system("cat ESC_H1.SRX080135.05.bed ESC_H1.SRX102975.05.bed ESC_H1.SRX150588.05.bed same_ESC_HUES6.SRX702150.05.bed | sort -k1,1 -k2,2n > combined.bed")
system("cut -f1-3 combined.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged.bed")
system("for file in *.05.bed; do bedtools intersect -wa -a merged.bed -b \"$file\" -c > \"${file%.bed}.count.bed\"; done")

#====================================
#count for each library
#R
path3=paste0(path_myc,"other_myc_chipseq/")
files=list.files(path=path3, pattern="count.bed")
files.names=gsub(".05.count.bed","",files)
files.names=gsub(".bed.peak.hg38.count.bed","",files.names)
data=read.delim(paste0(path3,files[1]),header=F, stringsAsFactors = F)
for (i in 2:length(files)){
  data1=read.delim(paste0(path3,files[i]),header=F, stringsAsFactors = F)
  data=cbind(data,data1[4])}
colnames(data)[c(4:10)]=files.names
data$n_sample=rowSums(data[,c(4:10)]>0)
data$n_sample6=rowSums(data[,c(4:6,8:10)]>0)
data$n_sample4=rowSums(data[,c(4:6,10)]>0)
data$width=data$V3-data$V2
data=data[which(nchar(data$V1)<=5),]

write.table(data,paste0(path3,"matrix.tsv"),row.names=F, col.names=T, sep="\t",quote=F)
data=read.delim(paste0(path3,"matrix.tsv"),header=T, stringsAsFactors = F, check.names = F)
write.table(data[which(data$n_sample4>0),c(1:3)], paste0(path3,"MYC4_union.bed"),col.names=F, row.names=F, sep="\t", quote=F)

#============================
#use the union from chipatlas for bedtools
#bash
setwd(paste0(path_myc,"other_myc_chipseq/"))
system("bedtools closest -a ../Neuron_THP1.S3.TES.table5.1bp.bed.gz -b MYC4_union.bed -D a> MYC4_union.bed_Neuron_THP1.S3.TES.table5.1bp.bed")
#system("for file in *.05.bed; do bedtools closest -a ../Neuron_THP1.S3.TES.table5.1bp.bed.gz -b \"$file\" -D a> \"${file%.05.bed}_Neuron_THP1.S3.TES.table5.1bp.bed\"; done")
#============================
#R
table5b=read.delim(paste0(path_TES,"CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
myc_result=read.delim(paste0(path3,"MYC4_union.bed_Neuron_THP1.S3.TES.table5.1bp.bed"), header=F, stringsAsFactors = F)
myc_result$up500="No"
myc_result$up500[which(myc_result$V6 == "+" & myc_result[,ncol(myc_result)-1]>= -500 & myc_result[,ncol(myc_result)-1]<= 0)]="Yes"
myc_result$up500[which(myc_result$V6 == "-" & myc_result[,ncol(myc_result)-1]<= 500 & myc_result[,ncol(myc_result)-1]>= 0)]="Yes"
length(which(myc_result$up500=="Yes")) #2599

myc_result=left_join(myc_result, table5b[,c(1:6)],by=c("V4"="TESID"),copy=F)
myc_result=myc_result[!is.na(myc_result$CpGTATA),]
myc_result=myc_result[!is.na(myc_result$V5),]
length(unique(myc_result$V4)) #84955

myc_result$label="Yes"
myc_resulta=spread(myc_result,key=16, value=17)
myc_resulta[is.na(myc_resulta)]="No"
myc_result1=unique(myc_resulta[,c(4,11:20)])%>%group_by(ex5cluster_class,downstream_CpG_island,up500)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
myc_result1$label=paste0(myc_result1$downstream_CpG_island,"_",myc_result1$up500)
myc_result2=spread(myc_result1[,c(1,6,4)], key=2, value=3)
myc_result2[is.na(myc_result2)]=0
for(i in 1:nrow(myc_result2)){
  GSEATasting <- matrix(c(myc_result2$No_No[i], myc_result2$No_Yes[i], myc_result2$Yes_No[i], myc_result2$Yes_Yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  myc_result2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  myc_result2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
write.table(myc_result2,gzfile(paste0(path_fig5_data,"union4_dCGI_cluster_myc_chipseq.FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

myc_result5=data.frame(matrix(nrow=0, ncol=8))
colnames(myc_result5)=c("ex5cluster_class","polyA","No_No","No_Yes","Yes_No","Yes_Yes","CpGTATA","Recursive")
need=c("TATA","Null","CGInap","CGIap")

for (j in 1:length(need)){
  myc_result3=unique(myc_resulta[,c(4,11:21)])%>%group_by(ex5cluster_class,polyA,!!sym(need[j]),up500)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
  myc_result3$label=paste0(myc_result3[[3]],"_",myc_result3$up500)
  myc_result4=spread(myc_result3[,c(1,2,7,5)], key=3, value=4)
  myc_result4[is.na(myc_result4)]=0
  myc_result4$CpGTATA=need[j]
  myc_result4$Recursive="Both"
  myc_result5=rbind(myc_result5,myc_result4)}

for (j in 1:length(need)){
  myc_result3=unique(myc_resulta[which(myc_resulta$V5 <3),c(4,11:21)])%>%group_by(ex5cluster_class,polyA,!!sym(need[j]),up500)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
  myc_result3$label=paste0(myc_result3[[3]],"_",myc_result3$up500)
  myc_result4=spread(myc_result3[,c(1,2,7,5)], key=3, value=4)
  myc_result4[is.na(myc_result4)]=0
  myc_result4$CpGTATA=need[j]
  myc_result4$Recursive="No"
  myc_result5=rbind(myc_result5,myc_result4)}

for (j in 1:length(need)){
  myc_result3=unique(myc_resulta[which(myc_resulta$V5 >=3),c(4,11:21)])%>%group_by(ex5cluster_class,polyA,!!sym(need[j]),up500)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
  myc_result3$label=paste0(myc_result3[[3]],"_",myc_result3$up500)
  myc_result4=spread(myc_result3[,c(1,2,7,5)], key=3, value=4)
  myc_result4[is.na(myc_result4)]=0
  myc_result4$CpGTATA=need[j]
  myc_result4$Recursive="Yes"
  myc_result5=rbind(myc_result5,myc_result4)}

for(i in 1:nrow(myc_result5)){
  GSEATasting <- matrix(c(myc_result5$No_No[i], myc_result5$No_Yes[i], myc_result5$Yes_No[i], myc_result5$Yes_Yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  myc_result5$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  myc_result5$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}

myc_result5$logOR=log(myc_result5$OR)
myc_result5$label="***"
myc_result5$label[which(myc_result5$p.val>=0.001)]="**"
myc_result5$label[which(myc_result5$p.val>=0.01)]="*"
myc_result5$label[which(myc_result5$p.val>=0.05)]="n.s."
myc_result5$label1=paste0("log(OR)=",signif(myc_result5$logOR,2),myc_result5$label)
write.table(myc_result5,gzfile(paste0(path_fig5_data,"union4_dCGI_cluster_polyA_myc_chipseq_up500.FE.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================



