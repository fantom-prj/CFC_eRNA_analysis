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
path_fig3_data=paste0(primary_folder,"fig3/data/")
SQANTI_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SQANTI/")

#===============================================================================
# RNA class by SQANTI3
setwd(SQANTI_path)
system("sh assemblies_vs_reference.sh")

#-> intermediate output files were removed
#-> gzip classification.txt

#===============================================================================
# ENST identification
TALON_t1=read.delim(paste0(primary_folder,"code_n_data/other_assemblers/TALON/gtf/TALON.table1.noPI.detectedENST.bed12.bed.gz"), header=F, stringsAsFactors = F)
iso_t1=read.delim(paste0(primary_folder,"code_n_data/other_assemblers/isoQuant/assemblies_stranded_corrected_internnal_priming_filtered/isoquant_sensitive_int_priming.bed12.bed.gz"), header=F, stringsAsFactors = F)
SALA_default=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_default/default/transcript/log/full_length_support_matrix.tsv.gz"), header=T, stringsAsFactors = F)
SALA_default=SALA_default[grep("ENS",SALA_default$model_ID),]
SALA_default=SALA_default[which(rowSums(SALA_default[,c(2:34)])>0),]

TALON_t1t=TALON_t1$V4[grep("ENST",TALON_t1$V4)]
SALA_t2t=SALA_default$model_ID
iso_t1t=iso_t1$V4[grep("ENST",iso_t1$V4)]
all=rbind(cbind(TALON_t1t,"TALON"), cbind(SALA_t2t,"SALA"), cbind(iso_t1t,"IsoQuant"))
colnames(all)=c("geneID","analysis")
write.table(all,gzfile(paste0(path_fig3_data,"venn_ENST_SALA_TALON_Iso.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext4e

#===============================================================================
# pick the raw SALAdefault for ENST comparison
# for -> Fig. Ext4e
transcript_info1a=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_default/default/transcript/log/Neuron_THP1.S3.table0_raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
partialENST=transcript_info1a$model_ID[which(transcript_info1a$ref_source=="partial_ref")]
raw=fread(paste0(SQANTI_path,"SALA_default_raw/SALA_default_raw_classification.txt.gz"), header=T, select=c(1:20))
all=read.delim(paste0(path_fig3_data,"venn_ENST_SALA_TALON_Iso.tsv.gz"), header=T, stringsAsFactors = F, check.names=F)
length(all$geneID[which(all$analysis == "SALA")])
length(all$geneID[which(all$analysis == "TALON")])
length(all$geneID[which(all$analysis == "IsoQuant")])

need=setdiff(all$geneID[which(all$analysis == "IsoQuant")], all$geneID[which(all$analysis == "SALA")])
need1=setdiff(need, all$geneID[which(all$analysis == "TALON")])
need2=setdiff(need, need1)
need3=setdiff(all$geneID[which(all$analysis == "TALON")], all$geneID[which(all$analysis == "SALA")])
need3=setdiff(need3,need2)

missENST=data.frame(rbind(cbind(need1,"IsoQuant alone"),cbind(need2,"IsoQuant & TALON"),cbind(need3,"TALON alone")))

raw1=raw[which(raw$associated_transcript %in% c(need,need3)),]
raw2=raw1[-grep("ENST",raw1$isoform),]
raw2$AbsSumEnd=abs(raw2$diff_to_TSS)+abs(raw2$diff_to_TTS)
raw3=raw2[which(raw2$structural_category == "full-splice_match"),]
raw3a=raw3%>%group_by(associated_transcript)%>%slice_min(AbsSumEnd)
raw4=raw2[-which(raw2$associated_transcript %in% unique(raw3$associated_transcript)),]
raw4a=raw4%>%group_by(associated_transcript)%>%slice_min(AbsSumEnd, n=1)

raw5=rbind(unique(raw3a[,c("associated_transcript","structural_category","AbsSumEnd")]),unique(raw4a[,c("associated_transcript","structural_category","AbsSumEnd")]))
missENST=left_join(missENST, raw5, by=c("need1"="associated_transcript"), copy=F)

colnames(missENST)[c(1,2)]=c("transcriptID","group")
missENST$SALA_partial="Absent"
missENST$SALA_partial[which(missENST$transcriptID %in% partialENST)]="Partial"
missENST%>%group_by(group, SALA_partial,structural_category)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
write.table(missENST,gzfile(paste0(path_fig3_data,"venn_ENST_SALA_missed.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
files=list.files(path=SQANTI_path, pattern="classification.txt.gz", recursive = T)
files.names=sapply(strsplit(files, "\\/"),"[",1)
data2=data.frame()
for (i in 1: length(files)){
  data=read.delim(paste0(SQANTI_path,files[i]), header=T, stringsAsFactors = F, check.names = F)
  data=data[-grep("ENS",data$isoform),]
  data1=data%>%group_by(structural_category)%>%dplyr::summarise(count=n())
  data1$label=files.names[i]
  data2=rbind(data2, data1)}
data3=spread(data2,key=3, value=2)
data3[is.na(data3)]=0
write.table(data3, gzfile(paste0(path_fig3_data,"SQANTI3.txt.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext4f

#===============================================================================


#transcript models intersect with CRE regions
#===============================================================================
#steps for intersect located in /code_n_data/SALA/SALA_Final.R, /code_n_data/SALA/SALA_Default.R, /code_n_data/other_assebmblers/talon.process.R, /code_n_data/other_assebmblers/isoquant.process.R
#===============================================================================
path1=paste0(primary_folder,"code_n_data/SALA/")
path2=paste0(primary_folder,"code_n_data/other_assemblers/TALON/gtf/")
path3=paste0(primary_folder,"code_n_data/other_assemblers/isoQuant/assemblies_stranded_corrected_internnal_priming_filtered/")

#====counting -> for all SALA 3 tables, TALON 2 tables, isoquant 2 tables
CRE.hit.rate=data.frame(matrix(nrow=48,ncol=7))
colnames(CRE.hit.rate)=c("platform","group","name","Tx.hit","Tx.total","novel.Tx.hit","novel.Tx.total")
CRE.hit.rate$group=c(rep("P",8),rep("PE",8),rep("cCRE",8),rep("ATAC",8),rep("F5_CAGE",8),rep("TSS_cluster",8))
group1=c("P","PE","cCRE","ATAC","F5_CAGE","TSS_cluster")
datab1=data.frame(matrix(nrow=0, ncol=3))
colnames(datab1)=c("V7","count", "name")

filesP=list.files(path=path1, pattern=".N5.Pcount.bed", recursive = T)
filesE=list.files(path=path1, pattern=".N5.Ecount.bed", recursive = T)
filesC=list.files(path=path1, pattern=".N5.Ccount.bed", recursive = T)
filesA=list.files(path=path1, pattern=".N5.ATACcount.bed", recursive = T)
filesCAGE=list.files(path=path1, pattern=".N5.CAGEcount.bed", recursive = T)
filesCluster=list.files(path=path1, pattern=".N5.SCAFEcount.bed", recursive = T)
names=sapply(strsplit(filesP,"\\/S"),"[", 2)
names=paste0("S",sapply(strsplit(names,"\\."),"[", 1))
for (i in 1:length(filesP)){
  dataP=read.delim(paste0(path1,filesP[i]), header=F, stringsAsFactors = F, check.names = F)
  dataE=read.delim(paste0(path1,filesE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataC=read.delim(paste0(path1,filesC[i]), header=F, stringsAsFactors = F, check.names = F)
  dataA=read.delim(paste0(path1,filesA[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCAGE=read.delim(paste0(path1,filesCAGE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCluster=read.delim(paste0(path1,filesCluster[i]), header=F, stringsAsFactors = F, check.names = F)
  
  data=cbind(dataP,dataE[,8],dataC[,8],dataA[,8],dataCAGE[,8],dataCluster[,8])
  datab=data%>%group_by(V7)%>%dplyr::summarise(count=n())
  datab$name=names[i]
  datab1=rbind(datab1, datab)
  colnames(data)[c(8:13)]=c("P","E","C","ATAC","F5_CAGE","TSS_cluster")
  data$cCRE=0
  data$cCRE[which(rowSums(data[,c(8:10)]>0)>0)]=1
  data$PE=0
  data$PE[which(rowSums(data[,c(8:9)]>0)>0)]=1
  data=data[which(data$V1 != "chrM"),]
  for (j in 1:length(group1)){
    data1=data
    colnames(data1)[which(colnames(data1)==group1[j])]="key"
    data2=data1[which(data1$V7 == "Transcript_from_novel_gene"),]
    CRE.hit.rate$platform[(j*8)-8+i] = "SALA"
    CRE.hit.rate$name[(j*8)-8+i] = names[i] 
    CRE.hit.rate$Tx.hit[(j*8)-8+i] = length(which(data1$key > 0))
    CRE.hit.rate$Tx.total[(j*8)-8+i] = nrow(data1)
    CRE.hit.rate$novel.Tx.hit[(j*8)-8+i] = length(which(data2$key > 0))
    CRE.hit.rate$novel.Tx.total[(j*8)-8+i] = nrow(data2)}}

filesP=list.files(path=path2, pattern=".n5.Pcount.bed")
filesE=list.files(path=path2, pattern=".n5.Ecount.bed")
filesC=list.files(path=path2, pattern=".n5.Ccount.bed")
filesA=list.files(path=path2, pattern=".n5.ATACcount.bed")
filesCAGE=list.files(path=path2, pattern=".n5.CAGEcount.bed")
filesCluster=list.files(path=path2, pattern=".n5.SCAFEcount.bed")
names=sapply(strsplit(filesP,"\\."),"[", 1)
for (i in 1:length(filesP)){
  dataP=read.delim(paste0(path2,filesP[i]), header=F, stringsAsFactors = F, check.names = F)
  dataE=read.delim(paste0(path2,filesE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataC=read.delim(paste0(path2,filesC[i]), header=F, stringsAsFactors = F, check.names = F)
  dataA=read.delim(paste0(path2,filesA[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCAGE=read.delim(paste0(path2,filesCAGE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCluster=read.delim(paste0(path2,filesCluster[i]), header=F, stringsAsFactors = F, check.names = F)
  
  data=cbind(dataP,dataE[,8],dataC[,8],dataA[,8],dataCAGE[,8],dataCluster[,8])
  datab=data%>%group_by(V7)%>%dplyr::summarise(count=n())
  datab$name=names[i]
  datab1=rbind(datab1, datab)
  colnames(data)[c(8:13)]=c("P","E","C","ATAC","F5_CAGE","TSS_cluster")
  data$cCRE=0
  data$cCRE[which(rowSums(data[,c(8:10)]>0)>0)]=1
  data$PE=0
  data$PE[which(rowSums(data[,c(8:9)]>0)>0)]=1
  data=data[which(data$V1 != "chrM"),]
  for (j in 1:length(group1)){
    data1=data
    colnames(data1)[which(colnames(data1)==group1[j])]="key"
    data2=data1[which(data1$V7 == "Transcript_from_novel_gene"),]
    CRE.hit.rate$platform[(j*8)-8+i+4] = "TALON" 
    CRE.hit.rate$name[(j*8)-8+i+4] = names[i] 
    CRE.hit.rate$Tx.hit[(j*8)-8+i+4] = length(which(data1$key > 0))
    CRE.hit.rate$Tx.total[(j*8)-8+i+4] = nrow(data1)
    CRE.hit.rate$novel.Tx.hit[(j*8)-8+i+4] = length(which(data2$key > 0))
    CRE.hit.rate$novel.Tx.total[(j*8)-8+i+4] = nrow(data2)}}

filesP=list.files(path=path3, pattern=".n5.Pcount.bed")
filesE=list.files(path=path3, pattern=".n5.Ecount.bed")
filesC=list.files(path=path3, pattern=".n5.Ccount.bed")
filesA=list.files(path=path3, pattern=".n5.ATACcount.bed")
filesCAGE=list.files(path=path3, pattern=".n5.CAGEcount.bed")
filesCluster=list.files(path=path3, pattern=".n5.SCAFEcount.bed")
names=sapply(strsplit(filesP,"\\."),"[", 1)
for (i in 1:length(filesP)){
  dataP=read.delim(paste0(path3,filesP[i]), header=F, stringsAsFactors = F, check.names = F)
  dataE=read.delim(paste0(path3,filesE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataC=read.delim(paste0(path3,filesC[i]), header=F, stringsAsFactors = F, check.names = F)
  dataA=read.delim(paste0(path3,filesA[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCAGE=read.delim(paste0(path3,filesCAGE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCluster=read.delim(paste0(path3,filesCluster[i]), header=F, stringsAsFactors = F, check.names = F)
  
  data=cbind(dataP,dataE[,8],dataC[,8],dataA[,8],dataCAGE[,8],dataCluster[,8])
  datab=data%>%group_by(V7)%>%dplyr::summarise(count=n())
  datab$name=names[i]
  datab1=rbind(datab1, datab)
  colnames(data)[c(8:13)]=c("P","E","C","ATAC","F5_CAGE","TSS_cluster")
  data$cCRE=0
  data$cCRE[which(rowSums(data[,c(8:10)]>0)>0)]=1
  data$PE=0
  data$PE[which(rowSums(data[,c(8:9)]>0)>0)]=1
  data=data[which(data$V1 != "chrM"),]
  for (j in 1:length(group1)){
    data1=data
    colnames(data1)[which(colnames(data1)==group1[j])]="key"
    data2=data1[which(data1$V7 == "Transcript_from_novel_gene"),]
    CRE.hit.rate$platform[(j*8)-8+i+6] = "IsoQuant"
    CRE.hit.rate$name[(j*8)-8+i+6] = names[i] 
    CRE.hit.rate$Tx.hit[(j*8)-8+i+6] = length(which(data1$key > 0))
    CRE.hit.rate$Tx.total[(j*8)-8+i+6] = nrow(data1)
    CRE.hit.rate$novel.Tx.hit[(j*8)-8+i+6] = length(which(data2$key > 0))
    CRE.hit.rate$novel.Tx.total[(j*8)-8+i+6] = nrow(data2)}}

CRE.hit.rate$rate=CRE.hit.rate$Tx.hit/CRE.hit.rate$Tx.total
CRE.hit.rate$novel.rate=CRE.hit.rate$novel.Tx.hit/CRE.hit.rate$novel.Tx.total
write.table(datab1, gzfile(paste0(path_fig3_data,"transcript.group.count.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(CRE.hit.rate, gzfile(paste0(path_fig3_data,"transcript.hit.rate.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext4g & h


