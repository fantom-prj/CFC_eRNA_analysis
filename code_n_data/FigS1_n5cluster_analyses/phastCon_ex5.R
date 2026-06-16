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
Phast_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/phastCon/")
Phastex5_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/phastCon/")

#===============================
##PhastCon intersect
#bash

#download from ucsc: 
# file size too big not included in the folder
# download 30way as an example:

rsync -avz --progress \ rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons30way/ /analysisdata/fantom6/Interactome/resources/hg38_phastcon/
  ~/bigWigToBedGraph hg38.phastCons30way.bw hg38.phastCons30way.bed
gzip hg38.phastCons30way.bed

# resource file not included, please download from UCSC
#=============
setwd(Phastex5_path)

cluster_summit=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), header=F, stringsAsFactors = F)
cluster_summit1=cluster_summit
cluster_summit1$V2=cluster_summit1$V2-2000
cluster_summit1$V3=cluster_summit1$V3+2000
cluster_summit1$V2[which(cluster_summit1$V2 <0)]=0
write.table(cluster_summit1[order(cluster_summit1$V1,cluster_summit1$V2),], gzfile("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#=============
#bash
#bedtools intersect

setwd(Phastex5_path)
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons4way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons4way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons7way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons7way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons17way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons17way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons20way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons20way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons30way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons30way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons100way.bed.gz -b Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  Neuron_THP1.S3.end5.summit.table5_2kb_extend.phastCons100way.bed.gz")

#===============
cluster_summit=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), header=F, stringsAsFactors = F)

#=========================================
files=list.files(path=Phastex5_path, pattern="Neuron_THP1.S3.end5.summit.table5_2kb_extend.gz")
files.names=sapply(strsplit(files,"phastCons"),"[",2)
files.names=gsub(".bed.gz","",files.names)
for (j in 1:length(files)){
  con100=fread(paste0(Phastex5_path,files[j]), header=F)
  need=unique(con100$V1)
  for (i in 1:length(need)){
    con100a=con100[which(con100$V1==need[i]),]
    con100a$length=con100a$V3-con100a$V2
    con100b=con100a[which(con100a$length>1),]
    con100b=con100b%>%group_by(V5,V2)%>%dplyr::mutate(V3=paste((V2+1):V3, collapse=";"))
    con100b=con100b%>%group_by(V5,V2)%>%dplyr::mutate(V4=paste(rep(V4,length), collapse=";"))
    con100bex=separate_rows(con100b[,c(3:6)],V3, V4, sep=";")
    con100a=rbind(con100a[which(con100a$length==1),c(3:6)],con100bex)
    con100a$V3=as.numeric(con100a$V3)
    con100a$V4=as.numeric(con100a$V4)
    con100a=left_join(con100a,cluster_summit[,c(4,3)], by=c("V5"="V4"), copy=F, suffix=c("","_summit"))
    con100a$positionV3=con100a$V3-con100a$V3_summit
    con100a$positionV3[which(con100a$V6=="-")]=-(con100a$V3[which(con100a$V6=="-")]-con100a$V3_summit[which(con100a$V6=="-")])
    con100a=con100a[which(abs(con100a$positionV3)<=2000),]
    colnames(con100a)[c(2,3)]=c("phastCon_score","n5_string")
    con100a=left_join(con100a, unique(table5[,c(62,64)]),by="n5_string",copy=F)
    con100a=left_join(con100a,CREanno[,c(1,32)],by="CREID",copy=F)
    
    con100a2=con100a%>%group_by(promoter_type,n5_string)%>%dplyr::summarise(mean_4001=mean(phastCon_score),max_4001=max(phastCon_score),count_4001=n())
    con100a3=con100a%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_4001 = mean(phastCon_score))
    con100a4=con100a%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_4001 = mean(phastCon_score))
    
    con100a5=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(n5_string)%>%dplyr::summarise(mean_2001=mean(phastCon_score),max_2001=max(phastCon_score),count_2001=n())
    con100a6=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_2001 = mean(phastCon_score))
    con100a7=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_2001 = mean(phastCon_score))
    
    con100a8=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(n5_string)%>%dplyr::summarise(mean_1001=mean(phastCon_score),max_1001=max(phastCon_score),count_1001=n())
    con100a9=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_1001 = mean(phastCon_score))
    con100a10=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_1001 = mean(phastCon_score))
    
    con100a11=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(n5_string)%>%dplyr::summarise(mean_up500=mean(phastCon_score),max_up500=max(phastCon_score),count_up500=n())
    con100a12=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_up500 = mean(phastCon_score))
    con100a13=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_up500 = mean(phastCon_score))
    
    con100a14=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(n5_string)%>%dplyr::summarise(mean_down500=mean(phastCon_score),max_down500=max(phastCon_score),count_down500=n())
    con100a15=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_down500 = mean(phastCon_score))
    con100a16=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(n5_string)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_down500 = mean(phastCon_score))
    
    con100a2=left_join(con100a2,con100a3, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a4, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a5, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a6, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a7, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a8, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a9, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a10, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a11, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a12, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a13, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a14, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a15, by="n5_string",copy=F)
    con100a2=left_join(con100a2,con100a16, by="n5_string",copy=F)
    
    CRE1=con100a2$n5_string[which(con100a2$count_4001 >=3950)]
    con100a1=con100a[which(con100a$n5_string %in% CRE1)]%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=mean(phastCon_score),count=n())
    
    write.table(con100a,gzfile(paste0(Phastex5_path,need[i],"_end5_cluster_phastcon",files.names[j],"_breakdown.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
    write.table(con100a1,gzfile(paste0(Phastex5_path,need[i],"_end5_cluster_phastcon",files.names[j],"_plot.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
    write.table(con100a2,gzfile(paste0(Phastex5_path,need[i],"_end5_cluster_phastcon",files.names[j],"_summary.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
  }}
#==
files=list.files(path=Phastex5_path, pattern="plot.tsv")
index=c("4way","7way","17way","20way","30way","100way")
files1=files[grep(index[1],files)]
con100a1=read.delim(paste0(Phastex5_path,files1[1]), header=T, stringsAsFactors = F, check.names = F)
for (i in 2:length(files1)){
  con100a1a=read.delim(paste0(Phastex5_path,files1[i]), header=T, stringsAsFactors = F, check.names = F)
  con100a1=rbind(con100a1,con100a1a)}
con100a1=con100a1%>%group_by(promoter_type,positionV3)%>%dplyr::mutate(total=sum(count))
con100a1a=con100a1%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=sum(score*count/total))
con100a1a$label=index[1]
for (j in 2:6){
  files1=files[grep(index[j],files)]
  con1000a1=read.delim(paste0(Phastex5_path,files1[1]), header=T, stringsAsFactors = F, check.names = F)
  for (i in 2:length(files1)){
    con1000a1a=read.delim(paste0(Phastex5_path,files1[i]), header=T, stringsAsFactors = F, check.names = F)
    con1000a1=rbind(con1000a1,con1000a1a)}
  con1000a1=con1000a1%>%group_by(promoter_type,positionV3)%>%dplyr::mutate(total=sum(count))
  con1000a1a=con1000a1%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=sum(score*count/total))
  con1000a1a$label=index[j]
  con100a1a=rbind(con100a1a,con1000a1a)}
write.table(con100a1a,gzfile(paste0(Phastex5_path,"end5cluster_phastconALL_plot.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)

files=list.files(path=Phastex5_path, pattern="summary.tsv.gz")
index=c("4way","7way","17way","20way","30way","100way")
files1=files[grep(index[1],files)]
con100a2=read.delim(paste0(Phastex5_path,files1[1]),header=T, stringsAsFactors = F, check.names = F)
for (i in 2:length(files1)){
  con100a2a=read.delim(paste0(Phastex5_path,files1[i]),header=T, stringsAsFactors = F, check.names = F)
  con100a2=rbind(con100a2,con100a2a)}
con100a2$label=index[1]
for (j in 2:6){
  files1=files[grep(index[j],files)]
  con1000a2=read.delim(paste0(Phastex5_path,files1[1]),header=T, stringsAsFactors = F, check.names = F)
  for (i in 2:length(files1)){
    con1000a2a=read.delim(paste0(Phastex5_path,files1[i]),header=T, stringsAsFactors = F, check.names = F)
    con1000a2=rbind(con1000a2,con1000a2a)}
  con1000a2$label=index[j]
  con100a2=rbind(con100a2,con1000a2)}
write.table(con100a2,gzfile(paste0(Phastex5_path,"end5cluster_phastconALL_summary.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)



#=====================================
#other plot - bidirectional
files=list.files(path=Phastex5_path, pattern="_breakdown.tsv.gz")
need=unique(sapply(strsplit(files,"_"),"[",4))
for (i in 1:length(need)){
  files1=files[grep(need[i],files)]
  con100a=fread(paste0(Phastex5_path,files1[1]), header=T)
  con100a=left_join(con100a, CRE[,c(1,35)], by="CREID",copy=F)
  con100a2=con100a%>%group_by(promoter_type,n5_string)%>%dplyr::summarise(count_4001=n())
  CRE1=con100a2$n5_string[which(con100a2$count_4001 >=3950)]
  con100a1=con100a[which(con100a$n5_string %in% CRE1)]%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=mean(phastCon_score),count=n())
  for (j in 2: length(files1)){
    con100aa=fread(paste0(Phastex5_path,files1[j]), header=T)
    con100aa=left_join(con100aa, CRE[,c(1,35)], by="CREID",copy=F)
    con100a2=con100aa%>%group_by(promoter_type,n5_string)%>%dplyr::summarise(count_4001=n())
    CRE1=con100a2$n5_string[which(con100a2$count_4001 >=3950)]
    con100aa1=con100aa[which(con100aa$n5_string %in% CRE1)]%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=mean(phastCon_score),count=n())
    con100a1=rbind(con100a1,con100aa1)}
  con100a1=con100a1%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::mutate(total=sum(count))
  con100a1a=con100a1%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=sum(score*count/total))
  write.table(con100a1a,gzfile(paste0(Phastex5_path,"end5cluster_bidirection_",need[i],"_plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)}

#bi-directional
files=list.files(path=Phastex5_path, pattern="bidirection")
names=c("100way","17way","20way","30way","4way","7way")
plotdata=read.delim(paste0(Phastex5_path,files[1]),header=T, stringsAsFactors = F, check.names = F)
plotdata$signalID=names[1]
for (i in 2: length(files)){
  plotdata1=read.delim(paste0(Phastex5_path,files[i]),header=T, stringsAsFactors = F, check.names = F)
  plotdata1$signalID=names[i]
  plotdata=rbind(plotdata,plotdata1)}
plotdata$group="phastCon_directionality"
plotdata=plotdata[which(plotdata$orientation != "Others"),]
write.table(plotdata,gzfile(paste0(path_fig2_data,"end5cluster_phastcon_directionality_all_way_plot.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#remove intermediate files
setwd(Phastex5_path)
system("rm Neuron_THP1.S3.end5.summit.table5_2kb_extend.*.bed.gz")
system("rm *_breakdown.tsv.gz")




