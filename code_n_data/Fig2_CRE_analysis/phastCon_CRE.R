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

#===============================
##PhastCon intersect
#bash

#download from ucsc: 
# file size too big not included in the folder
# download 30way as an example:

rsync -avz --progress \ rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons30way/ /analysisdata/fantom6/Interactome/resources/hg38_phastcon/
~/bigWigToBedGraph hg38.phastCons30way.bw hg38.phastCons30way.bed
gzip hg38.phastCons30way.bed

#=============
setwd(Phast_path)
CREsummit=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz"), header=F, stringsAsFactors = F)
CREsummit$V2=CREsummit$V2-2000
CREsummit$V3=CREsummit$V3+2000
CREsummit$V2[which(CREsummit$V2 <0)]=0
write.table(CREsummit, gzfile("ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#=============
#bash
#bedtools intersect

setwd(Phast_path)
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons4way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons4way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons7way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons7way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons17way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons17way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons20way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons20way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons30way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons30way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons100way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons100way.bed.gz")
system("bedtools intersect -wa -wb -a /analysisdata/fantom6/Interactome/resources/hg38_phastcon/hg38.phastCons470way.bed.gz -b ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.bed.gz | cut -f-4,8-8,10-10 | gzip >  ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.phastCons470way.bed.gz")
#===============

CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
CREsummit=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz"), header=F, stringsAsFactors = F)

#=========================================
files=list.files(path=Phast_path, pattern="ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend")
files.names=sapply(strsplit(files,"phastCons"),"[",2)
files.names=gsub(".bed.gz","",files.names)
for (j in 1:length(files)){
  con100=fread(paste0(Phast_path,files[j]), header=F)
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
    con100a=left_join(con100a,CREsummit[,c(4,3)], by=c("V5"="V4"), copy=F, suffix=c("","_summit"))
    con100a$positionV3=con100a$V3-con100a$V3_summit
    con100a$positionV3[which(con100a$V6=="-")]=-(con100a$V3[which(con100a$V6=="-")]-con100a$V3_summit[which(con100a$V6=="-")])
    con100a=con100a[which(abs(con100a$positionV3)<=2000),]
    con100a=left_join(con100a,CREanno[,c("CREID","promoter_type")],by=c("V5"="CREID"),copy=F)
    
    con100a2=con100a%>%group_by(promoter_type,V5)%>%dplyr::summarise(mean_4001=mean(phastCon_score),max_4001=max(phastCon_score),count_4001=n())
    con100a3=con100a%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_4001 = mean(phastCon_score))
    con100a4=con100a%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_4001 = mean(phastCon_score))
    
    con100a5=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(V5)%>%dplyr::summarise(mean_2001=mean(phastCon_score),max_2001=max(phastCon_score),count_2001=n())
    con100a6=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_2001 = mean(phastCon_score))
    con100a7=con100a[which(abs(con100a$positionV3)<=1000),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_2001 = mean(phastCon_score))
    
    con100a8=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(V5)%>%dplyr::summarise(mean_1001=mean(phastCon_score),max_1001=max(phastCon_score),count_1001=n())
    con100a9=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_1001 = mean(phastCon_score))
    con100a10=con100a[which(abs(con100a$positionV3)<=500),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_1001 = mean(phastCon_score))
    
    con100a11=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(V5)%>%dplyr::summarise(mean_up500=mean(phastCon_score),max_up500=max(phastCon_score),count_up500=n())
    con100a12=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_up500 = mean(phastCon_score))
    con100a13=con100a[which(con100a$positionV3<0 & con100a$positionV3>=(-500)),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_up500 = mean(phastCon_score))
    
    con100a14=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(V5)%>%dplyr::summarise(mean_down500=mean(phastCon_score),max_down500=max(phastCon_score),count_down500=n())
    con100a15=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 100)%>%summarise(mean_top_100_down500 = mean(phastCon_score))
    con100a16=con100a[which(con100a$positionV3>=0 & con100a$positionV3<500),]%>%group_by(V5)%>%arrange(desc(phastCon_score), .by_group = TRUE)%>%slice_head(n = 200)%>%summarise(mean_top_200_down500 = mean(phastCon_score))
    
    con100a2=left_join(con100a2,con100a3, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a4, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a5, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a6, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a7, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a8, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a9, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a10, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a11, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a12, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a13, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a14, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a15, by="V5",copy=F)
    con100a2=left_join(con100a2,con100a16, by="V5",copy=F)
    
    CRE1=con100a2$V5[which(con100a2$count_4001 >=3950)]
    con100a1=con100a[which(con100a$V5 %in% CRE1)]%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=mean(phastCon_score),count=n())
    
    write.table(con100a,gzfile(paste0(Phast_path,need[i],"_CRE_phastcon",files.names[j],"_breakdown.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
    write.table(con100a1,gzfile(paste0(Phast_path,need[i],"_CRE_phastcon",files.names[j],"_plot.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
    write.table(con100a2,gzfile(paste0(Phast_path,need[i],"_CRE_phastcon",files.names[j],"_summary.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)
  }}

#==
CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
CREanno=CREanno[which(CREanno$representative == "Yes"),] #take major strand alone

files=list.files(path=Phast_path, pattern="_breakdown.tsv.gz")
need=unique(sapply(strsplit(files,"_"),"[",3))
for (i in 1:length(need)){
  files1=files[grep(need[i],files)]
  con100a=fread(paste0(Phast_path,files1[1]), header=T)
  con100a=con100a[which(con100a$V5 %in% CREanno$CREID),]
  con100a2=con100a%>%group_by(promoter_type,V5)%>%dplyr::summarise(count_4001=n())
  CRE1=con100a2$V5[which(con100a2$count_4001 >=3950)]
  con100a1=con100a[which(con100a$V5 %in% CRE1)]%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=mean(V4),count=n())
  for (j in 2: length(files1)){
    con100aa=fread(paste0(Phast_path,files1[j]), header=T)
    con100aa=con100aa[which(con100aa$V5 %in% CREanno$CREID),]
    con100a2=con100aa%>%group_by(promoter_type,V5)%>%dplyr::summarise(count_4001=n())
    CRE1=con100a2$V5[which(con100a2$count_4001 >=3950)]
    con100aa1=con100aa[which(con100aa$V5 %in% CRE1)]%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=mean(V4),count=n())
    con100a1=rbind(con100a1,con100aa1)}
  con100a1=con100a1%>%group_by(promoter_type,positionV3)%>%dplyr::mutate(total=sum(count))
  con100a1a=con100a1%>%group_by(promoter_type,positionV3)%>%dplyr::summarise(score=sum(score*count/total))
  write.table(con100a1a,gzfile(paste0(Phast_path,"tCRE_nooverlap_",need[i],"_plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)}

files=list.files(path=Phast_path, pattern="tCRE_nooverlap_")
names=gsub("_plot.tsv","",files)
names=gsub("tCRE_nooverlap_","",names)
plotdata=read.delim(paste0(Phast_path,files[1]),header=T, stringsAsFactors = F, check.names = F)
colnames(plotdata)[1]="variable"
plotdata$signalID=names[1]
for (i in 2: length(files)){
  plotdata1=read.delim(paste0(Phast_path,files[i]),header=T, stringsAsFactors = F, check.names = F)
  colnames(plotdata1)[1]="variable"
  plotdata1$signalID=names[i]
  plotdata=rbind(plotdata,plotdata1)}
plotdata$group="phastCon"
plotdata$signalID=paste0(gsub("phastcon","",plotdata$signalID),"way")

write.table(plotdata,gzfile(paste0(path_fig2_data,"tCRE_nooverlap_phastcon_all_way_plot.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===
#other plot - bidirectional - non-overlap
CREanno=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv", header=T, stringsAsFactors = F, check.names = F)
CREanno=CREanno[which(CREanno$representative == "Yes"),]

files=list.files(path=Phast_path, pattern="_breakdown.tsv.gz")
need=unique(sapply(strsplit(files,"_"),"[",3))
for (i in 1:length(need)){
  files1=files[grep(need[i],files)]
  con100a=fread(paste0(Phast_path,files1[1]), header=T)
  con100a=left_join(con100a, CREanno[,c(1,35)], by=c("V5"="CREID"),copy=F)
  con100a2=con100a%>%group_by(promoter_type,V5)%>%dplyr::summarise(mean_4001=mean(V4),max_4001=max(V4),count_4001=n())
  CRE1=con100a2$V5[which(con100a2$count_4001 >=3950)]
  con100a1=con100a[which(con100a$V5 %in% CRE1)]%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=mean(V4),count=n())
  for (j in 2: length(files1)){
    con100aa=fread(paste0(Phast_path,files1[j]), header=T)
    con100aa=left_join(con100aa, CREanno[,c(1,35)], by=c("V5"="CREID"),copy=F)
    con100a2=con100aa%>%group_by(promoter_type,V5)%>%dplyr::summarise(mean_4001=mean(V4),max_4001=max(V4),count_4001=n())
    CRE1=con100a2$V5[which(con100a2$count_4001 >=3950)]
    con100aa1=con100aa[which(con100aa$V5 %in% CRE1)]%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=mean(V4),count=n())
    con100a1=rbind(con100a1,con100aa1)}
  con100a1=con100a1%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::mutate(total=sum(count))
  con100a1a=con100a1%>%group_by(promoter_type,orientation,positionV3)%>%dplyr::summarise(score=sum(score*count/total))
  write.table(con100a1a,gzfile(paste0(Phast_path,"CRE_nooverlap_bidirection_",need[i],"_plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)}

files=list.files(path=Phast_path, pattern="CRE_nooverlap_bidirection_")
names=gsub("_plot.tsv.gz","",files)
names=gsub("CRE_nooverlap_","",names)
plotdata=read.delim(paste0(Phast_path,files[1]),header=T, stringsAsFactors = F, check.names = F)
colnames(plotdata)[2]="variable"
plotdata$signalID=names[1]
for (i in 2: length(files)){
  plotdata1=read.delim(paste0(Phast_path,files[i]),header=T, stringsAsFactors = F, check.names = F)
  colnames(plotdata1)[2]="variable"
  plotdata1$signalID=names[i]
  plotdata=rbind(plotdata,plotdata1)}
plotdata$group="phastCon"
plotdata=plotdata[which(plotdata$variable != "Others"),]
plotdata=plotdata[which(!is.na(plotdata$variable)),]
plotdata$source=sapply(strsplit(plotdata$signalID,"_"),"[",1)
plotdata$signalID=sapply(strsplit(plotdata$signalID,"_"),"[",2)
plotdata$signalID=paste0(gsub("phastcon","",plotdata$signalID),"way")
write.table(plotdata,gzfile(paste0(path_fig2_data,"tCRE_nooverlap_phastcon_directionality_all_way_plot.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#remove intermediate files
setwd(Phast_path)
system("rm ontCAGE.Neuron_THP1.CRE.stranded_summit_2kb_extend.*.bed.gz")
system("rm *_breakdown.tsv.gz")




