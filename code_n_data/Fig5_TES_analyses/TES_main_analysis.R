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

#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig4_data=paste0(primary_folder,"fig4/data/")
path_fig5_data=paste0(primary_folder,"fig5/data/")
RNAfold_path=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/RNAfold/")
CGI_path=paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/")

#===============================================================================
#link TESID to ex3_cluster and ex5_cluster

# prepare extended TES bed file for mapping CGI, MYC
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
transcriptTES1=unique(table5[,c("n3_string","TESID","TEScount")])
transcriptTES1=transcriptTES1%>%group_by(n3_string)%>%slice_max(TEScount, n=1) #one TESID with highest read count from one ex3_cluster

transcriptTES1$chr=sapply(strsplit(transcriptTES1$TESID,"_"),"[",1)
transcriptTES1$start=as.numeric(sapply(strsplit(transcriptTES1$TESID,"_"),"[",2))-5000
transcriptTES1$end=as.numeric(sapply(strsplit(transcriptTES1$TESID,"_"),"[",3))+5000
transcriptTES1$strand=sapply(strsplit(transcriptTES1$TESID,"_"),"[",4)
write.table(transcriptTES1[order(transcriptTES1$start, transcriptTES1$start),c(4,5,6,2,3,7)],gzfile(paste0(CGI_path,"Neuron_THP1.S3.TES.table5.5000bp_extend.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

table5b=left_join(unique(table5[,c("n5_string","TESID","TEScount","polyA")]),data1[,c("n5_string","ex5cluster_class","downstream_CpG_island","CpGTATA")], by="n5_string", copy=F)
table5b=table5b[which(table5b$TESID %in% transcriptTES1$TESID),]
table5b=table5b[which(!is.na(table5b$ex5cluster_class)),]
table5b=table5b%>%group_by(TESID,TEScount,polyA)%>%dplyr::summarise(ex5cluster_class=paste(unique(ex5cluster_class),collapse=";"), downstream_CpG_island=paste(unique(downstream_CpG_island),collapse=";"),CpGTATA=paste(unique(CpGTATA),collapse=";"), n5_string=paste(unique(n5_string),collapse=";"),count=n())
table5b=table5b[-grep(";",table5b$ex5cluster_class),]
table5b=table5b[-grep(";",table5b$CpGTATA),]
table5b=table5b[-grep(";",table5b$downstream_CpG_island),]
#add back n3_string
table5b=left_join(table5b,unique(table5[,c("TESID","n3_string")]),by="TESID", copy=F)
write.table(table5b,gzfile(paste0(CGI_path,"TESID_restricted_to_n3_string.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

# remove 1105 zero count TES (they are detected known TES from known GENCODE transcript but the read TES located somewhere else)
table5b=table5b[which(!is.na(table5b$TEScount)),]
write.table(table5b,gzfile(paste0(path_TES,"CGI/TESID_restricted_to_n3_string.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#84955

#===============================================================================
#get cell type CTES count from transcript bed
path_t=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1_S3/bam_to_bed/bed/")
files=list.files(path=path_t, pattern=".bed.bgz")
files=files[-grep("tbi",files)]
files.names=gsub(".bed.bgz","",files)
for (i in 1:length(files)){
  data=read.delim(paste0(path_t,files[i]), header=F, stringsAsFactors = F)[,c(1:6)]
  data$V2[which(data$V6=="+")]=data$V3[which(data$V6=="+")]-1
  data$V3[which(data$V6=="-")]=data$V2[which(data$V6=="-")]+1
  data1=data%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(V5=n())
  data1$V4=paste0(data1$V1,"_",data1$V2,"_",data1$V3,"_",data1$V6)
  write.table(data1[order(data1$V1, data1$V2),c(1,2,3,6,5,4)],gzfile(paste0(path_t,"CTES/",files.names[i],".CTES.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)}

path_tes=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/bam_to_bed/bed/CTES/")
files=list.files(path=path_tes, pattern=".bed.gz")
need=c("iPSC","NSC","Neuron","DMSO","PMA")
for (i in 1:length(need)){
  files1=files[grep(need[i],files)]
  data=read.delim(paste0(path_tes,files1[1]), header=F, stringsAsFactors = F)
  for (j in 2: length(files1)){
    data1=read.delim(paste0(path_tes,files1[j]), header=F, stringsAsFactors = F)
    data=rbind(data,data1)
    data=data%>%group_by(V1,V2,V3,V4,V6)%>%dplyr::summarise(V5=sum(V5))
    data=data[,c(1,2,3,4,6,5)]}
  write.table(data[order(data$V1, data$V2),],gzfile(paste0(path_tes,need[i],".CTES.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)}

files=list.files(path=path_tes, pattern=".bed.gz")
files=files[-grep("rep",files)]
files.names=gsub(".CTES.bed.gz","",files)
data=read.delim(paste0(path_tes,files[1]), header=F, stringsAsFactors = F)
data=data[,c(4,5)]
for (i in 2: length(files)){
  data1=read.delim(paste0(path_tes,files[i]), header=F, stringsAsFactors = F)
  data=full_join(data, data1[,c(4,5)], by="V4",copy=F)}
colnames(data)[c(2:6)]=files.names
colnames(data)[c(2,6)]=c("THP1","dTHP1")
data[is.na(data)]=0

table5b=read.delim(paste0(CGI_path,"TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5b=left_join(table5b, data[,c(1,3,5,4,2,6)], by=c("TESID"="V4"), copy=F)

write.table(table5b,gzfile(paste0(CGI_path,"TESID_restricted_to_n3_string.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare TES region for RNAfold
table5b=read.delim(paste0(CGI_path,"TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

table5b=table5b[which(!is.na(table5b$TEScount)),]
table5b$V1=sapply(strsplit(table5b$TESID,"_"),"[",1)
table5b$V2=as.numeric(sapply(strsplit(table5b$TESID,"_"),"[",2))
table5b$V3=as.numeric(sapply(strsplit(table5b$TESID,"_"),"[",3))
table5b$V6=sapply(strsplit(table5b$TESID,"_"),"[",4)

#extend the 3' end upstream 200 downstream 200 bp
options(scipen=999)
table5b$V2=table5b$V2-200
table5b$V3=table5b$V3+200
table5b=table5b[order(table5b$V1, table5b$V2),]

table5.3nbedeRNAa=table5b[which(table5b$ex5cluster_class == "e_ncRNA" & table5b$polyA == "No" & table5b$TEScount>=3),]
table5.3nbedeRNAb=table5b[which(table5b$ex5cluster_class == "e_ncRNA" & table5b$polyA == "No" & table5b$TEScount<3),]
table5.3nbedeRNAc=table5b[which(table5b$ex5cluster_class == "e_ncRNA" & table5b$polyA == "Yes" & table5b$TEScount>=3),]
table5.3nbedeRNAd=table5b[which(table5b$ex5cluster_class == "e_ncRNA" & table5b$polyA == "Yes" & table5b$TEScount<3),]
table5.3nbedpRNAa=table5b[which(table5b$ex5cluster_class == "p_ncRNA" & table5b$polyA == "No" & table5b$TEScount>=3),]
table5.3nbedpRNAb=table5b[which(table5b$ex5cluster_class == "p_ncRNA" & table5b$polyA == "No" & table5b$TEScount<3),]
table5.3nbedpRNAc=table5b[which(table5b$ex5cluster_class == "p_ncRNA" & table5b$polyA == "Yes" & table5b$TEScount>=3),]
table5.3nbedpRNAd=table5b[which(table5b$ex5cluster_class == "p_ncRNA" & table5b$polyA == "Yes" & table5b$TEScount<3),]
table5.3nbedmRNAa=table5b[which(table5b$ex5cluster_class == "mRNA" & table5b$polyA == "No" & table5b$TEScount>=3),]
table5.3nbedmRNAb=table5b[which(table5b$ex5cluster_class == "mRNA" & table5b$polyA == "No" & table5b$TEScount<3),]
table5.3nbedmRNAc=table5b[which(table5b$ex5cluster_class == "mRNA" & table5b$polyA == "Yes" & table5b$TEScount>=3),]
table5.3nbedmRNAd=table5b[which(table5b$ex5cluster_class == "mRNA" & table5b$polyA == "Yes" & table5b$TEScount<3),]

setwd(RNAfold_path)
write.table(table5.3nbedeRNAa[,c(15,16,17,1,2,18)], gzfile("nonployA_eRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedeRNAc[,c(15,16,17,1,2,18)], gzfile("ployA_eRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedpRNAa[,c(15,16,17,1,2,18)], gzfile("nonployA_pRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedpRNAc[,c(15,16,17,1,2,18)], gzfile("ployA_pRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedmRNAa[,c(15,16,17,1,2,18)], gzfile("nonployA_mRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedmRNAc[,c(15,16,17,1,2,18)], gzfile("ployA_mRNA_3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

write.table(table5.3nbedeRNAb[,c(15,16,17,1,2,18)], gzfile("nonployA_eRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedeRNAd[,c(15,16,17,1,2,18)], gzfile("ployA_eRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedpRNAb[,c(15,16,17,1,2,18)], gzfile("nonployA_pRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedpRNAd[,c(15,16,17,1,2,18)], gzfile("ployA_pRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedmRNAb[,c(15,16,17,1,2,18)], gzfile("nonployA_mRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table5.3nbedmRNAd[,c(15,16,17,1,2,18)], gzfile("ployA_mRNA_less3read_up200_down200.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#getfasta for RNAfold
setwd(RNAfold_path)
system("for file in *.bed.gz; do bedtools getfasta -s -nameOnly -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed \"$file\" | gzip > \"${file%.bed.gz}.fasta.gz\"; done")
# fasta.gz transferred to RNAfold_path/output

#RNAfold
setwd(RNAfold_path)
system("Rscript RNAfold_run.R")
# this script run RNAfold, collapse the probability matrix (_dp.ps) into a single file, calculate depletion score and output the structure string

#===============================================================================
# process the RNAfold output file
# structure score per position in different class

gather1=data.frame(matrix(nrow=0,ncol=4))
colnames(gather1)=c("score","nt","event","group")

files=list.files(path=paste0(RNAfold_path,"output/"), pattern="prob_matrix.rds")
files.names=gsub("_full_prob_matrix.rds","",files)

for (j in 1: length(files)){
  RNAfold=readRDS(paste0(RNAfold_path,"output/", files[j]))
  RNAfold1=as.data.frame(colMeans(RNAfold))
  RNAfold1$nt=c(-200:200)
  RNAfold1$event=nrow(RNAfold)
  RNAfold1$group=files.names[j]
  colnames(RNAfold1)[1]="score"
  gather1=rbind(gather1,RNAfold1)}

gather1$group1=sapply(strsplit(gather1$group,"_"),"[",1)
gather1$group2=sapply(strsplit(gather1$group,"_"),"[",2)
gather1$group3=sapply(strsplit(gather1$group,"_"),"[",3)
gather1$group2[which(gather1$group2 != "mRNA")]=gsub("RNA","_ncRNA", gather1$group2[which(gather1$group2 != "mRNA")])
gather1$group1=gsub("non","non-",gather1$group1)
gather1$group1=gsub("ploy","poly",gather1$group1)
gather1$group3=gsub("less3read","< 3 reads", gather1$group3)
gather1$group3=gsub("3read",">= 3 reads", gather1$group3)
write.table(gather1,gzfile(paste0(path_fig5_data,"RNAfold_collapse_6groups.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
gather1=read.delim(paste0(path_fig5_data,"RNAfold_collapse_6groups.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#stored in [primary_folder]/fig5/data
#for -> Fig 5g

#===============================================================================
# depletion score for all the non-poly(A) RNAs
files=list.files(path=paste0(RNAfold_path,"output/"), pattern="structural_scoring.tsv.gz")
files.names=gsub("_structural_scoring.tsv.gz","",files)

both2=data.frame(matrix(nrow=0, ncol=7))
colnames(both2)=c("id","mean_body_P","mean_valley_P","depletion_score","group1","group2","group3")

for (i in 1:length(files)){
  both=read.delim(paste0(RNAfold_path,"output/",files[i]), header=T, stringsAsFactors = F, check.names = F)
  both$group1=sapply(strsplit(files.names[i],"_"),"[",1)
  both$group2=sapply(strsplit(files.names[i],"_"),"[",2)
  both$group3=sapply(strsplit(files.names[i],"_"),"[",3)
  colnames(both)=colnames(both2)
  both2=rbind(both2,both)}

both2$group2[which(both2$group2 != "mRNA")]=gsub("RNA","_ncRNA",both2$group2[which(both2$group2 != "mRNA")])
both2$body_05="Yes"
both2$body_05[which(both2$mean_body_P < 0.5)]="No"
both2$drop="Yes"
both2$drop[which(both2$depletion_score > 0.85)]="No"
both2$structual_depletion="No"
both2$structual_depletion[which(both2$body_05 == "Yes" & both2$drop =="Yes")]="Yes"
both2%>%group_by(group1,group2,group3)%>%dplyr::summarise(depleted=length(which(structual_depletion=="Yes")),count=n(), percent=length(which(structual_depletion=="Yes"))/n())
both2$id=sapply(strsplit(both2$id,"\\("),"[",1)

table5b=read.delim(paste0(CGI_path,"TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
both2=left_join(both2, table5b[,c(1,4,6,5,8:14)], by=c("id"="TESID"),copy=F)

write.table(both2,gzfile(paste0(path_fig5_data,"eRNA_structural_depleteion.raw.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
both2=read.delim(paste0(path_fig5_data,"eRNA_structural_depleteion.raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#stored in [primary_folder]/fig5/data
#for -> fig. Ext7g & table S14

#===============================================================================
# seq.composition per position in different class
files=list.files(path=paste0(RNAfold_path,"output/"), pattern=".fasta.gz")
files.names=gsub("_up200_down200.fasta.gz","",files)
aRNAfold4=data.frame(matrix(nrow=0, ncol=6))
colnames(aRNAfold4)=c("A","U","G","C","nt","group")
seq=data.frame(matrix(nrow=0, ncol=3))
colnames(seq)=c("id","last10","group")
aRNAfold4_r=aRNAfold4
for (j in 1: length(files)){
  aRNAfold=read.delim(paste0(RNAfold_path,"output/", files[j]), header=F, stringsAsFactors = F)
  aRNAfold2=data.frame(cbind(aRNAfold$V1[seq(1, nrow(aRNAfold), by=2)], aRNAfold$V1[seq(2, nrow(aRNAfold),by=2)]))
  aRNAfold2$X1=gsub(">","",aRNAfold2$X1)
  aRNAfold2$X1=sapply(strsplit(aRNAfold2$X1,"\\("),"[",1)
  for (i in 1:401){aRNAfold2[,i+2]=substr(aRNAfold2$X2, start = i, stop = i)}
  aRNAfold3A=data.frame(colSums(aRNAfold2[,c(3:403)]=="A")/nrow(aRNAfold2))
  aRNAfold3T=data.frame(colSums(aRNAfold2[,c(3:403)]=="T")/nrow(aRNAfold2))
  aRNAfold3G=data.frame(colSums(aRNAfold2[,c(3:403)]=="G")/nrow(aRNAfold2))
  aRNAfold3C=data.frame(colSums(aRNAfold2[,c(3:403)]=="C")/nrow(aRNAfold2))
  aRNAfold3=cbind(aRNAfold3A,aRNAfold3T,aRNAfold3G,aRNAfold3C)
  colnames(aRNAfold3)=c("A","U","G","C")
  aRNAfold3$nt=(-200 : 200)
  aRNAfold3$group=files.names[j]
  aRNAfold4=rbind(aRNAfold4,aRNAfold3)
  aRNAfold2$last10=substr(aRNAfold2$X2, start = 192, stop = 201)
  aRNAfold2=aRNAfold2[,c(1,404)]
  aRNAfold2$group=files.names[j]
  colnames(aRNAfold2)=colnames(seq)
  seq=rbind(seq,aRNAfold2)}

write.table(aRNAfold4,gzfile(paste0(path_fig5_data,"polyA_all_eRNA_ATGC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(seq,gzfile(paste0(path_fig5_data,"polyA_all_eRNA_last10nt.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig5/data
#for -> fig.5f

#===============================================================================
# break down into CpGTATA class
# structure score per position in different CpGTATA
table5b=read.delim(paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5ee=table5b[which(table5b$ex5cluster_class == "e_ncRNA"),]
k1=unique(table5ee[,c("TESID","CpGTATA")])

gather1=data.frame()

files=list.files(path=paste0(RNAfold_path,"output/"), pattern="prob_matrix.rds")
files=files[grep("nonployA_eRNA",files)]
files.names=gsub("_full_prob_matrix.rds","",files)

for (j in 1: length(files)){
  RNAfold=readRDS(paste0(RNAfold_path,"output/", files[j]))
  RNAfold=data.frame(RNAfold)
  RNAfold$TESID=sapply(strsplit(rownames(RNAfold),"\\("),"[",1)
  RNAfold=left_join(RNAfold,k1,by="TESID",copy=F)
  RNAfold1 <- RNAfold %>% group_by(CpGTATA) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), count=n())
  RNAfold2=reshape2::melt(RNAfold1[which(!is.na(RNAfold1$CpGTATA)),],id=c(1,403))
  RNAfold2$variable=gsub("X","",as.character(RNAfold2$variable))
  RNAfold2$group=files.names[j]
  colnames(RNAfold2)[c(3,4)]=c("nt","score")
  gather1=rbind(gather1,RNAfold2)}

gather1$group1=sapply(strsplit(gather1$group,"_"),"[",1)
gather1$group2=sapply(strsplit(gather1$group,"_"),"[",2)
gather1$group3=sapply(strsplit(gather1$group,"_"),"[",3)
gather1$group2[which(gather1$group2 != "mRNA")]=gsub("RNA","_ncRNA", gather1$group2[which(gather1$group2 != "mRNA")])
gather1$group1=gsub("non","non-",gather1$group1)
gather1$group1=gsub("ploy","poly",gather1$group1)
gather1$group3=gsub("less3read","< 3 reads", gather1$group3)
gather1$group3=gsub("3read",">= 3 reads", gather1$group3)
write.table(gather1,gzfile(paste0(path_fig5_data,"RNAfold_collapse_eRNACpGTATA.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# seq.composition per position in different CpGTATA
table5b=read.delim(paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5ee=table5b[which(table5b$ex5cluster_class == "e_ncRNA"),]
k1=unique(table5ee[,c("TESID","CpGTATA")])

aRNAfold=read.delim(paste0(RNAfold_path,"output/nonployA_eRNA_3read_up200_down200.fasta.gz"), header=F, stringsAsFactors = F)
aRNAfold2=data.frame(cbind(aRNAfold$V1[seq(1, nrow(aRNAfold), by=2)], aRNAfold$V1[seq(2, nrow(aRNAfold),by=2)]))
aRNAfold2$X1=gsub(">","",aRNAfold2$X1)
aRNAfold2$X1=sapply(strsplit(aRNAfold2$X1,"\\("),"[",1)
aRNAfold2=left_join(aRNAfold2, k1, by=c("X1"="TESID"),copy=F)
for (i in 1:401){aRNAfold2[,i+3]=substr(aRNAfold2$X2, start = i, stop = i)}
aRNAfold3A <- aRNAfold2[, c(3:404)] %>% group_by(CpGTATA) %>% summarise(across(everything(), ~ sum(.x == "A", na.rm = TRUE) / n()), count = n(), group="A")
aRNAfold3T <- aRNAfold2[, c(3:404)] %>% group_by(CpGTATA) %>% summarise(across(everything(), ~ sum(.x == "T", na.rm = TRUE) / n()), count = n(), group="U")
aRNAfold3G <- aRNAfold2[, c(3:404)] %>% group_by(CpGTATA) %>% summarise(across(everything(), ~ sum(.x == "G", na.rm = TRUE) / n()), count = n(), group="G")
aRNAfold3C <- aRNAfold2[, c(3:404)] %>% group_by(CpGTATA) %>% summarise(across(everything(), ~ sum(.x == "C", na.rm = TRUE) / n()), count = n(), group="C")
aRNAfold3=rbind(aRNAfold3A,aRNAfold3T,aRNAfold3G,aRNAfold3C)
colnames(aRNAfold3)[c(2:402)]=-200 : 200
aRNAfold4=reshape2::melt(aRNAfold3,id=c(1,403,404))
aRNAfold4=aRNAfold4[which(!is.na(aRNAfold4$CpGTATA)),]

write.table(aRNAfold4,gzfile(paste0(path_fig5_data,"polyA_all_eRNACpGTATA_ATGC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
# TES hit promoter 

table5b=read.delim(paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/TESID_restricted_to_n3_string.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5b$chr=sapply(strsplit(table5b$TESID,"_"),"[",1)
table5b$start=as.numeric(sapply(strsplit(table5b$TESID,"_"),"[",2))
table5b$end=as.numeric(sapply(strsplit(table5b$TESID,"_"),"[",3))
table5b$strand=sapply(strsplit(table5b$TESID,"_"),"[",4)
write.table(table5b[order(table5b$chr,table5b$start),c(15:17,1,2,18)],gzfile(paste0(primary_folder,"code_n_data/Fig5_TES_analyses/CGI/TESID_restricted_to_n3_string.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)

#========
#bash
setwd(CGI_path)
system(paste0("bedtools closest -a TESID_restricted_to_n3_string.bed.gz -b ",primary_folder,"code_n_data/n5_regions/GRCh38-PLS.all.promoter.sort.bed.gz -D a | gzip > TES_PLS_cCRE_distance.bed.gz"))

#========
setwd(CGI_path)
TES_pcluster=read.delim("TES_PLS_cCRE_distance.bed.gz",header=F, stringsAsFactors = F, check.names=F)
TES_pcluster=left_join(TES_pcluster, table5b[,c(1,3,4,5,6,7,9:14)], by=c("V4"="TESID"),copy=F)
TES_pcluster1=TES_pcluster[which(TES_pcluster$V5>=3),]
TES_pcluster1=unique(TES_pcluster1[, c(1:6, 13:24)])
TES_pcluster1$promoter_hinder="Not Reach"
TES_pcluster1$promoter_hinder[which(TES_pcluster1$V13>=0 & TES_pcluster1$V13< 500)]="End"
TES_pcluster1$promoter_hinder[which(TES_pcluster1$V13<0 )]="Through"

#========
# only allow downstream closest
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/"))
system(paste0("bedtools closest -a enhancer_ex5_cluster.bed.gz -b ",primary_folder,"code_n_data/n5_regions/GRCh38-PLS.all.promoter.sort.bed.gz -D a -iu | gzip > enhancer_downstream_promoter_distance_cCRE.bed.gz"))

#========
setwd(CGI_path)
pCREdistance=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/enhancer_downstream_promoter_distance_cCRE.bed.gz"),header=F, stringsAsFactors = F, check.names=F)
pCREdistance=pCREdistance%>%group_by(V4)%>%slice_max(V13)
pCREdistance1=pCREdistance[which(pCREdistance$V13>=0 & pCREdistance$V13<2000),]

TES_pcluster1$dpromoter="No"
TES_pcluster1$dpromoter[which(TES_pcluster1$n5_string %in% pCREdistance1$V4)]="Yes"
TES_pcluster1$dpromoter[grep(";",TES_pcluster1$n5_string)]=NA

TES_pcluster2=TES_pcluster1[which(TES_pcluster1$ex5cluster_class %in% c("e_ncRNA") & TES_pcluster1$CpGTATA %in% c("CGIap","CGInap")),]%>%group_by(ex5cluster_class,CpGTATA,polyA,dpromoter,promoter_hinder)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count))
TES_pcluster2$dPLS="others"
TES_pcluster2$dPLS[which(TES_pcluster2$dpromoter=="Yes")]="dPLS"
TES_pcluster2$label=paste0(TES_pcluster2$CpGTATA,"_",TES_pcluster2$dPLS)
TES_pcluster2$polyA[which(TES_pcluster2$polyA == "Yes")]="Poly(A)"
TES_pcluster2$polyA[which(TES_pcluster2$polyA == "No")]="Non-poly(A)"
write.table(TES_pcluster2,gzfile(paste0(path_fig5_data,"enhancer_downstream_promoter_distance_cCRE.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# for fig ex7i

#===============================================================================
# MYC and exosome sensitive
myc_result=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

quant_RLE.tran=read.delim(paste0(primary_folder,"code_n_data/Fig4_transcription_features/exosome_sensitivity/ExoKD.transcript.RLE.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
quant_RLE.tran$model_ID=rownames(quant_RLE.tran)
table5=left_join(table5,quant_RLE.tran[which(quant_RLE.tran$detection=="Yes"),c(10,7)], by="model_ID", copy=F)

quant_RLE.gene=read.delim(paste0(primary_folder,"code_n_data/Fig4_transcription_features/exosome_sensitivity/ExoKD.gene.RLE.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
quant_RLE.gene$T4_gene_ID=rownames(quant_RLE.gene)
table5=left_join(table5,quant_RLE.gene[which(quant_RLE.gene$detection=="Yes"),c(10,7)], by="T4_gene_ID", copy=F, suffix=c("_Tx","_Gene"))
write.table(table5, gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz")),col.names=T, row.names=F, sep="\t",quote=F)

Txn3string=table5[which(!is.na(table5$exo_sensitivity_Tx)),]%>%group_by(n3_string)%>%dplyr::summarise(exo_sensitivity_Tx=sum(exo_sensitivity_Tx * (full_qry_count+partial_qry_count))/sum(full_qry_count+partial_qry_count), count=n() )
myc_result=left_join(myc_result, Txn3string[,c(1,2)], by="n3_string", copy=F)
myc_result <- myc_result %>% group_by(ex5cluster_class,polyA) %>% mutate(
  Exo_sensitive = if_else(exo_sensitivity_Tx > median(exo_sensitivity_Tx, na.rm = TRUE), "Yes", "No")) %>% ungroup()

myc_result$TESrecur="No"
myc_result$TESrecur[which(myc_result$TEScount>=3)]="Yes"
write.table(myc_result,gzfile(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig5/data
#for -> fig.Ext7l


#combine different features
#===============================================================================
# add CpG island hit at TES
#==================
#bedtools
#bash
setwd(path_myc)
system(paste0("bedtools closest -a Neuron_THP1.S3.TES.table5.1bp.bed.gz -b ",CGI_path,"cpgIslandExt.hg38.main_chr.bed.gz -D a | gzip > CGI_Neuron_THP1.S3.TES.table5.1bp.bed.gz"))
#==================
myc_result=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CPG=read.delim("CGI_Neuron_THP1.S3.TES.table5.1bp.bed.gz", header=F, stringsAsFactors = F, check.names = F)
need=unique(CPG$V4[which(CPG$V17==0 & CPG$V12>=40)]) # at least 40 unit
myc_result$CPG_TES0="No"
myc_result$CPG_TES0[which(myc_result$TESID %in% need)]="Yes"
length(which(myc_result$CPG_TES0 == "Yes")) #7399
write.table(myc_result,gzfile(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# add structural depletion 
both2=read.delim(paste0(path_fig5_data,"eRNA_structural_depleteion.raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(both2)[9]="TES_structual_depletion"
myc_result=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
myc_result=left_join(myc_result, both2[,c("id","TES_structual_depletion")], by=c("TESID"="id"), copy=F)

#===============================================================================
# add sequence of the last 10 nt
myc_result=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
seq=read.delim(paste0(path_fig5_data,"polyA_all_eRNA_last10nt.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(seq)[2]="seq_last10"
myc_result=left_join(myc_result,seq[,c(1,2)],by=c("TESID"="id"),copy=F)
myc_result$C_rich=NA
myc_result$C_rich[which(str_count(myc_result$seq_last10, "C")<4)]="No"
myc_result$C_rich[which(str_count(myc_result$seq_last10, "C")>=4)]="Yes"
write.table(myc_result,gzfile(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================






