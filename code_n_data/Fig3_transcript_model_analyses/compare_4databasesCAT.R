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
library(shadowtext)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#===============================================================================
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig3_data=paste0(primary_folder,"fig3/data/")
CAT_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_FANTOMCAT/")
CAT_path_log=paste0(CAT_path,"sala/transcript/Neuron_THP1_FCAT/log/")
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/")
output_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/")

#=============================
#prepare files from FANTOMCAT
FCAT_gtf=read.delim(paste0(CAT_path,"lv2_permissive.hg38.gtf.gz"), header=F, stringsAsFactors = F)
FCAT_gtfg=FCAT_gtf[which(FCAT_gtf$V3 == "gene"),]
FCAT_gtfg$geneID=sapply(strsplit(FCAT_gtfg$V9, "; "),"[",1)
FCAT_gtfg$geneClass=sapply(strsplit(FCAT_gtfg$V9, "; "),"[",3)
FCAT_gtfg$geneID=gsub("gene_id ","",FCAT_gtfg$geneID)
FCAT_gtfg$geneClass=gsub("gene_class ","",FCAT_gtfg$geneClass)
FCAT_gtfg$geneClass=gsub(";","",FCAT_gtfg$geneClass)
FCAT_gtf=FCAT_gtf[which(FCAT_gtf$V3 == "transcript"),]
FCAT_gtf$geneID=sapply(strsplit(FCAT_gtf$V9,"; "),"[",2)
FCAT_gtf$geneID=gsub("gene_id ","",FCAT_gtf$geneID)
FCAT_gtf$geneID=gsub(";","",FCAT_gtf$geneID)
FCAT_gtf$transcriptID=sapply(strsplit(FCAT_gtf$V9,"; "),"[",1)
FCAT_gtf$transcriptID=gsub("transcript_id ","",FCAT_gtf$transcriptID)
FCAT_anno=left_join(FCAT_gtf[,c(11,10)],FCAT_gtfg[,c(10,11)], by="geneID", copy=F)
write.table(FCAT_anno,paste0(CAT_path,"transcript_gene_link.tsv"), col.names=T, row.names=F, sep="\t", quote=F)



#compare FANTOM CAT with SALA final
#==================================================================
#bash
#prepare reference
setwd(paste0(SALA_path,"all_gtf_file/"))
system("zcat table5pENST.bed12.bed.gz | bgzip > table5pENST.bed12.bed.bgz")
system("tabix -p bed table5pENST.bed12.bed.bgz")

#==================================================================
#run input preparation for SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_FANTOMCAT/input.sh
#==================================================================

#==================================================================
#run transcript annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_FANTOMCAT/transcript.sh
#==================================================================

read_info=read.delim(paste0(CAT_path_log,"Neuron_THP1_FCAT.trnscpt.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_info$group1=substr(read_info$trnscpt_ID, start = 1, stop = 4)
read_info$group2=substr(read_info$model_ID_str, start = 1, stop = 4)
read_info%>%group_by(group1, group2)%>%dplyr::summarise(count=n())
length(unique(read_info$trnscpt_ID))-length(unique(read_info$model_ID_str))
#get all the ONTT that are the same as CATT

read_info1=read_info[which(read_info$group1 %in% c("ENCT","FTMT","MICT","HBMT")),]
read_info2=read_info1[which(read_info1$group2 == "ONTT"),]%>%group_by(model_ID_str)%>%dplyr::summarise(FCAT_transcriptID=paste(trnscpt_ID,collapse=";"))
write.table(read_info2,gzfile(paste0(CAT_path_log,"ONTT_FCATT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(read_info2,gzfile(paste0(output_path,"ONTT_FCATT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare input file for gene annotation 

table5.info=read.delim(paste0(SALA_path,"all_gtf_file/build_gtf/table5pENST.transcript.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.info1=table5.info[,c(9,10,11,14,15)]
write.table(table5.info1,paste0(CAT_path,"sala/gene/table5pENST.info.tsv"), col.names=F, row.names=F, sep="\t", quote=F)
#===============================================================================
#run gene annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_FANTOMCAT/gene.sh
#===============================================================================

ONTCATgene.info=read.delim(paste0(CAT_path,"sala/gene/Neuron_THP1_FCAT_disable_yes_10percent/log/Neuron_THP1_FCAT_disable_yes_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONTCATgene.info$group1=substr(ONTCATgene.info$model_ID, start = 1, stop = 4)
ONTCATgene.info$group2=substr(ONTCATgene.info$gene_ID, start = 1, stop = 4)
ONTCATgene.info%>%group_by(group2, group1)%>%dplyr::summarise(count=n())
length(unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$group2 == "NEWG")])) #118620
length(unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$group2 == "ONTG")])) #39390
length(unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$group2 == "ENSG")])) #61476
ONTCATgene.info1=ONTCATgene.info[which(ONTCATgene.info$group1 == "NEWT"),]
#find how many ONTG gene without any NEWT
ONTGout1=unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$group2 == "ONTG" & ONTCATgene.info$group1 == "NEWT")])
ONTGout2=unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$model_ID %in% sala_same)])
ONTGout=union(ONTGout1,ONTGout2) #15718
length(unique(ONTCATgene.info$gene_ID[which(ONTCATgene.info$group2 == "ONTG")]))-length(ONTGout)#23672
write.table(ONTGout,gzfile(paste0(CAT_path,"/sala/gene/ONTG.CAT_out.tsv.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#show which ONTG contain CATT
ONTCATgene.info3=ONTCATgene.info[which(ONTCATgene.info$gene_ID %in% ONTGout),]
read_info3=read_info[which(read_info$model_ID_str %in% ONTCATgene.info3$model_ID),c(1,7,15,16)]
read_info3=read_info3[which(read_info3$group1 %in% c("ENCT","FTMT", "HBMT", "MICT")),]
read_info_collap=read_info3%>%group_by(model_ID_str)%>%dplyr::summarise(FCAT_transcriptID=paste(trnscpt_ID, collapse=";"))
ONTCATgene.info3=left_join(ONTCATgene.info3,read_info_collap, by=c("model_ID"="model_ID_str"), copy=F)
gene_list=ONTCATgene.info3[which(!is.na(ONTCATgene.info3$FCAT_transcriptID)),]%>%group_by(gene_ID)%>%dplyr::summarise(FCAT_transcriptID=paste(FCAT_transcriptID, collapse=";"))
write.table(gene_list,gzfile(paste0(CAT_path_log,"ENSG_ONTG_FCATT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(gene_list,gzfile(paste0(output_path,"ENSG_ONTG_FCATT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


