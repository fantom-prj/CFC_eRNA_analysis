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

################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig3_data=paste0(primary_folder,"fig3/data/")
gencode_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_gencodev47/")
gencode_path_log=paste0(gencode_path,"sala/transcript/gencodev47.Neuron_THP1/log/")
output_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/")


#=============================
#prepare files from GENCODE v47
gtf=read.delim(paste0(gencode_path,"gencode.v47.annotation.gtf.gz"), header=F, stringsAsFactors = F, skip=5)
gtft=gtf[which(gtf$V3 == "transcript"),]
gtft$geneID=sapply(strsplit(gtft$V9,"; "),"[",1)
gtft$transcriptID=sapply(strsplit(gtft$V9,"; "),"[",2)
gtft$transcriptType=sapply(strsplit(gtft$V9,"; "),"[",5)
gtft$geneID=gsub("gene_id ","",gtft$geneID)
gtft$transcriptID=gsub("transcript_id ","",gtft$transcriptID)
gtft$transcriptType=gsub("transcript_type ","",gtft$transcriptType)
write.table(gtft[,c(10:12)],paste0(gencode_path,"gene_transcript_link.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

#compare GENCODE v47 with SALA final
#==================================================================
#bash
#prepare reference
setwd(paste0(SALA_path,"all_gtf_file/"))
system("zcat table5pENST.bed12.bed.gz | bgzip > table5pENST.bed12.bed.bgz")
system("tabix -p bed table5pENST.bed12.bed.bgz")

#===============================================================================
#run transcript annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_gencodev47/sala/transcript/script.sh
#===============================================================================

read_info=read.delim(paste0(gencode_path_log,"gencodev47.Neuron_THP1.trnscpt.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_info$group1=substr(read_info$trnscpt_ID, start = 1, stop = 4)
read_info$group2=substr(read_info$model_ID_str, start = 1, stop = 4)
read_info%>%group_by(group1, group2)%>%dplyr::summarise(count=n())
length(unique(read_info$trnscpt_ID))-length(unique(read_info$model_ID_str))

read_info1=read_info[which(read_info$group1 %in% c("ENST")),]
read_info2=read_info1[which(read_info1$group2 == "ONTT"),]%>%group_by(model_ID_str)%>%dplyr::summarise(GENCODEv47_transcriptID=paste(trnscpt_ID,collapse=";"))
write.table(read_info2,gzfile(paste0(gencode_path_log,"ONTT_gencodev47T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(read_info2,gzfile(paste0(output_path,"ONTT_gencodev47T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
sala_same=unique(read_info1$model_ID_str[which(read_info1$group2 %in% c("ONTT"))])

#===============================================================================
#prepare input file for gene annotation 

table5.info=read.delim(paste0(SALA_path,"all_gtf_file/build_gtf/table5pENST.transcript.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.info1=table5.info[,c(9,10,11,14,15)]
write.table(table5.info1,paste0(CAT_path,"sala/gene/table5pENST.info.tsv"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#run gene annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_gencodev47/sala/gene/script.10percent.sh
#===============================================================================

ONTgencodev47gene.info=read.delim(paste0(gencode_path,"sala/gene/Neuron_THP1_gencodev47_disable_yes_10percent/log/Neuron_THP1_gencodev47_disable_yes_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONTgencodev47gene.info$group1=substr(ONTgencodev47gene.info$model_ID, start = 1, stop = 4)
ONTgencodev47gene.info$group2=substr(ONTgencodev47gene.info$gene_ID, start = 1, stop = 4)
ONTgencodev47gene.info%>%group_by(group2, group1)%>%dplyr::summarise(count=n())
length(unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$group2 == "NEWG")])) #13128
length(unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$group2 == "ONTG")])) #39371 #why not 39425?
length(unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$group2 == "ENSG")])) #61239
ONTgencodev47gene.info1=ONTgencodev47gene.info[which(ONTgencodev47gene.info$group1 == "NEWT"),]
#find how many ONTG gene without any NEWT
length(unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$group2 == "ONTG" & ONTgencodev47gene.info$group1 == "NEWT")]))
ONTGout1=unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$group2 == "ONTG" & ONTgencodev47gene.info$group1 == "NEWT")])
ONTGout2=unique(ONTgencodev47gene.info$gene_ID[which(ONTgencodev47gene.info$model_ID %in% sala_same)])
ONTGout=union(ONTGout1,ONTGout2)
write.table(ONTGout,gzfile(paste0(gencode_path,"sala/gene/ONTG.gencodev47_out.tsv.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#show which ONTG contain gencodev47
ONTgencodev47gene.info3=ONTgencodev47gene.info[which(ONTgencodev47gene.info$gene_ID %in% ONTGout),]
read_info3=read_info[which(read_info$model_ID_str %in% unique(ONTgencodev47gene.info3$model_ID)),c(1,7,15,16)]
read_info3=read_info3[which(read_info3$group1 == "ENST" & read_info3$group2 != "ENST" ),]

read_info_collap=read_info3%>%group_by(model_ID_str)%>%dplyr::summarise(gencodev47_transcriptID=paste(unique(trnscpt_ID), collapse=";"))
ONTgencodev47gene.info3=left_join(ONTgencodev47gene.info3,read_info_collap, by=c("model_ID"="model_ID_str"), copy=F)
gene_list=ONTgencodev47gene.info3[which(!is.na(ONTgencodev47gene.info3$gencodev47_transcriptID)),]%>%group_by(gene_ID)%>%dplyr::summarise(gencodev47_transcriptID=paste(gencodev47_transcriptID, collapse=";"))
write.table(gene_list,gzfile(paste0(gencode_path_log,"ENSG_ONTG_gencodev47T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(gene_list,gzfile(paste0(output_path,"ENSG_ONTG_gencodev47T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

