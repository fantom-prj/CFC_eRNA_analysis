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


################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig3_data=paste0(primary_folder,"fig3/data/")
lncbook_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_lncbook/")
lncbook_path_log=paste0(lncbook_path,"sala/transcript/lncbook.Neuron_THP1/log/")
output_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/")

#=============================================================================
#gtf dnowloaded from https://ngdc.cncb.ac.cn/lncbook/downloads
#make transcript bed12
setwd(lncbook_path)
system("zcat lncRNA_LncBookv2.0_GRCh38.gtf.gz | bedparse gtf2bed | sort -k1,1 -k2,2n | gzip > lncRNA_LncBookv2.0_GRCh38.bed12.bed.gz")
system("zcat lncRNA_LncBookv2.0_GRCh38.bed12.bed.gz | bgzip > lncRNA_LncBookv2.0_GRCh38.bed12.bed.bgz")
system("tabix -p bed lncRNA_LncBookv2.0_GRCh38.bed12.bed.bgz")

#===================================
#look at the gtf # remove 5 transcripts from the bed12
lncbk=read.delim(paste0(lncbook_path,"lncRNA_LncBookv2.0_GRCh38.gtf.gz"), header=F, stringsAsFactors = F)
lncbk1=lncbk[union(grep("HSALNT0289467",lncbk$V9),grep("HSALNT0289449",lncbk$V9)),]
lncbk$transcriptID=sapply(strsplit(lncbk$V9, "; "),"[",2)
lncbk$transcriptID=gsub("transcript_id ","",lncbk$transcriptID)
lncbke=lncbk[which(lncbk$V3 == "exon"),]
lncbke1=lncbke%>%group_by(transcriptID)%>%dplyr::summarise(start=min(V4), end=max(V5))
lncbkt=lncbk[which(lncbk$V3 == "transcript"),]
lncbkt=left_join(lncbkt,lncbke1, by="transcriptID", copy=F)
lncbkt1=lncbkt[which(lncbkt$V4 != lncbkt$start | lncbkt$V5 != lncbkt$end),]
lncbkt$gene_type=sapply(strsplit(lncbkt$V9,"gene_type"),"[",2)
lncbkt$gene_type=sapply(strsplit(lncbkt$gene_type,";"),"[",1)
lncbkt$transcript_type=sapply(strsplit(lncbkt$V9,"transcript_type "),"[",2)
lncbkt$transcript_type=sapply(strsplit(lncbkt$transcript_type,";"),"[",1)
lncbkt$geneID=sapply(strsplit(lncbkt$V9,"; "),"[",1)
lncbkt$geneID=gsub("gene_id","",lncbkt$geneID)
lncbkt=lncbkt[-which(lncbkt$transcriptID %in% lncbkt1$transcriptID),]
write.table(lncbkt[,c(15,10:14)],paste0(lncbook_path,"gene_transcript.tsv"), col.names=T, row.names=F, sep="\t", quote=F)
bed12=read.delim(paste0(lncbook_path,"lncRNA_LncBookv2.0_GRCh38.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
bed12a=bed12[-which(bed12$V4 %in% lncbkt1$transcriptID),]
write.table(bed12a[order(bed12a$V1,bed12a$V2),],gzfile(paste0(lncbook_path,"lncRNA_LncBookv2.0_GRCh38.removed5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#compare LncBook with SALA final
#==================================================================
#bash
#prepare reference
setwd(paste0(SALA_path,"all_gtf_file/"))
system("zcat table5pENST.bed12.bed.gz | bgzip > table5pENST.bed12.bed.bgz")
system("tabix -p bed table5pENST.bed12.bed.bgz")

#===============================================================================
#run input preparation for SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_lncbook/input.sh
#===============================================================================

#===============================================================================
#run transcript annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_lncbook/transcript.sh
#===============================================================================

read_info=read.delim(paste0(lncbook_path_log,"lncbook.Neuron_THP1.trnscpt.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_info$group1=substr(read_info$trnscpt_ID, start = 1, stop = 4)
read_info$group2=substr(read_info$model_ID_str, start = 1, stop = 4)
read_info%>%group_by(group1, group2)%>%dplyr::summarise(count=n())
length(unique(read_info$trnscpt_ID))-length(unique(read_info$model_ID_str))

read_info1=read_info[which(read_info$group1 %in% c("HSAL")),]
read_info2=read_info1[which(read_info1$group2 == "ONTT"),]%>%group_by(model_ID_str)%>%dplyr::summarise(lncBookT_transcriptID=paste(trnscpt_ID,collapse=";"))
write.table(read_info2,gzfile(paste0(lncbook_path_log,"ONTT_lncbookT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(read_info2,gzfile(paste0(output_path,"ONTT_lncbookT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
sala_same=unique(read_info1$model_ID_str[which(read_info1$group2 %in% c("ONTT"))])


#===============================================================================
#prepare input file for gene annotation 

table5.info=read.delim(paste0(SALA_path,"all_gtf_file/build_gtf/table5pENST.transcript.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.info1=table5.info[,c(9,10,11,14,15)]
write.table(table5.info1,paste0(CAT_path,"sala/gene/table5pENST.info.tsv"), col.names=F, row.names=F, sep="\t", quote=F)
#===============================================================================
#run gene annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_lncbook/gene.sh
#===============================================================================

ONTlncbookgene.info=read.delim(paste0(lncbook_path,"sala/gene/Neuron_THP1_lncbook_disable_yes_10percent/log/Neuron_THP1_lncbook_disable_yes_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONTlncbookgene.info$group1=substr(ONTlncbookgene.info$model_ID, start = 1, stop = 4)
ONTlncbookgene.info$group2=substr(ONTlncbookgene.info$gene_ID, start = 1, stop = 4)
ONTlncbookgene.info%>%group_by(group2, group1)%>%dplyr::summarise(count=n())
length(unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$group2 == "NEWG")])) #69982
length(unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$group2 == "ONTG")])) #39399
length(unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$group2 == "ENSG")])) #61446
ONTlncbookgene.info1=ONTlncbookgene.info[which(ONTlncbookgene.info$group1 == "NEWT"),]
#find how many ONTG gene without any NEWT
length(unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$group2 == "ONTG" & ONTlncbookgene.info$group1 == "NEWT")]))
ONTGout1=unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$group2 == "ONTG" & ONTlncbookgene.info$group1 == "NEWT")])
ONTGout2=unique(ONTlncbookgene.info$gene_ID[which(ONTlncbookgene.info$model_ID %in% sala_same)])
ONTGout=union(ONTGout1,ONTGout2)
write.table(ONTGout,gzfile(paste0(lncbook_path,"sala/gene/ONTG.lncbook_out.tsv.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#show which ONTG contain lncT
ONTlncbookgene.info3=ONTlncbookgene.info[which(ONTlncbookgene.info$gene_ID %in% ONTGout),]
read_info3=read_info[which(read_info$model_ID_str %in% ONTlncbookgene.info3$model_ID),c(1,7,15,16)]
read_info3=read_info3[which(read_info3$group1 %in% c("HSAL")),]
read_info_collap=read_info3%>%group_by(model_ID_str)%>%dplyr::summarise(lncBook_transcriptID=paste(trnscpt_ID, collapse=";"))
ONTlncbookgene.info3=left_join(ONTlncbookgene.info3,read_info_collap, by=c("model_ID"="model_ID_str"), copy=F)
gene_list=ONTlncbookgene.info3[which(!is.na(ONTlncbookgene.info3$lncBook_transcriptID)),]%>%group_by(gene_ID)%>%dplyr::summarise(lncBook_transcriptID=paste(lncBook_transcriptID, collapse=";"))
write.table(gene_list,gzfile(paste0(lncbook_path_log,"ENSG_ONTG_lncbookT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(gene_list,gzfile(paste0(output_path,"ENSG_ONTG_lncbookT_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

