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
refseq_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_refseq2024/")
refseq_path_log=paste0(refseq_path,"sala/transcript/refseq2024.Neuron_THP1/log/")
output_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/")

#============================
#prepare refseq 2024 transcriptome
#refseq2024 downloaded from https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/GCF_000001405.40-RS_2024_08/

#=============================
#gtf to bed12
setwd(refseq_path)
system("zcat GCF_000001405.40_GRCh38.p14_genomic.gtf.gz | bedparse gtf2bed | gzip > GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.gz")

#=============================

refseq=read.delim(paste0(refseq_path,"GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.gz"), header=F)
aa=refseq%>%group_by(V1)%>%dplyr::summarise(count=n())
chrom=read.delim(paste0(refseq_path,"sequence_report.tsv"), header=T, stringsAsFactors = F)
aa=left_join(aa, chrom[,c(10,4)], by=c("V1"="RefSeq.seq.accession"), copy=F)
aa=aa[c(1:25),]
aa$Chromosome.name=paste0("chr",aa$Chromosome.name)
refseq=left_join(refseq, aa[,c(1,3)], by="V1", copy=F)
refseq=refseq[which(!is.na(refseq$Chromosome.name)),]
write.table(refseq[order(refseq$Chromosome.name,refseq$V2),c(13,2:12)], gzfile(paste0(refseq_path,"GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=============================
#gtf to bed12
setwd(refseq_path)
system("zcat GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.gz | bgzip > GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.bgz")
system("tabix -p bed GCF_000001405.40_GRCh38.p14_genomic.bed12.bed.bgz")

#=============================

gtf=fread(paste0(refseq_path,"GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"), header=F)
gtf=left_join(gtf,aa[,c(1,3)], by="V1", copy=F)
gtf=gtf[which(!is.na(gtf$Chromosome.name)),]
write.table(gtf[,c(10,2:9)], gzfile(paste0(refseq_path,"GCF_000001405.40_GRCh38.p14_genomic.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

gtf=read.delim(paste0(refseq_path,"GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"), header=F, stringsAsFactors = F)
gtft=gtf[which(gtf$V3 == "transcript"),]
gtft$geneID=sapply(strsplit(gtft$V9,"; "),"[",1)
gtft$transcriptID=sapply(strsplit(gtft$V9,"; "),"[",2)
gtft$db_xref=sapply(strsplit(gtft$V9,"; "),"[",3)
gtft$gbkey=sapply(strsplit(gtft$V9,"gbkey "),"[",2)
gtft$gbkey=sapply(strsplit(gtft$gbkey,";"),"[",1)
gtft$geneID=gsub("gene_id ","",gtft$geneID)
gtft$transcriptID=gsub("transcript_id ","",gtft$transcriptID)
gtft$db_xref=sapply(strsplit(gtft$db_xref, ":"),"[",2)
gtft%>%group_by(gbkey)%>%dplyr::summarise(count=n())
write.table(gtft[,c(10:13)],paste0(refseq_path,"gene_transcript_link.tsv"), col.names=T, row.names=F, sep="\t", quote=F)


#==================================================================
#bash
#prepare reference
setwd(paste0(SALA_path,"all_gtf_file/"))
system("zcat table5pENST.bed12.bed.gz | bgzip > table5pENST.bed12.bed.bgz")
system("tabix -p bed table5pENST.bed12.bed.bgz")

#===============================================================================
#run transcript annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_refseq2024/sala/transcript/script.sh
#===============================================================================

read_info=read.delim(paste0(refseq_path_log,"refseq2024.Neuron_THP1.trnscpt.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_info$group1=substr(read_info$trnscpt_ID, start = 1, stop = 3)
read_info$group2=substr(read_info$model_ID_str, start = 1, stop = 3)
read_info%>%group_by(group1, group2)%>%dplyr::summarise(count=n())
length(unique(read_info$trnscpt_ID))-length(unique(read_info$model_ID_str))

read_info1=read_info[which(read_info$group1 != "ONT" & read_info$group1 !="ENS"),]
read_info2=read_info1[which(read_info1$group2 == "ONT"),]%>%group_by(model_ID_str)%>%dplyr::summarise(refseqT_transcriptID=paste(trnscpt_ID,collapse=";"))
write.table(read_info2,gzfile(paste0(refseq_path_log,"ONTT_refseq2024T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(read_info2,gzfile(paste0(output_path,"ONTT_refseq2024T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
sala_same=unique(read_info1$model_ID_str[which(read_info1$group2 %in% c("ONT"))])

#===============================================================================
#prepare input file for gene annotation 

table5.info=read.delim(paste0(SALA_path,"all_gtf_file/build_gtf/table5pENST.transcript.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.info1=table5.info[,c(9,10,11,14,15)]
write.table(table5.info1,paste0(CAT_path,"sala/gene/table5pENST.info.tsv"), col.names=F, row.names=F, sep="\t", quote=F)
#===============================================================================
#run gene annotation of SALA
#refer to [primary_folder]/code_n_data/transcript_model_analyses_Fig3/SALA_compare_exisiting_databases/Compare_refseq2024/sala/gene/script.10percent.sh
#===============================================================================

ONTrefseqgene.info=read.delim(paste0(refseq_path,"sala/gene/Neuron_THP1_refseq2024_disable_yes_10percent/log/Neuron_THP1_refseq2024_disable_yes_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ONTrefseqgene.info$group1=substr(ONTrefseqgene.info$model_ID, start = 1, stop = 4)
ONTrefseqgene.info$group2=substr(ONTrefseqgene.info$gene_ID, start = 1, stop = 4)
ONTrefseqgene.info%>%group_by(group2, group1)%>%dplyr::summarise(count=n())
length(unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$group2 == "NEWG")])) #7124
length(unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$group2 == "ONTG")])) #39401 #why not 39425?
length(unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$group2 == "ENSG")])) #61330
ONTrefseqgene.info1=ONTrefseqgene.info[which(ONTrefseqgene.info$group1 == "NEWT"),]
#find how many ONTG gene without any NEWT
length(unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$group2 == "ONTG" & ONTrefseqgene.info$group1 == "NEWT")]))
ONTGout1=unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$group2 == "ONTG" & ONTrefseqgene.info$group1 == "NEWT")])
ONTGout2=unique(ONTrefseqgene.info$gene_ID[which(ONTrefseqgene.info$model_ID %in% sala_same)])
ONTGout=union(ONTGout1,ONTGout2)
write.table(ONTGout,gzfile(paste0(refseq_path,"sala/gene/ONTG.refseq2024_out.tsv.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#show which ONTG contain refseq
ONTrefseqgene.info3=ONTrefseqgene.info[which(ONTrefseqgene.info$gene_ID %in% ONTGout),]
read_info3=read_info[which(read_info$model_ID_str %in% ONTrefseqgene.info3$model_ID),c(1,7,15,16)]
read_info3=read_info3[which(read_info3$group1 != "ONT" & read_info3$group1 != "ENS"),]
read_info_collap=read_info3%>%group_by(model_ID_str)%>%dplyr::summarise(refseq_transcriptID=paste(trnscpt_ID, collapse=";"))
ONTrefseqgene.info3=left_join(ONTrefseqgene.info3,read_info_collap, by=c("model_ID"="model_ID_str"), copy=F)
gene_list=ONTrefseqgene.info3[which(!is.na(ONTrefseqgene.info3$refseq_transcriptID)),]%>%group_by(gene_ID)%>%dplyr::summarise(refseq_transcriptID=paste(refseq_transcriptID, collapse=";"))
write.table(gene_list,gzfile(paste0(refseq_path_log,"ENSG_ONTG_refseq2024T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(gene_list,gzfile(paste0(output_path,"ENSG_ONTG_refseq2024T_match_transcript.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)




