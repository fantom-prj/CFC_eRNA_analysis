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
library(tidyverse)

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
SALA_path=paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/")

#===============================================================================
path1=paste0(SALA_path,"transcript/zenbu/")
path2=paste0(SALA_path,"transcript/log/")
path3=paste0(SALA_path,"transcript/gtf/build_gtf/")
path4=paste0(SALA_path,"transcript/gtf/")
gene0_path=paste0(SALA_path,"table0_gene/iPSchro.table4ref.disable_ref_chain_bound_gene_anno_10percent/")
gene4_path=paste0(SALA_path,"table4_gene/T4_10percent/")


#===============================================================================
# The following folders were removed after incorporation:
# table0_gene, table4_gene, CPAT, Input 
# only provide upon request
#===============================================================================


#============================================================
###make gtf from bed12, bed6 exon, gene bed, gene info and GENCODE gtf
#for table1 (not table 0 )
library(tidyverse)
setwd(path1)
gencode_ont.gtf=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.gtf.gz"), header=F, stringsAsFactors = F)
table0.bed12=read.delim("iPSchro.table0.bed12.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.bed6=read.delim("iPSchro.table0.bed6.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.genebed=read.delim(paste0(gene0_path,"bed/iPSchro.table4ref.disable_ref_chain_bound_gene_anno_10percent.gene.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
transcript_info1a=read.delim(paste0(path2,"table1.remove_undetect_gencode_and_internal_prime.573k.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

gencode_ont.gtf.g=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "gene"),]
gencode_ont.gtf.t=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "transcript"),]
gencode_ont.gtf.e=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "exon"),]

gencode_ont.gtf.g$gene_ID=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",2)
gencode_ont.gtf.g$gene_ID=gsub(";","",gencode_ont.gtf.g$gene_ID)
gencode_ont.gtf.g$gene_type=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",4)
gencode_ont.gtf.g$gene_type=gsub(";","",gencode_ont.gtf.g$gene_type)
gencode_ont.gtf.g$gene_name=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",6)
gencode_ont.gtf.g$gene_name=gsub(";","",gencode_ont.gtf.g$gene_name)
gencode_ont.gtf.g=gencode_ont.gtf.g %>% mutate(across(c(gene_ID,gene_type,gene_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ga=gencode_ont.gtf.g[-which(gencode_ont.gtf.g$gene_ID %in% unique(transcript_info1a$IN1_gene_ID)),] 
gencode_ont.gtf.ga=gencode_ont.gtf.ga[,c(1:8,10:12)] #component1
gencode_ont.gtf.ga$gene_novelty="GENCODE_ONT"

table0.genebed$gene_ID=sapply(strsplit(table0.genebed$V4,"\\|"),"[",1)
table1.genebed=table0.genebed[which(table0.genebed$gene_ID %in% unique(transcript_info1a$IN1_gene_ID)),]
table1.genebed$label="CHR"
table1.genebed$type="gene"
table1.genebed=table1.genebed[,c(1,14,15,2,3,5,6,5,13)]
table1.genebed$V2=table1.genebed$V2+1
table1.genebed$V5="."
table1.genebed$V5.1="."
table1a.genebed=table1.genebed
table1a.genebed=left_join(table1a.genebed,unique(gencode_ont.gtf.g[,c(10,11)]), by="gene_ID", copy=F)
table1a.genebed=left_join(table1a.genebed,unique(transcript_info1a[,c(22,52,23)]), by=c("gene_ID"="IN1_gene_ID"), copy=F)
table1a.genebed$gene_type[which(is.na(table1a.genebed$gene_type))]=table1a.genebed$Novel_geneClass[which(is.na(table1a.genebed$gene_type))]
table1a.genebed=table1a.genebed[,-11]
table1a.genebed$gene_novelty="GENCODE_ONT_detected"
table1a.genebed$gene_novelty[grep("IN1G",table1a.genebed$IN1_gene_name)]="novel"

colnames(table1a.genebed)=colnames(gencode_ont.gtf.ga)
table1.genebed=rbind(table1a.genebed,gencode_ont.gtf.ga)
table1.genebed$V9=paste0("gene_id \"",table1.genebed$gene_ID,"\"; gene_type \"",table1.genebed$gene_type,"\"; gene_name \"",table1.genebed$gene_name,"\"; gene_novelty \"",table1.genebed$gene_novelty,"\";")
#write.table(table1.genebed,gzfile(paste0(path3,"table1pENST.gene.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#

gencode_ont.gtf.t$transcript_ID=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",4)
gencode_ont.gtf.t$transcript_ID=gsub(";","",gencode_ont.gtf.t$transcript_ID)
gencode_ont.gtf.t$gene_ID=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",2)
gencode_ont.gtf.t$gene_ID=gsub(";","",gencode_ont.gtf.t$gene_ID)
gencode_ont.gtf.t$transcript_type=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",10)
gencode_ont.gtf.t$transcript_type=gsub(";","",gencode_ont.gtf.t$transcript_type)
gencode_ont.gtf.t$transcript_name=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",12)
gencode_ont.gtf.t$transcript_name=gsub(";","",gencode_ont.gtf.t$transcript_name)
gencode_ont.gtf.t=gencode_ont.gtf.t %>% mutate(across(c(gene_ID,transcript_ID,transcript_type,transcript_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ta=gencode_ont.gtf.t[-which(gencode_ont.gtf.t$transcript_ID %in% transcript_info1a$model_ID),] 
gencode_ont.gtf.ta=gencode_ont.gtf.ta[,c(1:8,10:13)] #component1
gencode_ont.gtf.ta$transcript_novelty="GENCODE_ONT"

table1.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1a$model_ID),]
table1.bed12$label="CHR"
table1.bed12$type="transcript"
table1.bed12=table1.bed12[,c(1,13,14,2,3,5,6,5,4)]
table1.bed12$V2=table1.bed12$V2+1
table1.bed12$V5="."
table1.bed12$V5.1="."
table1a.bed12=table1.bed12
table1a.bed12=left_join(table1a.bed12, gencode_ont.gtf.t[,c(10:13)], by=c("V4"="transcript_ID"),copy=F)
table1a.bed12=left_join(table1a.bed12, transcript_info1a[,c(1,22,50)], by=c("V4"="model_ID"),copy=F)
table1a.bed12$gene_ID[which(is.na(table1a.bed12$gene_ID))]=table1a.bed12$IN1_gene_ID[which(is.na(table1a.bed12$gene_ID))]
table1a.bed12$transcript_type[which(is.na(table1a.bed12$transcript_type))]=table1a.bed12$Novel_transcriptClass[which(is.na(table1a.bed12$transcript_type))]
table1a.bed12$transcript_name[which(is.na(table1a.bed12$transcript_name))]=sapply(strsplit(table1a.bed12$V4[which(is.na(table1a.bed12$transcript_name))], "\\."),"[",1)
table1a.bed12$transcript_novelty="GENCODE_ONT_detected"
table1a.bed12$transcript_novelty[grep("CHRT",table1a.bed12$V4)]="novel"
table1a.bed12=table1a.bed12[,-c(13,14)]

colnames(table1a.bed12)=colnames(gencode_ont.gtf.ta)

table1.bed12=rbind(table1a.bed12,gencode_ont.gtf.ta)
table1.bed12=left_join(table1.bed12,table1.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table1.bed12$V9=paste0("gene_id \"",table1.bed12$gene_ID,"\"; transcript_id \"",table1.bed12$transcript_ID,"\"; gene_type \"",table1.bed12$gene_type,"\"; gene_name \"",table1.bed12$gene_name,"\"; transcript_type \"",table1.bed12$transcript_type,"\"; transcript_name \"",table1.bed12$transcript_name,"\"; gene_novelty \"",table1.bed12$gene_novelty,"\"; transcript_novelty \"",table1.bed12$transcript_novelty,"\";")
#write.table(table1.bed12,gzfile(paste0(path3,"table1pENST.transcript.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#

gencode_ont.gtf.e$transcript_ID=sapply(strsplit(gencode_ont.gtf.e$V9," "),"[",4)
gencode_ont.gtf.e$transcript_ID=gsub(";","",gencode_ont.gtf.e$transcript_ID)
gencode_ont.gtf.e$exon_number=sapply(strsplit(gencode_ont.gtf.e$V9," "),"[",14)
gencode_ont.gtf.e$exon_number=gsub(";","",gencode_ont.gtf.e$exon_number)
gencode_ont.gtf.e$exon_number=as.numeric(gsub("\"","",gencode_ont.gtf.e$exon_number))
gencode_ont.gtf.e=gencode_ont.gtf.e %>% mutate(across(c(transcript_ID),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ea=gencode_ont.gtf.e[-which(gencode_ont.gtf.e$transcript_ID %in% transcript_info1a$model_ID),]  
gencode_ont.gtf.ea=gencode_ont.gtf.ea[,c(1:8,10:11)] #component1

table1.bed6=table0.bed6[which(table0.bed6$V4 %in% transcript_info1a$model_ID),]
table1.bed6$label="CHR"
table1.bed6$type="exon"
table1.bed6=table1.bed6[,c(1,7,8,2,3,5,6,5,4)]
table1.bed6$V2=table1.bed6$V2+1
table1.bed6$V5="."
table1.bed6$V5.1="."
table1.bed6p=table1.bed6[which(table1.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table1.bed6n=table1.bed6[which(table1.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table1.bed6=rbind(table1.bed6p,table1.bed6n)
table1a.bed6=table1.bed6
colnames(table1a.bed6)=colnames(gencode_ont.gtf.ea)

table1.bed6=rbind(table1a.bed6,gencode_ont.gtf.ea)

table1.bed6=left_join(table1.bed6, table1.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
table1.bed6$V9=paste0("gene_id \"",table1.bed6$gene_ID,"\"; transcript_id \"",table1.bed6$transcript_ID,"\"; gene_type \"",table1.bed6$gene_type,"\"; gene_name \"",table1.bed6$gene_name,"\"; transcript_type \"",table1.bed6$transcript_type,"\"; transcript_name \"",table1.bed6$transcript_name,"\"; exon_number \"",table1.bed6$exon_number,"\";")
write.table(table1.bed6,gzfile(paste0(path3,"table1pENST.exon.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

options(scipen=999)
final=rbind(table1.genebed[,c(1:8,13)],table1.bed12[,c(1:8,17)], table1.bed6[,c(1:8,16)])
final=final[order(final$V1,final$V4),]
write.table(final,gzfile(paste0(path4,"table1pENST_ONT.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#=========
#detectable alone
final1=final[which(final$V2=="CHR"),]
#======
#update gene range
gtf=final1
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)
gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]

gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.g=gtf.t%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.gene=left_join(gtf.gene,gtf.g,by="geneID", copy=F)
gtf.gene=gtf.gene[,c(1:3,11,12,6:9)]
colnames(gtf.gene)=colnames(gtf.nogene)[c(1:9)]

gtf=rbind(gtf.gene,gtf.nogene[,c(1:9)])
write.table(gtf[order(gtf$V1,gtf$V4),],gzfile(paste0(path4,"table1_partial_yes_detected.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#bash
#validate the gtf
#setwd(path4)
#system("zcat table1pENST_ONT.gtf.gz | bedparse gtf2bed | gzip > table1pENST_ONT.bed12.bed.gz")
#system("zcat table1_partial_yes_detected.gtf.gz | bedparse gtf2bed | gzip > table1_partial_yes_detected.bed12.bed.gz")

#================================================================================
#for chr table5
library(tidyverse)
gencode_ont.gtf=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.gtf.gz"), header=F, stringsAsFactors = F)
table0.bed12=read.delim("iPSchro.table0.bed12.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.bed6=read.delim("iPSchro.table0.bed6.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table4.genebed=read.delim(paste0(gene4_path,"bed/T4_10percent.gene.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
transcript_info1a=read.delim(paste0(path2,"table5.chimeric.76K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

gencode_ont.gtf.g=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "gene"),]
gencode_ont.gtf.t=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "transcript"),]
gencode_ont.gtf.e=gencode_ont.gtf[which(gencode_ont.gtf$V3 == "exon"),]

gencode_ont.gtf.g$gene_ID=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",2)
gencode_ont.gtf.g$gene_ID=gsub(";","",gencode_ont.gtf.g$gene_ID)
gencode_ont.gtf.g$gene_type=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",4)
gencode_ont.gtf.g$gene_type=gsub(";","",gencode_ont.gtf.g$gene_type)
gencode_ont.gtf.g$gene_name=sapply(strsplit(gencode_ont.gtf.g$V9," "),"[",6)
gencode_ont.gtf.g$gene_name=gsub(";","",gencode_ont.gtf.g$gene_name)
gencode_ont.gtf.g=gencode_ont.gtf.g %>% mutate(across(c(gene_ID,gene_type,gene_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ga=gencode_ont.gtf.g[-which(gencode_ont.gtf.g$gene_ID %in% unique(transcript_info1a$T4_gene_ID)),] 
gencode_ont.gtf.ga=gencode_ont.gtf.ga[,c(1:8,10:12)] #component1
gencode_ont.gtf.ga$gene_novelty="GENCODE_ONT"

table4.genebed$gene_ID=sapply(strsplit(table4.genebed$V4,"\\|"),"[",1)
table5.genebed=table4.genebed[which(table4.genebed$gene_ID %in% unique(transcript_info1a$T4_gene_ID)),]

table5.genebed$label="CHR"
table5.genebed$type="gene"
table5.genebed=table5.genebed[,c(1,14,15,2,3,5,6,5,13)]
table5.genebed$V2=table5.genebed$V2+1
table5.genebed$V5="."
table5.genebed$V5.1="."

table5a.genebed=table5.genebed
table5a.genebed=left_join(table5a.genebed,unique(gencode_ont.gtf.g[,c(10,11)]), by="gene_ID", copy=F)
table5a.genebed=left_join(table5a.genebed,unique(transcript_info1a[,c(53,59,54)]), by=c("gene_ID"="T4_gene_ID"), copy=F)
table5a.genebed$gene_type[which(is.na(table5a.genebed$gene_type))]=table5a.genebed$T4_Novel_geneClass[which(is.na(table5a.genebed$gene_type))]
table5a.genebed=table5a.genebed[,-11]

table5a.genebed$gene_novelty="GENCODE_ONT_detected"
table5a.genebed$gene_novelty[grep("CHRG",table5a.genebed$T4_gene_name)]="novel"

colnames(table5a.genebed)=colnames(gencode_ont.gtf.ga)
table5.genebed=rbind(table5a.genebed,gencode_ont.gtf.ga)
table5.genebed$V9=paste0("gene_id \"",table5.genebed$gene_ID,"\"; gene_type \"",table5.genebed$gene_type,"\"; gene_name \"",table5.genebed$gene_name,"\"; gene_novelty \"",table5.genebed$gene_novelty,"\";")
#write.table(table5.genebed,gzfile(paste0(path3,"table5pENST.gene.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#

gencode_ont.gtf.t$transcript_ID=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",4)
gencode_ont.gtf.t$transcript_ID=gsub(";","",gencode_ont.gtf.t$transcript_ID)
gencode_ont.gtf.t$gene_ID=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",2)
gencode_ont.gtf.t$gene_ID=gsub(";","",gencode_ont.gtf.t$gene_ID)
gencode_ont.gtf.t$transcript_type=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",10)
gencode_ont.gtf.t$transcript_type=gsub(";","",gencode_ont.gtf.t$transcript_type)
gencode_ont.gtf.t$transcript_name=sapply(strsplit(gencode_ont.gtf.t$V9," "),"[",12)
gencode_ont.gtf.t$transcript_name=gsub(";","",gencode_ont.gtf.t$transcript_name)
gencode_ont.gtf.t=gencode_ont.gtf.t %>% mutate(across(c(gene_ID,transcript_ID,transcript_type,transcript_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ta=gencode_ont.gtf.t[-which(gencode_ont.gtf.t$transcript_ID %in% transcript_info1a$model_ID),] 
gencode_ont.gtf.ta=gencode_ont.gtf.ta[,c(1:8,10:13)] #component1
gencode_ont.gtf.ta$transcript_novelty="GENCODE_ONT"

table5.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1a$model_ID),]
table5.bed12$label="CHR"
table5.bed12$type="transcript"
table5.bed12=table5.bed12[,c(1,13,14,2,3,5,6,5,4)]
table5.bed12$V2=table5.bed12$V2+1
table5.bed12$V5="."
table5.bed12$V5.1="."

table5a.bed12=table5.bed12
table5a.bed12=left_join(table5a.bed12, gencode_ont.gtf.t[,c(10:13)], by=c("V4"="transcript_ID"),copy=F)
table5a.bed12=left_join(table5a.bed12, transcript_info1a[,c(1,53,50)], by=c("V4"="model_ID"),copy=F)
table5a.bed12$gene_ID[which(is.na(table5a.bed12$gene_ID))]=table5a.bed12$T4_gene_ID[which(is.na(table5a.bed12$gene_ID))]
table5a.bed12$transcript_type[which(is.na(table5a.bed12$transcript_type))]=table5a.bed12$Novel_transcriptClass[which(is.na(table5a.bed12$transcript_type))]
table5a.bed12$transcript_name[which(is.na(table5a.bed12$transcript_name))]=sapply(strsplit(table5a.bed12$V4[which(is.na(table5a.bed12$transcript_name))], "\\."),"[",1)
table5a.bed12$transcript_novelty="GENCODE_ONT_detected"
table5a.bed12$transcript_novelty[grep("CHRT",table5a.bed12$V4)]="novel"
table5a.bed12=table5a.bed12[,-c(13,14)]

colnames(table5a.bed12)=colnames(gencode_ont.gtf.ta)

table5.bed12=rbind(table5a.bed12,gencode_ont.gtf.ta)
table5.bed12=left_join(table5.bed12,table5.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table5.bed12$V9=paste0("gene_id \"",table5.bed12$gene_ID,"\"; transcript_id \"",table5.bed12$transcript_ID,"\"; gene_type \"",table5.bed12$gene_type,"\"; gene_name \"",table5.bed12$gene_name,"\"; transcript_type \"",table5.bed12$transcript_type,"\"; transcript_name \"",table5.bed12$transcript_name,"\"; gene_novelty \"",table5.bed12$gene_novelty,"\"; transcript_novelty \"",table5.bed12$transcript_novelty,"\";")
#write.table(table5.bed12,gzfile(paste0(path3,"table5pENST.transcript.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#
gencode_ont.gtf.e$transcript_ID=sapply(strsplit(gencode_ont.gtf.e$V9," "),"[",4)
gencode_ont.gtf.e$transcript_ID=gsub(";","",gencode_ont.gtf.e$transcript_ID)
gencode_ont.gtf.e$exon_number=sapply(strsplit(gencode_ont.gtf.e$V9," "),"[",14)
gencode_ont.gtf.e$exon_number=gsub(";","",gencode_ont.gtf.e$exon_number)
gencode_ont.gtf.e$exon_number=as.numeric(gsub("\"","",gencode_ont.gtf.e$exon_number))
gencode_ont.gtf.e=gencode_ont.gtf.e %>% mutate(across(c(transcript_ID),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode_ont.gtf.ea=gencode_ont.gtf.e[-which(gencode_ont.gtf.e$transcript_ID %in% transcript_info1a$model_ID),]  
gencode_ont.gtf.ea=gencode_ont.gtf.ea[,c(1:8,10:11)] #component1

table5.bed6=table0.bed6[which(table0.bed6$V4 %in% transcript_info1a$model_ID),]
table5.bed6$label="CHR"
table5.bed6$type="exon"
table5.bed6=table5.bed6[,c(1,7,8,2,3,5,6,5,4)]
table5.bed6$V2=table5.bed6$V2+1
table5.bed6$V5="."
table5.bed6$V5.1="."
table5.bed6p=table5.bed6[which(table5.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table5.bed6n=table5.bed6[which(table5.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table5.bed6=rbind(table5.bed6p,table5.bed6n)
table5a.bed6=table5.bed6
colnames(table5a.bed6)=colnames(gencode_ont.gtf.ea)

table5.bed6=rbind(table5a.bed6,gencode_ont.gtf.ea)

table5.bed6=left_join(table5.bed6, table5.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
table5.bed6$V9=paste0("gene_id \"",table5.bed6$gene_ID,"\"; transcript_id \"",table5.bed6$transcript_ID,"\"; gene_type \"",table5.bed6$gene_type,"\"; gene_name \"",table5.bed6$gene_name,"\"; transcript_type \"",table5.bed6$transcript_type,"\"; transcript_name \"",table5.bed6$transcript_name,"\"; exon_number \"",table5.bed6$exon_number,"\";")
#write.table(table5.bed6,gzfile(paste0(path3,"table5pENST.exon.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

options(scipen=999)
final=rbind(table5.genebed[,c(1:8,13)],table5.bed12[,c(1:8,17)], table5.bed6[,c(1:8,16)])
final=final[order(final$V1,final$V4),]
write.table(gtf[order(gtf$V1,gtf$V4),],gzfile(paste0(path4,"table5pENST_ONT.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#detectable alone
final1=final[which(final$V2=="CHR"),]

#======
#update gene range
gtf=final1
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)
gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]

gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.g=gtf.t%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.gene=left_join(gtf.gene,gtf.g,by="geneID", copy=F)
gtf.gene=gtf.gene[,c(1:3,11,12,6:9)]
colnames(gtf.gene)=colnames(gtf.nogene)[c(1:9)]

gtf=rbind(gtf.gene,gtf.nogene[,c(1:9)])
write.table(gtf[order(gtf$V1,gtf$V4),],gzfile(paste0(path4,"table5detected.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#bash
#validate the gtf
#setwd(path4)
#system("zcat table5pENST_ONT.gtf.gz | bedparse gtf2bed | gzip > table5pENST_ONT.bed12.bed.gz")
#system("zcat table5detected.gtf.gz | bedparse gtf2bed | gzip > table5detected.bed12.bed.gz")

#===============================================================================
#build final table 5 chromatin-bound with all full-length and partial ENST -> for long-read and short-read quantification

setwd(path1)

chr.table1=read.delim(paste0(path2,"table1.remove_undetect_and_partial_gencode_and_internal_prime.518k.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
chr.partial=read.delim(paste0(path2,"table1.only_partial_gencode_ONT.55K.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
trans1=chr.table1$model_ID[which(chr.table1$source %in% c("fulllength_gencode", "fulllength_ONT"))]
trans2=chr.partial$model_ID
gene1=unique(chr.table1$IN1_gene_ID[which(chr.table1$source %in% c("fulllength_gencode", "fulllength_ONT"))])
gene2=unique(chr.partial$IN1_gene_ID)

gtf.ori=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.gtf.gz"), header=F, stringsAsFactors = F)
gtf.ori.t=gtf.ori[which(gtf.ori$V3=="transcript"),]
gtf.ori.t$geneID=sapply(strsplit(gtf.ori.t$V9,";"),"[",1)
gtf.ori.t$geneID=gsub("gene_id \"","",gtf.ori.t$geneID)
gtf.ori.t$geneID=gsub("\"","",gtf.ori.t$geneID)
gtf.ori.t$transcriptID=sapply(strsplit(gtf.ori.t$V9,";"),"[",2)
gtf.ori.t$transcriptID=gsub(" transcript_id \"","",gtf.ori.t$transcriptID)
gtf.ori.t$transcriptID=gsub("\"","",gtf.ori.t$transcriptID)

t5transcript1=gtf.ori.t$transcriptID
t5gene1=unique(gtf.ori.t$geneID)

chr.table5=read.delim(paste0(path2,"table5.chimeric.76K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5transcript=unique(chr.table5$model_ID)
t5gene=unique(chr.table5$T4_gene_ID)

t5transcript=unique(c(t5transcript,t5transcript1,trans1,trans2))
t5gene=unique(c(t5gene,t5gene1,gene1,gene2))

gtf=fread(paste0(path4,"table5pENST_ONT.gtf.gz"), header=F, stringsAsFactors = F)
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)

gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]
gtf.gene5=gtf.gene[which(gtf.gene$geneID %in% t5gene),]
gtf.nogene5=gtf.nogene[which(gtf.nogene$transcriptID %in% t5transcript),]

gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.t5=gtf.t[which(gtf.t$transcriptID %in% t5transcript),]
gtf.g5=gtf.t5%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene5=left_join(gtf.gene5,gtf.g5,by="geneID", copy=F)
gtf.gene5=gtf.gene5[,c(1:3,11,12,6:9)]
colnames(gtf.gene5)=colnames(gtf.nogene5)[c(1:9)]

gtf5=rbind(gtf.gene5,gtf.nogene5[,c(1:9)])

write.table(gtf5[order(gtf5$V1,gtf5$V4),],gzfile(paste0(path4,"table5_partial_yes_detected.alone_allNeuron_THP1t5.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#bash
#validate the gtf
setwd(path4)
system("zcat table5_partial_yes_detected.alone_allNeuron_THP1t5.gtf.gz | bedparse gtf2bed | gzip > table5_partial_yes_detected.alone_allNeuron_THP1t5.bed12.bed.gz")
