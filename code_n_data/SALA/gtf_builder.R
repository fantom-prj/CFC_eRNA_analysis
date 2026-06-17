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
library(purrr)

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
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/")
GENCODE_path=paste0(primary_folder,"code_n_data/GENCODEv39/")

#===============================================================================
path1=paste0(SALA_path,"transcript/zenbu/")
path2=paste0(SALA_path,"transcript/log/")
path3=paste0(SALA_path,"transcript/all_gtf_file/build_gtf/")
path4=paste0(SALA_path,"transcript/all_gtf_file/")
gene0_path=paste0(SALA_path,"table0_gene/iPSC_NSC_Neuron.S3.disable_ref_chain_bound_gene_anno_10percent/")
gene4_path=paste0(SALA_path,"table4_gene/Neuron_THP1_T4_10percent/")

setwd(path1)

#============================================================
###make gtf from bed12, bed6 exon, gene bed, gene info and GENCODE gtf
#for table1 (raw)

gencode.gtf=fread(paste0(GENCODE_path,"gencode.v39.annotation.gtf.gz"), header=F, stringsAsFactors = F)

table0.bed12=read.delim("Neuron_THP1.table0.bed12.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.bed6=read.delim("Neuron_THP1.table0.bed6.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.genebed=read.delim(paste0(gene0_path,"bed/iPSC_NSC_Neuron.S3.disable_ref_chain_bound_gene_anno_10percent.gene.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
transcript_info1a=read.delim(paste0(path2,"table1.remove_undetect_and_partial_gencode_and_internal_prime.2.47M.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

transcript_info1aa=transcript_info1a[which(transcript_info1a$ENST_adjust != "No" & !is.na(transcript_info1a$ENST_adjust)),]
transcript_info1bb=transcript_info1a[which(transcript_info1a$ENSG_adjust != "No" & !is.na(transcript_info1a$ENSG_adjust)),]

gencode.gtf.g=gencode.gtf[which(gencode.gtf$V3 == "gene"),]
gencode.gtf.t=gencode.gtf[which(gencode.gtf$V3 == "transcript"),]
gencode.gtf.e=gencode.gtf[which(gencode.gtf$V3 == "exon"),]

gencode.gtf.g$gene_ID=sapply(strsplit(gencode.gtf.g$V9," "),"[",2)
gencode.gtf.g$gene_ID=gsub(";","",gencode.gtf.g$gene_ID)
gencode.gtf.g$gene_type=sapply(strsplit(gencode.gtf.g$V9," "),"[",4)
gencode.gtf.g$gene_type=gsub(";","",gencode.gtf.g$gene_type)
gencode.gtf.g$gene_name=sapply(strsplit(gencode.gtf.g$V9," "),"[",6)
gencode.gtf.g$gene_name=gsub(";","",gencode.gtf.g$gene_name)
gencode.gtf.g=gencode.gtf.g %>% mutate(across(c(gene_ID,gene_type,gene_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ga=gencode.gtf.g[-which(gencode.gtf.g$gene_ID %in% unique(transcript_info1bb$IN1_gene_ID)),] 
gencode.gtf.ga=gencode.gtf.ga[,c(1:8,10:12)] #component1
gencode.gtf.ga$V2="GENCODE"
gencode.gtf.ga$gene_novelty="GENCODEv39"

table0.genebed$gene_ID=sapply(strsplit(table0.genebed$V4,"\\|"),"[",1)
table1.genebed=table0.genebed[which(table0.genebed$gene_ID %in% unique(transcript_info1a$IN1_gene_ID)),]
table1.genebed$label="ONTCAGE"
table1.genebed$type="gene"
table1.genebed=table1.genebed[,c(1,14,15,2,3,5,6,5,13)]
table1.genebed$V2=table1.genebed$V2+1
table1.genebed$V5="."
table1.genebed$V5.1="."
table1a.genebed=table1.genebed[-grep("ENSG",table1.genebed$gene_ID),]
table1a.genebed=left_join(table1a.genebed,unique(transcript_info1a[,c(42,82)]), by=c("gene_ID"="IN1_gene_ID"), copy=F)
table1a.genebed$gene_name=sapply(strsplit(table1a.genebed$gene_ID,"\\."),"[",1)
table1a.genebed$gene_novelty="novel"
table1b.genebed=table1.genebed[which(table1.genebed$gene_ID %in% unique(transcript_info1bb$IN1_gene_ID)),]
table1b.genebed=left_join(table1b.genebed,gencode.gtf.g[,c(10:12)],by="gene_ID",copy=F)
table1b.genebed$novelty="GENCODE_updated"
table1b.genebed$label="GENCODE"
colnames(table1a.genebed)=colnames(gencode.gtf.ga)
colnames(table1b.genebed)=colnames(gencode.gtf.ga)
table1.genebed=rbind(table1a.genebed,table1b.genebed,gencode.gtf.ga)
table1.genebed$V9=paste0("gene_id \"",table1.genebed$gene_ID,"\"; gene_type \"",table1.genebed$gene_type,"\"; gene_name \"",table1.genebed$gene_name,"\"; gene_novelty \"",table1.genebed$gene_novelty,"\";")
write.table(table1.genebed,gzfile(paste0(path3,"table1pENST.gene.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#

gencode.gtf.t$transcript_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",4)
gencode.gtf.t$transcript_ID=gsub(";","",gencode.gtf.t$transcript_ID)
gencode.gtf.t$gene_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",2)
gencode.gtf.t$gene_ID=gsub(";","",gencode.gtf.t$gene_ID)
gencode.gtf.t$transcript_type=sapply(strsplit(gencode.gtf.t$V9," "),"[",10)
gencode.gtf.t$transcript_type=gsub(";","",gencode.gtf.t$transcript_type)
gencode.gtf.t$transcript_name=sapply(strsplit(gencode.gtf.t$V9," "),"[",12)
gencode.gtf.t$transcript_name=gsub(";","",gencode.gtf.t$transcript_name)
gencode.gtf.t=gencode.gtf.t %>% mutate(across(c(gene_ID,transcript_ID,transcript_type,transcript_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ta=gencode.gtf.t[-which(gencode.gtf.t$transcript_ID %in% transcript_info1aa$model_ID),] 
gencode.gtf.ta=gencode.gtf.ta[,c(1:8,10:13)] #component1
gencode.gtf.ta$V2="GENCODE"
gencode.gtf.ta$transcript_novelty="GENCODEv39"

table1.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1a$model_ID),]
table1.bed12$label="ONTCAGE"
table1.bed12$type="transcript"
table1.bed12=table1.bed12[,c(1,13,14,2,3,5,6,5,4)]
table1.bed12$V2=table1.bed12$V2+1
table1.bed12$V5="."
table1.bed12$V5.1="."
table1a.bed12=table1.bed12[-grep("ENST",table1.bed12$V4),]
table1a.bed12=left_join(table1a.bed12, transcript_info1a[,c(1,42,78)], by=c("V4"="model_ID"),copy=F)
table1a.bed12$transcript_name=sapply(strsplit(table1a.bed12$V4,"\\."),"[",1)
table1a.bed12$transcript_novelty="novel"
table1b.bed12=table1.bed12[which(table1.bed12$V4 %in% transcript_info1aa$model_ID),]
table1b.bed12=left_join(table1b.bed12,gencode.gtf.t[,c(10:13)],by=c("V4"="transcript_ID"),copy=F) #take geneID, transcript_type, transcript_name
table1b.bed12$transcript_novelty="GENCODE_updated"
table1b.bed12$label="GENCODE"
colnames(table1a.bed12)=colnames(gencode.gtf.ta)
colnames(table1b.bed12)=colnames(gencode.gtf.ta)
table1.bed12=rbind(table1a.bed12,table1b.bed12,gencode.gtf.ta)
table1.bed12=left_join(table1.bed12,table1.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table1.bed12$V9=paste0("gene_id \"",table1.bed12$gene_ID,"\"; transcript_id \"",table1.bed12$transcript_ID,"\"; gene_type \"",table1.bed12$gene_type,"\"; gene_name \"",table1.bed12$gene_name,"\"; transcript_type \"",table1.bed12$transcript_type,"\"; transcript_name \"",table1.bed12$transcript_name,"\"; gene_novelty \"",table1.bed12$gene_novelty,"\"; transcript_novelty \"",table1.bed12$transcript_novelty,"\";")
write.table(table1.bed12,gzfile(paste0(path3,"table1pENST.transcript.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

gencode.gtf.e$transcript_ID=sapply(strsplit(gencode.gtf.e$V9," "),"[",4)
gencode.gtf.e$transcript_ID=gsub(";","",gencode.gtf.e$transcript_ID)
gencode.gtf.e$exon_number=sapply(strsplit(gencode.gtf.e$V9," "),"[",14)
gencode.gtf.e$exon_number=as.numeric(gsub(";","",gencode.gtf.e$exon_number))
gencode.gtf.e$exon_id=sapply(strsplit(gencode.gtf.e$V9," "),"[",16)
gencode.gtf.e$exon_id=gsub(";","",gencode.gtf.e$exon_id)
gencode.gtf.e=gencode.gtf.e %>% mutate(across(c(transcript_ID,exon_id),~ map_chr(.x, ~ gsub("\"", "", .x))))
exon.list=unique(gencode.gtf.e[,c(10:12)])
gencode.gtf.ea=gencode.gtf.e[-which(gencode.gtf.e$transcript_ID %in% transcript_info1aa$model_ID),]  
gencode.gtf.ea=gencode.gtf.ea[,c(1:8,10:11)] #component1

table1.bed6=table0.bed6[which(table0.bed6$V4 %in% transcript_info1a$model_ID),]
table1.bed6$label="ONTCAGE"
table1.bed6$type="exon"
table1.bed6=table1.bed6[,c(1,7,8,2,3,5,6,5,4)]
table1.bed6$V2=table1.bed6$V2+1
table1.bed6$V5="."
table1.bed6$V5.1="."
table1.bed6p=table1.bed6[which(table1.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table1.bed6n=table1.bed6[which(table1.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table1.bed6=rbind(table1.bed6p,table1.bed6n)

table1a.bed6=table1.bed6[-grep("ENST",table1.bed6$V4),]
table1b.bed6=table1.bed6[which(table1.bed6$V4 %in% transcript_info1aa$model_ID),]
colnames(table1a.bed6)=colnames(gencode.gtf.ea)
colnames(table1b.bed6)=colnames(gencode.gtf.ea)
table1.bed6=rbind(table1a.bed6,table1b.bed6,gencode.gtf.ea)

table1.bed6=left_join(table1.bed6, table1.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
#table1.bed6=left_join(table1.bed6, exon.list, by=c("transcript_ID"="transcript_ID","exon_number"="exon_number"), copy=F)
table1.bed6$V9=paste0("gene_id \"",table1.bed6$gene_ID,"\"; transcript_id \"",table1.bed6$transcript_ID,"\"; gene_type \"",table1.bed6$gene_type,"\"; gene_name \"",table1.bed6$gene_name,"\"; transcript_type \"",table1.bed6$transcript_type,"\"; transcript_name \"",table1.bed6$transcript_name,"\"; exon_number \"",table1.bed6$exon_number,"\";")
write.table(table1.bed6,gzfile(paste0(path3,"table1pENST.exon.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

options(scipen=999)
final=rbind(table1.genebed[,c(1:8,13)],table1.bed12[,c(1:8,17)], table1.bed6[,c(1:8,16)])
final=final[order(final$V1,final$V4),]
write.table(final,gzfile(paste0(path4,"table1pENST.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#bash
#validate the gtf
#setwd(path4)
#system("zcat table1pENST.gtf.gz | bedparse gtf2bed | gzip > table1pENST.bed12.bed.gz")

#===============================================================================
#table1 gtf -> table2 gtf (read-filtered)

gtf=final
transcript_info2=read.delim(paste0(path2,"table2.standard.153K.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

t2transcript=unique(transcript_info2$model_ID)
t2gene=unique(transcript_info2$IN1_gene_ID)

#gtf trimmer
#===split the gtf into table1.noIP and table2 and revise gene range for all (all contain all ENSG and ENST)===
gtf=fread("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/homemade_full/build_gtf/table1pENST.gtf", header=F, stringsAsFactors = F)
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)
gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]
gtf.gene2=gtf.gene[union(which(gtf.gene$geneID %in% t2gene),grep("ENSG",gtf.gene$geneID)),]
gtf.nogene2=gtf.nogene[union(which(gtf.nogene$transcriptID %in% t2transcript), grep("ENST",gtf.nogene$transcriptID)),]
#revise gene region start end, due to loss of transcript model
gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.t2=gtf.t[union(which(gtf.t$transcriptID %in% t2transcript),grep("ENST",gtf.t$transcriptID)),]
gtf.g=gtf.t%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.g2=gtf.t2%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene=left_join(gtf.gene,gtf.g,by="geneID", copy=F)
gtf.gene2=left_join(gtf.gene2,gtf.g2,by="geneID", copy=F)
gtf.gene=gtf.gene[,c(1:3,11,12,6:9)]
gtf.gene2=gtf.gene2[,c(1:3,11,12,6:9)]
colnames(gtf.gene)=colnames(gtf.nogene)[c(1:9)]
colnames(gtf.gene2)=colnames(gtf.nogene)[c(1:9)]

gtf=rbind(gtf.gene,gtf.nogene[,c(1:9)])
gtf2=rbind(gtf.gene2,gtf.nogene2[,c(1:9)])
write.table(gtf2[order(gtf2$V1,gtf2$V4),],gzfile(paste0(path4,"table2pENST.5read.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#============================================================
#bash
#validate the gtf
#setwd(path4)
#system("zcat table2pENST.5read.gtf.gz | bedparse gtf2bed | gzip > table2pENST.5read.bed12.bed.gz")


#===============================================================================
#build gtf for table5 (finalized transcriptome)
#make gtf from bed12, bed6 exon, gene bed, gene info and GENCODE gtf

setwd(path1)
gencode.gtf=fread(paste0(GENCODE_path,"gencode.v39.annotation.gtf.gz"), header=F, stringsAsFactors = F)
table0.bed12=read.delim("Neuron_THP1.table0.bed12.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.bed6=read.delim("Neuron_THP1.table0.bed6.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table4.genebed=read.delim(paste0(gene4_path,"bed/Neuron_THP1_T4_10percent.gene.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
transcript_info5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

transcript_info1aa=transcript_info5[which(transcript_info5$ENST_adjust != "No" & !is.na(transcript_info5$ENST_adjust)),]
transcript_info1bb=transcript_info5[which(transcript_info5$T4_ENSG_adjust != "No" & !is.na(transcript_info5$T4_ENSG_adjust)),]
gencode.gtf.g=gencode.gtf[which(gencode.gtf$V3 == "gene"),]
gencode.gtf.t=gencode.gtf[which(gencode.gtf$V3 == "transcript"),]
gencode.gtf.e=gencode.gtf[which(gencode.gtf$V3 == "exon"),]

gencode.gtf.g$gene_ID=sapply(strsplit(gencode.gtf.g$V9," "),"[",2)
gencode.gtf.g$gene_ID=gsub(";","",gencode.gtf.g$gene_ID)
gencode.gtf.g$gene_type=sapply(strsplit(gencode.gtf.g$V9," "),"[",4)
gencode.gtf.g$gene_type=gsub(";","",gencode.gtf.g$gene_type)
gencode.gtf.g$gene_name=sapply(strsplit(gencode.gtf.g$V9," "),"[",6)
gencode.gtf.g$gene_name=gsub(";","",gencode.gtf.g$gene_name)
gencode.gtf.g=gencode.gtf.g %>% mutate(across(c(gene_ID,gene_type,gene_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ga=gencode.gtf.g[-which(gencode.gtf.g$gene_ID %in% unique(transcript_info1bb$T4_gene_ID)),] 
gencode.gtf.ga=gencode.gtf.ga[,c(1:8,10:12)] #component1
gencode.gtf.ga$V2="GENCODE"
gencode.gtf.ga$gene_novelty="GENCODEv39"

table4.genebed$gene_ID=sapply(strsplit(table4.genebed$V4,"\\|"),"[",1)
table5.genebed=table4.genebed[which(table4.genebed$gene_ID %in% unique(transcript_info5$T4_gene_ID)),]
table5.genebed$label="ONTCAGE"
table5.genebed$type="gene"
table5.genebed=table5.genebed[,c(1,14,15,2,3,5,6,5,13)]
table5.genebed$V2=table5.genebed$V2+1
table5.genebed$V5="."
table5.genebed$V5.1="."
table5a.genebed=table5.genebed[-grep("ENSG",table5.genebed$gene_ID),]
table5a.genebed=left_join(table5a.genebed,unique(transcript_info5[,c(86,93)]), by=c("gene_ID"="T4_gene_ID"), copy=F)
table5a.genebed$gene_name=sapply(strsplit(table5a.genebed$gene_ID,"\\."),"[",1)
table5a.genebed$gene_novelty="novel"
table5b.genebed=table5.genebed[which(table5.genebed$gene_ID %in% unique(transcript_info1bb$T4_gene_ID)),]
table5b.genebed=left_join(table5b.genebed,gencode.gtf.g[,c(10:12)],by="gene_ID",copy=F)
table5b.genebed$novelty="GENCODE_updated"
table5b.genebed$label="GENCODE"
colnames(table5a.genebed)=colnames(gencode.gtf.ga)
colnames(table5b.genebed)=colnames(gencode.gtf.ga)
table5.genebed=rbind(table5a.genebed,table5b.genebed,gencode.gtf.ga)
table5.genebed$V9=paste0("gene_id \"",table5.genebed$gene_ID,"\"; gene_type \"",table5.genebed$gene_type,"\"; gene_name \"",table5.genebed$gene_name,"\"; gene_novelty \"",table5.genebed$gene_novelty,"\";")
#write.table(table5.genebed,paste0(path3,"table5pENST.gene.gtf"), col.names=F, row.names=F, sep="\t", quote=F)
#

gencode.gtf.t$transcript_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",4)
gencode.gtf.t$transcript_ID=gsub(";","",gencode.gtf.t$transcript_ID)
gencode.gtf.t$gene_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",2)
gencode.gtf.t$gene_ID=gsub(";","",gencode.gtf.t$gene_ID)
gencode.gtf.t$transcript_type=sapply(strsplit(gencode.gtf.t$V9," "),"[",10)
gencode.gtf.t$transcript_type=gsub(";","",gencode.gtf.t$transcript_type)
gencode.gtf.t$transcript_name=sapply(strsplit(gencode.gtf.t$V9," "),"[",12)
gencode.gtf.t$transcript_name=gsub(";","",gencode.gtf.t$transcript_name)
gencode.gtf.t=gencode.gtf.t %>% mutate(across(c(gene_ID,transcript_ID,transcript_type,transcript_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ta=gencode.gtf.t[-which(gencode.gtf.t$transcript_ID %in% transcript_info1aa$model_ID),] 
gencode.gtf.ta=gencode.gtf.ta[,c(1:8,10:13)] #component1
gencode.gtf.ta$V2="GENCODE"
gencode.gtf.ta$transcript_novelty="GENCODEv39"

table5.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info5$model_ID),]
table5.bed12$label="ONTCAGE"
table5.bed12$type="transcript"
table5.bed12=table5.bed12[,c(1,13,14,2,3,5,6,5,4)]
table5.bed12$V2=table5.bed12$V2+1
table5.bed12$V5="."
table5.bed12$V5.1="."
table5a.bed12=table5.bed12[-grep("ENST",table5.bed12$V4),]
table5a.bed12=left_join(table5a.bed12, transcript_info5[,c(1,86,78)], by=c("V4"="model_ID"),copy=F)
table5a.bed12$transcript_name=sapply(strsplit(table5a.bed12$V4,"\\."),"[",1)
table5a.bed12$transcript_novelty="novel"
table5b.bed12=table5.bed12[which(table5.bed12$V4 %in% transcript_info1aa$model_ID),]
table5b.bed12=left_join(table5b.bed12,gencode.gtf.t[,c(10:13)],by=c("V4"="transcript_ID"),copy=F) #take geneID, transcript_type, transcript_name
table5b.bed12$transcript_novelty="GENCODE_updated"
table5b.bed12$label="GENCODE"
colnames(table5a.bed12)=colnames(gencode.gtf.ta)
colnames(table5b.bed12)=colnames(gencode.gtf.ta)
table5.bed12=rbind(table5a.bed12,table5b.bed12,gencode.gtf.ta)
table5.bed12=left_join(table5.bed12,table5.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table5.bed12$V9=paste0("gene_id \"",table5.bed12$gene_ID,"\"; transcript_id \"",table5.bed12$transcript_ID,"\"; gene_type \"",table5.bed12$gene_type,"\"; gene_name \"",table5.bed12$gene_name,"\"; transcript_type \"",table5.bed12$transcript_type,"\"; transcript_name \"",table5.bed12$transcript_name,"\"; gene_novelty \"",table5.bed12$gene_novelty,"\"; transcript_novelty \"",table5.bed12$transcript_novelty,"\";")
#write.table(table5.bed12,paste0(path3,"table5pENST.transcript.gtf"), col.names=F, row.names=F, sep="\t", quote=F)

gencode.gtf.e$transcript_ID=sapply(strsplit(gencode.gtf.e$V9," "),"[",4)
gencode.gtf.e$transcript_ID=gsub(";","",gencode.gtf.e$transcript_ID)
gencode.gtf.e$exon_number=sapply(strsplit(gencode.gtf.e$V9," "),"[",14)
gencode.gtf.e$exon_number=as.numeric(gsub(";","",gencode.gtf.e$exon_number))
gencode.gtf.e$exon_id=sapply(strsplit(gencode.gtf.e$V9," "),"[",16)
gencode.gtf.e$exon_id=gsub(";","",gencode.gtf.e$exon_id)
gencode.gtf.e=gencode.gtf.e %>% mutate(across(c(transcript_ID,exon_id),~ map_chr(.x, ~ gsub("\"", "", .x))))
exon.list=unique(gencode.gtf.e[,c(10:12)])
gencode.gtf.ea=gencode.gtf.e[-which(gencode.gtf.e$transcript_ID %in% transcript_info1aa$model_ID),]  
gencode.gtf.ea=gencode.gtf.ea[,c(1:8,10:11)] #component1

table5.bed6=table0.bed6[which(table0.bed6$V4 %in% transcript_info5$model_ID),]
table5.bed6$label="ONTCAGE"
table5.bed6$type="exon"
table5.bed6=table5.bed6[,c(1,7,8,2,3,5,6,5,4)]
table5.bed6$V2=table5.bed6$V2+1
table5.bed6$V5="."
table5.bed6$V5.1="."
table5.bed6p=table5.bed6[which(table5.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table5.bed6n=table5.bed6[which(table5.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table5.bed6=rbind(table5.bed6p,table5.bed6n)

table5a.bed6=table5.bed6[-grep("ENST",table5.bed6$V4),]
table5b.bed6=table5.bed6[which(table5.bed6$V4 %in% transcript_info1aa$model_ID),]
colnames(table5a.bed6)=colnames(gencode.gtf.ea)
colnames(table5b.bed6)=colnames(gencode.gtf.ea)
table5.bed6=rbind(table5a.bed6,table5b.bed6,gencode.gtf.ea)

table5.bed6=left_join(table5.bed6, table5.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
table5.bed6$V9=paste0("gene_id \"",table5.bed6$gene_ID,"\"; transcript_id \"",table5.bed6$transcript_ID,"\"; gene_type \"",table5.bed6$gene_type,"\"; gene_name \"",table5.bed6$gene_name,"\"; transcript_type \"",table5.bed6$transcript_type,"\"; transcript_name \"",table5.bed6$transcript_name,"\"; exon_number \"",table5.bed6$exon_number,"\";")
#write.table(table5.bed6,gzfile(paste0(path3,"table5pENST.exon.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

options(scipen=999)
final=rbind(table5.genebed[,c(1:8,13)],table5.bed12[,c(1:8,17)], table5.bed6[,c(1:8,16)])
final=final[order(final$V1,final$V4),]

#revise gene region start end
gtf=final
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
write.table(gtf[order(gtf$V1,gtf$V4),],gzfile(paste0(path4,"table5pENST.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#=========
#bash
#validate the gtf
#setwd(path4)
#system("zcat table5pENST.gtf.gz | bedparse gtf2bed | gzip > table5pENST.bed12.bed.gz")
#system("gzip -dc table5pENST.bed12.bed.gz | bgzip -c > table5pENST.bed12.bed.bgz")
#system("tabix -p bed table5pENST.bed12.bed.bgz")


#===============================================================================
#gtf trimmer: table1 remove un-detectable ENST but keep partially detected ENST

transcript_info1c=read.delim(paste0(path2,"table1.remove_undetect_gencode_and_internal_prime.2.51M.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t1transcript=unique(transcript_info1c$model_ID)
t1gene=unique(transcript_info1c$IN1_gene_ID)
#===split the gtf into table1.noIP and table2 and revise gene range for all (all contain all ENSG and ENST)===
gtf=fread(paste0(path4,"table1pENST.gtf.gz"), header=F, stringsAsFactors = F)
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)
gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]
gtf.gene1=gtf.gene[which(gtf.gene$geneID %in% t1gene),]
gtf.nogene1=gtf.nogene[which(gtf.nogene$transcriptID %in% t1transcript),]
#revise gene region start end
gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.t1=gtf.t[which(gtf.t$transcriptID %in% t1transcript),]
gtf.g1=gtf.t1%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene1=left_join(gtf.gene1,gtf.g1,by="geneID", copy=F)
gtf.gene1=gtf.gene1[,c(1:3,11,12,6:9)]
colnames(gtf.gene1)=colnames(gtf.nogene)[c(1:9)]

gtf1=rbind(gtf.gene1,gtf.nogene1[,c(1:9)])
write.table(gtf1[order(gtf1$V1,gtf1$V4),],gzfile(paste0(path4,"table1.noIP.yesPartial.detected.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#============================================================
#gtf trimmer: table5 remove un-detectable ENST and partially detected ENST. --> for zenbu visualization
transcript_info5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5transcript=unique(transcript_info5$model_ID)
t5gene=unique(transcript_info5$T4_gene_ID)

#===split the gtf into table5 and revise gene range for all (all contain all ENSG and ENST)===
gtf=fread(paste0(path4,"table5pENST.gtf.gz"), header=F, stringsAsFactors = F)
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

#revise gene region start end
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
write.table(gtf5[order(gtf5$V1,gtf5$V4),],gzfile(paste0(path4,"table5.final.fulllength_detected.alone.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#========================================
#bash
#validate the gtf
#setwd(path4)
#system("zcat table5.final.fulllength_detected.alone.gtf.gz | bedparse gtf2bed | gzip > table5.final.fulllength_detected.alone.bed12.bed.gz")


#===============================================================================
#build final table 5 with all full-length and partial ENST -> for long-read and short-read quantification

table1=read.delim(paste0(path2,"table1.remove_undetect_gencode_and_internal_prime.2.51M.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5transcript1=table1$model_ID[which(table1$source %in% c("partial_gencode", "fulllength_gencode"))]
t5gene1=unique(table1$IN1_gene_ID[which(table1$source %in% c("partial_gencode", "fulllength_gencode"))])

transcript_info5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5transcript=unique(transcript_info5$model_ID)
t5gene=unique(transcript_info5$T4_gene_ID)

t5transcript=union(t5transcript,t5transcript1)
t5gene=union(t5gene,t5gene1)
#===split the gtf into table5 and revise gene range for all (all contain all ENSG and ENST)===
gtf=fread(paste0(path4,"table5pENST.gtf.gz"), header=F, stringsAsFactors = F)
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

#revise gene region start end
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
write.table(gtf5[order(gtf5$V1,gtf5$V4),],gzfile(paste0(path4,"table5.final.partial_yes_detected.alone.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#===============================================================================
#bash
#validate the gtf
#setwd(path4)
#system("zcat table5.final.partial_yes_detected.alone.gtf.gz | bedparse gtf2bed | gzip > table5.final.partial_yes_detected.alone.bed12.bed.gz")






