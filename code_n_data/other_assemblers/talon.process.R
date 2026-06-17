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
TALON_path=paste0(primary_folder,"code_n_data/other_assemblers/TALON/")

#===============================================================================
#Neuron+THP1
#collapse per chromosome data into whole transcriptome

path3=paste0(TALON_path,"permissive/")
path4=paste0(TALON_path,"gtf/")
files=list.files(path=path3, pattern="_permissive_genomic_talon.gtf")

files.names=gsub("F6_interactome_","",files)
files.names=gsub("_permissive_genomic_talon.gtf.gz","",files.names)
data0=fread(paste0(path3,files[1]), header=F, stringsAsFactors = F)

data0a=data0[grep("TALON",data0$V9),]
write.table(data0a, gzfile(paste0(path4,files.names[1],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
for (i in 2: 10){
  data=fread(paste0(path3,files[i]), header=F, stringsAsFactors = F)
  dataa=data[which(data$V1 == files.names[i]),]
  write.table(dataa, gzfile(paste0(path4,files.names[i],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
  data0a=rbind(data0a,dataa)}

data=fread(paste0(path3,files[11]), header=F, stringsAsFactors = F)
dataa=data[which(data$V1 == "chr19" & data$V4 < 26000000),]
write.table(dataa, gzfile(paste0(path4,files.names[11],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
data0a=rbind(data0a,dataa)

data=fread(paste0(path3,files[12]), header=F, stringsAsFactors = F)
dataa=data[which(data$V1 == "chr19" & data$V4 > 26000000),]
write.table(dataa, gzfile(paste0(path4,files.names[12],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
data0a=rbind(data0a,dataa)

data=fread(paste0(path3,files[13]), header=F, stringsAsFactors = F)
dataa=data[which(data$V1 == "chr1" & data$V4 < 122050000),]
write.table(dataa, gzfile(paste0(path4,files.names[13],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
data0a=rbind(data0a,dataa)

data=fread(paste0(path3,files[14]), header=F, stringsAsFactors = F)
dataa=data[which(data$V1 == "chr1" & data$V4 > 122050000),]
write.table(dataa, gzfile(paste0(path4,files.names[14],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
data0a=rbind(data0a,dataa)

for (i in 15: length(files)){
  data=fread(paste0(path3,files[i]), header=F, stringsAsFactors = F)
  dataa=data[which(data$V1 == files.names[i]),]
  write.table(dataa, gzfile(paste0(path4,files.names[i],".gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
  data0a=rbind(data0a,dataa)}
write.table(data0a, gzfile(paste0(path4,"genecode_plus_table0.gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
#
data0a=fread(paste0(path4,"genecode_plus_table0.gtf.gz"), header=F, stringsAsFactors = F)
data0at1=data0a[-which(data0a$V3 == "gene"),]
data0at1$geneid=sapply(strsplit(data0at1$V9,";"),"[",1)
data0at1$transcriptid=sapply(strsplit(data0at1$V9,";"),"[",2)
data0at1$geneid=gsub("gene_id ","",data0at1$geneid)
data0at1$transcriptid=gsub(" transcript_id ","",data0at1$transcriptid)
data0at2=data0at1%>%group_by(transcriptid,V3,V7)%>%dplyr::summarise(count=n())
data0at2=data0at2%>%group_by(transcriptid)%>%dplyr::mutate(count2=n())
#find TALON011T000267817 having a starting exon with only 1 bp and annotated as "-" strand (while the transcript is positive)
#remove this transcript model from the gtf
data0b=data0a[-grep("TALON011T000267817",data0a$V9),]
                     
write.table(data0b, gzfile(paste0(path4,"genecode_plus_table0.gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)
data0b=data0b[which(nchar(data0b$V1)<=5),]
write.table(data0b, gzfile(paste0(path4,"genecode_plus_table1.main.chr.gtf.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#===================================
#bash
setwd(path4)
system("zcat genecode_plus_table1.main.chr.gtf.gz | bedparse gtf2bed | gzip > genecode_plus_table1.main.chr.bed12.bed.gz")
system("bed12ToBed6 -i genecode_plus_table1.main.chr.bed12.bed.gz | gzip > genecode_plus_table1.main.chr.bed16.bed.gz")

#===================================
#supporting read for each transcript model
summary=data.frame(matrix(nrow=28, ncol=3))
colnames(summary)=c("chr","before","after")
path2=paste0(TALON_path,"talon_read_annot/")
files=list.files(path=path2, pattern=".tsv.gz")
files.names=gsub(".tsv.gz","",files)
data=fread(paste0(path2,files[1]), header=T, stringsAsFactors = F)
summary$chr[1]=files.names[1]
summary$before[1]=nrow(data)

summary$after[1]=nrow(data)
data1=data%>%group_by(annot_gene_id,annot_transcript_id, gene_novelty, transcript_novelty, ISM_subtype, dataset)%>%dplyr::summarise(count=n())
data2=spread(data1, key=6, value=7)
data2[is.na(data2)]=0
data3=data[,c(1,2,13)]
for (i in 2:length(files)){
  data=fread(paste0(path2,files[i]), header=T, stringsAsFactors = F)
  summary$chr[i]=files.names[i]
  summary$before[i]=nrow(data)
  summary$after[i]=nrow(data)
  data1a=data%>%group_by(annot_gene_id,annot_transcript_id, gene_novelty, transcript_novelty, ISM_subtype, dataset)%>%dplyr::summarise(count=n())
  data2a=spread(data1a, key=6, value=7)
  data2a[is.na(data2a)]=0
  data2=rbind(data2,data2a)
  data3a=data[,c(1,2,13)]
  data3=rbind(data3,data3a)}
data4=data3[-which(data3$annot_transcript_id == "TALON011T000267817"),]
write.table(data4, paste0(TALON_path,"read_anno/combined_id_match.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
sum(summary$before)-sum(summary$after) #0

transcript.bed=read.delim(paste0(path4,"genecode_plus_table1.main.chr.bed12.bed.gz"), header=F, stringsAsFactors = F)
#number of detectable transcripts: 1,849,706
#number of non-detectable ENST added to gtf: 164,438
#1849706+172232

data5=data2[-which(data2$annot_transcript_id == "TALON011T000267817"),]
write.table(data5,gzfile(paste0(TALON_path,"read_anno/supporting_read.matrix.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#generate table1 and 2 from table0
table0=read.delim(paste0(TALON_path,"read_anno/supporting_read.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table1.bed12=read.delim(paste0(path4,"genecode_plus_table1.main.chr.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table1=table0[which(table0$annot_transcript_id %in% table1.bed12$V4),]
options(scipen=999)
table1.bed12$n5_start=table1.bed12$V2
table1.bed12$n5_end=table1.bed12$V3
table1.bed12$n5_start[which(table1.bed12$V6=="-")]=table1.bed12$V3[which(table1.bed12$V6=="-")]-1
table1.bed12$n5_end[which(table1.bed12$V6=="+")]=table1.bed12$V2[which(table1.bed12$V6=="+")]+1
table1.bed12$n3_start=table1.bed12$V2
table1.bed12$n3_end=table1.bed12$V3
table1.bed12$n3_start[which(table1.bed12$V6=="+")]=table1.bed12$V3[which(table1.bed12$V6=="+")]-1
table1.bed12$n3_end[which(table1.bed12$V6=="-")]=table1.bed12$V2[which(table1.bed12$V6=="-")]+1
colnames(table1.bed12)[c(1,2,3,6)]=c("chr","start","end","strand")
table1=left_join(table1,table1.bed12[,c(1,2,3,4,6,13:16)], by=c("annot_transcript_id"="V4"),copy=F)

bed6=read.delim(paste0(path4,"genecode_plus_table1.main.chr.bed16.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
bed6$length=bed6$V3-bed6$V2
bed6t=bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(),transcript_length=sum(length))
table1=left_join(table1,bed6t, by=c("annot_transcript_id"="V4"),copy=F)

#=========================
#internal priming for transcript 3n
setwd(paste0(TALON_path,"bed/"))
trans_n3=read.delim("transcript.n3.sorted.bed.gz", header=F, stringsAsFactors = F)
options(scipen=999)
trans_n3$V4=paste0(trans_n3$V1,"_",trans_n3$V2,"_",trans_n3$V3,"_",trans_n3$V6)
trans_n3=trans_n3%>%group_by(V1,V2,V3,V4,V6)%>%dplyr::summarise(V5=n())
trans_n3$V2=trans_n3$V2-20
trans_n3$V3=trans_n3$V3+20
trans_n3$V2[which(trans_n3$V2<0)]=0 #1 n3 from chroM, minus strand, end at the first bp, can be ignored
write.table(trans_n3[order(trans_n3$V1,trans_n3$V2),c(1:4,6,5)], gzfile("transcript.n3collapsed.CES41.sorted.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#==========================
#bash
#getfasta
setwd(paste0(TALON_path,"bed/"))
system("bedtools getfasta -s -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed transcript.n3collapsed.CES41.sorted.bed.gz | gzip > transcript.n3collapsed.CES41.sorted.FASTA.gz")
#please download the fasta file separately:
#wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

#======================================
setwd(paste0(TALON_path,"bed/"))
trans_n3=read.delim("transcript.n3collapsed.CES41.sorted.FASTA.gz", header=F, stringsAsFactors = F, check.names = F)
for (i in 1:41){trans_n3[,i+7]=substr(trans_n3$V7, start = i, stop = i)}
trans_n3$fracA08=(rowSums(trans_n3[,c(29:36)] == "A"))/8
trans_n3$fracA16=(rowSums(trans_n3[,c(29:44)] == "A"))/16
trans_n3$fracA20=(rowSums(trans_n3[,c(29:48)] == "A"))/20
trans_n3$internal_prime="no"
trans_n3$internal_prime[which(trans_n3$fracA08 > 0.75)]="fracA08"
trans_n3$internal_prime[which(trans_n3$fracA16 > 0.5)]="fracA16"
sum(trans_n3$V5[which(trans_n3$internal_prime == "no")])/sum(trans_n3$V5) #0.9623914
trans_n3=trans_n3[,-c(8:48)]
write.table(trans_n3, gzfile("transcript.n3collapsed.CES41.sorted.internal_prime.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

transcript=read.delim("transcript.n3.sorted.bed.gz", header=F, stringsAsFactors = F)
transcript$ID=paste0(transcript$V1,"_",transcript$V2,"_",transcript$V3,"_",transcript$V6)
internal_prime_gencode=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/n3_Ivano/n3_Gencode_internal_priming.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
transcript=left_join(transcript,trans_n3[,c(4,11)], by=c("ID"="V4"), copy=F)
transcript=left_join(transcript,internal_prime_gencode[,c(7,8)], by=c("ID"="label"), copy=F)
transcript$internal_prime[which(transcript$internal_prime == "no")]="No"
transcript$internal_prime[which(transcript$internal_prime != "No")]="Yes"
transcript$GENCODE[which(is.na(transcript$GENCODE))]="No"

table1=left_join(table1,transcript[,c(4,8,9)],by=c("annot_transcript_id"="V4"), copy=F)
table1%>%group_by(GENCODE,internal_prime)%>%dplyr::summarise(count=n())

table1$source="permissive_output"
table1$source[which(rowSums(table1[,c(6,7)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(8,9)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(10,11)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(12,13)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(14,15)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(16,17)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(18,19)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(20,21)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(22,23)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(24,25)] >=5)>1)]="standard_output"
table1$source[which(rowSums(table1[,c(26,27)] >=5)>1)]="standard_output"
table1$source[which(table1$transcript_novelty == "Known")]="GENCODE"
table1$source[which(rowSums(table1[,c(6:27)]) ==0)]="undetected_GENCODE"
table1%>%group_by(source)%>%dplyr::summarise(count=n())

setwd(TALON_path)
write.table(table1,gzfile("TALON_Neuron_THP1.table0.2M.no_alt.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
table1=read.delim("TALON_Neuron_THP1.table0.2M.no_alt.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
table1a=table1[-which(table1$internal_prime == "Yes" & table1$GENCODE == "No"),]
write.table(table1a,gzfile("TALON_Neuron_THP1.table1.2M.no_IP.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
table1a=read.delim("TALON_Neuron_THP1.table1.2M.no_IP.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
table2=table1a[which(table1a$source == "standard_output" | table1a$source == "GENCODE"),]
write.table(table2,gzfile("TALON_Neuron_THP1.table2.174k.5reads.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)
length(unique(table2$annot_gene_id))

# raw intermediate files, read-based files were removed.
# Transcript info table in [primary_folder]/code_n_data/other_assemblers/TALON

#===============================================================================
#gtf trimmer
#split the gtf into table1.noIP and table2 and revise gene range for all (all contain all ENSG and ENST)===
gtf=fread(paste0(path4,"genecode_plus_table1.main.chr.gtf.gz"), header=F, stringsAsFactors = F)
gtf.gene=gtf[which(gtf$V3=="gene"),]
gtf.nogene=gtf[which(gtf$V3!="gene"),]
gtf.gene$geneID=sapply(strsplit(gtf.gene$V9,";"),"[",1)
gtf.gene$geneID=gsub("gene_id \"","",gtf.gene$geneID)
gtf.gene$geneID=gsub("\"","",gtf.gene$geneID)
gtf.nogene$transcriptID=sapply(strsplit(gtf.nogene$V9,";"),"[",2)
gtf.nogene$transcriptID=gsub(" transcript_id \"","",gtf.nogene$transcriptID)
gtf.nogene$transcriptID=gsub("\"","",gtf.nogene$transcriptID)
gtf.t=gtf.nogene[which(gtf.nogene$V3=="transcript"),]
gtf.gene1=gtf.gene[union(which(gtf.gene$geneID %in% t1gene),grep("ENSG",gtf.gene$geneID)),]
gtf.gene2=gtf.gene[union(which(gtf.gene$geneID %in% t2gene),grep("ENSG",gtf.gene$geneID)),]
gtf.nogene1=gtf.nogene[union(which(gtf.nogene$transcriptID %in% t1transcript), grep("ENST",gtf.nogene$transcriptID)),]
gtf.nogene2=gtf.nogene[union(which(gtf.nogene$transcriptID %in% t2transcript), grep("ENST",gtf.nogene$transcriptID)),]
#revise gene region start end
gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.t1=gtf.t[union(which(gtf.t$transcriptID %in% t1transcript),grep("ENST",gtf.t$transcriptID)),]
gtf.t2=gtf.t[union(which(gtf.t$transcriptID %in% t2transcript),grep("ENST",gtf.t$transcriptID)),]
gtf.g=gtf.t%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.g1=gtf.t1%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.g2=gtf.t2%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene=left_join(gtf.gene,gtf.g,by="geneID", copy=F)
gtf.gene1=left_join(gtf.gene1,gtf.g1,by="geneID", copy=F)
gtf.gene2=left_join(gtf.gene2,gtf.g2,by="geneID", copy=F)
gtf.gene=gtf.gene[,c(1:3,11,12,6:9)]
gtf.gene1=gtf.gene1[,c(1:3,11,12,6:9)]
gtf.gene2=gtf.gene2[,c(1:3,11,12,6:9)]
colnames(gtf.gene)=colnames(gtf.nogene)[c(1:9)]
colnames(gtf.gene1)=colnames(gtf.nogene)[c(1:9)]
colnames(gtf.gene2)=colnames(gtf.nogene)[c(1:9)]

gtf=rbind(gtf.gene,gtf.nogene[,c(1:9)])
gtf1=rbind(gtf.gene1,gtf.nogene1[,c(1:9)])
gtf2=rbind(gtf.gene2,gtf.nogene2[,c(1:9)])
write.table(gtf[order(gtf$V1,gtf$V4),],gzfile(paste0(path4,"genecode_plus_table1.main.chr.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(gtf1[order(gtf1$V1,gtf1$V4),],gzfile(paste0(path4,"genecode_plus_table1.noPI.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(gtf2[order(gtf2$V1,gtf2$V4),],gzfile(paste0(path4,"genecode_plus_table2.noPI.5read.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#============================================================
#bash
#validate the gtf
setwd(path4)
system("zcat genecode_plus_table1.noPI.gtf.gz | bedparse gtf2bed | gzip > gencode_plus_table1.noIP.bed12.bed.gz")
system("zcat genecode_plus_table2.noPI.5read.gtf.gz | bedparse gtf2bed | gzip > gencode_plus_table2.bed12.bed.gz")

#===================================
#gtf trimmer2, without undetected ENST
t1transcript=unique(table1$annot_transcript_id)
t2transcript=unique(table2$annot_transcript_id)
t1gene=unique(table1$annot_gene_id)
t2gene=unique(table2$annot_gene_id)

gtf=fread(paste0(path4,"genecode_plus_table1.noPI.gtf.gz"), header=F, stringsAsFactors = F)
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
gtf.gene2=gtf.gene[which(gtf.gene$geneID %in% t2gene),]
gtf.nogene1=gtf.nogene[which(gtf.nogene$transcriptID %in% t1transcript),]
gtf.nogene2=gtf.nogene[which(gtf.nogene$transcriptID %in% t2transcript),]
#revise gene region start end
gtf.t$geneID=sapply(strsplit(gtf.t$V9,";"),"[",1)
gtf.t$geneID=gsub("gene_id \"","",gtf.t$geneID)
gtf.t$geneID=gsub("\"","",gtf.t$geneID)
gtf.t=gtf.t[,c(4,5,10,11)]
gtf.t1=gtf.t[which(gtf.t$transcriptID %in% t1transcript),]
gtf.t2=gtf.t[which(gtf.t$transcriptID %in% t2transcript),]
gtf.g1=gtf.t1%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))
gtf.g2=gtf.t2%>%group_by(geneID)%>%dplyr::summarise(gene_start=min(V4), gene_end=max(V5))

gtf.gene1=left_join(gtf.gene1,gtf.g1,by="geneID", copy=F)
gtf.gene2=left_join(gtf.gene2,gtf.g2,by="geneID", copy=F)
gtf.gene1=gtf.gene1[,c(1:3,11,12,6:9)]
gtf.gene2=gtf.gene2[,c(1:3,11,12,6:9)]
colnames(gtf.gene1)=colnames(gtf.nogene)[c(1:9)]
colnames(gtf.gene2)=colnames(gtf.nogene)[c(1:9)]

gtf1=rbind(gtf.gene1,gtf.nogene1[,c(1:9)])
gtf2=rbind(gtf.gene2,gtf.nogene2[,c(1:9)])
write.table(gtf1[order(gtf1$V1,gtf1$V4),],gzfile(paste0(path4,"TALON.table1.noPI.detectedENST.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(gtf2[order(gtf2$V1,gtf2$V4),],gzfile(paste0(path4,"TALON.table2.noPI.5read.detectedENST.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

# gtf files in [primary_folder]/code_n_data/other_assemblers/TALON/gtf
#============================================================
#bash
#validate the gtf
setwd(path4)
system("zcat TALON.table1.noPI.detectedENST.gtf.gz | bedparse gtf2bed | gzip > TALON.table1.noPI.detectedENST.bed12.bed.gz")
system("zcat TALON.table2.noPI.5read.detectedENST.gtf.gz | bedparse gtf2bed | gzip > TALON.table2.noPI.5read.detectedENST.bed12.bed.gz")
#system("bed12ToBed6 -i TALON.table1.noPI.detectedENST.bed12.bed.gz | gzip > TALON.table1.noPI.detectedENST.exon.bed6.bed.gz")
#system("zcat TALON.table1.noPI.detectedENST.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > TALON.table1.noPI.detectedENST.intron.bed6.bed.gz")
#system("bed12ToBed6 -i TALON.table2.noPI.5read.detectedENST.bed12.bed.gz | gzip > TALON.table2.noPI.5read.exon.detectedENST.bed6.bed.gz")
#system("zcat TALON.table2.noPI.5read.detectedENST.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > TALON.table2.noPI.5read.detectedENST.intron.bed6.bed.gz")

#===================================
TALONt1=read.delim(paste0(TALON_path,"TALON_Neuron_THP1.table0.2M.no_alt.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
TALONt1$transcript_group="ENST"
TALONt1$transcript_group[which(TALONt1$transcript_novelty != "Known" & TALONt1$gene_novelty != "Known")]="Transcript_from_novel_gene"
TALONt1$transcript_group[which(TALONt1$transcript_novelty != "Known" & TALONt1$gene_novelty == "Known")]="Novel_isoform"
TALONt1%>%group_by(transcript_group)%>%dplyr::summarise(n())

#prepare n5 feature
options(scipen=999)
bed12=read.delim(paste0(path4,"TALON.table1.noPI.detectedENST.bed12.bed.gz"), header=F, stringsAsFactors = F)
bed12$V3[which(bed12$V6=="+")]=bed12$V2[which(bed12$V6=="+")]+1
bed12$V2[which(bed12$V6=="-")]=bed12$V3[which(bed12$V6=="-")]-1
bed12=bed12[which(nchar(bed12$V1)<6),]
bed12=left_join(bed12, TALONt1[,c(2,41)], by=c("V4"="annot_transcript_id"), copy=F)
write.table(bed12[order(bed12$V1,bed12$V2),c(1:6,13)], gzfile(paste0(path4,"TALON_Raw.table1.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
bed12=read.delim(paste0(path4,"TALON.table2.noPI.5read.detectedENST.bed12.bed.gz"), header=F, stringsAsFactors = F)
bed12$V3[which(bed12$V6=="+")]=bed12$V2[which(bed12$V6=="+")]+1
bed12$V2[which(bed12$V6=="-")]=bed12$V3[which(bed12$V6=="-")]-1
bed12=bed12[which(nchar(bed12$V1)<6),]
bed12=left_join(bed12, TALONt1[,c(2,41)], by=c("V4"="annot_transcript_id"), copy=F)
write.table(bed12[order(bed12$V1,bed12$V2),c(1:6,13)], gzfile(paste0(path4,"TALON_Read_filtered.table2.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#intersect with external and data-driven 5' features
n5_path=paste0(primary_folder,"code_n_data/n5_regions/")
cluster_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")
setwd(path4)
cmd <- paste0(
  "for file in *n5.bed.gz; do ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-ELS.all.enhancer.sort.bed.gz | gzip > \"${file%.bed.gz}.Ecount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-PLS.all.promoter.sort.bed.gz | gzip > \"${file%.bed.gz}.Pcount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-CTCF.sort.bed.gz | gzip > \"${file%.bed.gz}.Ccount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"F6_CAT.promoter.bed.gz | gzip > \"${file%.bed.gz}.CAGEcount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"peaks.merged.all.main_chr.bed.gz | gzip > \"${file%.bed.gz}.ATACcount.bed.gz\" ;",
  "bedtools intersect -c -s -a \"$file\" -b ",cluster_path,"ontCAGE.Neuron_THP1.cluster.coord.bed.gz | gzip > \"${file%.bed.gz}.SCAFEcount.bed.gz\" ;",
  "done")
system(cmd)

# bed files in [primary_folder]/code_n_data/other_assemblers/TALON/gtf
#===================================



