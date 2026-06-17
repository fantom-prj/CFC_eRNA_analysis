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
path_fig3_data=paste0(primary_folder,"fig3/data/")
isoquant_path=paste0(primary_folder,"code_n_data/other_assemblers/isoQuant/assemblies/")
isoquant_path2=paste0(primary_folder,"code_n_data/other_assemblers/isoQuant/assemblies_stranded_corrected_internnal_priming_filtered/")

#==================================================
#isoQuant was run according to the "pipeline.txt" in [primary_folder]/code_n_data/, 
#using the all the libraries of Neuron-series and THP-1-series after transcriptClean.

#==================================================
setwd(isoquant_path)

sensitive.gtf=fread(paste0(isoquant_path, "isoquant_sensitive/Neuron_Series_THP1.transcript_models.gtf.gz"), header=F, skip=3, stringsAsFactors = F)
standard.gtf=fread(paste0(isoquant_path, "isoquant_standard/Neuron_Series_THP1.transcript_models.gtf.gz"), header=F, skip=3, stringsAsFactors = F)

length(which(sensitive.gtf$V7 == "+")) + length(which(sensitive.gtf$V7 == "-")) #827620
length(which(standard.gtf$V7 == "+")) + length(which(standard.gtf$V7 == "-")) #635815
#excluded two gene without strand
sensitive.gtf1=sensitive.gtf[which(sensitive.gtf$V7 != "+" & sensitive.gtf$V7 != "-"),]
sensitive.gtf=sensitive.gtf[which(sensitive.gtf$V7 == "+" | sensitive.gtf$V7 == "-"),]
write.table(sensitive.gtf,paste0(path1, "isoquant_sensitive/Neuron_Series_THP1.transcript_models_remove.no.strand.gtf.gz"), col.names=F, row.name=F, sep="\t", quote=F)
standard.gtf1=standard.gtf[which(standard.gtf$V7 != "+" & standard.gtf$V7 != "-"),]
standard.gtf=standard.gtf[which(standard.gtf$V7 == "+" | standard.gtf$V7 == "-"),]
write.table(standard.gtf,paste0(path1, "isoquant_standard/Neuron_Series_THP1.transcript_models_remove.no.strand.gtf.gz"), col.names=F, row.name=F, sep="\t", quote=F)

#novel gene and transcript ID is portable between sensitive and standard

#=========================================
#bash
#get bed 12
setwd(paste0(isoquant_path,"isoquant_sensitive/"))
system("zact Neuron_Series_THP1.transcript_models_remove.no.strand.gtf.gz | bedparse gtf2bed | gzip > Neuron_Series_THP1.transcript_models_remove.no.strand.bed12.bed.gz")
setwd(paste0(isoquant_path,"isoquant_standard/"))
system("zcat Neuron_Series_THP1.transcript_models_remove.no.strand.gtf.gz | bedparse gtf2bed | gzip > Neuron_Series_THP1.transcript_models_remove.no.strand.bed12.bed.gz")

setwd(isoquant_path2)
system("zcat isoquant_sensitive_int_priming.gtf.gz | bedparse gtf2bed | gzip > isoquant_sensitive_int_priming.bed12.bed.gz")
system("bed12ToBed6 -i isoquant_sensitive_int_priming.bed12.bed.gz | gzip > isoquant_sensitive_int_priming.exon.bed6.bed.gz")
system("zcat isoquant_sensitive_int_priming.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > isoquant_sensitive_int_priming.intron.bed6.bed.gz")
system("zcat isoquant_standard_int_priming.gtf.gz | bedparse gtf2bed | gzip > isoquant_standard_int_priming.bed12.bed.gz")
system("bed12ToBed6 -i isoquant_standard_int_priming.bed12.bed.gz | gzip > isoquant_standard_int_priming.exon.bed6.bed.gz")
system("zcat isoquant_standard_int_priming.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > isoquant_standard_int_priming.intron.bed6.bed.gz")

#===============================================================================
sensitive.gtf=read.delim(paste0(isoquant_path, "isoquant_sensitive/Neuron_Series_THP1.transcript_models_remove.no.strand.gtf"), header=F, stringsAsFactors = F)
standard.gtf=read.delim(paste0(isoquant_path, "isoquant_standard/Neuron_Series_THP1.transcript_models_remove.no.strand.gtf"), header=F, stringsAsFactors = F)

sensitive.gtf$geneID=sapply(strsplit(sensitive.gtf$V9,"; "),"[",1)
sensitive.gtf$geneID=gsub("gene_id ","",sensitive.gtf$geneID)
sensitive.gtfg=sensitive.gtf[which(sensitive.gtf$V3 == "gene"),]
sensitive.gtfng=sensitive.gtf[which(sensitive.gtf$V3 != "gene"),]
sensitive.gtfng$transcriptID=sapply(strsplit(sensitive.gtfng$V9,"; "),"[",2)
sensitive.gtfng$transcriptID=gsub("transcript_id ","",sensitive.gtfng$transcriptID)
sensitive.gtft=sensitive.gtfng[which(sensitive.gtfng$V3 == "transcript"),]
sensitive.gtfe=sensitive.gtfng[which(sensitive.gtfng$V3 == "exon"),]
sensitive.gtft$transcript_novelty="novel"
sensitive.gtft$transcript_novelty[grep("ENST",sensitive.gtft$transcriptID)]="GENCODE"
sensitive.gtft%>%group_by(transcript_novelty)%>%dplyr::summarise(count=n())

#==========================================
setwd(isoquant_path2)
sensitive.gtf2=read.delim("isoquant_sensitive_int_priming.gtf.gz", header=F, skip=3, stringsAsFactors = F)
sensitive.gtf2$geneID=sapply(strsplit(sensitive.gtf2$V9,"; "),"[",1)
sensitive.gtf2$geneID=gsub("gene_id ","",sensitive.gtf2$geneID)
sensitive.gtf2g=sensitive.gtf2[which(sensitive.gtf2$V3 == "gene"),]
sensitive.gtf2ng=sensitive.gtf2[which(sensitive.gtf2$V3 != "gene"),]
sensitive.gtf2ng$transcriptID=sapply(strsplit(sensitive.gtf2ng$V9,"; "),"[",2)
sensitive.gtf2ng$transcriptID=gsub("transcript_id ","",sensitive.gtf2ng$transcriptID)
sensitive.gtf2t=sensitive.gtf2ng[which(sensitive.gtf2ng$V3 == "transcript"),]
sensitive.gtf2e=sensitive.gtf2ng[which(sensitive.gtf2ng$V3 == "exon"),]
sensitive.gtf2t$transcript_novelty="novel"
sensitive.gtf2t$transcript_novelty[grep("ENST",sensitive.gtf2t$transcriptID)]="GENCODE"
sensitive.gtf2t%>%group_by(transcript_novelty)%>%dplyr::summarise(count=n())
sensitive.gtf2t$gene_novelty="novel"
sensitive.gtf2t$gene_novelty[grep("ENSG",sensitive.gtf2t$geneID)]="GENCODE"
unique(sensitive.gtf2t[,c(10,13)])%>%group_by(gene_novelty)%>%dplyr::summarise(count=n())

sensitive.gtf2t$n5_1base=sensitive.gtf2t$V4
sensitive.gtf2t$n5_1base[which(sensitive.gtf2t$V7=="-")]=sensitive.gtf2t$V5[which(sensitive.gtf2t$V7=="-")]
sensitive.gtf2t$n3_1base=sensitive.gtf2t$V5
sensitive.gtf2t$n3_1base[which(sensitive.gtf2t$V7=="-")]=sensitive.gtf2t$V4[which(sensitive.gtf2t$V7=="-")]

nrow(unique(sensitive.gtf2t[,c(1,14)])) #73613
nrow(unique(sensitive.gtf2t[,c(1,15)])) #65409

#===================================
standard.gtf$geneID=sapply(strsplit(standard.gtf$V9,"; "),"[",1)
standard.gtf$geneID=gsub("gene_id ","",standard.gtf$geneID)
standard.gtfg=standard.gtf[which(standard.gtf$V3 == "gene"),]
standard.gtfng=standard.gtf[which(standard.gtf$V3 != "gene"),]
standard.gtfng$transcriptID=sapply(strsplit(standard.gtfng$V9,"; "),"[",2)
standard.gtfng$transcriptID=gsub("transcript_id ","",standard.gtfng$transcriptID)
standard.gtft=standard.gtfng[which(standard.gtfng$V3 == "transcript"),]
standard.gtfe=standard.gtfng[which(standard.gtfng$V3 == "exon"),]
standard.gtft$transcript_novelty="novel"
standard.gtft$transcript_novelty[grep("ENST",standard.gtft$transcriptID)]="GENCODE"
standard.gtft%>%group_by(transcript_novelty)%>%dplyr::summarise(count=n())

setwd(isoquant_path2)
standard.gtf2=read.delim("isoquant_standard_int_priming.gtf.gz", header=F, skip=3, stringsAsFactors = F)
standard.gtf2$geneID=sapply(strsplit(standard.gtf2$V9,"; "),"[",1)
standard.gtf2$geneID=gsub("gene_id ","",standard.gtf2$geneID)
standard.gtf2g=standard.gtf2[which(standard.gtf2$V3 == "gene"),]
standard.gtf2ng=standard.gtf2[which(standard.gtf2$V3 != "gene"),]
standard.gtf2ng$transcriptID=sapply(strsplit(standard.gtf2ng$V9,"; "),"[",2)
standard.gtf2ng$transcriptID=gsub("transcript_id ","",standard.gtf2ng$transcriptID)
standard.gtf2t=standard.gtf2ng[which(standard.gtf2ng$V3 == "transcript"),]
standard.gtf2e=standard.gtf2ng[which(standard.gtf2ng$V3 == "exon"),]
standard.gtf2t$transcript_novelty="novel"
standard.gtf2t$transcript_novelty[grep("ENST",standard.gtf2t$transcriptID)]="GENCODE"
standard.gtf2t%>%group_by(transcript_novelty)%>%dplyr::summarise(count=n())
standard.gtf2t$gene_novelty="novel"
standard.gtf2t$gene_novelty[grep("ENSG",standard.gtf2t$geneID)]="GENCODE"
unique(standard.gtf2t[,c(10,13)])%>%group_by(gene_novelty)%>%dplyr::summarise(count=n())

standard.gtf2t$n5_1base=standard.gtf2t$V4
standard.gtf2t$n5_1base[which(standard.gtf2t$V7=="-")]=standard.gtf2t$V5[which(standard.gtf2t$V7=="-")]
standard.gtf2t$n3_1base=standard.gtf2t$V5
standard.gtf2t$n3_1base[which(standard.gtf2t$V7=="-")]=standard.gtf2t$V4[which(standard.gtf2t$V7=="-")]

nrow(unique(standard.gtf2t[,c(1,14)])) #63100
nrow(unique(standard.gtf2t[,c(1,15)])) #56501

bed6=read.delim("isoquant_sensitive_int_priming.exon.bed6.bed.gz", header=F, stringsAsFactors = F)
bed6$length=bed6$V3-bed6$V2
bed6a=bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(), length=sum(length))
sensitive.gtf2t=left_join(sensitive.gtf2t,bed6a, by=c("transcriptID"="V4"),copy=F)
sensitive.gtf2t$transcript_group="ENST"
sensitive.gtf2t$transcript_group[which(sensitive.gtf2t$transcript_novelty == "novel" & sensitive.gtf2t$gene_novelty == "novel")]="Transcript_from_novel_gene"
sensitive.gtf2t$transcript_group[which(sensitive.gtf2t$transcript_novelty == "novel" & sensitive.gtf2t$gene_novelty == "GENCODE")]="Novel_isoform"
write.table(sensitive.gtf2t,gzfile("sensitive_info_table.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

bed6=read.delim("isoquant_standard_int_priming.exon.bed6.bed.gz", header=F, stringsAsFactors = F)
bed6$length=bed6$V3-bed6$V2
bed6a=bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(), length=sum(length))
standard.gtf2t=left_join(standard.gtf2t,bed6a, by=c("transcriptID"="V4"),copy=F)
standard.gtf2t$transcript_group="ENST"
standard.gtf2t$transcript_group[which(standard.gtf2t$transcript_novelty == "novel" & standard.gtf2t$gene_novelty == "novel")]="Transcript_from_novel_gene"
standard.gtf2t$transcript_group[which(standard.gtf2t$transcript_novelty == "novel" & standard.gtf2t$gene_novelty == "GENCODE")]="Novel_isoform"
write.table(standard.gtf2t,gzfile("standard_info_table.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

sensitive=read.delim("sensitive_info_table.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
sensitive$n5_0base=sensitive$n5_1base-1
sensitive$n3_0base=sensitive$n3_1base-1
write.table(sensitive[order(sensitive$V1,sensitive$n5_0base),c(1,19,14,11,6,7,18)], gzfile("IsoQuant_Sensitive.n5.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(sensitive[order(sensitive$V1,sensitive$n3_0base),c(1,20,15,11,6,7,18)], gzfile("IsoQuant_Sensitive.n3.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

standard=read.delim("standard_info_table.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
standard$n5_0base=standard$n5_1base-1
standard$n3_0base=standard$n3_1base-1
write.table(standard[order(standard$V1,standard$n5_0base),c(1,19,14,11,6,7,18)], gzfile("IsoQuant_Default.n5.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(standard[order(standard$V1,standard$n3_0base),c(1,20,15,11,6,7,18)], gzfile("IsoQuant_Default.n3.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#intersect with external and data-driven 5' features
n5_path=paste0(primary_folder,"code_n_data/n5_regions/")
cluster_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")
setwd(isoquant_path2)
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


