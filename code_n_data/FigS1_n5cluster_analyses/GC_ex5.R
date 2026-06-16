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
GC_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/GCcontent/")
GCex5_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/GCcontent/")

#===============================================================================
#bash
#prepare UCSC GC content bed file
#cd /analysisdata/fantom6/Interactome/resources/UCSC
#~/bigWigToBedGraph hg38.gc5Base.bw | gzip > hg38.gc5Base.bed.gz

setwd("/analysisdata/fantom6/Interactome/resources/UCSC")
system(paste0("cp hg38.gc5Base.bed.gz ",GCex5_path))

#===============================
setwd(GCex5_path)
cluster_summit=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/Neuron_THP1.S3.end5.summit.table5.bed"), header=F, stringsAsFactors = F)
cluster_summit1=cluster_summit
cluster_summit1$V2=cluster_summit1$V2-2000
cluster_summit1$V3=cluster_summit1$V3+2000
cluster_summit1$V2[which(cluster_summit1$V2 <0)]=0
write.table(cluster_summit1[order(cluster_summit1$V1,cluster_summit1$V2),], gzfile("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz"), row.names=F, col.names=F, sep="\t", quote=F)

#===============================
#bash
#bedtools intersect
setwd(GCex5_path)
system("bedtools intersect -wb -a Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz -b hg38.gc5Base.bed.gz | gzip > Neuron_THP1.S3.end5.summit.table5.2000bp_extend.gc.bed.gz")

#generate random seq -> exclude the tCRE 4001 bp region, exclude all exons from gencode
setwd(GCex5_path)
system("zcat Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz gencode.v39.annotation.all_exon.bed.gz | sort -k1,1 -k2,2n | gzip > end5cluster_gencode_exon.bed.gz")
system("bedtools shuffle -i Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz -g hg38.chrom.sizes.bed.gz -noOverlapping -excl end5cluster_gencode_exon.bed.gz | gzip > random_2kb_extend.bed.gz")
system("bedtools intersect -wa -wb -a random_2kb_extend.bed.gz -b hg38.gc5Base.bed.gz | gzip > random_2kb_extend.gc.bed.gz")

system("rm hg38.gc5Base.bed.gz")

#===============================================================================
#R
setwd(GCex5_path)
options(scipen=999)
end5cluster=read.delim("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.bed.gz",header=F, stringsAsFactors = F)
gc=fread("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.gc.bed.gz", header=F, stringsAsFactors = F)
gc=left_join(gc, end5cluster[,c(4,2)], by="V4", suffix=c("",".ori"))

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

gc$locS=gc$V2-gc$V2.ori-2000
gc$locE=gc$V3-gc$V2.ori-2000
gc1=gc[which(gc$V6 == "+" & gc$V2==gc$V8),c(1,12,13,4,10)]
gc2=gc[which(gc$V6 == "-" & gc$V3==gc$V9),c(1,13,12,4,10)]
gc2$locE=gc2$locE * (-1)
gc2$locS=gc2$locS * (-1)
gc2$V1="chr1"
gc1$V1="chr1"
colnames(gc2)=colnames(gc1)
gc1$locS=gc1$locS+2000
gc1$locE=gc1$locE+2000
gc2$locS=gc2$locS+2001
gc2$locE=gc2$locE+2001
gc=rbind(gc1,gc2)
colnames(gc)[c(4,5)]=c("n5_string","percent")
gc=left_join(gc, unique(table5[,c(62,64)]),by="n5_string",copy=F)
gc=gc[order(gc$locS),]
write.table(gc,gzfile("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.gc.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)

gcp=gc[which(gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]),]
gce=gc[which(gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]),]
gcu=gc[which(gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")]),]
gcc=gc[which(gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]),]

#add random seq GC content
random_gc=fread("random_2kb_extend.gc.bed.gz", header=F, stringsAsFactors = F)

random_gc$locS=random_gc$V8-random_gc$V2-2000
random_gc$locE=random_gc$V9-random_gc$V2-2000
random_gc1=random_gc[which(random_gc$V6 == "+" ),c(1,11,12,4,10)]
random_gc2=random_gc[which(random_gc$V6 == "-" ),c(1,12,11,4,10)]
random_gc2$locE=random_gc2$locE * (-1)
random_gc2$locS=random_gc2$locS * (-1)
random_gc2$V1="chr1"
random_gc1$V1="chr1"
colnames(random_gc2)=colnames(random_gc1)
random_gc1$locS=random_gc1$locS+2000
random_gc1$locE=random_gc1$locE+2000
random_gc2$locS=random_gc2$locS+2001
random_gc2$locE=random_gc2$locE+2001
random_gc=rbind(random_gc1,random_gc2)
colnames(random_gc)[c(4,5)]=c("n5_string","percent")
random_gc=left_join(random_gc, unique(table5[,c(62,64)]),by="n5_string",copy=F)
random_gc=random_gc[order(random_gc$locS),]
write.table(random_gc,gzfile("Neuron_THP1.S3.end5.summit.table5.2000bp_extend.random_gc.fakebed.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

random_gc$locS[which(random_gc$locS<0)]=0
random_gcp=random_gc[which(random_gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "promoter-like")]),]
random_gce=random_gc[which(random_gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "enhancer-like")]),]
random_gcu=random_gc[which(random_gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "unclassed")]),]
random_gcc=random_gc[which(random_gc$CREID %in% CREanno$CREID[which(CREanno$promoter_type == "CTCF-alone")]),]

fake2=data.frame(cbind("chr1",c(0:4000),c(1:4001),c(-2000:2000)))
fake2$GCpercent=0
colnames(fake2)=c("chr","start","end","position","GCprecent")

gcpCount=fake2
for (i in 1:4001){gcpCount$GCprecent[i]=mean(gcp$percent[which(gcp$locS<i & gcp$locE>=i)])}
gceCount=fake2
for (i in 1:4001){gceCount$GCprecent[i]=mean(gce$percent[which(gce$locS<i & gce$locE>=i)])}
gcuCount=fake2
for (i in 1:4001){gcuCount$GCprecent[i]=mean(gcu$percent[which(gcu$locS<i & gcu$locE>=i)])}
gccCount=fake2
for (i in 1:4001){gccCount$GCprecent[i]=mean(gcc$percent[which(gcc$locS<i & gcc$locE>=i)],na.rm=T)}
gcpCount$promoter_type="promoter-like"
gceCount$promoter_type="enhancer-like"
gcuCount$promoter_type="unclassed"
gccCount$promoter_type="CTCF-alone"
gcCount=rbind(gcpCount,gceCount,gcuCount,gccCount)
gcCount$group="anno_region"
write.table(gcCount,gzfile("end5_cluster_GCcontent_result.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)

random_gcpCount=fake2
for (i in 1:4001){random_gcpCount$GCprecent[i]=mean(random_gcp$percent[which(random_gcp$locS<i & random_gcp$locE>=i)])}
random_gceCount=fake2
for (i in 1:4001){random_gceCount$GCprecent[i]=mean(random_gce$percent[which(random_gce$locS<i & random_gce$locE>=i)])}
random_gcuCount=fake2
for (i in 1:4001){random_gcuCount$GCprecent[i]=mean(random_gcu$percent[which(random_gcu$locS<i & random_gcu$locE>=i)])}
random_gccCount=fake2
for (i in 1:4001){random_gccCount$GCprecent[i]=mean(random_gcc$percent[which(random_gcc$locS<i & random_gcc$locE>=i)],na.rm=T)}
random_gcpCount$promoter_type="promoter-like"
random_gceCount$promoter_type="enhancer-like"
random_gcuCount$promoter_type="unclassed"
random_gccCount$promoter_type="CTCF-alone"
random_gcCount=rbind(random_gcpCount,random_gceCount,random_gcuCount,random_gccCount)
random_gcCount$group="random"
write.table(random_gcCount,gzfile("end5_cluster_random_GCcontent_result.tsv.gz"),col.names=T, row.names=F, sep="\t", quote=F)

#===============================
data1=read.delim("end5_cluster_GCcontent_result.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
data2=read.delim("end5_cluster_random_GCcontent_result.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
data1=rbind(data1, data2)
write.table(data1,gzfile(paste0(path_fig2_data,"end5_cluster_both_GCcontent_result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# remove intermediate files
system("rm Neuron_THP1.S3.end5.summit.table5.2000bp_extend.gc.bed.gz")
system("rm random_2kb_extend.gc.bed.gz")
system("rm Neuron_THP1.S3.end5.summit.table5.2000bp_extend.gc.fakebed.tsv.gz")
system("rm Neuron_THP1.S3.end5.summit.table5.2000bp_extend.random_gc.fakebed.tsv.gz")

