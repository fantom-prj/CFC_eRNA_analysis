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
intersect_folder=paste0(primary_folder,"code_n_data/Fig6_GWAS_eQTL_SNP/intersect/")
path_fig6_data=paste0(primary_folder,"fig6/data/")


#GWAS with fine mapping, get the coordination from all SNP 151 from hg38
#========================================================================================
#use fine mapping GWAS downloaded from http://www.mulinlab.org/causaldb/index.html, it is from hg19
#use all_snp151_USCS_hg38 downloaded from UCSC, cross over the rsid 
#downloadable files are not provided
options(scipen=999)
fine=fread("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/gwas_eqtl/chung_gwas/v2.0/credible_set.txt", header=T, stringsAsFactors = F, check.names = F)
meta=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/gwas_eqtl/chung_gwas/v2.0/meta.txt", header=T, stringsAsFactors = F, check.names = F)
hg38snp=fread("/home/hon-chun/resources/genome/human/inUse/hg38/snp/all_snp151_USCS_hg38_chrom_chromEnd_Name.txt.gz", header=F, nrows = 100)
hg38snp1=hg38snp[which(hg38snp$V3 %in% unique(fine$rsid)),]
hg38snp1=unique(hg38snp1)
hg38snp1=hg38snp1[which(nchar(hg38snp1$V1) <=5),]
length(unique(hg38snp1$V3))
fine=left_join(fine, hg38snp1, by=c("rsid"="V3"),copy=F)
fine=fine[which(!is.na(fine$V2)),]
#give up 16685 SNP
rm(hg38snp)
colnames(fine)[c(25,26)]=c("hg38_chr","hg38_V3")
fine$hg38_V2=fine$hg38_V3-1

fine$ID=paste0(fine$hg38_chr,"_",fine$hg38_V2,"_",fine$hg38_V3)
fine1=unique(fine[order(fine$hg38_chr,fine$hg38_V2),c(25,27,26,28)])
write.table(fine1, paste0(intersect_folder,"credible_set.hg38.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
write.table(fine, gzfile(paste0(intersect_folder,"credible_set.hg38.txt.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#trimmed GWAS
meta_brain=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/gwas_eqtl/chung_gwas/v2.0/brain_neuron_related_traits.tsv", header=T, stringsAsFactors = F, check.names = F)
fine_trim=fine[which(fine$meta_id %in% unique(meta_brain$meta_id)),]
fine_trim1=unique(fine_trim[order(fine_trim$hg38_chr,fine_trim$hg38_V2),c(25,27,26,28)])
write.table(fine_trim1, gzfile(paste0(intersect_folder,"credible_set.brain_hg38.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(fine_trim, gzfile(paste0(intersect_folder,"credible_set.brain_hg38.txt.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#process eQTL SNP
#data from https://www.ebi.ac.uk/eqtl/
#downloadable files are not provided

metadata=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/tabix_ftp_paths.tsv", header=T, stringsAsFactors = F, check.names = F)
metadata1=metadata[union(grep("brain",metadata$tissue_label),grep("iPSC",metadata$tissue_label)),]

path3=("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/download/")
files=list.files(path=path3, pattern="tsv.gz")
metadata2=data.frame(files)
metadata2$dataset_id=sapply(strsplit(metadata2$files,"\\."),"[",5)
metadata2$quant_method=sapply(strsplit(metadata2$files,"\\."),"[",3)
metadata2=left_join(metadata2, metadata, by=c("dataset_id","quant_method"), copy=F)
metadata2a=metadata2
metadata2=metadata2[which(metadata2$dataset_id %in% metadata1$dataset_id),]

data=read.delim(paste0(path3,metadata2$files[1]), header=T, stringsAsFactors = F, check.names = F)
data=data[which(data$pip>0.3),]
data$dataset_id=metadata2$dataset_id[1]
data$quant_method=metadata2$quant_method[1]
for (i in 2: nrow(metadata2)){
  data1=read.delim(paste0(path3,metadata2$files[i]), header=T, stringsAsFactors = F, check.names = F)
  data1=data1[which(data1$pip>0.3),]
  data1$dataset_id=metadata2$dataset_id[i]
  data1$quant_method=metadata2$quant_method[i]
  data=rbind(data,data1)}
data$chr=sapply(strsplit(data$variant,"_"),"[",1)
data$end=sapply(strsplit(data$variant,"_"),"[",2)
data$start=as.numeric(data$end)-1
data$ID=paste0(data$chr,"_",data$start,"_",data$end)
data=left_join(data, metadata2[,c(2,6)], by="dataset_id", copy=F)
datak=data%>%group_by(ID,variant, pip, gene_id,sample_group,quant_method)%>%dplyr::summarise(rsid=paste(rsid,collapse=";"))
datab=unique(datak[,c(1,4)])%>%group_by(ID)%>%dplyr::summarise(no_gene=n())
datac=unique(datak[,c(1,5)])%>%group_by(ID)%>%dplyr::summarise(no_sample=n())
dataa=unique(datak[,c(1,2)])%>%group_by(ID)%>%dplyr::summarise(no_variant=n())
datad=unique(datak[,c(1,6)])%>%group_by(ID)%>%dplyr::summarise(no_method=n())
data2=datak%>%group_by(ID)%>%dplyr::summarise(variant=paste(variant, collapse=";"),pip=paste(pip, collapse=";"),gene_id=paste(gene_id,collapse=";"),sample_group=paste(sample_group, collapse=";"),quant_method=paste(quant_method, collapse=";"))
data2=left_join(data2, dataa,by="ID",copy=F)
data2=left_join(data2, datab,by="ID",copy=F)
data2=left_join(data2, datac,by="ID",copy=F)
data2=left_join(data2, datad,by="ID",copy=F)
data3=unique(data[,c(16,18,17,19)])

write.table(data, gzfile(paste0(intersect_folder,"eqtl_LD_pip3.tsv.gz")), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data3, gzfile(paste0(intersect_folder,"eqtl_LD_pip3.collapsed.bed.gz")), col.names=F, row.names=F, quote=F, sep="\t")

#===============================================================================
#generate eQTL file for all cell-types
data=read.delim(paste0(path3,metadata2a$files[1]), header=T, stringsAsFactors = F, check.names = F)
data=data[which(data$pip>0.3),]
data$dataset_id=metadata2a$dataset_id[1]
data$quant_method=metadata2a$quant_method[1]
for (i in 2: nrow(metadata2a)){
  data1=read.delim(paste0(path3,metadata2a$files[i]), header=T, stringsAsFactors = F, check.names = F)
  data1=data1[which(data1$pip>0.3),]
  data1$dataset_id=metadata2a$dataset_id[i]
  data1$quant_method=metadata2a$quant_method[i]
  data=rbind(data,data1)}
data$chr=sapply(strsplit(data$variant,"_"),"[",1)
data$end=sapply(strsplit(data$variant,"_"),"[",2)
data$start=as.numeric(data$end)-1
data$ID=paste0(data$chr,"_",data$start,"_",data$end)
data=left_join(data, metadata2a[,c(2,6)], by="dataset_id", copy=F)
datak=data%>%group_by(ID,variant, pip, gene_id,sample_group,quant_method)%>%dplyr::summarise(rsid=paste(rsid,collapse=";"))
datab=unique(datak[,c(1,4)])%>%group_by(ID)%>%dplyr::summarise(no_gene=n())
datac=unique(datak[,c(1,5)])%>%group_by(ID)%>%dplyr::summarise(no_sample=n())
dataa=unique(datak[,c(1,2)])%>%group_by(ID)%>%dplyr::summarise(no_variant=n())
datad=unique(datak[,c(1,6)])%>%group_by(ID)%>%dplyr::summarise(no_method=n())
data2=datak%>%group_by(ID)%>%dplyr::summarise(variant=paste(variant, collapse=";"),pip=paste(pip, collapse=";"),gene_id=paste(gene_id,collapse=";"),sample_group=paste(sample_group, collapse=";"),quant_method=paste(quant_method, collapse=";"))
data2=left_join(data2, dataa,by="ID",copy=F)
data2=left_join(data2, datab,by="ID",copy=F)
data2=left_join(data2, datac,by="ID",copy=F)
data2=left_join(data2, datad,by="ID",copy=F)
data3=unique(data[,c(16,18,17,19)])

write.table(data, gzfile(paste0(intersect_folder,"allsamples_eqtl_LD_pip3.tsv.gz")), col.names=T, row.names=F, quote=F, sep="\t")
write.table(data3, gzfile(paste0(intersect_folder,"allsamples_eqtl_LD_pip3.collapsed.bed.gz")), col.names=F, row.names=F, quote=F, sep="\t")

#===============================================================================
#bash
#bedtools intersect the bed6 exon from table5 plus all ENST (v39)

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/"))
system(paste0("bed12ToBed6 -i table5pENST.bed12.bed.gz | gzip > ",intersect_folder,"table5pENST.bed6.bed.gz"))

setwd(intersect_folder)
#eqtl, brain or iPSC
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.bed6.bed.gz -b ",intersect_folder,"eqtl_LD_pip3.collapsed.bed.gz | gzip > ",intersect_folder,"table5pENST.bed6.eqtl_finemap.bed.gz"))
#eqtl3, all celltypes
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.bed6.bed.gz -b ",intersect_folder,"allsamples_eqtl_LD_pip3.collapsed.bed.gz | gzip > ",intersect_folder,"table5pENST.bed6.allsamples_eqtl_finemap.bed.gz"))
#GWAS
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.bed6.bed.gz -b ",intersect_folder,"credible_set.hg38.bed.gz | gzip > ",intersect_folder,"table5pENST.bed6.fine.bed.gz"))
#trimmed GWAS
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.bed6.bed.gz -b ",intersect_folder,"credible_set.brain_hg38.bed.gz | gzip > ",intersect_folder,"table5pENST.bed6.fine_brain.bed.gz"))

#===============================================================================
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

t5ENSTbed12eqtl=read.delim(paste0(intersect_folder,"table5pENST.bed6.eqtl_finemap.bed.gz"), header=F, stringsAsFactors = F)
t5ENSTbed12eqtl=unique(t5ENSTbed12eqtl[,c(4,10)])
t5ENSTbed12eqtl$group=substr(t5ENSTbed12eqtl$V4, 1, 4)

t5ENSTbed12eqtl=left_join(t5ENSTbed12eqtl,table5[,c("model_ID","transcript_novelty")], by=c("V4"="model_ID"),copy=F)
t5ENSTbed12eqtl$transcript_novelty[which(is.na(t5ENSTbed12eqtl$transcript_novelty))]="ENST"
write.table(t5ENSTbed12eqtl,gzfile(paste0(path_fig6_data,"eqtl_table5pENST.intersect.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#exon-eqtl file located in [primary_folder]/fig6/data/

#===============================================================================
t5ENSTbed12eqtl_all=read.delim(paste0(intersect_folder,"table5pENST.bed6.allsamples_eqtl_finemap.bed.gz"), header=F, stringsAsFactors = F)
t5ENSTbed12eqtl_all=unique(t5ENSTbed12eqtl_all[,c(4,10)])
t5ENSTbed12eqtl_all$group=substr(t5ENSTbed12eqtl_all$V4, 1, 4)

t5ENSTbed12eqtl_all=left_join(t5ENSTbed12eqtl_all,table5[,c("model_ID","transcript_novelty")], by=c("V4"="model_ID"),copy=F)
t5ENSTbed12eqtl_all$transcript_novelty[which(is.na(t5ENSTbed12eqtl_all$transcript_novelty))]="ENST"
write.table(t5ENSTbed12eqtl_all,gzfile(paste0(path_fig6_data,"all_eqtl_table5pENST.intersect.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#exon-eqtl_allsample file located in [primary_folder]/fig6/data/

#===============================================================================
t5ENSTbed12fine=read.delim(paste0(intersect_folder,"table5pENST.bed6.fine.bed.gz"), header=F, stringsAsFactors = F)
t5ENSTbed12fine=unique(t5ENSTbed12fine[,c(4,10)])
t5ENSTbed12fine$group=substr(t5ENSTbed12fine$V4, 1, 4)

t5ENSTbed12fine=left_join(t5ENSTbed12fine,table5[,c("model_ID","transcript_novelty")], by=c("V4"="model_ID"),copy=F)
t5ENSTbed12fine$transcript_novelty[which(is.na(t5ENSTbed12fine$transcript_novelty))]="ENST"
write.table(t5ENSTbed12fine,gzfile(paste0(path_fig6_data,"GWAS_table5pENST.intersect.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#exon-GWAS file located in [primary_folder]/fig6/data/

#===============================================================================
t5ENSTbed12fine_trim=read.delim(paste0(intersect_folder,"table5pENST.bed6.fine_brain.bed.gz"), header=F, stringsAsFactors = F)
t5ENSTbed12fine_trim=unique(t5ENSTbed12fine_trim[,c(4,10)])
t5ENSTbed12fine_trim$group=substr(t5ENSTbed12fine_trim$V4, 1, 4)

t5ENSTbed12fine_trim=left_join(t5ENSTbed12fine_trim,table5[,c("model_ID","transcript_novelty")], by=c("V4"="model_ID"),copy=F)
t5ENSTbed12fine_trim$transcript_novelty[which(is.na(t5ENSTbed12fine_trim$transcript_novelty))]="ENST"
write.table(t5ENSTbed12fine_trim,gzfile(paste0(path_fig6_data,"GWAS_trim_table5pENST.intersect.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#exon-trimmed GWAS file located in [primary_folder]/fig6/data/

#===============================================================================
#===============================================================================
#extract intersect result and put together
#eQTL and GWAS SNP at SJ by spliceAI (generated from [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp.R)
t5SJ_both=read.delim(paste0(path_fig6_data,"spliceAI_eQTL_GWAS_acceprot_donor_gain_loss.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5SJ_both=t5SJ_both[which(t5SJ_both$region %in% c("Donor","Acceptor")),]%>%group_by(region, group,SNP_pos) %>% dplyr::summarise(transcript_novelty=paste(unique(transcript_novelty),collapse=";"))
t5SJ_both$transcript_novelty[grep("ENST",t5SJ_both$transcript_novelty)]="ENST"
t5SJ_both$transcript_novelty[grep("Novel isoform",t5SJ_both$transcript_novelty)]="Novel isoform"

#GWAS SNP in exon by intersect
t5ENST_gwas=read.delim(paste0(path_fig6_data,"GWAS_table5pENST.intersect.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
t5ENST_gwas=t5ENST_gwas%>%group_by(V10) %>% dplyr::summarise(transcript_novelty=paste(unique(transcript_novelty),collapse=";"))
t5ENST_gwas$transcript_novelty[grep("ENST",t5ENST_gwas$transcript_novelty)]="ENST"
t5ENST_gwas$transcript_novelty[grep("Novel isoform",t5ENST_gwas$transcript_novelty)]="Novel isoform"
t5ENST_gwas$group="GWAS"
t5ENST_gwas$region="Exon"
t5ENST_gwas=t5ENST_gwas[,c(4,3,1,2)]
colnames(t5ENST_gwas)=colnames(t5SJ_both)

#eQTL SNP in exon by intersect
t5ENST_eqtl=read.delim(paste0(path_fig6_data,"eqtl_table5pENST.intersect.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
t5ENST_eqtl=t5ENST_eqtl%>%group_by(V10) %>% dplyr::summarise(transcript_novelty=paste(unique(transcript_novelty),collapse=";"))
t5ENST_eqtl$transcript_novelty[grep("ENST",t5ENST_eqtl$transcript_novelty)]="ENST"
t5ENST_eqtl$transcript_novelty[grep("Novel isoform",t5ENST_eqtl$transcript_novelty)]="Novel isoform"
t5ENST_eqtl$group="eQTL"
t5ENST_eqtl$region="Exon"
t5ENST_eqtl=t5ENST_eqtl[,c(4,3,1,2)]
colnames(t5ENST_eqtl)=colnames(t5SJ_both)

#eQTL_allsample SNP in exon by intersect
t5ENST_eqtlall=read.delim(paste0(path_fig6_data,"all_eqtl_table5pENST.intersect.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
t5ENST_eqtlall=t5ENST_eqtlall%>%group_by(V10) %>% dplyr::summarise(transcript_novelty=paste(unique(transcript_novelty),collapse=";"))
t5ENST_eqtlall$transcript_novelty[grep("ENST",t5ENST_eqtlall$transcript_novelty)]="ENST"
t5ENST_eqtlall$transcript_novelty[grep("Novel isoform",t5ENST_eqtlall$transcript_novelty)]="Novel isoform"
t5ENST_eqtlall$group="eQTL_all"
t5ENST_eqtlall$region="Exon"
t5ENST_eqtlall=t5ENST_eqtlall[,c(4,3,1,2)]
colnames(t5ENST_eqtlall)=colnames(t5SJ_both)

data=rbind(t5SJ_both, t5ENST_gwas, t5ENST_eqtl, t5ENST_eqtlall)
write.table(data, gzfile(paste0(path_fig6_data,"eQTL_GWAS_table5_exon_acceptor_donor.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data/


#===============================================================================
#gwas eqtl only common snp version - for comparison
#for Fig. 6c & ex8a
#generate bed file from vcf
#downloadable files are not provided

commonall=fread("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/SNP_background/common_all_20180418.vcf.gz", header=T)
commonall$info1=sapply(strsplit(commonall$INFO,"CAF="),"[",2)
commonall$info1=sapply(strsplit(commonall$info1,";"),"[",1)
commonall$info1=sapply(strsplit(commonall$info1,","),"[",1)
options(scipen=999)
commonall$start=commonall$POS-1
commonall$chr=paste0("chr",commonall$`#CHROM`)
write.table(unique(commonall[order(commonall$chr,commonall$start),c(11,10,2)]), gzfile(paste0(intersect_folder, "common_all_20180418.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
commonall$V4=paste0(commonall$chr,"_",commonall$start,"_",commonall$POS)
commonall=unique(commonall[,c(1:3,12)])

#common all unique bed located in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/intersect
#===============================================================================


#===============================================================================
#prepared merged exon grouped by ex5_cluster
table5a=table5[which(table5$ex5cluster_class != "others"),]
table5b=table5a%>%group_by(n5_string)%>%dplyr::summarise(max.length=max(transcript_length))
table5c=table5[which(table5$ex5cluster_class != "others"),c(1,62)]
t5bed6=read.delim(paste0(intersect_folder,"table5pENST.bed6.bed.gz"), header=F, stringsAsFactors = F)
t5bed6=t5bed6[which(t5bed6$V4 %in% table5c$model_ID),]
t5bed6=left_join(t5bed6,table5c, by=c("V4"="model_ID"), copy=F)
t5bed6=t5bed6%>%group_by(n5_string)%>%dplyr::mutate(count=n())
write.table(t5bed6[order(t5bed6$V1,t5bed6$V2),c(1,2,3,7,5,6)],gzfile(paste0(intersect_folder,"table5pENST.for_merge.bed6.bed.gz")),row.name =F, col.names = F, sep="\t", quote=F)

#===============================================================================
#bash
#merge the exons if they overlap
system(paste0("sh ",intersect_folder,"merge_exon.sh"))
system("gzip table5pENST.merged_by_ex5cluster.bed6.bed")
# this takes time

#===============================================================================

merged=read.delim(paste0(intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.bed.gz"), header=F, stringsAsFactors = F)
merged$length=merged$V3-merged$V2

#===============================================================================
#bash
#intersect merged exon grouped by ex5, and ex5_cluster with background SNP (all_common)
#for ex5_cluster, intersect with all here but restricted to those from Final SALA -> paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")

system(paste0("bedtools intersect -wa -wb -a ",intersect_folder, "table5pENST.merged_by_ex5cluster.bed6.bed.gz -b ",intersect_folder, "common_all_20180418.bed.gz | gzip > ",intersect_folder, "table5pENST.merged_by_ex5cluster.eqtl_common20180418_background.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder, "Neuron_THP1.S3.end5.bed.gz -b ",intersect_folder, "common_all_20180418.bed.gz | gzip > ",intersect_folder, "Neuron_THP1.S3.end5.eqtl_common20180418_background.bed.gz"))
#===============================================================================


#===============================================================================
#add background from common to the data
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

cluster1=fread(paste0(intersect_folder, "Neuron_THP1.S3.end5.eqtl_common20180418_background.bed.gz"), header=F, stringsAsFactors = F)
cluster1$V16=paste0(cluster1$V13,"_",cluster1$V14,"_",cluster1$V15)
cluster11=unique(cluster1[,c(4,16)])%>%group_by(V4)%>%dplyr::summarise(n_common_SNP=n())
data1=left_join(data1, cluster11,by=c("n5_string"="V4"),copy=F)
data1$n_common_SNP[which(is.na(data1$n_common_SNP))]=0

exon1=fread(paste0(intersect_folder, "table5pENST.merged_by_ex5cluster.eqtl_common20180418_background.bed.gz"), header=F, stringsAsFactors = F)
exon1$V9=paste0(exon1$V6,"_",exon1$V7,"_",exon1$V8)
exon11=unique(exon1[,c(4,9)])%>%group_by(V4)%>%dplyr::summarise(n_common_SNP=n())
data1=left_join(data1, exon11, by=c("n5_string"="V4"),copy=F, suffix=c("_ex5","_exon"))
data1$n_common_SNP_exon[which(is.na(data1$n_common_SNP_exon))]=0

#===============================================================================
#get the common only gwas and eqtl
options(scipen=999)
fine=fread(paste0(intersect_folder,"credible_set.hg38.txt.gz"), header=T, stringsAsFactors = F)
fine1=fine[which(fine$rsid %in% commonall$ID),]
fine1=fine1[which(fine1$ID %in% commonall$V4),]
fine1=unique(fine1[order(fine1$hg38_chr,fine1$hg38_V2),c(25,27,26,28)])
length(unique(fine$ID)) #917107
length(unique(fine1$ID)) #855819
write.table(fine1, gzfile(paste0(intersect_folder,"credible_set.commononly.hg38.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

fine_trim=fread(paste0(intersect_folder,"credible_set.brain_hg38.txt.gz"), header=T, stringsAsFactors = F)
fine_trim1=fine_trim[which(fine_trim$rsid %in% commonall$ID),]
fine_trim1=fine_trim1[which(fine_trim1$ID %in% commonall$V4),]
fine_trim1=unique(fine_trim1[order(fine_trim1$hg38_chr,fine_trim1$hg38_V2),c(25,27,26,28)])
length(unique(fine_trim$ID)) #36197
length(unique(fine_trim1$ID)) #34404
write.table(fine_trim1, gzfile(paste0(intersect_folder,"credible_set.brain.commononly.hg38.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

data=read.delim(paste0(intersect_folder,"eqtl_LD_pip3.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
dataa=data[which(data$rsid %in% commonall$ID),]
dataa=dataa[which(dataa$ID %in% commonall$V4),]
dataa=unique(dataa[order(dataa$chr,dataa$start),c(16,18,17,19)])
length(unique(data$ID)) #92918
length(unique(dataa$ID)) #81782
write.table(dataa, gzfile(paste0(intersect_folder,"eqtl_LD_pip3.collapsed.commononly.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

data=read.delim(paste0(intersect_folder,"allsamples_eqtl_LD_pip3.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
dataa=data[which(data$rsid %in% commonall$ID),]
dataa=dataa[which(dataa$ID %in% commonall$V4),]
dataa=unique(dataa[order(dataa$chr,dataa$start),c(16,18,17,19)])
length(unique(data$ID)) #367989
length(unique(dataa$ID)) #320993
write.table(dataa, gzfile(paste0(intersect_folder,"allsamples_eqtl_LD_pip3.collapsed.commononly.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#intersect with ex5_cluster
#for ex5_cluster, intersect with all here but restricted to those from Final SALA -> paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")

setwd(intersect_folder)
#eqtl3_ipsc_brain
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"Neuron_THP1.S3.end5.bed.gz -b ",intersect_folder,"eqtl_LD_pip3.collapsed.commononly.bed.gz | gzip > ",intersect_folder,"Neuron_THP1.S3.end5.eqtl_finemap.commononly.bed.gz"))
#eqtl3_allcelltype
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"Neuron_THP1.S3.end5.bed.gz -b ",intersect_folder,"allsamples_eqtl_LD_pip3.collapsed.commononly.bed.gz | gzip > ",intersect_folder,"Neuron_THP1.S3.end5.allsamples_eqtl_finemap.commononly.bed.gz"))
#GWAS_trimmed
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"Neuron_THP1.S3.end5.bed.gz -b ",intersect_folder,"credible_set.brain.commononly.hg38.bed.gz | gzip > ",intersect_folder,"Neuron_THP1.S3.end5.fine_brain.commononly.bed.gz"))
#GWAS_all
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"Neuron_THP1.S3.end5.bed.gz -b ",intersect_folder,"credible_set.commononly.hg38.bed.gz | gzip > ",intersect_folder,"Neuron_THP1.S3.end5.fine.commononly.bed.gz"))

#intersect with exon grouped by ex5_cluster
setwd(intersect_folder)
#eqtl3_ipsc_brain
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.bed.gz -b ",intersect_folder,"eqtl_LD_pip3.collapsed.commononly.bed.gz | gzip > ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.eqtl_finemap.commononly.bed.gz"))
#eqtl3_allcelltype
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.bed.gz -b ",intersect_folder,"allsamples_eqtl_LD_pip3.collapsed.commononly.bed.gz | gzip > ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.allsamples_eqtl_finemap.commononly.bed.gz"))
#GWAS_trimmed
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.bed.gz -b ",intersect_folder,"credible_set.brain.commononly.hg38.bed.gz | gzip > ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.fine_brain.commononly.bed.gz"))
#GWAS_all
system(paste0("bedtools intersect -wa -wb -a ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.bed.gz -b ",intersect_folder,"credible_set.commononly.hg38.bed.gz | gzip > ",intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.fine.commononly.bed.gz"))


#===============================================================================
# extract intersect data (common only) -> eQTL use iPSC and brain samples alone

eqtln5=read.delim(paste0(intersect_folder,"Neuron_THP1.S3.end5.eqtl_finemap.commononly.bed.gz"),header=F, stringsAsFactors = F)
gwasn5_trim=read.delim(paste0(intersect_folder,"Neuron_THP1.S3.end5.fine_brain.commononly.bed.gz"),header=F, stringsAsFactors = F)
eqtln51=unique(eqtln5[,c(4,16)])%>%group_by(V4)%>%dplyr::summarise(n_eQTL_common_ex5=n())
gwasn5_trim1=unique(gwasn5_trim[,c(4,16)])%>%group_by(V4)%>%dplyr::summarise(n_GWAS_common_ex5=n())

path_fig4_data=paste0(primary_folder,"fig4/data/")
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=left_join(data1, eqtln51, by=c("n5_string"="V4"), copy=F)
data1=left_join(data1, gwasn5_trim1, by=c("n5_string"="V4"), copy=F)
#
t5ENSTbed12eqtl=read.delim(paste0(intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.eqtl_finemap.commononly.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
t5ENSTbed12fine_trim=read.delim(paste0(intersect_folder,"table5pENST.merged_by_ex5cluster.bed6.fine_brain.commononly.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
t5ENSTbed12eqtl1=unique(t5ENSTbed12eqtl[,c(4,9)])%>%group_by(V4)%>%dplyr::summarise(n_eQTL_common_exon=n())
t5ENSTbed12fine_trim1=unique(t5ENSTbed12fine_trim[,c(4,9)])%>%group_by(V4)%>%dplyr::summarise(n_GWAS_common_exon=n())
data1=left_join(data1, t5ENSTbed12eqtl1, by=c("n5_string"="V4"), copy=F)
data1=left_join(data1, t5ENSTbed12fine_trim1, by=c("n5_string"="V4"), copy=F)
write.table(data1,paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#==============
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2=data1[,c("n5_string","ex5cluster_class","ex5_length","exon_coverage","n_eQTL_common_ex5","n_GWAS_common_ex5","n_common_SNP_ex5","n_eQTL_common_exon","n_GWAS_common_exon","n_common_SNP_exon","CpGTATA")]
colnames(data2)=c("n5_string","ex5cluster_class","ex5_length","exon_coverage","n_eQTL_ex5","n_GWAS_ex5","n_commonSNP_ex5","n_eQTL_exon","n_GWAS_exon","n_commonSNP_exon","CpGTATA")
data2[is.na(data2)]=0
data2$n_eQTL_ex5_rate=data2$n_eQTL_ex5/data2$n_commonSNP_ex5
data2$n_GWAS_ex5_rate=data2$n_GWAS_ex5/data2$n_commonSNP_ex5
data2$n_eQTL_exon_rate=data2$n_eQTL_exon/data2$n_commonSNP_exon
data2$n_GWAS_exon_rate=data2$n_GWAS_exon/data2$n_commonSNP_exon
data2$n_eQTL_ex5_rate[which(data2$n_eQTL_ex5_rate == "NaN")]=NA
data2$n_GWAS_ex5_rate[which(data2$n_GWAS_ex5_rate == "NaN")]=NA
data2$n_eQTL_exon_rate[which(data2$n_eQTL_exon_rate == "NaN")]=NA
data2$n_GWAS_exon_rate[which(data2$n_GWAS_exon_rate == "NaN")]=NA

write.table(data2,gzfile(paste0(path_fig6_data,"end5_cluster_SNP_number_commononly.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#number of SNP by ex5_cluster located in [primary_folder]/fig6/data

#===============================================================================




