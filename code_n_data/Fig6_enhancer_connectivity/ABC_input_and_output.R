library(dplyr) 
library(magrittr)
library(edgeR)
library(knitr)
library(ggplot2)
library(stringr)
library(ggthemes)
library(ggrepel)
library(reshape2)
library(ggsignif)
library(gridExtra)
library(grid)
library(data.table)
library(tidyr)
library(ggrastr)

#===============================================================================
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig6_data=paste0(primary_folder,"fig6/data/")
ABC_folder=paste0(primary_folder,"code_n_data/Fig6_enhancer_connectivity/ABC/")
SCAFE_folder=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")

#===============================================================================
# merge tCRE region from opposite strands
setwd(SCAFE_folder)
system("bedtools merge -i ontCAGE.Neuron_THP1.CRE.coord.bed.gz | gzip > merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz")
system("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz | gzip > label_merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz")

#================
merge=read.delim("merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz", header=F, stringsAsFactors = F)
merge$V4=paste0(merge$V1,"_",merge$V2,"_",merge$V3)
merge$V5="."
merge$V6="."
write.table(merge,gzfile("merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz"),col.names=F, row.names=F, sep="\t", quote=F)

#=================
labeling=read.delim("label_merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz", header=F, stringsAsFactors = F)

CRE=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE=left_join(CRE,labeling[,c(4,16)], by=c("CREID"="V4"), copy=F)
colnames(CRE)[86]="merge_CREID"
write.table(CRE,gzfile(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#=================
# incorporate merged tCRE region as enhancer/promoter region

peak_a=read.delim("merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz", header=F)
summary(peak_a$V3-peak_a$V2)

peak_a <- peak_a %>%
  mutate(chr_num = as.numeric(str_extract(V1, "\\d+")),
         chr_num = ifelse(is.na(chr_num), 99, chr_num)) %>% 
  arrange(chr_num, V1, V2) %>%
  select(-chr_num)

write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/iPSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/NSC/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
write.table(peak_a[,c(1:3)],"/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/Neuron/Peaks/macs2_peaks.narrowPeak.sorted.candidateRegions.bed", row.names=F, col.names=F, sep="\t", quote=F)
# put them into ABC pipeline

#===============================================================================
# prepare input: scATAC to tagalign
# sc BAM file from DRA019608 (DRR618513- DRR618515)
# intermediate files too big, not provided.
setwd(paste0(ABC_folder,"input"))
system("sbatch bam_to_tagalign.sh")

#=====
#keep only main chr and resort according to ABC_model requirement
setwd(paste0(ABC_folder,"input/tagalign"))
files=list.files()
files=files[-grep("tbi", files)]

for (i in 1:3){
data=fread(files[i], header=F)
data=data[which(nchar(data$V1)<6),]

data <- data %>%
  mutate(chr_num = as.numeric(str_extract(V1, "\\d+")),
         chr_num = ifelse(is.na(chr_num), 99, chr_num)) %>% 
    arrange(chr_num, V1, V2) %>%
  select(-chr_num)

write.table(data,files[i],row.names=F, col.names=F, sep="\t", quote=F)
}

#===============================================================================
# prepare Hi-C data to bedpe
# HiC data from DRA019572 (DRR614942- DRR614953)
# intermediate files too big, not provided.
setwd(paste0(ABC_folder,"input"))
system("sbatch HiCPE.sbatch")

# split Hi-C files into chromasome base
setwd(paste0(ABC_folder,"input"))
system("sh split_hic_bedpe.sh")

#===============================================================================
# H3K27ac BAM files used as they are
# CUT&Tag data from DRA019568 (DRR614873- DRR614905)

#===============================================================================
# run ABC model separtely on iPSC, NSC and Neuron
# run according to the "config.yaml" and "config_biosamples_NeuronSeries.tsv" sopied from the config folder
setwd(ABC_folder)
system("sh ABC.sh")

#===============================================================================
# select representative merged tCRE region without redundant features
setwd(ABC_folder)
path_fig4_data="/osc-fs_home/yip/CFC_seq_paper_fig_data/fig4/data/"
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

CRE=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1=left_join(data1, CRE[,c(1,86)], by="CREID", copy=F)
data1$CpGTATA[which(data1$CpGTATA=="Null")]=NA
data2=data1%>%group_by(merge_CREID)%>%dplyr::summarise(CREID=paste(CREID,collapse=";"),n5_string=paste(n5_string,collapse=";"),ex5cluster_class=paste(unique(ex5cluster_class),collapse=";"),CpGTATA=paste(unique(na.omit(CpGTATA)),collapse=";"))
data3=data2[-grep(";",data2$ex5cluster_class),]
data3=data3[-grep(";",data3$CpGTATA),]
data3=data3[which(data3$ex5cluster_class == "e_ncRNA"),]
data3$CpGTATA[which(data3$CpGTATA=="")]="Null"
data3%>%group_by(CpGTATA)%>%dplyr::summarise(count=n())

#===============================================================================
# gather ABC results
setwd(ABC_folder)
peak_a=read.delim(paste0(SCAFE_folder,"merge_ontCAGE.Neuron_THP1.CRE.coord.bed.gz"), header=F)

path1="/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/"
files=list.files(path=path1, pattern="EnhancerPredictionsFull_threshold0.02_self_promoter.tsv", recursive =T)
cells=sapply(strsplit(files,"\\/"),"[",1)

data0=data.frame()
for (i in 1:3){
data=fread(paste0(path1,files[i]), header=T)
data=left_join(data, peak_a[,c(1:4)], by=c("chr"="V1", "start"="V2", "end"="V3"),copy=F)
data=left_join(data, data3, by=c("V4"="merge_CREID"),copy=F)
data=data[which(!is.na(data$ex5cluster_class)),]
data=data[which(data$class != "promoter"),] #remove self
data$cell=cells[i]
data0=rbind(data0, data)}
write.table(data0,gzfile("EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

data00=data0%>%group_by(cell,CpGTATA, V4)%>%dplyr::summarise(count=n(), ABC.Score=mean(ABC.Score))
write.table(data00,gzfile(paste0(path_fig6_data,"EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone_summary.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
# -> for fig ex8 g-h

#=====================
#get the overall matrix
setwd(ABC_folder)
path1="/home/yip/tool2026/ABC-Enhancer-Gene-Prediction-main/results/"
files=list.files(path=path1, pattern="EnhancerPredictionsAllPutative.tsv", recursive =T)
cells=sapply(strsplit(files,"\\/"),"[",1)
data0=data.frame()
for (i in 1:3){
data=fread(paste0(path1,files[i]), header=T, select=c(1,2,3,11,28))
data$celltype=cells[i]
data0=rbind(data0, data)}

data0=left_join(data0, peak_a[,c(1:4)], by=c("chr"="V1", "start"="V2", "end"="V3"),copy=F)
data0$pairID=paste0(data0$V4,"::",data0$TargetGene)
data0a=spread(data0[,c(8,6,5)], key=2, value=3)
data0a$merge_CREID=sapply(strsplit(data0a$pairID,"::"),"[",1)
data0a$TargetGene=sapply(strsplit(data0a$pairID,"::"),"[",2)
write.table(data0a,gzfile("ABCscore_matrix_3cells.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
# supplementary table 
data0=read.delim("EnhancerPredictionsFull_threshold0.02_3cell_eRNAalone.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
data0a=spread(data0[,c("V4","TargetGene","ex5cluster_class","CpGTATA","n5_string","cell","ABC.Score")], key=6, value=7)
colnames(data0a)[1]="merged_tCRE_ID"
data0a=data0a[which(data0a$CpGTATA != "Others"),]

#====
# add eRNA geneID to it
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5a=unique(table5[which(table5$ex5cluster_class == "e_ncRNA"),c("n5_string","T4_gene_ID","T4_gene_name")])

data1=read.delim(paste0(primary_folder,"fig4/data/features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data0b=separate_rows(data0a, n5_string)
data0b=left_join(data0b, table5a, by=c("n5_string"), copy=F)
data0b=left_join(data0b,data1[,c("n5_string","SE_source","repeat_ex5cluster","NFY_bind_iPSC","TEAD_bind_iPSC")], by="n5_string", copy=F)

data0c=data0b%>%group_by(merged_tCRE_ID)%>%dplyr::summarise(SE_source=paste(unique(SE_source),collapse=";"), repeat_ex5cluster=paste(unique(repeat_ex5cluster),collapse=";"),
                                                            NFY_bind_iPSC=paste(unique(NFY_bind_iPSC),collapse=";"),TEAD_bind_iPSC=paste(unique(TEAD_bind_iPSC),collapse=";"),
                                                            eRNA_ID=paste(unique(T4_gene_ID),collapse=";"), eRNA_name=paste(unique(T4_gene_name),collapse=";"))

data0c$NFY_bind_iPSC[grep(";",data0c$NFY_bind_iPSC)]="Yes"
data0c$TEAD_bind_iPSC[grep(";",data0c$TEAD_bind_iPSC)]="Yes"

data0a=left_join(data0a,data0c, by="merged_tCRE_ID", copy=F)
write.table(data0a,paste0(primary_folder,"supplementary_table/TableS16_ABCmodel.tsv"), col.names=T, row.names=F, sep="\t", quote=F)



