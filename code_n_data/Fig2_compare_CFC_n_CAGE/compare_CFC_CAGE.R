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
library(ggrastr)

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
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/")
compare_path=paste0(primary_folder, "code_n_data/Fig2_compare_CFC_n_CAGE/compare/")
repeat_path=paste0(primary_folder, "code_n_data/Fig2_compare_CFC_n_CAGE/repeat/")
long_read_bambu=paste0(primary_folder,"code_n_data/transcript_model_analyses_Fig3/bambu_long_t5_partialYes.ENST/")

#===============================================================================
# Mapping of ssCAGE data using both PE and SE
# FASTQ and BAM not included, please find them from DRA019567 (DRR614867- DRR614872)
setwd(paste0(primary_folder, "code_n_data/Fig2_compare_CFC_n_CAGE/"))
system("sh ssCAGE.STAR.map.sh")
# -> map to hg38 and extract bed file with MAPQ from read1

# sorted bam files with mapq 0 were applied to SCAFE for independent run, joint SCAFE combine CFC-seq data was also performed
# code and results of independent ssCAGE SCAFE refer to [primary_folder]/code_n_data/SCAFE/CAGE_Neuronalone_mapq0/
# code and result of joint SCAFE run refer to [primary_folder]/code_n_data/SCAFE/joint_CAGE_CFC_Neuronalone/

#===============================================================================
#compare TSS location -> use ssCAGE SCAFE independent run here
path8=paste0(SCAFE_path,"CFC_Neuron_THP1/ontCAGE/bam_to_ctss/")

need=c("iPS","NSC","Neuron")
files=list.files(pattern="collapse.ctss.bed.gz", path=path8, recursive = T)
files=files[-grep("tbi",files)]
for (i in 1 : length(need)){
  files1=files[grep(need[i],files)]
  files1a=files1[grep("unencoded_G",files1)]
  files1b=files1[-grep("unencoded_G",files1)]
  dataa=read.delim(paste0(path8,files1a[1]), header=F, stringsAsFactors = F)
  datab=read.delim(paste0(path8,files1b[1]), header=F, stringsAsFactors = F)
  for (j in 2: length(files1a)){
    data1a=read.delim(paste0(path8,files1a[j]), header=F, stringsAsFactors = F)
    dataa=rbind(dataa,data1a)
    dataa=dataa%>%group_by(V1,V2,V3,V6)%>%summarise(V4=sum(V4),V5=sum(V5))
    dataa=dataa[order(dataa$V1,dataa$V2),c(1,2,3,5,6,4)]
    data1b=read.delim(paste0(path8,files1b[j]), header=F, stringsAsFactors = F)
    datab=rbind(datab,data1b)
    datab=datab%>%group_by(V1,V2,V3,V6)%>%summarise(V4=sum(V4),V5=sum(V5))
    datab=datab[order(datab$V1,datab$V2),c(1,2,3,5,6,4)]}
  write.table(dataa,gzfile(paste0(compare_path,need[i],".CFC.unG.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
  write.table(datab,gzfile(paste0(compare_path,need[i],".CFC.CTSS.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)}

path9=paste0(SCAFE_path,"CAGE_Neuronalone_mapq0/bam_to_ctss/")
files=list.files(pattern="collapse.ctss.bed.gz", path=path9, recursive = T)
files=files[-grep("tbi",files)]
for (i in 1 : length(need)){
  files1=files[grep(need[i],files)]
  files1a=files1[grep("unencoded_G",files1)]
  files1b=files1[-grep("unencoded_G",files1)]
  dataa=read.delim(paste0(path9,files1a[1]), header=F, stringsAsFactors = F)
  datab=read.delim(paste0(path9,files1b[1]), header=F, stringsAsFactors = F)
  data1a=read.delim(paste0(path9,files1a[2]), header=F, stringsAsFactors = F)
  dataa=rbind(dataa,data1a)
  dataa=dataa%>%group_by(V1,V2,V3,V6)%>%summarise(V4=sum(V4),V5=sum(V5))
  dataa=dataa[order(dataa$V1,dataa$V2),c(1,2,3,5,6,4)]
  data1b=read.delim(paste0(path9,files1b[2]), header=F, stringsAsFactors = F)
  datab=rbind(datab,data1b)
  datab=datab%>%group_by(V1,V2,V3,V6)%>%summarise(V4=sum(V4),V5=sum(V5))
  datab=datab[order(datab$V1,datab$V2),c(1,2,3,5,6,4)]
  write.table(dataa,gzfile(paste0(compare_path,need[i],".CAGE.unG.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
  write.table(datab,gzfile(paste0(compare_path,need[i],".CAGE.CTSS.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)}


files=list.files(path=compare_path, pattern="bed.gz")
need=c("iPS","NSC","Neuron")
class=c("CTSS","unG")
summary=data.frame(matrix(nrow=0, ncol=5))
for (i in 1: length(need)){
  files1=files[grep(need[i], files)]
  for (j in 1:2){
    files1a=files1[grep(class[j], files1)]
    ID=sapply(strsplit(files1a,"\\."),"[",2)
    data=read.delim(paste0(path,files1a[1]), header=F, stringsAsFactors = F)
    data1=read.delim(paste0(path,files1a[2]), header=F, stringsAsFactors = F)
    data$label=paste0(data$V1,"_",data$V2,"_",data$V3,"_",data$V6)
    data1$label=paste0(data1$V1,"_",data1$V2,"_",data1$V3,"_",data1$V6)
    colnames(data)[c(5)]=ID[1]
    colnames(data1)[c(5)]=ID[2]
    data=full_join(data[,c(7,5)],data1[,c(7,5)], by="label",copy=F)
    data[is.na(data)]=0
    data$group=need[i]
    data$class=class[j]
    colnames(summary)=colnames(data)
    summary=rbind(summary, data)}}
write.table(summary,gzfile(paste0(compare_path,"CTSS_1nt_CAGE_CFC_compare.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

summary$group2="both"
summary$group2[which(summary$CAGE == 0)]="CFC-\nspecific"
summary$group2[which(summary$CFC == 0)]="CAGE-\nspecific"
summary$count=summary$CAGE+summary$CFC
summary0=summary%>%group_by(group, class, group2)%>%dplyr::summarise(count=sum(count))%>%dplyr::mutate(percent=count/sum(count))
summary1=summary%>%group_by(class, group2)%>%dplyr::summarise(count=sum(count))%>%dplyr::mutate(percent=count/sum(count), group="Neuron-series")
summary0$group[which(summary0$group == "iPS")]="iPSC"
summary0=rbind(summary0, summary1[,c(5,1:4)])

write.table(summary0,gzfile(paste0(path_fig2_data,"CTSS_1nt_CAGE_CFC_compare.summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#===========================================
#tCRE base CFC CAGE compare
##plot venn diagram compare ss and ont-cage from SCAFE using both CFC and CAGE

count_m=read.delim(paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/count/output/count_matrix_cluster/ont_ss.Neuronalone.count.txt"), header=T, check.names = F, stringsAsFactors = F)
anno=read.delim(paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/aggregate/run_full/out/annotate/ont_ss.Neuronalone/log/ont_ss.Neuronalone.cluster.info.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
colnames(count_m)[1]="clusterID"
count_m=left_join(count_m, anno[,c(1,16)], by="clusterID",copy=F)
count_m=count_m[,-1]
count_m=count_m %>% group_by(CREID) %>% summarise_each(list(sum))
count_m1=data.frame(colSums(count_m[,c(2:13)]))
count_m1$lib=sapply(strsplit(rownames(count_m1),"\\.rep"),"[",1)
colnames(count_m1)[1]="qualified_read_cluster"
count_m1=count_m1%>%group_by(lib)%>%dplyr::summarise(qualified_read_cluster=sum(qualified_read_cluster))
count_m$CFC1=0
count_m$CFC1[which(rowSums(count_m[,c(2:7)])>=1)]=1
count_m$CAGE1=0
count_m$CAGE1[which(rowSums(count_m[,c(8:13)])>=1)]=1
count_m$CFC2=0
count_m$CFC2[which(rowSums(count_m[,c(2:7)])>=2)]=1
count_m$CAGE2=0
count_m$CAGE2[which(rowSums(count_m[,c(8:13)])>=2)]=1
count_m$CFC3=0
count_m$CFC3[which(rowSums(count_m[,c(2:7)])>=3)]=1
count_m$CAGE3=0
count_m$CAGE3[which(rowSums(count_m[,c(8:13)])>=3)]=1
write.table(count_m,gzfile(paste0(path_fig2_data,"tCRE.identification.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#===============================================================================
#tCRE quantification between CAGE and CFC -> use joint SCAFE run here
count_m=read.delim(paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/count/output/count_matrix_cluster/ont_ss.Neuronalone.count.txt"), header=T, check.names = F, stringsAsFactors = F)
anno=read.delim(paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/aggregate/run_full/out/annotate/ont_ss.Neuronalone/log/ont_ss.Neuronalone.cluster.info.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
colnames(count_m)[1]="clusterID"
count_m=left_join(count_m, anno[,c(1,16)], by="clusterID",copy=F)
count_m=count_m[,-1]
count_m=count_m %>% group_by(CREID) %>% summarise(across(everything(), sum))
count_m$CFC_iPSC=count_m$ontCAGE.iPSC.rep1+count_m$ontCAGE.iPSC.rep2
count_m$CAGE_iPSC=count_m$ssCAGE.iPS.rep1+count_m$ssCAGE.iPS.rep2
count_m$CFC_NSC=count_m$ontCAGE.NSC.rep1+count_m$ontCAGE.NSC.rep2
count_m$CAGE_NSC=count_m$ssCAGE.NSC.rep1+count_m$ssCAGE.NSC.rep2
count_m$CFC_Neuron=count_m$ontCAGE.Neuron.rep1+count_m$ontCAGE.Neuron.rep2
count_m$CAGE_Neuron=count_m$ssCAGE.Neuron.rep1+count_m$ssCAGE.Neuron.rep2
count_m=data.frame(count_m)
rownames(count_m)=count_m$CREID
count_m=count_m[,c(14:19)]
summary1=data.frame(matrix(nrow=0, ncol=4))
need=c("iPSC","NSC","Neuron")
for (i in 1: length(need)){
  count_m1=count_m[,grep(need[i],colnames(count_m))]
  count_m1=count_m1[which(rowSums(count_m1)>0),]
  colnames(count_m1)=sapply(strsplit(colnames(count_m1),"_"),"[",1)
  d <- DGEList(counts=count_m1)
  RLE <- calcNormFactors(d, method="RLE")
  RLE.CRE=cpm(RLE, normalized.lib.sizes=TRUE)
  RLE.CRE=as.data.frame(RLE.CRE)
  RLE.CRE=log10(RLE.CRE+0.01)
  label1=cor.test(RLE.CRE[,1],RLE.CRE[,2])$estimate
  RLE.CRE$CREID=rownames(RLE.CRE)
  RLE.CRE$cell=need[i]
  colnames(summary1)=colnames(RLE.CRE)
  summary1=rbind(summary1, RLE.CRE)}
write.table(summary1, gzfile(paste0(path_fig2_data,"tCRE_quantification_CFC_CAGE_RLE_CPM.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data
#===============================================================================


# repetitive elements
#===============================================================================
# bash
# take MAPQ and NM from CFC-seq bam file
setwd(paste0(primary_folder, "code_n_data/Fig2_compare_CFC_n_CAGE/"))
system("sh ont.bamto_aln_info.sh")

# bam to bed to extract readname and CTSS
system("sh ont.bamtobed.sh")

#===============================================================================
# combine all the CTSS files into ont, PE and SE
# bash
# bedtools intersect read CTSS with joint SCAFE TSS cluster
# cluster base, use the summit to the 150nt downstream, for comparing CAGE and CFC-seq

SCAFE_TSScluster=paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/bed/ont_ss.Neuronalone.cluster.coord.bed.gz")
setwd(repeat_path)
system(paste0("bedtools intersect -wa -wb -s -a ss_all.PE.CTSS.tsv.gz -b ",SCAFE_TSScluster," | gzip > ss_all.PE.CTSS.scafe_cluster.bed.gz"))
system(paste0("bedtools intersect -wa -wb -s -a ss_all.SE.CTSS.tsv.gz -b ",SCAFE_TSScluster," | gzip > ss_all.SE.CTSS.scafe_cluster.bed.gz"))
system(paste0("bedtools intersect -wa -wb -s -a ont_all.CTSS.tsv.gz -b ",SCAFE_TSScluster," | gzip > ont_all.pri.CTSS.scafe_cluster.bed.gz"))

#bedtools intersect with repeat
#bash
setwd(repeat_path)
system("bedtools intersect -wa -wb -s -a ont_ss.Neuronalone.cluster.coord.submit_down150.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/UCSC.hg38.repeat.bed.gz | gzip > ont_ss.Neuronalone.cluster.coord.submit_down150.rmsk.bed.gz")

#=================================================
#R

ontcluster=read.delim(paste0(repeat_path,"ont_all.pri.CTSS.scafe_cluster.bed.gz"), header=F, stringsAsFactors = F)
ontcluster1=ontcluster%>%group_by(V13)%>%dplyr::summarise(min_NM=min(V7), max_MAPQ=max(V8), rRNA= max(V9), total_count=sum(V5) )
colnames(ontcluster1)[c(2:5)]=paste0("ontCAGE_",colnames(ontcluster1)[c(2:5)])
ssPEcluster=read.delim(paste0(repeat_path,"ss_all.PE.CTSS.scafe_cluster.bed.gz"), header=F, stringsAsFactors = F)
ssPEcluster1=ssPEcluster%>%group_by(V11)%>%dplyr::summarise(max_MAPQ=max(V7), total_count=sum(V5) )
colnames(ssPEcluster1)[c(2,3)]=paste0("PEssCAGE_",colnames(ssPEcluster1)[c(2,3)])
ssSEcluster=read.delim(paste0(repeat_path,"ss_all.SE.CTSS.scafe_cluster.bed.gz" ), header=F, stringsAsFactors = F)
ssSEcluster1=ssSEcluster%>%group_by(V11)%>%dplyr::summarise(max_MAPQ=max(V7), total_count=sum(V5) )
colnames(ssSEcluster1)[c(2,3)]=paste0("SEssCAGE_",colnames(ssSEcluster1)[c(2,3)])

allcluster=read.delim(paste0(SCAFE_path,"joint_CAGE_CFC_Neuronalone/aggregate/run_full/out/annotate/ont_ss.Neuronalone/bed/ont_ss.Neuronalone.cluster.coord.bed.gz"), header=F, stringsAsFactors = F)
allcluster=left_join(allcluster,ontcluster1,by=c("V4"="V13"),copy=F)
allcluster=left_join(allcluster,ssPEcluster1,by=c("V4"="V11"),copy=F)
allcluster=left_join(allcluster,ssSEcluster1,by=c("V4"="V11"),copy=F)

write.table(allcluster,gzfile(paste0(repeat_path,"ont_ss.Neuronalone.cluster.final.result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
allcluster=read.delim(paste0(repeat_path,"ont_ss.Neuronalone.cluster.final.result.tsv.gz"), header=T, stringsAsFactors = T, check.names = F)

repeat_anno=read.delim("/analysisdata/fantom6/Interactome/resources/UCSC/UCSC.hg38.repeat.tsv", header=T, stringsAsFactors = F, check.names = F)
repeat_anno$ID=paste0(repeat_anno$genoName,"_",repeat_anno$genoStart,"_",repeat_anno$genoEnd,"_",repeat_anno$strand)

clusterrep=read.delim(paste0(repeat_path,"ont_ss.Neuronalone.cluster.coord.submit_down150.rmsk.bed.gz"), header=F, stringsAsFactors = F)
clusterrep$region=clusterrep$V15-clusterrep$V14
clusterrep=left_join(clusterrep,repeat_anno[,c(18,11:13,3)], by=c("V16"="ID"),copy=F)
clusterrep=clusterrep[,-c(19:21)]
write.table(clusterrep,gzfile(paste0(repeat_path,"ont_ss.Neuronalone.cluster.hit.rep.info.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
                              
aa=clusterrep%>%group_by(repClass)%>%dplyr::summarise(count=n(), average_region=mean(region))
aa1=aa[which(aa$count>240),] #select 7 major classes

clusterrep1=clusterrep%>%group_by(V4)%>%dplyr::summarise(repName=paste(repName,collapse=";"),repClass=paste(repClass,collapse=";"), repFamily=paste(repFamily,collapse=";"), repSize=paste(region,collapse=";"), div=paste(milliDiv,collapse=";"))
allcluster=left_join(allcluster,clusterrep1,by="V4",copy=F)
length(which(allcluster$ontCAGE_max_MAPQ <=3 & !is.na(allcluster$repFamily)))
length(which(allcluster$PEssCAGE_max_MAPQ <=3 & !is.na(allcluster$repFamily)))
length(which(allcluster$SEssCAGE_max_MAPQ <=3 & !is.na(allcluster$repFamily)))
write.table(allcluster,gzfile(paste0(repeat_path,"ont_ss.Neuronalone.cluster.final.result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#=========================
##include only TSS cluster that were detected by all 3 library types
allcluster1=allcluster[which(!is.na(allcluster$ontCAGE_max_MAPQ)),]
allcluster1=allcluster1[which(!is.na(allcluster1$PEssCAGE_max_MAPQ)),]
allcluster1=allcluster1[which(!is.na(allcluster1$SEssCAGE_max_MAPQ)),]
write.table(allcluster1,gzfile(paste0(repeat_path,"ont_ss.Neuronalone.cluster.final.rm_na.result.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

allcluster0=allcluster1[which( !is.na(allcluster1$repFamily)),]

asummary1=data.frame(matrix(nrow=7*4, ncol=7))
colnames(asummary1)=c("repClass","ont","PE","SE","total","average_region","max_MQ")
score=c(0,5,10,20)
for (k in c(1:4)){
  for (i in 1:nrow(aa1)){
    allcluster0a=allcluster0[grep(aa1$repClass[i],allcluster0$repClass),]
    allcluster0a$PEssCAGE_max_MAPQ[which(allcluster0a$PEssCAGE_max_MAPQ==225)]=20
    allcluster0a$SEssCAGE_max_MAPQ[which(allcluster0a$SEssCAGE_max_MAPQ==225)]=20
    allcluster0a$PEssCAGE_max_MAPQ[which(allcluster0a$PEssCAGE_max_MAPQ==2)]=5
    allcluster0a$SEssCAGE_max_MAPQ[which(allcluster0a$SEssCAGE_max_MAPQ==2)]=5
    allcluster0a$PEssCAGE_max_MAPQ[which(allcluster0a$PEssCAGE_max_MAPQ==3)]=10
    allcluster0a$SEssCAGE_max_MAPQ[which(allcluster0a$SEssCAGE_max_MAPQ==3)]=10
    asummary1$repClass[7*(k-1)+i]=aa1$repClass[i]
    asummary1$ont[7*(k-1)+i]=length(which(allcluster0a$ontCAGE_max_MAPQ >= score[k]))
    asummary1$PE[7*(k-1)+i]=length(which(allcluster0a$PEssCAGE_max_MAPQ >= score[k]))
    asummary1$SE[7*(k-1)+i]=length(which(allcluster0a$SEssCAGE_max_MAPQ >= score[k]))
    asummary1$total[7*(k-1)+i]=nrow(allcluster0a)
    asummary1$average_region[7*(k-1)+i]=aa1$average_region[i]
    asummary1$max_MQ[7*(k-1)+i]=k
  }}
asummary1$percent_ont=asummary1$ont/asummary1$total
asummary1$percent_PE=asummary1$PE/asummary1$total
asummary1$percent_SE=asummary1$SE/asummary1$total

write.table(asummary1,gzfile(paste0(path_fig2_data,"RE_MAPQ_compare.NA_removed.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#===============================================================================
#see the evolutionary age

clusterrep=read.delim(paste0(repeat_path,"ont_ss.Neuronalone.cluster.hit.rep.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
clusterrep=clusterrep[which(!is.na(clusterrep$ontCAGE_max_MAPQ)),]
clusterrep=clusterrep[which(!is.na(clusterrep$PEssCAGE_max_MAPQ)),]
clusterrep=clusterrep[which(!is.na(clusterrep$SEssCAGE_max_MAPQ)),]
clusterrep=clusterrep[which(clusterrep$repClass %in% c("DNA","LTR","LINE","SINE")),]
clusterrep$group5="others"
clusterrep$group5[which(clusterrep$ontCAGE_max_MAPQ>=10 & clusterrep$SEssCAGE_max_MAPQ <225)]="CFC"
clusterrep$repClass=factor(clusterrep$repClass, levels=c("DNA","LTR","LINE","SINE"))
repcluster4=clusterrep%>%group_by(repClass)%>%dplyr::summarise(p=wilcox.test(milliDiv ~group5, alternative = "two.sided")$p.value)
repcluster4$label="***"
repcluster4$label[which(repcluster4$p>=0.001)]="**"
repcluster4$label[which(repcluster4$p>=0.01)]="*"
repcluster4$label[which(repcluster4$p>=0.05)]="n.s."
clusterrep$group5=factor(clusterrep$group5, levels=c("CFC","others"))


#==========================================================================================
#repetitive elements in tCRE (summit -400 + 100 region)
#for fig.ex3k&l

#bedtools intersect
#bash
setwd(repeat_path)
system("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/UCSC.hg38.repeat.bed.gz | gzip > Neuron_THP1_S3.CRE.coord.rmsk.bed.gz")
system("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.stranded_summit_extend_100_400.bed.gz -b /analysisdata/fantom6/Interactome/resources/UCSC/UCSC.hg38.repeat.bed.gz | gzip > Neuron_THP1_S3.CRE.coord.submit_extended.rmsk.bed.gz")

#=================================================
#R

CREanno=read.delim(paste0(SCAFE_path,"CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno1=unique(CREanno[which(CREanno$representative == "Yes"),c(1,32,33,34,35,40)])

CRE_re=read.delim(paste0(repeat_path, "Neuron_THP1_S3.CRE.coord.rmsk.bed.gz"), header=F, stringsAsFactors = F)

CRE_re$V23=1
CRE_re1=CRE_re%>%group_by(V4,V20)%>%dplyr::summarise(region=max(V23))
CRE_re1=spread(CRE_re1,key=2, value=3)
CRE_re1=CRE_re1[,c(1,2,4,5,6,13,15,16)]

t51a=left_join(CREanno1,CRE_re1, by=c("CREID"="V4"),copy=F)
t51a[is.na(t51a)]=0

RE=colnames(t51a)[c(7:13)]
pro_type=unique(t51a$promoter_type)
direction=c("1D","2D")
summaryRE=data.frame(matrix(nrow=56,ncol=6))
colnames(summaryRE)=c("promoter_type","RE","yesyes","yesno","noyes","nono")

for (i in 1: length(pro_type)){
  for (k in 1: length(RE)){
    summaryRE$promoter_type[(i-1)*7+k]=pro_type[i]
    summaryRE$RE[(i-1)*7+k]=RE[k]
    summaryRE$yesyes[(i-1)*7+k]=length(which(t51a$promoter_type == pro_type[i] & t51a[,k+6] == 1))
    summaryRE$yesno[(i-1)*7+k]=length(which(t51a$promoter_type == pro_type[i] & t51a[,k+6] == 0))
    summaryRE$noyes[(i-1)*7+k]=length(which(t51a$promoter_type != pro_type[i] & t51a[,k+6] == 1))
    summaryRE$nono[(i-1)*7+k]=length(which(t51a$promoter_type != pro_type[i] & t51a[,k+6] == 0))}}
t51b=t51a[which(t51a$promoter_type=="enhancer-like" & t51a$orientation %in% c("1D","2D")),]
for (i in 1: length(direction)){
  for (k in 1: length(RE)){
    summaryRE$promoter_type[(i-1)*7+k+28]=paste0("enhancer_",direction[i])
    summaryRE$RE[(i-1)*7+k+28]=RE[k]
    summaryRE$yesyes[(i-1)*7+k+28]=length(which(t51b$orientation == direction[i] & t51b[,k+6] == 1))
    summaryRE$yesno[(i-1)*7+k+28]=length(which(t51b$orientation == direction[i] & t51b[,k+6] == 0))
    summaryRE$noyes[(i-1)*7+k+28]=length(which(t51b$orientation != direction[i] & t51b[,k+6] == 1))
    summaryRE$nono[(i-1)*7+k+28]=length(which(t51b$orientation != direction[i] & t51b[,k+6] == 0))}}
t51c=t51a[which(t51a$promoter_type=="promoter-like"& t51a$orientation %in% c("1D","2D")),]
for (i in 1: length(direction)){
  for (k in 1: length(RE)){
    summaryRE$promoter_type[(i-1)*7+k+42]=paste0("promoter_",direction[i])
    summaryRE$RE[(i-1)*7+k+42]=RE[k]
    summaryRE$yesyes[(i-1)*7+k+42]=length(which(t51b$orientation == direction[i] & t51b[,k+6] == 1))
    summaryRE$yesno[(i-1)*7+k+42]=length(which(t51b$orientation == direction[i] & t51b[,k+6] == 0))
    summaryRE$noyes[(i-1)*7+k+42]=length(which(t51b$orientation != direction[i] & t51b[,k+6] == 1))
    summaryRE$nono[(i-1)*7+k+42]=length(which(t51b$orientation != direction[i] & t51b[,k+6] == 0))}}

for(i in 1:56){
  GSEATasting <- matrix(c(summaryRE$yesyes[i], summaryRE$yesno[i], summaryRE$noyes[i], summaryRE$nono[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  summaryRE$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  summaryRE$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
write.table(summaryRE, gzfile(paste0(path_fig2_data,"promoter_RE_FEresult.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#=====================================
#number of repetitive elements in tCRE
CREanno=read.delim(paste0(SCAFE_path,"CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno1=unique(CREanno[which(CREanno$representative == "Yes"),c(1,32,34,35)])
CRE_re=read.delim(paste0(repeat_path, "Neuron_THP1_S3.CRE.coord.rmsk.bed.gz"), header=F, stringsAsFactors = F)
RE=c("DNA", "LTR", "LINE","SINE","Low_complexity","Satellite","Simple_repeat")
CRE_re1=unique(CRE_re[which(CRE_re$V20 %in% RE),c(4,20)])
CREanno1=left_join(CREanno1,CRE_re1, by=c("CREID"="V4"),copy=F)
CREanno2=CREanno1%>%group_by(CREID,promoter_type,orientation)%>%dplyr::summarise(RE=paste(unique(V20),collapse=";"))
CREanno2$RE[which(CREanno2$RE == "NA")]="Null"
CREanno2$RE[grep(";", CREanno2$RE)]="Multiple"
write.table(CREanno2,gzfile(paste0(path_fig2_data,"tCRE_RE.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#=====================================
#HipSTR at tCRE summit

#=======================================
#if tCRE summit initiate from simple repeat
#hipstr reference was downloaded from https://github.com/HipSTR-Tool/HipSTR-references

#=================================
#bash
#intersect this hipstr with table5 transcript model 5' end and gencode 5' end
setwd(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",long_read_bambu,"hg38.hipstr_reference.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.stranded_summit_hipstr.bed.gz"))

# -> hipstr reference file was downloaded from https://github.com/HipSTR-Tool/HipSTR-references/blob/master/human/hg38.hipstr_reference.bed.gz,
# and removed
#================================
setwd(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/"))
summit_hip=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_hipstr.bed.gz", header=F, stringsAsFactors = F, check.names = F)
summit_hip=left_join(summit_hip, CREanno[,c(1,32)], by=c("V4"="CREID"),copy=F)
aa=summit_hip%>%group_by(V13,promoter_type)%>%dplyr::summarise(count=n())
aa1=spread(aa, key=2, value=3)
aa1$all=rowSums(aa1[,c(2:5)],na.rm=T)
length(intersect(table5_hip$V4,gencode_hip$V4))
table0_hip=left_join(table0_hip, table5[,c(1,62,64,65,86,87,90,110)], by=c("V4"="model_ID"),copy=F)
table5_hip=table0_hip[which(!is.na(table0_hip$T4_gene_ID)),]
k0=table0_hip%>%group_by(V13,V6)%>%dplyr::summarise(count=n())
k5=table5_hip%>%group_by(V13,V6)%>%dplyr::summarise(count=n())

CREanno=read.delim(paste0(SCAFE_path,"CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
CREanno1=unique(CREanno[which(CREanno$representative == "Yes"),c(1,32,33,34,35,40)])
summit_hip=read.delim("ontCAGE.Neuron_THP1.CRE.stranded_summit_hipstr.bed.gz", header=F, stringsAsFactors = F, check.names = F)

summit_hip$V14=1
summit_hip1=summit_hip%>%group_by(V13)%>%dplyr::summarise(count=n())
need=summit_hip1$V13[which(summit_hip1$count>30)]
summit_hip2=spread(summit_hip[which(summit_hip$V13 %in% need),c(4,13,14)],key=2, value=3)

t51b=left_join(CREanno1,summit_hip2, by=c("CREID"="V4"),copy=F)
t51b[is.na(t51b)]=0

STR=colnames(t51b)[c(7:14)]
pro_type=unique(t51b$promoter_type)
summarySTR=data.frame(matrix(nrow=29,ncol=6))
colnames(summarySTR)=c("promoter_type","STR","yesyes","yesno","noyes","nono")
for (i in 1: length(pro_type)){
  for (k in 1: length(STR)){
    summarySTR$promoter_type[(i-1)*7+k]=pro_type[i]
    summarySTR$STR[(i-1)*7+k]=STR[k]
    summarySTR$yesyes[(i-1)*7+k]=length(which(t51b$promoter_type == pro_type[i] & t51b[,k+6] == 1))
    summarySTR$yesno[(i-1)*7+k]=length(which(t51b$promoter_type == pro_type[i] & t51b[,k+6] == 0))
    summarySTR$noyes[(i-1)*7+k]=length(which(t51b$promoter_type != pro_type[i] & t51b[,k+6] == 1))
    summarySTR$nono[(i-1)*7+k]=length(which(t51b$promoter_type != pro_type[i] & t51b[,k+6] == 0))}}

for(i in 1:29){
  GSEATasting <- matrix(c(summarySTR$yesyes[i], summarySTR$yesno[i], summarySTR$noyes[i], summarySTR$nono[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  summarySTR$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  summarySTR$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
write.table(summarySTR, gzfile(paste0(path_fig2_data,"promoter_STR_FEresult.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


