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
path_fig1_data=paste0(primary_folder,"fig1/data/")
read_path=paste0(primary_folder,"code_n_data/Fig1_read_analyses/")
CTSS_path=paste0(read_path,"WTC11.CTSS/")
bed_path=paste0(read_path,"bamtobed/")
TES_path=paste0(read_path,"TES/")
n5_path=paste0(primary_folder,"code_n_data/n5_regions/")
GENCODE_path=paste0(primary_folder,"GENCODEv39/")
intersect_path=paste0(read_path,"n3_intersect/")

#===============================================================================
#long-read data of WTC11 was downloaded from ENCODE
setwd(read_path)
system("sh longread.minimap2.sh")
system("sh longread.bamtobed.intersect.sh")
system("sh longread.intersect.mRNAbed6.sh")
system("sh longread.intersect.mRNAbed6.myWTC11.sh")

#-> These steps generated bed12.bed showing all the reads that mapped to mRNA from GENCODE v39 from each library.
#-> These bam files were not maintained

#===============================================================================
#process by R
files=list.files(path=bed_path, pattern=".mRNA.bed", recursive = T)
groups=sapply(strsplit(files,"\\/"),"[", 1)
groups=sapply(strsplit(groups,"_EN"),"[", 1)
names=sapply(strsplit(files,"\\/"),"[", 2)
names=gsub(".mRNA.bed.gz","",names)
options(scipen=999)
for (i in 1:length(files)){
  data=read.delim(paste0(bed_path,files[i]), header=F, stringsAsFactors = F, check.names = F)
  data=data[which(data$V13 > 0),c(1:6)]
  data$V3[which(data$V6 == "+")]=data$V2[which(data$V6 == "+")]+1
  data$V2[which(data$V6 == "-")]=data$V3[which(data$V6 == "-")]-1
  data=data%>%group_by(V1,V2,V3,V6)%>%dplyr::summarise(count=n())
  data$label=paste0("n",1:nrow(data))
  write.table(data[,c(1,2,3,6,5,4)],gzfile(paste0(CTSS_path,names[i],".CTSS.mRNA.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)}

# CTSS data stored in [primary_folder]/code_n_data/Fig1_read_analyses/data/WTC11.CTSS
# after this step the bed12 files were also removed.

#============================
#bash
#intersect with cCRE p/e/ctcf with counting
setwd(CTSS_path)
system(paste0(
  "for file in *.CTSS.mRNA.bed.gz; do ",
  "bedtools intersect -c -a \"$file\" -b ", n5_path, "GRCh38-ELS.all.enhancer.sort.bed.gz | gzip > \"${file%.bed.gz}.Pcount.bed\"; ",
  "bedtools intersect -c -a \"$file\" -b ", n5_path, "GRCh38-PLS.all.promoter.sort.bed.gz | gzip > \"${file%.bed.gz}.Ecount.bed\"; ",
  "bedtools intersect -c -a \"$file\" -b ", n5_path, "GRCh38-CTCF.sort.bed.gz | gzip > \"${file%.bed.gz}.Ccount.bed\"; ",
  "done"
))
#intersect iPSC ATAC with counting
setwd(CTSS_path)
system(paste0("for file in *.CTSS.mRNA.bed.gz; do bedtools intersect -c -a \"$file\" -b ",n5_path, "peaks.iPS.bed | gzip > \"${file%.bed.gz}.ATACcount.bed\" ; done"))

#intersect with FANTOM5 CAGE cluster with counting
setwd(CTSS_path)
system(paste0("for file in *.CTSS.mRNA.bed.gz; do bedtools intersect -c -a \"$file\" -b ", n5_path, "F6_CAT.promoter.bed | gzip > \"${file%.bed.gz}.CAGEcount.bed\" ; done"))

#===========================
#after promoter/enhancer/ctcf counting

CRE.hit.rate=data.frame(matrix(nrow=105,ncol=6))
colnames(CRE.hit.rate)=c("group","name","read.hit","read.total","read.hit.smallRNAd","read.total.smallRNAd")
smallRNA=read.delim(paste0(CTSS_path,"gencode.v39.annotation.smallRNA.bed12.5n.sort.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
small_deplete=paste0(smallRNA$V1,"_",smallRNA$V2,"_",smallRNA$V3,"_",smallRNA$V6)
CRE.hit.rate$group=c(rep("P",21),rep("PE",21),rep("cCRE",21),rep("ATAC",21),rep("F5_CAGE",21))
group1=c("P","PE","cCRE","ATAC","F5_CAGE")

filesP=list.files(path=CTSS_path, pattern=".CTSS.mRNA.Pcount.bed.gz")
filesE=list.files(path=CTSS_path, pattern=".CTSS.mRNA.Ecount.bed.gz")
filesC=list.files(path=CTSS_path, pattern=".CTSS.mRNA.Ccount.bed.gz")
filesA=list.files(path=CTSS_path, pattern=".CTSS.mRNA.ATACcount.bed.gz")
filesCAGE=list.files(path=CTSS_path, pattern=".CTSS.mRNA.CAGEcount.bed.gz")
names=sapply(strsplit(filesP,"\\."),"[", 1)
for (i in 1:length(filesP)){
  dataP=read.delim(paste0(CTSS_path,filesP[i]), header=F, stringsAsFactors = F, check.names = F)
  dataE=read.delim(paste0(CTSS_path,filesE[i]), header=F, stringsAsFactors = F, check.names = F)
  dataC=read.delim(paste0(CTSS_path,filesC[i]), header=F, stringsAsFactors = F, check.names = F)
  dataA=read.delim(paste0(CTSS_path,filesA[i]), header=F, stringsAsFactors = F, check.names = F)
  dataCAGE=read.delim(paste0(CTSS_path,filesCAGE[i]), header=F, stringsAsFactors = F, check.names = F)
  
  data=cbind(dataP,dataE[,7],dataC[,7],dataA[,7],dataCAGE[,7])
  colnames(data)[c(7:11)]=c("P","E","C","ATAC","F5_CAGE")
  data$cCRE=0
  data$cCRE[which(rowSums(data[,c(7:9)]>0)>0)]=1
  data$PE=0
  data$PE[which(rowSums(data[,c(7:8)]>0)>0)]=1
  data$name=paste0(data$V1,"_",data$V2,"_",data$V3,"_",data$V6)
  data$small=0
  data$small[which(data$name %in% small_deplete)]=1
  data=data[which(data$V1 != "chrM"),]
  for (j in 1:length(group1)){
    data1=data
    colnames(data1)[which(colnames(data1)==group1[j])]="key"
    CRE.hit.rate$name[(j*21)-21+i] = names[i] 
    CRE.hit.rate$read.hit[(j*21)-21+i] = sum(data1$V5[which(data1$key > 0)])
    CRE.hit.rate$read.total[(j*21)-21+i] = sum(data1$V5)
    CRE.hit.rate$read.hit.smallRNAd[(j*21)-21+i] = sum(data1$V5[which(data1$key > 0 & data1$small == 0)])
    CRE.hit.rate$read.total.smallRNAd[(j*21)-21+i] = sum(data1$V5[which(data1$small == 0)])}}

CRE.hit.rate$rate=CRE.hit.rate$read.hit/CRE.hit.rate$read.total
CRE.hit.rate$rate2=CRE.hit.rate$read.hit.smallRNAd/CRE.hit.rate$read.total.smallRNAd
write.table(CRE.hit.rate,gzfile(paste0(path_fig1_data,"CRE.hit.rate.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#for PAS motif location around the 3'end
#THP-1 get all the TES
path7=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/")
files=list.files(path=path7, pattern="read.info.tsv.gz")
files=files[grep("PAP", files)]
i=1
data=fread(paste0(path7,files[i]), stringsAsFactors = F, select=c(1,9))
for (i in 2:length(files)){
  data1=fread(paste0(path7,files[i]), stringsAsFactors = F, select=c(1,9))
  data=rbind(data,data1)}

all.thp1=fread(paste0(TES_path,"all_thp1.tsv.gz"), header=T, stringsAsFactors =F, check.names=F, select=c(1:8,11,12,17))
all.thp1$readname=gsub("_rep","",all.thp1$V4)
all.thp1$readname=gsub("DMSO_","D",all.thp1$readname)
all.thp1$readname=gsub("PMA_","P",all.thp1$readname)
all.thp1$readname=gsub("woPAP","N",all.thp1$readname)
all.thp1$readname=gsub("PAP","P",all.thp1$readname)
all.thp1$readname=gsub("4h_","",all.thp1$readname)
all.thp1$readname=gsub("6h_","",all.thp1$readname)

all.thp1=all.thp1[-grep("D9P2",all.thp1$readname),]
all.thp1=left_join(all.thp1,data, by=c("readname"="V4"),copy=F)
all.thp2=all.thp1%>%group_by(PAtailing,model_ID_str,label,internal_prime2)%>%dplyr::summarise(TES_count=n())

all3n=fread(paste0(TES_path,"all3n_up25.PASstandard.bed.gz"), header=F, stringsAsFactors = F)
all3n$V10[which(all3n$V7 == ".")]=max(all3n$V10)
all3n$V10=all3n$V10-25
colnames(all3n)[c(7,8,10)]=c("PAS_motif","motif_score","PAS_distance")
all3n=all3n[which(all3n$V4 %in% unique(all.thp2$label)),]
all3n=all3n%>%group_by(V4)%>%dplyr::mutate(n_V4=n())
all3na=all3n[which(all3n$n_V4 ==1),]
all3nb=all3n[which(all3n$n_V4 >1),]
all3nb=all3nb%>%group_by(V4)%>%dplyr::slice(which.max(motif_score))
all3n=rbind(all3na,all3nb)
write.table(all3n, gzfile(paste0(TES_path,"THP1_PAS_position.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
all.thp3=left_join(all.thp2, all3n[,c(4,7,8,10)], by=c("label"="V4"),copy=F)

#=take back transcript class from table 1 and table 5
T1=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table1.remove_undetect_and_partial_gencode_and_internal_prime.2.47M.tsv.gz"), header=T, check.names=F, stringsAsFactors = F)
T5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors = F)

T1$geneClass="others"
T1$geneClass[which(T1$Gencode_geneClass2 == "protein_coding")]="mRNA"
T1$geneClass[which(T1$Novel_geneClass == "lncRNA" | T1$Novel_geneClass == "ncRNA" | T1$Gencode_geneClass2 == "lncRNA")]="lncRNA"
T1$geneClass[which(T1$promoter_type == "promoter-like" & T1$geneClass == "lncRNA")]="p_lncRNA"
T1$geneClass[which(T1$promoter_type == "enhancer-like" & T1$geneClass == "lncRNA")]="e_lncRNA"
T1$transcriptClass="others"
T1$transcriptClass[which(T1$Gencode_transcriptClass2 == "protein_coding")]="mRNA"
T1$transcriptClass[which(T1$Novel_transcriptClass == "lncRNA" | T1$Novel_transcriptClass == "ncRNA" | T1$Gencode_transcriptClass2 == "lncRNA")]="lncRNA"
T1$transcriptClass[which(T1$promoter_type == "promoter-like" & T1$transcriptClass == "lncRNA")]="p_lncRNA"
T1$transcriptClass[which(T1$promoter_type == "enhancer-like" & T1$transcriptClass == "lncRNA")]="e_lncRNA"

T5$geneClass="others"
T5$geneClass[which(T5$T4_Gencode_geneCalss2 == "protein_coding")]="mRNA"
T5$geneClass[which(T5$T4_Novel_geneClass == "lncRNA" | T5$T4_Novel_geneClass == "ncRNA" | T5$T4_Gencode_geneCalss2 == "lncRNA")]="lncRNA"
T5$geneClass[which(T5$promoter_type == "promoter-like" & T5$geneClass == "lncRNA")]="p_lncRNA"
T5$geneClass[which(T5$promoter_type == "enhancer-like" & T5$geneClass == "lncRNA")]="e_lncRNA"
T5$transcriptClass="others"
T5$transcriptClass[which(T5$Gencode_transcriptClass2 == "protein_coding")]="mRNA"
T5$transcriptClass[which(T5$Novel_transcriptClass == "lncRNA" | T5$Novel_transcriptClass == "ncRNA" | T5$Gencode_transcriptClass2 == "lncRNA")]="lncRNA"
T5$transcriptClass[which(T5$promoter_type == "promoter-like" & T5$transcriptClass == "lncRNA")]="p_lncRNA"
T5$transcriptClass[which(T5$promoter_type == "enhancer-like" & T5$transcriptClass == "lncRNA")]="e_lncRNA"

colnames(T1)[c(1,96,97)] #"model_ID"  "geneClass" "transcriptClass"
colnames(T5)[c(1,121,122)] #"model_ID"  "geneClass" "transcriptClass"
all.thp3=left_join(all.thp3[,c(1:8)], T1[,c("model_ID",  "geneClass", "transcriptClass","IN1_gene_ID")],by=c("model_ID_str"="model_ID"),copy=F)
all.thp3=left_join(all.thp3, T5[,c("model_ID",  "geneClass", "transcriptClass", "T4_gene_ID")],by=c("model_ID_str"="model_ID"),copy=F, suffix=c("_T1","_T5"))

matrixT5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/partial_length_support_matrix.tsv.gz"), header=T, check.names=F, stringsAsFactors = F)
matrixT5=matrixT5[,-c(2:7)]
matrixT5$PAT=rowSums(matrixT5[,c(4,5,8,9,12,13,16,17)])
matrixT5$noPAT=rowSums(matrixT5[,c(2,3,6,7,10,11,14,15)])
matrixT5z=reshape2::melt(matrixT5[,c(1,18,19)],id=1)
matrixT5z=left_join(matrixT5z, T5[,c(1,121,122)], by="model_ID", copy=F)
matrixT5z=matrixT5z[which(!is.na(matrixT5z$geneClass)),]
matrixT5z$geneClass[grep("ncRNA",matrixT5z$geneClass)]="lncRNA"
matrixT5z1=matrixT5z%>%group_by(variable,geneClass)%>%dplyr::summarise(count=sum(value))%>%mutate(percent=count/sum(count))

#===============================================================================
#add back gencode 3n to filter internal priming here:
all.thp1=fread(paste0(TES_path,"all_thp1.tsv.gz"), header=T, stringsAsFactors =F, check.names=F, select=c(7,8))
all.thp1=unique(all.thp1)
all.thp3=left_join(all.thp3,all.thp1, by="label", copy=F)
all.thp3%>%group_by(internal_prime2,gencode3n)%>%summarise(count=n())
all.thp3$internal_prime2[which(all.thp3$gencode3n == "yes")]="no"

#===============================================================================
#locate the small RNA and histone gene
#source: https://www.genenames.org/data/genegroup/#!/group/864

sRNA=read.delim(paste0(TES_path,"GENCODE.rRNA.list.bed"), header=F, check.names=F, stringsAsFactors = F)
sRNA_transcript=unique(sRNA$V4)
histonelist=read.delim(paste0(TES_path,"histone.list.HGNC.118.txt"), header=T, stringsAsFactors = F, check.names = F)
histoneGene=histonelist$`Ensembl gene ID`

#remove the rRNA, snoRNA, histone gene
all.thp3=left_join(all.thp3, T1[,c( "model_ID","Gencode_geneClass")],by=c("model_ID_str"="model_ID"), copy=F)
all.thp3$class="others"
all.thp3$class[which(all.thp3$model_ID_str %in% sRNA_transcript)]="rRNA"
all.thp3$IN1_gene_ID=sapply(strsplit(all.thp3$IN1_gene_ID,"\\."),"[",1)
all.thp3$class[which(all.thp3$IN1_gene_ID %in% histoneGene)]="histone"
all.thp3$class[grep("snoRNA",all.thp3$Gencode_geneClass)]="snoRNA"
write.table(all.thp3, gzfile(paste0(TES_path,"THP1_PAT_modelID_TES_geneClass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# the file keeping the original read [all_thp1.tsv.gz] was removed
#===============================================================================
#Neuron series
path7=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/")
files=list.files(path=path7, pattern="read.info.tsv.gz")
files=files[-grep("PAP", files)]
i=1
data=fread(paste0(path7,files[i]), stringsAsFactors = F, select=c(1,9))
for (i in 2:length(files)){
  data1=fread(paste0(path7,files[i]), stringsAsFactors = F, select=c(1,9))
  data=rbind(data,data1)}

all.neuron1=fread(paste0(TES_path,"ontAll.read.0base.tsv.gz"), header=T, select=c(1,7,8,4,6,12,13))
all.neuron1$readname=gsub("_rep","",all.neuron1$V4)
all.neuron1$readname=gsub("_run","",all.neuron1$readname)
all.neuron1$readname=gsub("iPSC","I",all.neuron1$readname)
all.neuron1$readname=gsub("NSC","S",all.neuron1$readname)
all.neuron1$readname=gsub("Neuron","N",all.neuron1$readname)
all.neuron1=left_join(all.neuron1,data, by=c("readname"="V4"),copy=F)

all.neuron1$label=paste0(all.neuron1$V1,"_",all.neuron1$n3V2, "_", all.neuron1$n3V3, "_",all.neuron1$V6)
all.neuron1=all.neuron1[which(!is.na(all.neuron1$internal_prime)),]
all.neuron1$internal_prime2=all.neuron1$internal_prime
all.neuron1$internal_prime2[which(all.neuron1$internal_prime2 != "no")]="yes"
all.neuron2=all.neuron1%>%group_by(model_ID_str,label,internal_prime2)%>%dplyr::summarise(TES_count=n())

remove(all.neuron1)
remove(data)

all3n=fread(paste0(TES_path,"all3n_up25.PASstandard.bed.gz"), header=F, stringsAsFactors = F)
all3n$V10[which(all3n$V7 == ".")]=max(all3n$V10)
all3n$V10=all3n$V10-25
colnames(all3n)[c(7,8,10)]=c("PAS_motif","motif_score","PAS_distance")
all3n=all3n[which(all3n$V4 %in% unique(all.neuron2$label)),]
all3n=all3n%>%group_by(V4)%>%dplyr::mutate(n_V4=n())
all3na=all3n[which(all3n$n_V4 ==1),]
all3nb=all3n[which(all3n$n_V4 >1),]
all3nb=all3nb%>%group_by(V4)%>%dplyr::slice(which.max(motif_score))
all3n=rbind(all3na,all3nb)
write.table(all3n, gzfile(paste0(TES_path,"Neuron_PAS_position.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
all.neuron3=left_join(all.neuron2, all3n[,c(4,7,8,10)], by=c("label"="V4"),copy=F)

all.neuron3=left_join(all.neuron3, T1[,c(1,93,94)],by=c("model_ID_str"="model_ID"),copy=F)
all.neuron3=left_join(all.neuron3, T5[,c(1,113,114)],by=c("model_ID_str"="model_ID"),copy=F, suffix=c("_T1","_T5"))

#add back Gencode 3n to filter internal priming here:
all.neuron1=fread(paste0(TES_path,"ontAll.read.0base.tsv.gz"), header=T, stringsAsFactors =F, check.names=F, select=c(9,13))
all.neuron1=unique(all.neuron1)
all.neuron3=left_join(all.neuron3,all.neuron1, by=c("label"="n3ID"), copy=F)
all.neuron3%>%group_by(internal_prime2,gencode3n)%>%summarise(count=n())
all.neuron3$internal_prime2[which(all.neuron3$gencode3n == "yes")]="no"

#remove the rRNA, snoRNA, histone gene
all.neuron3=left_join(all.neuron3, T1[,c(1,42, 74)],by=c("model_ID_str"="model_ID"), copy=F)
all.neuron3$class="others"
all.neuron3$class[which(all.neuron3$model_ID_str %in% sRNA_transcript)]="rRNA"
all.neuron3$IN1_gene_ID=sapply(strsplit(all.neuron3$IN1_gene_ID,"\\."),"[",1)
all.neuron3$class[which(all.neuron3$IN1_gene_ID %in% histoneGene)]="histone"
all.neuron3$class[grep("snoRNA",all.neuron3$Gencode_geneClass)]="snoRNA"
write.table(all.neuron3, gzfile(paste0(TES_path,"Neuron_PAT_modelID_TES_geneClass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

# the file keeping the original read [ontAll.read.0base.tsv.gz] was removed
#==================================================
# plot file for ex1e
all.thp3=read.delim(paste0(TES_path,"THP1_PAT_modelID_TES_geneClass.tsv.gz"), header=T, stringsAsFactor=F, check.names=F)
all.neuron3=read.delim(paste0(TES_path,"Neuron_PAT_modelID_TES_geneClass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all.neuron6=all.neuron3%>%group_by(PAtailing,class)%>%dplyr::summarise(count=sum(TES_count))%>%dplyr::mutate(percent=count/sum(count))
all.thp6=all.thp3%>%group_by(PAtailing,class)%>%dplyr::summarise(count=sum(TES_count))%>%dplyr::mutate(percent=count/sum(count))
all6=rbind(all.neuron6, all.thp6)
all6$group=c(rep("Neuron.series_PAT",4),rep("THP-1.series_PAT",4),rep("THP-1.series_noPAT",4))
write.table(all6, gzfile(paste0(path_fig1_data,"non_polyA_gene_group.tsv.gz")), col.names=T, row.names = F, sep="\t", quote=F)


#==================================================
#TSS and TES intersect with different properties
options(scipen=999)
gencodegtf=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.gtf.gz"), skip=5, header=F, stringsAsFactors = F)
unique(gencodegtf$V3)
gencodegtf$geneID=sapply(strsplit(gencodegtf$V9,"; "),"[",1)
gencodegtf$transcriptID=sapply(strsplit(gencodegtf$V9,"; "),"[",2)
gencodegtf$geneClass=sapply(strsplit(gencodegtf$V9,"; "),"[",3)
gencodegtf$transcriptClass=sapply(strsplit(gencodegtf$V9,"; "),"[",5)
gencodegtf$exonNumber=sapply(strsplit(gencodegtf$V9,"; "),"[",7)

gencodegtfu=gencodegtf[which(gencodegtf$V3 == "UTR"),]
gencodegtfustart=unique(gencodegtf[which(gencodegtf$V3=="start_codon"),c(4,5,11)])
gencodegtfustart=gencodegtfustart%>%group_by(transcriptID)%>%dplyr::slice_min(V4)
gencodegtfuend=unique(gencodegtf[which(gencodegtf$V3=="stop_codon"),c(4,5,11)])
gencodegtfuend=gencodegtfuend%>%group_by(transcriptID)%>%dplyr::slice_min(V4)
gencodegtfu=left_join(gencodegtfu, gencodegtfustart, by="transcriptID", copy=F, suffix=c("","_start"))
gencodegtfu=left_join(gencodegtfu, gencodegtfuend, by="transcriptID", copy=F, suffix=c("","_end"))

gencodegtfup=gencodegtfu[which(gencodegtfu$V7 == "+"),]
gencodegtfun=gencodegtfu[which(gencodegtfu$V7 == "-"),]
gencodegtfup$ID=NA
gencodegtfup$ID[which(gencodegtfup$V4 < gencodegtfup$V4_start)]="5UTR"
gencodegtfup$ID[which(gencodegtfup$V5 > gencodegtfup$V5_end)]="3UTR"
gencodegtfup$size=gencodegtfup$V5-gencodegtfup$V4
gencodegtfun$ID=NA
gencodegtfun$ID[which(gencodegtfun$V5 > gencodegtfun$V5_start)]="5UTR"
gencodegtfun$ID[which(gencodegtfun$V4 < gencodegtfun$V4_end)]="3UTR"
gencodegtfun$size=gencodegtfun$V5-gencodegtfun$V4
gencodegtfu=rbind(gencodegtfup,gencodegtfun)
gencodegtfu1=gencodegtfu[which(gencodegtfu$geneClass=="gene_type protein_coding"),]
gencodegtfu1$V4=gencodegtfu1$V4-1
gencodegtfu1$transcriptID=gsub("transcript_id ","",gencodegtfu1$transcriptID)
gencodegtfu1$exonNumber=gsub("exon_number ","",gencodegtfu1$exonNumber)
gencodegtfu1$transcriptID=paste0(gencodegtfu1$transcriptID,"_",gencodegtfu1$exonNumber)
UTR3=gencodegtfu1[which(gencodegtfu1$ID=="3UTR"),c(1,4,5,11,6,7)]
UTR5=gencodegtfu1[which(gencodegtfu1$ID=="5UTR"),c(1,4,5,11,6,7)]

write.table(UTR3[order(UTR3$V1,UTR3$V4),],gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.3UTR.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(UTR5[order(UTR5$V1,UTR5$V4),],gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.5UTR.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#make gencode 3' end +/- 5 bp bed
gencodegtft=gencodegtf[which(gencodegtf$V3 == "transcript"),]
gencodegtft$V4=gencodegtft$V4-1
gencodegtft$V43n=gencodegtft$V5-1
gencodegtft$V53n=gencodegtft$V5
gencodegtft$V43n[which(gencodegtft$V7 == "-")]=gencodegtft$V4[which(gencodegtft$V7 == "-")]
gencodegtft$V53n[which(gencodegtft$V7 == "-")]=gencodegtft$V4[which(gencodegtft$V7 == "-")]+1
gencodegtft$V43n=gencodegtft$V43n-5
gencodegtft$V53n=gencodegtft$V53n+5
gencodegtft1=unique(gencodegtft[,c(1,15,16,6,7)])
gencodegtft1$name=paste0("n",1:nrow(gencodegtft1))
write.table(gencodegtft1[order(gencodegtft1$V1, gencodegtft1$V43n),c(1:3,6,4,5)],gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.n3_pm5bp.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
gencodegtft1=read.delim(gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.n3_pm5bp.bed.gz")), header=F, stringsAsFactors = F)

#make gencode junction bed
intron=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.bed6.intron.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
intronp=intron[which(intron$V6 == "+"),]
intronn=intron[which(intron$V6 == "-"),]
intronp$donorV2=intronp$V2-1
intronp$donorV3=intronp$V2+1
intronp$acceptorV2=intronp$V3-1
intronp$acceptorV3=intronp$V3+1

intronn$donorV2=intronn$V3-1
intronn$donorV3=intronn$V3+1
intronn$acceptorV2=intronn$V2-1
intronn$acceptorV3=intronn$V2+1
intron=rbind(intronp,intronn)
write.table(intron[order(intron$V1,intron$donorV2),c(1,7,8,4,5,6)], gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.donor2bp.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(intron[order(intron$V1,intron$acceptorV2),c(1,9,10,4,5,6)], gzfile(paste0(GENCODE_path,"gencode.v39.annotation.bed6.acceptor2bp.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

n3bed=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed"), header=F, stringsAsFactors = F, check.names = F)
n5bed=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end5.bed"), header=F, stringsAsFactors = F, check.names = F)
#===============================================================================

#========================
#bedtools intersect
#bash
setwd(intersect_path)
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.donor2bp.bed.gz -s | gzip > n3_donor.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.acceptor2bp.bed.gz -s | gzip > n3_acceptor.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.intron.bed.gz -s | gzip > n3_intron.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.exon.bed.gz -s | gzip > n3_exon.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.3UTR.bed.gz -s | gzip > n3_3UTR.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.5UTR.bed.gz -s | gzip > n3_5UTR.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed -b ",GENCODE_path,"gencode.v39.annotation.bed6.n3_pm5bp.bed.gz -s | gzip > n3_n3_pm5bp.bed.gz"))

#===============================================================================
n3bed=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed"), header=F, stringsAsFactors = F, check.names = F)
n3bed$label=paste0(n3bed$V1,"_",n3bed$V2,"_",n3bed$V3,"_",n3bed$V6)
files=list.files(path=intersect_path, pattern=".bed.gz")
files.names=gsub(".bed.gz","",files)
for (i in 1:length(files)){
  data=read.delim(paste0(intersect_path, files[i]),header=F, stringsAsFactors=F)
  data$label=paste0(data$V1,"_",data$V2,"_",data$V3,"_",data$V6)
  data$V7[which(data$V7>0)]=1
  n3bed=left_join(n3bed, data[,c(7,8)], by="label", copy=F)}
colnames(n3bed)[c(8:14)]=files.names
write.table(n3bed, gzfile(paste0(TES_path, "n3bed_all_info.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=============================
# plot table for ex1h
all.thp3=all.thp3%>%group_by(IN1_gene_ID)%>%dplyr::mutate(percent_transcript=TES_count/sum(TES_count))
all.thp5=all.thp3[which(!is.na(all.thp3$IN1_gene_ID)),]
all.thp5=all.thp5[which(all.thp5$PAS_distance < (-35) | all.thp5$PAS_distance > (-5)),]
all.thp5=all.thp5[which(all.thp5$internal_prime2 == "no"),]
all.thp5$group="THP1"
all.neuron3=all.neuron3%>%group_by(PAtailing,IN1_gene_ID)%>%dplyr::mutate(percent_transcript=TES_count/sum(TES_count))
all.neuron5=all.neuron3[which(!is.na(all.neuron3$IN1_gene_ID)),]
all.neuron5=all.neuron5[which(all.neuron5$PAS_distance < (-35) | all.neuron5$PAS_distance > (-5)),]
all.neuron5=all.neuron5[which(all.neuron5$internal_prime2 == "no"),]
all.neuron5$group="Neuron"
all5=rbind(all.thp5,all.neuron5)
all5=all5[which(all5$geneClass_T1 != "others"),]
all5=all5[which(all5$class == "others"),] #remove histon, snoRNA etc
n3bed=read.delim(paste0(TES_path, "n3bed_all_info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all5=left_join(all5, n3bed[,c(7:14)], by="label", copy=F)

all5$class2="Intergenic/others"
all5$class2[which(all5$n3_intron >0)]="Intron"
all5$class2[which(all5$n3_exon >0)]="Exon"
all5$class2[which(all5$n3_3UTR >0)]="3'UTR"
all5$class2[which(all5$n3_donor >0)]="Donor_site"
all5$class2[which(all5$n3_n3_pm5bp > 0)]="GENCODE_3n"

write.table(all5, gzfile(paste0(TES_path,"neuron_THP1_TES_genomic_location.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

allsummary=all5[which(all5$geneClass_T1 %in% c("mRNA","p_lncRNA","e_lncRNA")),]%>%group_by(PAtailing, geneClass_T1, class2)%>%dplyr::summarise(count=sum(TES_count))%>%dplyr::mutate(percent=count/sum(count))
allsummary1=all.thp3[which(all.thp3$geneClass_T1 %in% c("mRNA","p_lncRNA","e_lncRNA") & all.thp3$internal_prime2 == "no" & all.thp3$class == "others"),]%>%group_by(PAtailing, geneClass_T1)%>%dplyr::summarise(total=sum(TES_count))
allsummary2=all.neuron3[which(all.neuron3$geneClass_T1 %in% c("mRNA","p_lncRNA","e_lncRNA")  & all.neuron3$internal_prime2 == "no" & all.neuron3$class == "others"),]%>%group_by(PAtailing, geneClass_T1)%>%dplyr::summarise(total=sum(TES_count))
allsummary3=rbind(allsummary1,allsummary2)
allsummary=left_join(allsummary,allsummary3, by=c("PAtailing","geneClass_T1"),copy=F)
allsummary$percent2=allsummary$count/allsummary$total
allsummary$PAtailing[which(allsummary$PAtailing == "PAP")]="THP-1.series_PAT"
allsummary$PAtailing[which(allsummary$PAtailing == "woPAP")]="THP-1.series_noPAT"
allsummary$PAtailing=factor(allsummary$PAtailing, levels=c("Neuron.series_PAT","THP-1.series_PAT","THP-1.series_noPAT"))
allsummary$class2=factor(allsummary$class2, levels=c("GENCODE_3n","Intron","Donor_site","Exon","3'UTR","Intergenic/others"))
allsummary$scale=allsummary$percent/allsummary$percent2
allsummary$label=paste0(signif(allsummary$percent,2)*100,"%")

write.table(allsummary,gzfile(paste0(path_fig1_data,"neuron_THP1_TES_genomic_location.summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#find PAT-specific TES and gene with PAS end >80%
all.thp4=all.thp3[which(all.thp3$PAS_distance >= (-35) & all.thp3$PAS_distance <= (-5)),] #real polyA
all.thp4a=all.thp4%>%group_by(IN1_gene_ID,geneClass_T1)%>%dplyr::summarise(percent_transcript=sum(percent_transcript))
all.thp4b=all.thp4a[which(all.thp4a$percent_transcript >0.8),]
all.thp4a%>%group_by(geneClass_T1)%>%dplyr::summarise(count=n())
all.thp4b%>%group_by(geneClass_T1)%>%dplyr::summarise(count=n())

#exc=unique(all5$label[which(all5$PAtailing == "woPAP" & all5$TES_count >=3)])
exc1=unique(all5$label[which(all5$PAtailing == "woPAP")])
all7=all5[which(all5$PAtailing=="PAP" & all5$group == "THP1"),]
all7a=all7[-which(all7$label %in% exc1),]
all7aNascent=all7a[which(all7a$IN1_gene_ID %in% all.thp4b$IN1_gene_ID),]
all7aOther=all7a[-which(all7a$IN1_gene_ID %in% all.thp4b$IN1_gene_ID),]
all7aNascent$group2="PAT-specific: Partial end"
all7aOther$group2="PAT-specific: Others"
all7b=all5[which(all5$PAtailing=="woPAP" & all5$group == "THP1"),]
all7b$group2="no PAT"
all7c=rbind(all7aNascent,all7aOther,all7b)
write.table(all7c, gzfile(paste0(path_fig1_data,"noPAT_PATspecific_TES_location.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#============================
# plot table for ex1k
all.thp3=read.delim(paste0(TES_path,"THP1_PAT_modelID_TES_geneClass.tsv.gz"), header=T, stringsAsFactor=F, check.names=F)
all.thp7=all.thp3[which(all.thp3$internal_prime2 == "no" & all.thp3$class == "others" & !is.na(all.thp3$geneClass_T5)),]

PATspecific=all.thp7[which(all.thp7$PAtailing == "PAP" & all.thp7$T4_gene_ID %in% PATspecificRNA),]
n1=unique(PATspecific[,c("T4_gene_ID","geneClass_T5")])%>%group_by(geneClass_T5)%>%dplyr::summarise(count=n())
noPATspecific=all.thp7[which(all.thp7$PAtailing == "woPAP" & all.thp7$T4_gene_ID %in% noPATspecificRNA),]
n3=unique(noPATspecific[,c("T4_gene_ID","geneClass_T5")])%>%group_by(geneClass_T5)%>%dplyr::summarise(count=n())
common=all.thp7[which(all.thp7$T4_gene_ID %in% both),]
n2=unique(common[,c("T4_gene_ID","geneClass_T5")])%>%group_by(geneClass_T5)%>%dplyr::summarise(count=n())
vennsummary=rbind(n1,n2,n3)
vennsummary$group=c(rep("PAT-specific",5),rep("common",5),rep("noPAT-specific",5))
vennsummary$geneClass_T5[grep("lncRNA",vennsummary$geneClass_T5)]="lncRNA"
vennsummary=vennsummary%>%group_by(group,geneClass_T5)%>%dplyr::summarise(count=sum(count))%>%dplyr::mutate(percent=count/sum(count))
write.table(vennsummary, gzfile(paste0(path_fig1_data,"noPAT_PATspecific_venn.summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


