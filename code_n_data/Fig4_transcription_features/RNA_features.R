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
path_fig4_data=paste0(primary_folder,"fig4/data/")
exo_path=paste0(primary_folder,"code_n_data/Fig4_transcription_features/exosome_sensitivity/")
chromatin_path=paste0(primary_folder,"code_n_data/Fig4_transcription_features/quantification_chr_total/")
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/")

#===============================================================================
#collapse reads into ex5_cluster with length, genomic region and number of exon 
#take the read info from SALA output in the /log folder
path1=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/")

files=list.files(path=path1, pattern="read.info.tsv.gz")
i=1
data=fread(paste0(path1,files[i]), stringsAsFactors = F, select=c(1:7,9))
data=data[which(data$model_ID_str %in% table5$model_ID),]
data$V4=substr(data$V4, start = 1, stop = 1)
data$genomic_range=data$V3-data$V2
for (i in 2:length(files)){
  data1=fread(paste0(path1,files[i]), stringsAsFactors = F, select=c(1:7,9))
  data1=data1[which(data1$model_ID_str %in% table5$model_ID),]
  data1$V4=substr(data1$V4, start = 1, stop = 1)
  data1$genomic_range=data1$V3-data1$V2
  data=rbind(data,data1)}

data=left_join(data,table5[,c("model_ID","n5_string")],by=c("model_ID_str"="model_ID"),copy=F)
data1c=data%>%group_by(n5_string)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
data1e=data%>%group_by(V4,n5_string)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())

write.table(data1c,gzfile(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
# for -> table S9
write.table(data1e,gzfile(paste0(path_fig4_data,"length.exon.info.end5_cluster_celltype.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#collapse the read length and exon number into transcript model base
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

data2=data%>%group_by(model_ID_str)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), n_read=n())
data2a=data%>%group_by(V4,model_ID_str)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), n_read=n())

data2=left_join(data2,table5[,c(1,101:105,63:68,2,11,56,57,86,90,93,14)], by=c("model_ID_str"="model_ID"),copy=F)
data2=left_join(data2,table5[,c(1,107,108,113)], by=c("model_ID_str"="model_ID"),copy=F) #add poly(A)
data2$loc=sapply(strsplit(data2$loc,":"),"[",2)
data2$transcript_range=as.numeric(sapply(strsplit(data2$loc,"-"),"[",2))-as.numeric(sapply(strsplit(data2$loc,"-"),"[",1))
write.table(data2,gzfile(paste0(path_fig4_data,"features_by_transcriptModel.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# for -> table S10
write.table(data2a,gzfile(paste0(path_fig4_data,"length.exon.info.withcelltype.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

rm(data)
gc()

#===============================================================================
#define ex5_cluster_class from SALA Final 
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

table5=table5%>%group_by(n5_string)%>%dplyr::mutate(ex5cluster_ratio_overall=full_qry_count/sum(full_qry_count))
table5enst=table5[which(table5$Gencode_transcriptClass2 %in% c("lncRNA") ),]
table5novel=table5[which(table5$Novel_transcriptClass %in% c("ncRNA","lncRNA") ),]
table5both=rbind(table5enst,table5novel)
table5both1=table5both%>%group_by(promoter_type,ATAC,n5_string)%>%dplyr::summarise(ncRNA_rate=sum(ex5cluster_ratio_overall))
table5both1=table5both1[which(table5both1$ncRNA_rate > 0.5),]
table5both1$ex5cluster_class=paste0(substr(table5both1$promoter_type,1,1),"_ncRNA")
table5both1$ex5cluster_class=gsub("C_ncRNA","CTCF_ncRNA",table5both1$ex5cluster_class)
table5both1$ex5cluster_class=gsub("u_ncRNA","other_ncRNA",table5both1$ex5cluster_class)
table5both1$ex5cluster_class[which(table5both1$ex5cluster_class %in% c("CTCF_ncRNA","other_ncRNA") & table5both1$ATAC == "noATAC")]=NA

table5enst=unique(table5$n5_string[which(table5$Gencode_transcriptClass2 == "protein_coding" )])
table5enstnc=unique(table5$n5_string[which(table5$Gencode_transcriptClass2 == "lncRNA" )])
table5enstot=unique(table5$n5_string[which(table5$Gencode_transcriptClass2 == "others" )])
table5enstfinal=setdiff(table5enst,table5enstnc)
table5enstfinal=setdiff(table5enstfinal,table5enstot)

table5=left_join(table5,table5both1[,c(3,5)], by = "n5_string", copy=F)
table5$ex5cluster_class[which(table5$n5_string %in% table5enstfinal)]="mRNA"

kk=table5[which(table5$ex5cluster_class == "mRNA" & table5$Novel_transcriptClass == "lncRNA"),]
kk1=table5[which(table5$ex5cluster_class == "mRNA" & table5$Novel_transcriptClass == "short_ncRNA"),]
kk2=table5[which(table5$ex5cluster_class == "p_ncRNA" & table5$Gencode_transcriptClass2 == "protein_coding"),]
kk3=table5[which(table5$ex5cluster_class == "e_ncRNA" & table5$Gencode_transcriptClass2 == "protein_coding"),]

data12=table5%>%group_by(n5_string,ex5cluster_class)%>%dplyr::summarise(Gencode_geneClass2=paste(unique(Gencode_geneClass2),collapse=";"),Gencode_transcriptClass2=paste(unique(Gencode_transcriptClass2), collapse=";"),Novel_transcriptClass=paste(unique(Novel_transcriptClass), collapse=";"), promoter_type=paste(unique(promoter_type), collapse=";"))
data12a=data12[which(data12$n5_string %in% c(kk$n5_string, kk1$n5_string,kk2$n5_string, kk3$n5_string)),]
data12c=data12[-which(data12$n5_string %in% c(kk$n5_string, kk1$n5_string,kk2$n5_string, kk3$n5_string)),]
table5$ex5cluster_class[which(table5$n5_string %in% data12a$n5_string)]=NA

data1a=data12c[which(!is.na(data12c$ex5cluster_class)),]

#also remove ex5_cluster that have gene level conflict
data1k=data1a[-which(data1a$ex5cluster_class != "mRNA" & data1a$Gencode_geneClass2=="protein_coding"),]
data1k=data1k[-which(data1k$ex5cluster_class == "mRNA" & data1k$Gencode_geneClass2!="protein_coding"),]

#revise data1
data1=right_join(data1,data1k[,c(1,3,4,5)], by="n5_string",copy=F)
write.table(data1,gzfile(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
#number of ex5_cluster -> 58109
table5$ex5cluster_class[-which(table5$n5_string %in% data1$n5_string)]=NA

#===============================================================================
#ex5_cluster-based analyses
n5cluster_info=read.delim(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(n5cluster_info)[4]="n5_string"
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

#join the ex5_cluster feature (from Fig.S1) with the ex5_cluster base RNA features

data1=left_join(data1, n5cluster_info, by="n5_string",copy=F)
data1=data1[which(!is.na(data1$ex5cluster_class)),]
data1%>%group_by(ex5cluster_class,ATAC)%>%dplyr::summarise(count=n())
data1%>%group_by(promoter_type,ATAC)%>%dplyr::summarise(count=n())

data1$ex5cluster_class=factor(data1$ex5cluster_class,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA", "other_ncRNA"))
data1%>%group_by(ex5cluster_class)%>%dplyr::summarise(rlength=median(read_median_length),exon=median(median_exon),range=median(median_range), count=n())


data1$CpGTATA[which(data1$any_CpG_island == "Yes" & data1$TATA_box == "No" & data1$distance_cCRE_PLS < 2000)]="CGIap"
data1$CpGTATA[which(data1$any_CpG_island == "Yes" & data1$TATA_box == "No" & data1$distance_cCRE_PLS >= 2000)]="CGInap"
data1$CpGTATA[which(data1$CpGTATA %in% c("uCGI","dCGI"))]="CGI"

data1$CpGTATA[which(data1$CpGTATA %in% c("CGInap","CGIap") & data1$ex5cluster_class=="mRNA")]="CGI"
data1$CGIap="No"
data1$CGIap[which(data1$CpGTATA == "CGIap")]="Yes"
data1$CGInap="No"
data1$CGInap[which(data1$CpGTATA == "CGInap")]="Yes"
data1$CGI="No"
data1$CGI[which(data1$CpGTATA == "CGI")]="Yes"
data1$TATA="No"
data1$TATA[which(data1$CpGTATA == "TATA")]="Yes"

#=====
#update bidirectionality from Andersson study
data1$Andersson_permissive[which(data1$Andersson_permissive==0)]="others"
data1$Andersson_permissive[which(data1$Andersson_permissive==1)]="bidirectional"
data1$Andersson_robust[which(data1$Andersson_robust==0)]="others"
data1$Andersson_robust[which(data1$Andersson_robust==1)]="bidirectional"

#=====
#update conservation
data1$conserve="Others"
data1$conserve[which(data1$mean_up500 <= median(data1$mean_up500, na.rm=T))]="No"
data1$conserve[which(data1$mean_up500 > median(data1$mean_up500, na.rm=T))]="Yes"
write.table(data1,gzfile(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
#stored as "features_by_ex5cluster.tsv.gz" in [primary_folder]/fig4/data
# for -> Fig. 4a-g, i,j & table S9



#enrichment across features
#===============================================================================
#enrichment with CGIap ex5_cluster
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$conserve_up="No"
data1$conserve_up[which(data1$phastCon17_mean_up500 > median(data1$phastCon17_mean_up500, na.rm=T))]="Yes"
data1$conserve_down="No"
data1$conserve_down[which(data1$phastCon17_mean_down500 > median(data1$phastCon17_mean_down500, na.rm=T))]="Yes"

a3=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGIap,orientation)%>%dplyr::summarise(count=n(), group="2D")
a4=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGIap,SE_all)%>%dplyr::summarise(count=n(), group="SE")
a5=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGIap,ubiquitous)%>%dplyr::summarise(count=n(), group="ubiquitous")
a6=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGIap,conserve_up)%>%dplyr::summarise(count=n(), group="conserve_up")
a7=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGIap,conserve_down)%>%dplyr::summarise(count=n(), group="conserve_down")

colnames(a4)=colnames(a3)
colnames(a5)=colnames(a3)
colnames(a6)=colnames(a3)
colnames(a7)=colnames(a3)
a1=rbind(a3,a4,a5,a6,a7)
a1=a1[which(a1$orientation != "Others"),]
a1=a1[which(!is.na(a1$orientation)),]
a1$orientation[which(a1$orientation == "1D")]="No"
a1$orientation[which(a1$orientation == "2D")]="Yes"
a1$orientation[which(a1$orientation == "TE")]="No"
a1$orientation[which(a1$orientation == "SE")]="Yes"
a1$label=paste0(a1$CGIap,"_",a1$orientation)
a2=spread(a1[,c(1,5,6,4)], key=3,value=4)
for(i in 1:nrow(a2)){
  GSEATasting <- matrix(c(a2$Yes_Yes[i], a2$Yes_No[i], a2$No_Yes[i], a2$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  a2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  a2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
a2$FE_logOR=log(a2$OR)
a2$sig_level="ns"
a2$sig_level[which(a2$p.val<0.05)]="*"
a2$sig_level[which(a2$p.val<0.01)]="**"
a2$sig_level[which(a2$p.val<0.001)]="***"

a2$sig_level=factor(a2$sig_level, levels=c("ns","*","**","***"))

#enrichment with CGInap ex5_cluster
b3=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGInap,orientation)%>%dplyr::summarise(count=n(), group="2D")
b4=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGInap,SE_all)%>%dplyr::summarise(count=n(), group="SE")
b5=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGInap,ubiquitous)%>%dplyr::summarise(count=n(), group="ubiquitous")
b6=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGInap,conserve_up)%>%dplyr::summarise(count=n(), group="conserve_up")
b7=data1[which(data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CGInap,conserve_down)%>%dplyr::summarise(count=n(), group="conserve_down")

colnames(b4)=colnames(b3)
colnames(b5)=colnames(b3)
colnames(b6)=colnames(b3)
colnames(b7)=colnames(b3)
b1=rbind(b3,b4,b5,b6,b7)
b1=b1[which(b1$orientation != "Others"),]
b1=b1[which(!is.na(b1$orientation)),]
b1$orientation[which(b1$orientation == "1D")]="No"
b1$orientation[which(b1$orientation == "2D")]="Yes"
b1$orientation[which(b1$orientation == "TE")]="No"
b1$orientation[which(b1$orientation == "SE")]="Yes"
b1$label=paste0(b1$CGInap,"_",b1$orientation)
b2=spread(b1[,c(1,5,6,4)], key=3,value=4)
for(i in 1:nrow(b2)){
  GSEATasting <- matrix(c(b2$Yes_Yes[i], b2$Yes_No[i], b2$No_Yes[i], b2$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  b2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  b2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
b2$FE_logOR=log(b2$OR)
b2$sig_level="ns"
b2$sig_level[which(b2$p.val<0.05)]="*"
b2$sig_level[which(b2$p.val<0.01)]="**"
b2$sig_level[which(b2$p.val<0.001)]="***"

b2$sig_level=factor(b2$sig_level, levels=c("ns","*","**","***"))

#======
#enrichment with CGI ex5_cluster
c3=data1[which(data1$ex5cluster_class %in% c("mRNA","p_ncRNA","other_ncRNA")),]%>%group_by(ex5cluster_class,CGI,orientation)%>%dplyr::summarise(count=n(), group="2D")
c4=data1[which(data1$ex5cluster_class %in% c("mRNA","p_ncRNA","other_ncRNA")),]%>%group_by(ex5cluster_class,CGI,SE_all)%>%dplyr::summarise(count=n(), group="SE")
c5=data1[which(data1$ex5cluster_class %in% c("mRNA","p_ncRNA","other_ncRNA")),]%>%group_by(ex5cluster_class,CGI,ubiquitous)%>%dplyr::summarise(count=n(), group="ubiquitous")
c6=data1[which(data1$ex5cluster_class %in% c("mRNA","p_ncRNA","other_ncRNA")),]%>%group_by(ex5cluster_class,CGI,conserve_up)%>%dplyr::summarise(count=n(), group="conserve_up")
c7=data1[which(data1$ex5cluster_class %in% c("mRNA","p_ncRNA","other_ncRNA")),]%>%group_by(ex5cluster_class,CGI,conserve_down)%>%dplyr::summarise(count=n(), group="conserve_down")

colnames(c4)=colnames(c3)
colnames(c5)=colnames(c3)
colnames(c6)=colnames(c3)
colnames(c7)=colnames(c3)
c1=rbind(c3,c4,c5,c6,c7)
c1=c1[which(c1$orientation != "Others"),]
c1=c1[which(!is.na(c1$orientation)),]
c1$orientation[which(c1$orientation == "1D")]="No"
c1$orientation[which(c1$orientation == "2D")]="Yes"
c1$orientation[which(c1$orientation == "TE")]="No"
c1$orientation[which(c1$orientation == "SE")]="Yes"
c1$label=paste0(c1$CGI,"_",c1$orientation)
c2=spread(c1[,c(1,5,6,4)], key=3,value=4)
c2[is.na(c2)]=0
for(i in 1:nrow(c2)){
  GSEATasting <- matrix(c(c2$Yes_Yes[i], c2$Yes_No[i], c2$No_Yes[i], c2$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  c2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  c2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
c2$FE_logOR=log(c2$OR)
c2$sig_level="ns"
c2$sig_level[which(c2$p.val<0.05)]="*"
c2$sig_level[which(c2$p.val<0.01)]="**"
c2$sig_level[which(c2$p.val<0.001)]="***"
c2=c2[-which(c2$group == "SE" & c2$ex5cluster_class %in% c("mRNA","p_ncRNA")),]

##TATA
d3=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,TATA,orientation)%>%dplyr::summarise(count=n(), group="2D")
d4=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,TATA,SE_all)%>%dplyr::summarise(count=n(), group="SE")
d5=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,TATA,ubiquitous)%>%dplyr::summarise(count=n(), group="ubiquitous")
d6=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,TATA,conserve_up)%>%dplyr::summarise(count=n(), group="conserve_up")
d7=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,TATA,conserve_down)%>%dplyr::summarise(count=n(), group="conserve_down")

colnames(d4)=colnames(d3)
colnames(d5)=colnames(d3)
colnames(d6)=colnames(d3)
colnames(d7)=colnames(d3)
d1=rbind(d3,d4,d5,d6,d7)
d1=d1[which(d1$orientation != "Others"),]
d1=d1[which(!is.na(d1$orientation)),]
d1$orientation[which(d1$orientation == "1D")]="No"
d1$orientation[which(d1$orientation == "2D")]="Yes"
d1$orientation[which(d1$orientation == "TE")]="No"
d1$orientation[which(d1$orientation == "SE")]="Yes"
d1$label=paste0(d1$TATA,"_",d1$orientation)
d2=spread(d1[,c(1,5,6,4)], key=3,value=4)
d2[is.na(d2)]=0
for(i in 1:nrow(d2)){
  GSEATasting <- matrix(c(d2$Yes_Yes[i], d2$Yes_No[i], d2$No_Yes[i], d2$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  d2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  d2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
d2$FE_logOR=log(d2$OR)
d2$sig_level="ns"
d2$sig_level[which(d2$p.val<0.05)]="*"
d2$sig_level[which(d2$p.val<0.01)]="**"
d2$sig_level[which(d2$p.val<0.001)]="***"
d2=d2[-which(d2$group == "SE" & d2$ex5cluster_class %in% c("mRNA","p_ncRNA")),]

##Null
data1$Null="No"
data1$Null[which(data1$CpGTATA == "Null")]="Yes"
e3=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,Null,orientation)%>%dplyr::summarise(count=n(), group="2D")
e4=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,Null,SE_all)%>%dplyr::summarise(count=n(), group="SE")
e5=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,Null,ubiquitous)%>%dplyr::summarise(count=n(), group="ubiquitous")
e6=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,Null,conserve_up)%>%dplyr::summarise(count=n(), group="conserve_up")
e7=data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,Null,conserve_down)%>%dplyr::summarise(count=n(), group="conserve_down")

colnames(e4)=colnames(e3)
colnames(e5)=colnames(e3)
colnames(e6)=colnames(e3)
colnames(e7)=colnames(e3)
e1=rbind(e3,e4,e5,e6,e7)
e1=e1[which(e1$orientation != "Others"),]
e1=e1[which(!is.na(e1$orientation)),]
e1$orientation[which(e1$orientation == "1D")]="No"
e1$orientation[which(e1$orientation == "2D")]="Yes"
e1$orientation[which(e1$orientation == "TE")]="No"
e1$orientation[which(e1$orientation == "SE")]="Yes"
e1$label=paste0(e1$Null,"_",e1$orientation)
e2=spread(e1[,c(1,5,6,4)], key=3,value=4)
for(i in 1:nrow(e2)){
  GSEATasting <- matrix(c(e2$Yes_Yes[i], e2$Yes_No[i], e2$No_Yes[i], e2$No_No[i]), nrow = 2, dimnames = list(oligo1 = c("yes", "no"), oligo2 = c("yes", "no")))
  e2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  e2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
e2$FE_logOR=log(e2$OR)
e2$sig_level="ns"
e2$sig_level[which(e2$p.val<0.05)]="*"
e2$sig_level[which(e2$p.val<0.01)]="**"
e2$sig_level[which(e2$p.val<0.001)]="***"
e2=e2[-which(e2$group == "SE" & e2$ex5cluster_class %in% c("mRNA","p_ncRNA")),]

a2$enrichment="CGIap"
b2$enrichment="CGInap"
c2$enrichment="CGI"
d2$enrichment="TATA"
e2$enrichment="Null"
e2=rbind(a2,b2,c2,d2,e2)

write.table(e2,gzfile(paste0(path_fig4_data,"FE_CGI_dCGI_TATA.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig4/data
# for -> Fig. 4g

#===============================================================================
#transcript model base: polyA and splicing
data1=read.delim(paste0(path_fig4_data,"features_by_transcriptModel.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$group=factor(data1$group,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA","other_ncRNA"))

#splicing & polyA  separately
data2=data1[which(data1$group=="e_ncRNA" | data1$group=="p_ncRNA"), c("group","model_ID_str","orientation","read_median_length","n_exon")]
data2$orientation[which(data2$n_exon>1)]="Yes"
data2$orientation[which(data2$n_exon==1)]="No"
data3=data1[which(data1$group=="e_ncRNA" | data1$group=="p_ncRNA"),]
data3=data3[which(data3$TES_recur == "Yes"), c("group","model_ID_str","polyA","read_median_length","n_exon")]
data3$polyA[which(data3$polyA=="poly(A)")]="Yes"
data3$polyA[which(data3$polyA=="non-poly(A)")]="No"

colnames(data3)[3]="orientation"
data2$group2="Spliced"
data3$group2="poly(A)"

data2=rbind(data2,data3)
colnames(data2)[3]="value"
write.table(data2, gzfile(paste0(path_fig4_data,"transcript.splice.and.polyA.RNA.length.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig4/data
# for -> Fig. 4h

#===============================================================================
#relative MAD of different ex5_cluster class
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$rel_MAD=data1$MAD/data1$read_median_length
data1$group3="ex5_cluster-based"
data1=data1[,c("ex5cluster_class","rel_MAD","polyArate","group3","n_read")]
colnames(data1)[1]="group"
data1$polyA="non-poly(A)"
data1$polyA[which(data1$polyArate>0.8)]="poly(A)"
data1$polyArate_bin=round(data1$polyArate*10)/10

write.table(data1,gzfile(paste0(path_fig4_data,"ex5_cluster.n.transcript.base.length.relative.MAD.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig4/data
# for -> Fig. ex6j
