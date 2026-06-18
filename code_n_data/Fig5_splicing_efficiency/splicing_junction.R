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
path_fig5_data=paste0(primary_folder,"fig5/data/")
SJ_path=paste0(primary_folder,"code_n_data/Fig5_splicing_efficiency/SJ/")
spliceAI_path=paste0(primary_folder,"code_n_data/Fig5_splicing_efficiency/spliceAI_score_output/")

#============================================================================================
#splice junctions from external data

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
system("zcat Neuron_THP1.S3.trnscpt.ref.bed.bgz | bedparse introns | bed12ToBed6 -i stdin | gzip > Neuron_THP1.S3.trnscpt.ref.intron.bed6.bed.gz")

GENCODE_SJ=read.delim("Neuron_THP1.S3.trnscpt.ref.intron.bed6.bed.gz", header=F, stringsAsFactors = F)
ref_intronID=unique(paste0(GENCODE_SJ$V1,"_",GENCODE_SJ$V2,"_",GENCODE_SJ$V3,"_",GENCODE_SJ$V6))


#====
short=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/valeria_short_read_junction/neuron_series.SJ.out.unique1.gencode_back.tab", header=F, stringsAsFactors = F, check.names = F)
short$V4[which(short$V4 == 1)]="+"
short$V4[which(short$V4 == 2)]="-"
short$V2=short$V2-1
short$site_label="0_noncanonical"
short$site_label[which(short$V5 == 1 & short$V4=="+")]="1_GT/AG"
short$site_label[which(short$V5 == 2 & short$V4=="-")]="1_GT/AG"
short$site_label[which(short$V5 == 3 & short$V4=="+")]="3_GC/AG"
short$site_label[which(short$V5 == 4 & short$V4=="-")]="3_GC/AG"
short$site_label[which(short$V5 == 5 & short$V4=="+")]="5_AT/AC"
short$site_label[which(short$V5 == 6 & short$V4=="-")]="5_AT/AC"
short$group="Short_read"
short$ID=paste0(short$V1,"_",short$V2,"_",short$V3,"_",short$V4)
length(unique(short$ID[which(short$V6 == 1)]))

ref_intronID1=setdiff(ref_intronID,unique(short$ID[which(short$V6 == 1)]))
GENCODE_SJ=data.frame(cbind("","GENCODEv39",ref_intronID1))
colnames(GENCODE_SJ)=c("site_label","group","ID")

GENCODE_SJ2=short[which(short$V6 ==1), c(10:12)]
GENCODE_SJ2$group="GENCODEv39"
GENCODE_SJ=rbind(GENCODE_SJ,GENCODE_SJ2)

short=short[which(!is.na(short$V7)),]

#short$group[which(short$ID %in% ref_intronID)]="GENCODEv39"

length(unique(short$ID))

intro=fread("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/read_bed/junction/intropolis.v1.hg19_with_liftover_to_hg38.tsv.gz", header=F, select=c(9:11,5,6,12,8))
intro$V10=intro$V10-1
options(scipen=999)
intro$site_label="0_noncanonical"
intro$site_label[which(intro$V5 == "GT" & intro$V6 == "AG")]="1_GT/AG"
intro$site_label[which(intro$V5 == "GC" & intro$V6 == "AG")]="3_GC/AG"
intro$site_label[which(intro$V5 == "AT" & intro$V6 == "AC")]="5_AT/AC"
intro$ID=paste0(intro$V9,"_",intro$V10,"_",intro$V11,"_",intro$V12)
intro%>%group_by(site_label)%>%dplyr::summarise(count=n())
intro$group="Intropolis"


support_intro=rbind(GENCODE_SJ,short[,c(10:12)],intro[,c(8,10,9)])
support_intro=unique(support_intro)

write.table(support_intro,gzfile(paste0(SJ_path,"all_support_intro.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

support_intro1=support_intro%>%group_by(ID, site_label)%>%dplyr::summarise(Supported = list(group), .groups = "drop")
support_sets1=support_intro1%>%group_by(Supported,site_label)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
support_sets2=support_sets1%>%group_by(Supported)%>%dplyr::summarise(count=sum(count))

saveRDS(support_sets1, paste0(path_fig5_data,"all_support_intro.summary.RDS"))
#stored in [primary_folder]/fig5/data

rm(intro)
rm(short)
rm(support_intro)
#===============================================================================
#calculate splicing efficiency
#need to use the read files to count the frequency of donor and acceptor sites

#convert bed12 to intron bed6 and exon bed6

system("for file in *.bed12.bed.gz; do zcat \"$file\" | bedparse introns | bed12ToBed6 -i stdin | gzip > \"${file%.bed12.bed.gz}.intron.bed6.bed.gz\"; done")
system("for file in *.bed12.bed.gz; do zcat \"$file\" | bed12ToBed6 -i stdin | gzip > \"${file%.bed12.bed.gz}.exon.bed6.bed.gz\"; done")


path1="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/"
#path2="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/junction/"
files=list.files(path=path1, pattern="intron.bed6.bed.gz")
intronbed6=fread(paste0(path1,files[1]), header=F, stringsAsFactors = F)
##set junction as 2 basepair##
options(scipen=999)
intronbed6$intronID=paste0(intronbed6$V1,"_",intronbed6$V2,"_",intronbed6$V3,"_",intronbed6$V6)
intronbed6$V2a=intronbed6$V2+1
intronbed6$V3a=intronbed6$V3-1
intronbed6$V2=intronbed6$V2-1
intronbed6$V3=intronbed6$V3+1
SS3p=intronbed6[which(intronbed6$V6 == "+"),c(1,9,3,4,5,6)]
SS5p=intronbed6[which(intronbed6$V6 == "+"),c(1,2,8,4,5,6)]
SS3n=intronbed6[which(intronbed6$V6 == "-"),c(1,2,8,4,5,6)]
SS5n=intronbed6[which(intronbed6$V6 == "-"),c(1,9,3,4,5,6)]
colnames(SS3n)=colnames(SS3p)
colnames(SS5n)=colnames(SS5p)
SS3p=rbind(SS3p,SS3n)
SS5p=rbind(SS5p,SS5n)
SS3p=SS3p%>%group_by(V1,V3a,V3,V6)%>%dplyr::summarise(count=n())
SS5p=SS5p%>%group_by(V1,V2,V2a,V6)%>%dplyr::summarise(count=n())
label=unique(intronbed6[,c(1:3,6:9)])
intron=intronbed6[,c(7,5)]
label$small=paste0(label$V1,"_",label$V2,"_",label$V2a,"_",label$V6)
label$big=paste0(label$V1,"_",label$V3a,"_",label$V3,"_",label$V6)
labelp=unique(label[which(label$V6=="+"),c(5,8,9)])
labeln=unique(label[which(label$V6=="-"),c(5,9,8)])
colnames(labelp)[c(2,3)]=c("SS5J_ID","SS3J_ID")
colnames(labeln)=colnames(labelp)
label=rbind(labelp,labeln)
intron=intron%>%group_by(intronID)%>%dplyr::summarise(V5=n())
options(scipen=999)
for (i in 2: length(files)){
  intronbed6=fread(paste0(path1,files[i]), header=F, stringsAsFactors = F)
  intronbed6$intronID=paste0(intronbed6$V1,"_",intronbed6$V2,"_",intronbed6$V3,"_",intronbed6$V6)
  intronbed6$V2a=intronbed6$V2+1
  intronbed6$V3a=intronbed6$V3-1
  intronbed6$V2=intronbed6$V2-1
  intronbed6$V3=intronbed6$V3+1
  aSS3p=intronbed6[which(intronbed6$V6 == "+"),c(1,9,3,4,5,6)]
  aSS5p=intronbed6[which(intronbed6$V6 == "+"),c(1,2,8,4,5,6)]
  aSS3n=intronbed6[which(intronbed6$V6 == "-"),c(1,2,8,4,5,6)]
  aSS5n=intronbed6[which(intronbed6$V6 == "-"),c(1,9,3,4,5,6)]
  colnames(aSS3n)=colnames(aSS3p)
  colnames(aSS5n)=colnames(aSS5p)
  aSS3p=rbind(aSS3p,aSS3n)
  aSS5p=rbind(aSS5p,aSS5n)
  aSS3p=aSS3p%>%group_by(V1,V3a,V3,V6)%>%dplyr::summarise(count=n())
  aSS5p=aSS5p%>%group_by(V1,V2,V2a,V6)%>%dplyr::summarise(count=n())
  alabel=unique(intronbed6[,c(1:3,6:9)])
  aintron=intronbed6[,c(7,5)]
  alabel$small=paste0(alabel$V1,"_",alabel$V2,"_",alabel$V2a,"_",alabel$V6)
  alabel$big=paste0(alabel$V1,"_",alabel$V3a,"_",alabel$V3,"_",alabel$V6)
  alabelp=unique(alabel[which(alabel$V6=="+"),c(5,8,9)])
  alabeln=unique(alabel[which(alabel$V6=="-"),c(5,9,8)])
  colnames(alabelp)[c(2,3)]=c("SS5J_ID","SS3J_ID")
  colnames(alabeln)=colnames(alabelp)
  alabel=rbind(alabelp,alabeln)
  SS3p=rbind(SS3p,aSS3p)
  SS5p=rbind(SS5p,aSS5p)
  label=rbind(label,alabel)
  aintron=aintron%>%group_by(intronID)%>%dplyr::summarise(V5=n())
  intron=rbind(intron, aintron)}

options(scipen=999)
SS3p=SS3p%>%group_by(V1,V3a,V3,V6)%>%dplyr::summarise(count=sum(count))
SS5p=SS5p%>%group_by(V1,V2,V2a,V6)%>%dplyr::summarise(count=sum(count))
SS3p=SS3p[order(SS3p$V1, SS3p$V3a),]
SS5p=SS5p[order(SS5p$V1, SS5p$V2),]
label=unique(label)
label=label[order(label$intronID),]
intron=intron%>%group_by(intronID)%>%dplyr::summarise(V5=sum(V5))
intron$V1=sapply(strsplit(intron$intronID,"_"),"[",1)
intron$V2=sapply(strsplit(intron$intronID,"_"),"[",2)
intron$V3=sapply(strsplit(intron$intronID,"_"),"[",3)
intron$V6=sapply(strsplit(intron$intronID,"_"),"[",4)
intron=intron[-grep("random", intron$intronID),]
intron=intron[-grep("chrUn", intron$intronID),]
intron=intron[order(intron$V1,as.numeric(intron$V2)),]

write.table(SS3p[,c(1,2,3,5,5,4)],gzfile(paste0(SJ_path,"all_ontCAGE.SS3J.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS5p[,c(1,2,3,5,5,4)],gzfile(paste0(SJ_path,"all_ontCAGE.SS5J.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(label,gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(intron[,c(3,4,5,1,2,6)], gzfile(paste0(SJ_path,"all_ontCAGE.intron.for.zenbu.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
label=read.delim(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#======================================================
#bash
#bedtools intersect
setwd("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/")
system("for file in *exon.bed6.bed.gz; do bedtools intersect -a /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/junction/all_ontCAGE.SS3J.count.bed.gz -b \"$file\" -c -s -f 1 -wa | gzip > \"${file%.exon.bed6.bed.gz}.SS3J.count.span.bed.gz\" ; done")
system("for file in *exon.bed6.bed.gz; do bedtools intersect -a /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/junction/all_ontCAGE.SS5J.count.bed.gz -b \"$file\" -c -s -f 1 -wa | gzip > \"${file%.exon.bed6.bed.gz}.SS5J.count.span.bed.gz\" ; done")
#======================================================

path1="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/"
#path2="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/read_bed_dorado/all_read_TC/junction/"
files=list.files(path=path1, pattern="SS3J.count.span.bed")
eff=fread(paste0(path1,files[1]), header=F, stringsAsFactors = F)
for (i in 2: length(files)){
  eff1=fread(paste0(path1,files[i]), header=F, stringsAsFactors = F)
  eff=left_join(eff,eff1[,c(1,2,3,6,7)], by=c("V1","V2","V3","V6"),copy=F)}
eff$V5=rowSums(eff[,c(7:39)])
eff=eff[,c(1:6)]
colnames(eff)[c(4,5)]=c("splice","span")
eff$eff=eff$splice/(eff$span+eff$splice)
eff$ID=paste0(eff$V1,"_",eff$V2,"_",eff$V3,"_",eff$V6)
length(which(eff$eff == 1))/nrow(eff) #0.3005347
write.table(eff,paste0(SJ_path,"all_ontCAGE.SS3J.splice.span.eff.tsv"),col.names=T, row.names=F, sep="\t", quote=F)

files=list.files(path=path1, pattern="SS5J.count.span.bed")
eff=fread(paste0(path1,files[1]), header=F, stringsAsFactors = F)
for (i in 2: length(files)){
  eff1=fread(paste0(path1,files[i]), header=F, stringsAsFactors = F)
  eff=left_join(eff,eff1[,c(1,2,3,6,7)], by=c("V1","V2","V3","V6"),copy=F)}
eff$V5=rowSums(eff[,c(7:39)])
eff=eff[,c(1:6)]
colnames(eff)[c(4,5)]=c("splice","span")
eff$eff=eff$splice/(eff$span+eff$splice)
eff$ID=paste0(eff$V1,"_",eff$V2,"_",eff$V3,"_",eff$V6)
length(which(eff$eff == 1))/nrow(eff) #0.2182673
write.table(eff,gzfile(paste0(SJ_path,"all_ontCAGE.SS5J.splice.span.eff.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)


#=================
junction_info_TC=read.delim(paste0(primary_folder,"/code_n_data/SALA/Neuron_THP1_full/transcript/log/Neuron_THP1.S3.junct.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
junction_info_TC$intron_ID=paste0(junction_info_TC$chrom,"_",junction_info_TC$intron_start,"_",junction_info_TC$intron_end,"_",junction_info_TC$strand)

junction_info_TC$V2a=junction_info_TC$intron_start+1
junction_info_TC$V3a=junction_info_TC$intron_end-1
junction_info_TC$V2=junction_info_TC$intron_start-1
junction_info_TC$V3=junction_info_TC$intron_end+1
SS3p=junction_info_TC[which(junction_info_TC$strand == "+"),c(4,12,14,10,9,7)]
SS5p=junction_info_TC[which(junction_info_TC$strand == "+"),c(4,13,11,10,9,7)]
SS3n=junction_info_TC[which(junction_info_TC$strand == "-"),c(4,13,11,10,9,7)]
SS5n=junction_info_TC[which(junction_info_TC$strand == "-"),c(4,12,14,10,9,7)]
colnames(SS3n)=colnames(SS3p)
colnames(SS5n)=colnames(SS5p)
SS3p=rbind(SS3p,SS3n)
SS5p=rbind(SS5p,SS5n)

label=junction_info_TC[,c(4,13,14,7,10,11,12)]
label$small=paste0(label$chrom,"_",label$V2,"_",label$V2a,"_",label$strand)
label$big=paste0(label$chrom,"_",label$V3a,"_",label$V3,"_",label$strand)
labelp=unique(label[which(label$strand=="+"),c(5,8,9)])
labeln=unique(label[which(label$strand=="-"),c(5,9,8)])
colnames(labelp)[c(2,3)]=c("SS5J_ID","SS3J_ID")
colnames(labeln)=colnames(labelp)
label=rbind(labelp,labeln)

options(scipen=999)
SS3p=SS3p[order(SS3p$chrom, SS3p$V3a),]
SS5p=SS5p[order(SS5p$chrom, SS5p$V2),]
write.table(SS3p,gzfile(paste0(SJ_path,"all_ontCAGE.SS3J.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS5p,gzfile(paste0(SJ_path,"all_ontCAGE.SS5J.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(label,gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

junction_info_TC=left_join(junction_info_TC[,c(1:10)], label, by="intron_ID")
junction_info_TC=left_join(junction_info_TC,junction_info[,c(1,8:13)],by=c("intron_ID"="junct_ID"),copy=F)
junction_info_TC$site_label="0_noncanonical"
junction_info_TC$site_label[which(junction_info_TC$splicing_site == "GT-AG")]="1_GT/AG"
junction_info_TC$site_label[which(junction_info_TC$splicing_site ==  "GC-AG")]="3_GC/AG"
junction_info_TC$site_label[which(junction_info_TC$splicing_site ==  "AT-AC")]="5_AT/AC"
junction_info_TC$source="detected"
junction_info_TC$source[which(junction_info_TC$qry_count == 0)]="undetected"

write.table(junction_info_TC, gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
junction_info_TC=read.delim(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz") , header=T, stringsAsFactors = F, check.names = F)

SS5J=read.delim(paste0(SJ_path,"all_ontCAGE.SS5J.splice.span.eff.tsv"), header=T, stringsAsFactors = F, check.names = F)
SS3J=read.delim(paste0(SJ_path,"all_ontCAGE.SS3J.splice.span.eff.tsv"), header=T, stringsAsFactors = F, check.names = F)
junction_info_TC=left_join(junction_info_TC, SS5J[,c(8,4,5,7)], by=c("SS5J_ID"="ID"),copy=F)
junction_info_TC=left_join(junction_info_TC, SS3J[,c(8,4,5,7)], by=c("SS3J_ID"="ID"),copy=F, suffix=c("_donor","_acceptor"))
colnames(intron)[2]="intron_count"
junction_info_TC=left_join(junction_info_TC, intron[,c(1,2)], by=c("intron_ID"="intronID"), copy=F)
junction_info_TC$eff_junction=junction_info_TC$intron_count/((junction_info_TC$span_donor+junction_info_TC$span_acceptor)/2+junction_info_TC$intron_count)
write.table(junction_info_TC, gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/code_n_data/Fig5_splicing_efficiency/SJ
# for -> Fig. 5c,d,e, supplementary table S12

#===========================================


#spliceAI
#===============================================================================
###prepare SS5 for getfasta for spliceAI###

SS3p=read.delim(paste0(SJ_path,"all_ontCAGE.SS3J.count.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
SS5p=read.delim(paste0(SJ_path,"all_ontCAGE.SS5J.count.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
SS3_20=SS3p
SS5_20=SS5p
SS3_20$V4=paste0(SS3_20$V1,"_",SS3_20$V2,"_",SS3_20$V3,"_",SS3_20$V6)
SS5_20$V4=paste0(SS5_20$V1,"_",SS5_20$V2,"_",SS5_20$V3,"_",SS5_20$V6)
SS3_20$V2[which(SS3_20$V6 == "+")]=SS3_20$V2[which(SS3_20$V6 == "+")]-20
SS3_20$V2[which(SS3_20$V6 == "-")]=SS3_20$V2[which(SS3_20$V6 == "-")]-19
SS3_20$V3[which(SS3_20$V6 == "+")]=SS3_20$V3[which(SS3_20$V6 == "+")]+19
SS3_20$V3[which(SS3_20$V6 == "-")]=SS3_20$V3[which(SS3_20$V6 == "-")]+20
SS3_30=SS3_20
SS3_30$V2=SS3_30$V2-10
SS3_30$V3=SS3_30$V3+10
SS3_40=SS3_20
SS3_40$V2=SS3_40$V2-20
SS3_40$V3=SS3_40$V3+20
write.table(SS3_20,gzfile(paste0(SJ_path,"all_ontCAGE.SS3J41.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS3_30,gzfile(paste0(SJ_path,"all_ontCAGE.SS3J61.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS3_40,gzfile(paste0(SJ_path,"all_ontCAGE.SS3J81.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
SS5_20$V2[which(SS5_20$V6 == "+")]=SS5_20$V2[which(SS5_20$V6 == "+")]-19
SS5_20$V2[which(SS5_20$V6 == "-")]=SS5_20$V2[which(SS5_20$V6 == "-")]-20
SS5_20$V3[which(SS5_20$V6 == "+")]=SS5_20$V3[which(SS5_20$V6 == "+")]+20
SS5_20$V3[which(SS5_20$V6 == "-")]=SS5_20$V3[which(SS5_20$V6 == "-")]+19
SS5_30=SS5_20
SS5_30$V2=SS5_30$V2-10
SS5_30$V3=SS5_30$V3+10
SS5_40=SS5_20
SS5_40$V2=SS5_40$V2-20
SS5_40$V3=SS5_40$V3+20
write.table(SS5_20,gzfile(paste0(SJ_path,"all_ontCAGE.SS5J41.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS5_30,gzfile(paste0(SJ_path,"all_ontCAGE.SS5J61.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(SS5_40,gzfile(paste0(SJ_path,"all_ontCAGE.SS5J81.count.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#getfasta for spliceAI
setwd (SJ_path)
system("for file in *1.count.bed.gz; do bedtools getfasta -s -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed \"$file\" |cut -f4,7 > \"./spliceAI/${file%.count.bed.gz}.tsv\" ; done")
#input of spliceAI (fasta in tsv format) located in [primary_folder]/code_n_data/Fig5_splicing_efficiency/spliceAI_score_input

#run spliceAI
system("nohup python3 junction_run.py 1 all_ontCAGE.SS3J81.tsv > run1.log 2>&1 &")
system("nohup python3 junction_run.py 1 all_ontCAGE.SS5J81.tsv > run1.log 2>&1 &")
#junction_run.py located in [primary_folder]/code_n_data/Fig5_splicing_efficiency
#output of spliceAI was zipped and placed in [primary_folder]/code_n_data/Fig5_splicing_efficiency/spliceAI_score_output

#===============================================================================
AI3=read.delim(paste0(spliceAI_path,"all_ontCAGE.SS3J81_5SpliceAI.tsv.gz"), header=T, stringsAsFactors = F)
AI3a=separate_rows(AI3[,c(1,4,5)], acceptor,donor, sep=",")
AI3a=AI3a%>%group_by(ID)%>%dplyr::mutate(position=1:n())
AI3a$position = AI3a$position-41
AI3b=AI3a[which(AI3a$position == 1),]
AI3a1=AI3a%>%group_by(position)%>%dplyr::summarise(mean_acceptor_score=mean(as.numeric(acceptor)),mean_donor_score=mean(as.numeric(donor)))
AI3a2=melt(AI3a1, id=1)
AI3a2$group1="Acceptor_site"
AI5=read.delim(paste0(spliceAI_path,"all_ontCAGE.SS5J81_5SpliceAI.tsv.gz"), header=T, stringsAsFactors = F)
AI5a=separate_rows(AI5[,c(1,4,5)], acceptor,donor, sep=",")
AI5a=AI5a%>%group_by(ID)%>%dplyr::mutate(position=1:n())
AI5a$position = AI5a$position-41
AI5b=AI5a[which(AI5a$position == (-1)),]
AI5a1=AI5a%>%group_by(position)%>%dplyr::summarise(mean_acceptor_score=mean(as.numeric(acceptor)),mean_donor_score=mean(as.numeric(donor)))
AI5a2=melt(AI5a1, id=1)
AI5a2$group1="Donor_site"
final.plot=rbind(AI3a2,AI5a2)
write.table(final.plot, gzfile(paste0(spliceAI_path,"spliceAIdonor.acceptor.n81.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=========================
junction_info_TC=read.delim(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
junction_info_TC=junction_info_TC[,c(1:28, 32:37)]
colnames(AI5b)[3]="spliceAI_donor"
colnames(AI3b)[2]="spliceAI_acceptor"
junction_info_TC=left_join(junction_info_TC,AI5b[,c(1,3)],by=c("SS5J_ID"="ID"),copy=F)
junction_info_TC=left_join(junction_info_TC,AI3b[,c(1,2)],by=c("SS3J_ID"="ID"),copy=F)
junction_info_TC$spliceAI_donor=as.numeric(junction_info_TC$spliceAI_donor)
junction_info_TC$spliceAI_acceptor=as.numeric(junction_info_TC$spliceAI_acceptor)
junction_info_TC$spliceAI_junction= apply(junction_info_TC[29:30], MARGIN =  1, FUN = min, na.rm = T)
junction_info_TC=junction_info_TC[,c(1:28,35:37,29:34)]
write.table(junction_info_TC, gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#===============================================================================
junction_info_TC=read.delim(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
junction_info_TC$GENCODEv39 = "not_supported"
junction_info_TC$GENCODEv39[which(junction_info_TC$intron_ID %in% unique(support_intro$ID[which(support_intro$group == "GENCODEv39")]))] = "supported"
junction_info_TC$short_read = "not_supported"
junction_info_TC$short_read[which(junction_info_TC$intron_ID %in% unique(support_intro$ID[which(support_intro$group == "Short_read")]))] = "supported"
junction_info_TC$intropolis = "not_supported"
junction_info_TC$intropolis[which(junction_info_TC$intron_ID %in% unique(support_intro$ID[which(support_intro$group == "Intropolis")]))] = "supported"
junction_info_TC%>%group_by(T5,GENCODEv39,short_read,intropolis)%>%dplyr::summarise(count=n())

support_intro2=support_intro%>%arrange(ID,group)%>%group_by(ID)%>%dplyr::summarise(Supported=paste(unique(group),collapse=";"))

junction_info_TC=left_join(junction_info_TC, support_intro2, by=c("intron_ID"="ID"),copy=F)
junction_info_TC$Supported[which(is.na(as.character(junction_info_TC$Supported)))]="not_supported"
junction_info_TC[which(junction_info_TC$qry_count>0),]%>%group_by(Supported,TC_introduced)%>%dplyr::summarise(count=n())

write.table(junction_info_TC, gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
#stored in [primary_folder]/code_n_data/Fig5_splicing_efficiency/SJ
#for -> Fig. Ext7e,f, supplementary table S12
rm(support_intro)

#===============================================================================
#splicing junctions before and after transcriptClean
junction_info_TC=read.delim(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
junction_info_raw=fread(paste0(primary_folder,"/code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/junction_extractor/pool/Neuron_THP1.full.junct.info.tsv.gz"), header=T, stringsAsFactors = F)

junction_info_raw$TC_removed="Yes"
junction_info_raw$TC_removed[which(junction_info_raw$junct_ID %in% junction_info_TC$intron_ID[which(junction_info_TC$qry_count>0)])]="No"
junction_info_raw%>%group_by(TC_removed, canonical)%>%dplyr::summarise(locus=n(), count=sum(total_count),avg_score=median(avg_score_per_nt))%>%dplyr::mutate(percent_locus=locus/sum(locus),percent_count=count/sum(count))
junction_info_raw%>%group_by(TC_removed)%>%dplyr::summarise(locus=n(), count=sum(total_count),avg_score=median(avg_score_per_nt))%>%dplyr::mutate(percent_locus=locus/sum(locus),percent_count=count/sum(count))
#TC_removed   locus     count avg_score percent_locus percent_count
#<chr>        <int>     <int>     <dbl>         <dbl>         <dbl>
#1 No         1654272 844423685      19.2        0.976      1.00     
#2 Yes          40024     67871      16.7        0.0236     0.0000804

junction_info_TC$TC_introduced="Yes"
junction_info_TC$TC_introduced[which(junction_info_TC$intron_ID %in% junction_info_raw$junct_ID)]="No"
junction_info_TC[which(junction_info_TC$qry_count>0),]%>%group_by(TC_introduced, canonical)%>%dplyr::summarise(locus=n(), count=sum(qry_count))%>%dplyr::mutate(percent_locus=locus/sum(locus),percent_count=count/sum(count))
junction_info_TC[which(junction_info_TC$qry_count>0),]%>%group_by(TC_introduced)%>%dplyr::summarise(locus=n(), count=sum(qry_count))%>%dplyr::mutate(percent_locus=locus/sum(locus),percent_count=count/sum(count))
#TC_introduced   locus     count percent_locus percent_count
#<chr>           <int>     <int>         <dbl>         <dbl>
#1 No            1654272 843643343       0.995       1.00     
#2 Yes              8107     15706       0.00488     0.0000186

1654272/(1654272+40024+8107) #97.17%


#===============================================================================
#define SJ that were found in finalized transcriptome
#bash
#get intron from SALA Final

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file"))
system("zcat table5.final.fulllength_detected.alone.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > table5.final.fulllength_detected.alone.intron.bed6.bed.gz")

#===============================================================================
t5_intron=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.fulllength_detected.alone.intron.bed6.bed.gz"), header=F)
t5_intronID=unique(paste0(t5_intron$V1,"_",t5_intron$V2,"_",t5_intron$V3,"_",t5_intron$V6))

junction_info_TC$T5 = "No"
junction_info_TC$T5[which(junction_info_TC$intron_ID %in% t5_intronID)] = "Yes"
junction_info_TC_t5=junction_info_TC[which(junction_info_TC$T5 == "Yes"),]
write.table(junction_info_TC, gzfile(paste0(SJ_path,"all_ontCAGE.intronID_SS5J_SS3J.label.eff.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
# stored in [primary_folder]/code_n_data/Fig5_splicing_efficiency/SJ

junction_info_TC1=junction_info_TC%>%group_by(GENCODEv39, site_label, Supported)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
junction_info_TC1$Supported[which(junction_info_TC1$GENCODEv39 == "supported")]="GENCODE"
junction_info_TC1$Supported[which(junction_info_TC1$Supported == "Intropolis;Short_read")]="Both"

junction_info_TC2=junction_info_TC_t5%>%group_by(GENCODEv39, site_label, Supported)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
junction_info_TC2$Supported[which(junction_info_TC2$GENCODEv39 == "supported")]="GENCODE"
junction_info_TC2$Supported[which(junction_info_TC2$Supported == "Intropolis;Short_read")]="Both"
junction_info_TC1$group="All read"
junction_info_TC2$group="Final transcroptome"
junction_info_TC1=rbind(junction_info_TC1,junction_info_TC2)
write.table(junction_info_TC1,gzfile(paste0(path_fig5_data,"splicing.junction.Neuron_n_THP1.from.read.n.table5.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
# stored in [primary_folder]/fig5/data
# summary for -> fig5a
#================================================

junction_info_TC$Support="Supported"
junction_info_TC$Support[which(junction_info_TC$Supported == "not_supported")]="Novel"

junction_info_TC$final_confidence="No"
junction_info_TC$final_confidence[which(junction_info_TC$Support == "Supported")]="Yes"
junction_info_TC$final_confidence[which(junction_info_TC$max_score_per_nt >= 10 & junction_info_TC$total_count >=3)]="Yes"

junction_info_TC_t5=junction_info_TC[which(junction_info_TC$T5 == "Yes"),]
junction_info_TC_t5%>%group_by(final_confidence,Support)%>%dplyr::summarise(count=n())%>%dplyr::mutate(count/sum(count))

junction_info_TC_t5$canonical[which(junction_info_TC_t5$canonical=="Y")]="Canonical"
junction_info_TC_t5$canonical[which(junction_info_TC_t5$canonical=="N")]="Non-canonical"
write.table(junction_info_TC_t5, gzfile(paste0(path_fig5_data,"All_splicing.junction.Neuron_n_THP1.from.read.n.table5.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
#stored in [primary_folder]/fig5/data
# for -> Fig. Ext7b

#===============================================================================
#geneID to ex5_cluster link
table5b=unique(table5[,c(62,86,114)])%>%group_by(n5_string)%>%dplyr::summarise(T4_gene_ID=paste(T4_gene_ID, collapse=";"), gene_group=paste(unique(gene_group),collapse=";"), count=n())
data1=read.delim(paste0(primary_folder,"fig4/data/features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
length(which(table5b$count==1))/nrow(table5b) #98.9%
table5b=table5b[which(table5b$count==1),]
table5c=table5b%>%group_by(T4_gene_ID)%>%dplyr::summarise(n5_string=paste(unique(n5_string),collapse=";"), count=n())
length(which(table5c$count==1))/nrow(table5c)  #79%
table5c=table5c[which(table5c$count==1),]
table5c=left_join(table5c, data1[,c(1,9,21:24,84,85)], by="n5_string", copy=F)
table5c=table5c[which(!is.na(table5c$CpGTATA)),]
write.table(table5c,gzfile(paste0(path_fig5_data,"geneID_to_ex5cluster.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
###separate to RNA class
t5_intron=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.fulllength_detected.alone.intron.bed6.bed.gz"), header=F)
t5_intron$intron_ID=paste0(t5_intron$V1,"_",t5_intron$V2,"_",t5_intron$V3,"_",t5_intron$V6)
t5_intron1=t5_intron[which(t5_intron$V4 %in% table5$model_ID[which(table5$group %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))]),]
t5_intron1=left_join(t5_intron1, table5[,c("model_ID","group")], by=c("V4"="model_ID"),copy=F)
t5_intron1=left_join(t5_intron1, junction_info_TC,by="intron_ID", copy=F)
write.table(t5_intron1, gzfile(paste0(path_fig5_data,"splicing.junction.RNAclass.plot.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
#stored in [primary_folder]/fig5/data
# for -> Fig. Ext7d

#===============================================================================

t5_intron1=read.delim(paste0(path_fig5_data,"splicing.junction.RNAclass.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t5_intron_t=t5_intron1%>%group_by(group,V4)%>%dplyr::summarise(All_SJ_support_confidence=paste(unique(final_confidence),collapse=";"))
t5_intron_t$All_SJ_support_confidence[which(t5_intron_t$All_SJ_support_confidence == "Yes")]="All SJ supported"
t5_intron_t$All_SJ_support_confidence[which(t5_intron_t$All_SJ_support_confidence == "No")]="Partial/not supported"
t5_intron_t$All_SJ_support_confidence[grep(";",t5_intron_t$All_SJ_support_confidence)]="Partial/not supported"
t5_intron_t%>%group_by(All_SJ_support_confidence)%>%dplyr::summarise(count=n())
table5a=table5[which(table5$group %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),]
table5a=left_join(table5a,t5_intron_t,by=c("model_ID" = "V4"),copy=F)
table5a$Exon="Spliced"
table5a$Exon[which(table5a$n_exon == 1)]="Un-spliced"
colnames(table5a)[102]="group"
write.table(table5a[,c(1,102,122:123)], gzfile(paste0(path_fig5_data,"splicing_support_group.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig5/data
# for -> Fig. Ext7c



#gene-base splicing efficiency
#===============================================================================
#take the gene region of above and locate all the junction inside
#grouping at gene-level
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5$gene_group=table5$T4_Novel_geneClass
table5$gene_group[which(is.na(table5$gene_group) & table5$T4_Gencode_geneCalss2 == "protein_coding")]="mRNA"
table5$gene_group[which(is.na(table5$gene_group) & table5$T4_Gencode_geneCalss2 == "lncRNA")]="ncRNA"
table5$gene_group[which(is.na(table5$gene_group))]="others"
table5$gene_group[which(table5$gene_group=="lncRNA")]="ncRNA"
table5$gene_group=paste0(table5$T4_gene_promoter_type,"_",table5$gene_group)
table5$gene_group[grep("mRNA",table5$gene_group)]="mRNA"
table5$gene_group[grep("others",table5$gene_group)]="excluded"
table5$gene_group[grep("excluded",table5$gene_group)]="excluded"
table5$gene_group=gsub("nhancer-like","",table5$gene_group)
table5$gene_group=gsub("romoter-like","",table5$gene_group)
table5$gene_group=gsub("-alone","",table5$gene_group)
table5$gene_group=gsub("unclassed","other",table5$gene_group)
unique(table5[,c(86,114)])%>%group_by(gene_group)%>%dplyr::summarise(count=n())
write.table(table5, gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz")),col.names=T, row.names = F, sep="\t", quote=F)
#as in Table S4

options(scipen=999)
table5g=table5[which(table5$gene_group %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),c("T4_gene_ID","gene_group","n_exon")]%>%group_by(gene_group,T4_gene_ID)%>%dplyr::summarise(max_exon=max(n_exon))
table5t=table5[which(table5$T4_gene_ID %in% table5g$T4_gene_ID),]
t5_intron2=right_join(t5_intron,table5t[,c("model_ID","T4_gene_ID")], by=c("V4"="model_ID"),copy=F)
t5_intron2=unique(t5_intron2[,c(1:3,8,5:7)])
t5_intron3=t5_intron2[,c(4,7)]%>%group_by(T4_gene_ID)%>%dplyr::mutate(count=n())
t5_intron3=right_join(t5_intron3, table5g, by="T4_gene_ID",copy=F)
t5_intron3=t5_intron3[which(!is.na(t5_intron3$intron_ID)),] #contain gene with junctions only
t5_intron3=left_join(t5_intron3, junction_info_TC[,c(10,9,19:37)], by="intron_ID",copy=F)
intron_gene1=t5_intron3

#==============================
#add expression filter
filter1=which(intron_gene1$splice_donor+intron_gene1$span_donor >=1 & intron_gene1$splice_acceptor+intron_gene1$span_acceptor >=1)
filter2=which(intron_gene1$splice_donor+intron_gene1$span_donor >=3 & intron_gene1$splice_acceptor+intron_gene1$span_acceptor >=3)
filter3=which(intron_gene1$splice_donor+intron_gene1$span_donor >=5 & intron_gene1$splice_acceptor+intron_gene1$span_acceptor >=5)
intron_gene1a=intron_gene1[filter1,]
intron_gene1b=intron_gene1[filter2,]
intron_gene1c=intron_gene1[filter3,]
intron_gene1a$occurence=">=1"
intron_gene1b$occurence=">=3"
intron_gene1c$occurence=">=5"
intron_gene1=rbind(intron_gene1a,intron_gene1b,intron_gene1c)

intron_gene1=intron_gene1%>%group_by(occurence,T4_gene_ID)%>%dplyr::mutate(weight1=splice_acceptor/sum(splice_acceptor), weight2=splice_donor/sum(splice_donor), weight3=qry_count/sum(qry_count))
intron_gene2=intron_gene1%>%group_by(occurence,gene_group, T4_gene_ID)%>%dplyr::summarise(eff_acceptor=sum(eff_acceptor*weight1), spliceAI_acceptor=sum(as.numeric(spliceAI_acceptor)*weight1),
                                                                                  eff_donor=sum(eff_donor*weight2), spliceAI_donor=sum(as.numeric(spliceAI_donor)*weight2),
                                                                                  eff_junction=sum(eff_junction*weight3), spliceAI_junction=sum(as.numeric(spliceAI_junction)*weight3))

write.table(intron_gene2,gzfile(paste0(path_fig5_data,"gene_base_SJ_eff_AI.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig5/data
# for -> Fig. 5d&e

#===============================================================================
#update finalized transcriptome log table with min splicing efficiency to each transcript model

t5_intron=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.fulllength_detected.alone.intron.bed6.bed.gz"), header=F)
t5_intron$intron_ID=paste0(t5_intron$V1,"_",t5_intron$V2,"_",t5_intron$V3,"_",t5_intron$V6)

t5_intron=left_join(t5_intron, junction_info_TC[,c(10,28,31,39)], by="intron_ID", copy=F)
t5_intron_t=t5_intron%>%group_by(V4)%>%dplyr::summarise(min_sj_eff=min(eff_junction,na.rm = T), min_sj_AI=min(spliceAI_junction), All_SJ_support=paste(unique(final_confidence),collapse=";"))
t5_intron_t$All_SJ_support[grep(";",t5_intron_t$All_SJ_support)]="No"
t5_intron_t$min_sj_eff[which(t5_intron_t$min_sj_eff == Inf)]=NA
t5_intron_t$min_sj_AI[which(t5_intron_t$min_sj_AI == Inf)]=NA

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5=left_join(table5, t5_intron_t, by=c("model_ID"="V4"),copy=F)
table5$min_sj_eff[which(table5$n_exon ==1)]="single exon"
table5$min_sj_AI[which(table5$n_exon ==1)]="single exon"
table5$All_SJ_support[which(table5$n_exon ==1)]="single exon"
write.table(table5, gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#for -> supplementary table S4
#===============================================================================



