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
SE_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/SE_identification/")
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/")

#===============================================================================
# promoter-typing and super enhancer identification
#===============================================================================
# ROSE
#ROSE_main.py -g hg38 -i [significant H3K27Ac peaks ] -r [H3K27ac bam file] –c [IgG bam file] -o [output folder] -s 10000 -t 2500
#CUT&Tag analyses from mapping to rose please refer to [primary_folder]/code_n_data/Fig2_CRE_analysis/mapping_n_SE.txt
#BAM files of CUT&Tag please download from DDBJ

#==========================
#super enhancer by 27ac#
setwd(SE_path)
ips=read.delim(paste0(SE_path,"iPS_k27ac5/ips_27ac_narrowpeak_AllStitched.table.txt"), header=F, stringsAsFactors = F, check.names = F, skip=6)
nsc=read.delim(paste0(SE_path,"NSC_k27ac5/nsc_27ac_narrowpeak_AllStitched.table.txt"), header=F, stringsAsFactors = F, check.names = F, skip=6)
neuron=read.delim(paste0(SE_path,"NRN_k27ac5/nrn_27ac_narrowpeak_AllStitched.table.txt"), header=F, stringsAsFactors = F, check.names = F, skip=6)
ips$group="iPSC"
nsc$group="NSC"
neuron$group="Neuron"

ips=ips[order(ips$V7),]
ips$rank=1:nrow(ips)
ips$cutoff=16537.506
nsc=nsc[order(nsc$V7),]
nsc$rank=1:nrow(nsc)
nsc$cutoff=8480.589
neuron=neuron[order(neuron$V7),]
neuron$rank=1:nrow(neuron)
neuron$cutoff=13252.344
combine=rbind(ips, nsc, neuron)
colnames(combine)[7]="total_count"
write.table(combine,gzfile(paste0(path_fig2_data,"rose_SE_info.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data
#for -> figExt3.b & c


####SE bed file prepare####
options(scipen=999)
ips.27ac1=read.delim(paste0(SE_path,"iPS_k27ac5/ips_27ac_narrowpeak_SuperStitched.table.txt"), header=T, stringsAsFactors = F, check.names = F, skip=5)
write.table(ips.27ac1[,c(2,3,4,1)], paste0(SE_path,"iPS_k27ac5/iPS_27ac_SE.bed"), row.names=F, col.names=F, sep="\t", quote=F)
nsc.27ac1=read.delim(paste0(SE_path,"NSC_k27ac5/nsc_27ac_narrowpeak_SuperStitched.table.txt"), header=T, stringsAsFactors = F, check.names = F, skip=5)
write.table(nsc.27ac1[,c(2,3,4,1)], paste0(SE_path,"NSC_k27ac5/NSC_27ac_SE.bed"), row.names=F, col.names=F, sep="\t", quote=F)
nrn.27ac1=read.delim(paste0(SE_path,"NRN_k27ac5/nrn_27ac_narrowpeak_SuperStitched.table.txt"), header=T, stringsAsFactors = F, check.names = F, skip=5)
write.table(nrn.27ac1[,c(2,3,4,1)], paste0(SE_path,"NRN_k27ac5/NRN_27ac_SE.bed"), row.names=F, col.names=F, sep="\t", quote=F)

#promoter-typing & SE cluster
#======================================================================
#bash
#intersect with SCREEN cCRE and rose SE
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
system(paste0("bedtools closest -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",primary_folder,"code_n_data/n5_regions/GRCh38-ELS.all.enhancer.sort.bed.gz -D a | gzip > ontCAGE.Neuron_THP1.CRE.coord.e.bed.gz"))
system(paste0("bedtools closest -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",primary_folder,"code_n_data/n5_regions/GRCh38-PLS.all.promoter.sort.bed.gz -D a | gzip > ontCAGE.Neuron_THP1.CRE.coord.p.bed.gz"))
system(paste0("bedtools closest -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",primary_folder,"code_n_data/n5_regions/GRCh38-CTCF.sort.bed.gz -D a | gzip > ontCAGE.Neuron_THP1.CRE.coord.ctcf.bed.gz"))

system(paste0("bedtools intersect -c -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",SE_path,"iPS_k27ac5/iPS_27ac_SE.bed | gzip > ontCAGE.Neuron_THP1.CRE.coord.SE_ips27ac.bed.gz"))
system(paste0("bedtools intersect -c -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",SE_path,"NSC_k27ac5/NSC_27ac_SE.bed | gzip > ontCAGE.Neuron_THP1.CRE.coord.SE_nsc27ac.bed.gz"))
system(paste0("bedtools intersect -c -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",SE_path,"NRN_k27ac5/NRN_27ac_SE.bed | gzip > ontCAGE.Neuron_THP1.CRE.coord.SE_neur27ac.bed.gz"))
#======================================================================

CRE=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/log/ontCAGE.Neuron_THP1.CRE.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE.p=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.p.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
CRE.e=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.e.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
CRE.ctcf=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.ctcf.bed.gz"), header=F, stringsAsFactors = F, check.names = F)

CRE.se.ips=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.SE_ips27ac.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
CRE.se.nsc=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.SE_nsc27ac.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
CRE.se.neur=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.SE_neur27ac.bed.gz"), header=F, stringsAsFactors = F, check.names = F)

se.ips=unique(CRE.se.ips$V4[which(CRE.se.ips$V13>0)])
se.nsc=unique(CRE.se.nsc$V4[which(CRE.se.nsc$V13>0)])
se.neur=unique(CRE.se.neur$V4[which(CRE.se.neur$V13>0)])
se.inter=intersect(intersect(se.ips,se.nsc),se.neur)
se.all2=union(se.ips,union(se.nsc,se.neur))

CRE.e$SE2=0
CRE.e$SE2[which(CRE.e$V4 %in% se.all2)]=1
CRE.e$SE.ips=0
CRE.e$SE.ips[which(CRE.e$V4 %in% se.ips)]=1
CRE.e$SE.nsc=0
CRE.e$SE.nsc[which(CRE.e$V4 %in% se.nsc)]=1
CRE.e$SE.neuron=0
CRE.e$SE.neuron[which(CRE.e$V4 %in% se.neur)]=1

CREp=unique(CRE.p$V4[which(CRE.p$V19 ==0)])
CREe=unique(CRE.e$V4[which(CRE.e$V19 ==0)])
CREc=unique(CRE.ctcf$V4[which(CRE.ctcf$V19 ==0)])
CREseips=unique(CRE.e$V4[which(CRE.e$V19 ==0 & CRE.e$SE.ips == 1)])
CREsensc=unique(CRE.e$V4[which(CRE.e$V19 ==0 & CRE.e$SE.nsc == 1)])
CREseneuron=unique(CRE.e$V4[which(CRE.e$V19 ==0 & CRE.e$SE.neuron == 1)])

CRE$promoter=0
CRE$promoter[which(CRE$CREID %in% CREp)]=1
CRE$enhancer=0
CRE$enhancer[which(CRE$CREID %in% CREe)]=1
CRE$CTCF=0
CRE$CTCF[which(CRE$CREID %in% CREc)]=1
CRE$SE.ips=0
CRE$SE.ips[which(CRE$enhancer == 1 & CRE$CREID %in% CREseips)] =1
CRE$SE.nsc=0
CRE$SE.nsc[which(CRE$enhancer == 1 & CRE$CREID %in% CREsensc)] =1
CRE$SE.neuron=0
CRE$SE.neuron[which(CRE$enhancer == 1 & CRE$CREID %in% CREseneuron)] =1

CRE$promoter_type="unclassed"
CRE$promoter_type[which(CRE$CTCF ==1)]="CTCF-alone"
CRE$promoter_type[which(CRE$enhancer ==1)]="enhancer-like"
CRE$promoter_type[which(CRE$promoter ==1)]="promoter-like"
CRE$SE_source=NA
CRE$SE_source[which(CRE$SE.ips==1 & CRE$SE.nsc==0 & CRE$SE.neuron==0)]="iPSC"
CRE$SE_source[which(CRE$SE.ips==0 & CRE$SE.nsc==1 & CRE$SE.neuron==0)]="NSC"
CRE$SE_source[which(CRE$SE.ips==0 & CRE$SE.nsc==0 & CRE$SE.neuron==1)]="Neuron"
CRE$SE_source[which(CRE$SE.ips==1 & CRE$SE.nsc==1 & CRE$SE.neuron==0)]="iPSC & NSC"
CRE$SE_source[which(CRE$SE.ips==1 & CRE$SE.nsc==0 & CRE$SE.neuron==1)]="iPSC & Neuron"
CRE$SE_source[which(CRE$SE.ips==0 & CRE$SE.nsc==1 & CRE$SE.neuron==1)]="NSC & Neuron"
CRE$SE_source[which(CRE$SE.ips==1 & CRE$SE.nsc==1 & CRE$SE.neuron==1)]="iPSC & NSC & Neuron"
CRE$SE_source[which(CRE$promoter_type=="promoter-like" | CRE$promoter_type=="unclassed")]=NA 

CRE%>%group_by(promoter_type)%>%dplyr::summarise(count=n())
venn4=CRE%>%group_by(SE_source)%>%dplyr::summarise(count=n())
###sample-specific bi-directional enhancer###
CREdir=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/directionality/ontCAGE.Neuron_THP1/log/ontCAGE.Neuron_THP1.directionality.log.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREdir$fwd_rev_count=CREdir$fwd_count+CREdir$rev_count
CREdir$orientation[which(abs(CREdir$directionality)<=0.2)]="unidirectional"
CREdir$orientation[which(CREdir$orientation == "divergent")]="bidirectional"
CRE=left_join(CRE,CREdir[,c(1,4,5,7,8,13,10)],by="CREID",copy=F)
CRE%>%group_by(promoter_type,orientation)%>%dplyr::summarise(count=n())

#Andersson enhancer
#=====================bedtools closest=================================
#bash
#intersect with enhancer region defined by Andersson et al.
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
system("bedtools closest -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/Andersson.enhancer/hg38_robust_enhancers.sort.bed -D a> ontCAGE.Neuron_THP1.CRE.coord.andersson.robust.e.bed")
system("bedtools closest -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/Andersson.enhancer/hg38_permissive_enhancers.sort.bed -D a> ontCAGE.Neuron_THP1.CRE.coord.andersson.permissive.e.bed")
#======================================================================

re=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.andersson.robust.e.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
pe=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.andersson.permissive.e.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
re_e=unique(re$V4[which(re$V25 ==0)])
pe_e=unique(pe$V4[which(pe$V25 ==0)])
CRE$Andersson_robust=0
CRE$Andersson_robust[which(CRE$CREID %in% re_e)]=1
CRE$Andersson_permissive=0
CRE$Andersson_permissive[which(CRE$CREID %in% pe_e)]=1
CRE%>%group_by(promoter_type,Andersson_permissive)%>%dplyr::summarise(count=n())

write.table(CRE, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#==================================
# add ATAC support
# only consider tCRE from Neuron series, exclude THP-1 specific tCRE
count_qualified_read=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count/output/count_matrix/ontCAGE.Neuron_THP1.count.txt"), header=T, stringsAsFactors = F, check.names = F)
count_qualified_read1=count_qualified_read[,-grep("THP",colnames(count_qualified_read))]
count_qualified_read1=count_qualified_read1[rowSums(count_qualified_read1[2:7])>0,]
count_qualified_read1$iPS=0
count_qualified_read1$iPS[which(count_qualified_read1$ontCAGE.iPSC.rep1+count_qualified_read1$ontCAGE.iPSC.rep2 >0)]=1
count_qualified_read1$NSC=0
count_qualified_read1$NSC[which(count_qualified_read1$ontCAGE.NSC.rep1+count_qualified_read1$ontCAGE.NSC.rep2 >0)]=1
count_qualified_read1$NRN=0
count_qualified_read1$NRN[which(count_qualified_read1$ontCAGE.Neuron.rep1+count_qualified_read1$ontCAGE.Neuron.rep2 >0)]=1
  
tCRE_bed=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.bed.gz"), header=F, stringsAsFactors = F)
tCRE_bed1=right_join(tCRE_bed[,c(1:6)],count_qualified_read1[,c(1,8:10)], by=c("V4"="CREID"),copy=F)
write.table(tCRE_bed1[order(tCRE_bed1$V1,tCRE_bed1$V2),c(1:6)],gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/neuron_series_tCRE.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(tCRE_bed1[which(tCRE_bed1$iPS ==1),c(1:6)],gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/iPS_tCRE.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(tCRE_bed1[which(tCRE_bed1$NSC ==1),c(1:6)],gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/NSC_tCRE.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(tCRE_bed1[which(tCRE_bed1$NRN ==1),c(1:6)],gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/NRN_tCRE.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

path4=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/ATAC/")
aCRE_count=read.delim(paste0(path4,"intersect.iPS.bed"), header=F, stringsAsFactors = F)
aCRE_count1=read.delim(paste0(path4,"intersect.NSC.bed"), header=F, stringsAsFactors = F)
aCRE_count2=read.delim(paste0(path4,"intersect.NRN.bed"), header=F, stringsAsFactors = F)

write.table(aCRE_count[,c(1:6)], gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/aCRE.neuron_series.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(aCRE_count[which(aCRE_count$V7 !=0),c(1:6)], gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/aCRE.iPS.filtered.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(aCRE_count1[which(aCRE_count1$V7 !=0),c(1:6)], gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/aCRE.NSC.filtered.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(aCRE_count2[which(aCRE_count2$V7 !=0),c(1:6)], gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/aCRE.NRN.filtered.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#================================================================
# bash
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
system("bedtools intersect -c -a neuron_series_tCRE.bed.gz -b aCRE.neuron_series.bed.gz | gzip > neuron_series_tCRE.ATAC1.bed.gz")
system("bedtools intersect -c -a iPS_tCRE.bed.gz -b aCRE.iPS.filtered.bed.gz | gzip > iPS_tCRE.ATAC1.bed.gz")
system("bedtools intersect -c -a NSC_tCRE.bed.gz -b aCRE.NSC.filtered.bed.gz | gzip > NSC_tCRE.ATAC1.bed.gz")
system("bedtools intersect -c -a NRN_tCRE.bed.gz -b aCRE.NRN.filtered.bed.gz | gzip > NRN_tCRE.ATAC1.bed.gz")

#================================================================

atac1=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/neuron_series_tCRE.ATAC1.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
atac2=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/iPS_tCRE.ATAC1.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
atac3=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/NSC_tCRE.ATAC1.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
atac4=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/NRN_tCRE.ATAC1.bed.gz"), header=F, stringsAsFactors = F, check.names = F)

atac1=left_join(atac1,atac2[,c(4,7)],by="V4",copy=F)
atac1=left_join(atac1,atac3[,c(4,7)],by="V4",copy=F)
atac1=left_join(atac1,atac4[,c(4,7)],by="V4",copy=F)
colnames(atac1)[c(7:10)]=c("Neuron-series","iPS","NSC","Neuron")
atac1[is.na(atac1)]=0

atac2=reshape2::melt(atac1[,c(4,7:10)], id=1)
atac2$value[which(atac2$value!=0)]="ATAC"
atac2$value[which(atac2$value==0)]="No_ATAC"
colnames(atac2)[3]="ATAC_supported"
CREanno=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno=left_join(CREanno,atac1[,c(4,7)], by=c("CREID"="V4"),copy=F)
CREanno$Neuron_series = "No"
CREanno$Neuron_series[!is.na(CREanno$'Neuron-series')] = "Yes"
CREanno$ATAC="noATAC"
CREanno$ATAC[which(CREanno$'Neuron-series'==1)]="withATAC"
CREanno%>%group_by(Neuron_series,ATAC)%>%dplyr::summarise(count=n())

atac2=left_join(atac2,CREanno[,c(1,32)], by=c("V4"="CREID"),copy=F)
atac3=atac2%>%group_by(variable, promoter_type,ATAC_supported)%>%summarise(count=n())%>%mutate(percent=count/sum(count)*100)

atac3=atac3[which(atac3$promoter_type != "CTCF-alone"),]
atac3$label=paste0(signif(atac3$percent,3),"%")
write.table(atac3,gzfile(paste0(path_fig2_data,"Neuron_tCRE_ATAC.plot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
# use cell type agnostic ATAC support in the table
# stored in [primary_folder]/fig2/data
# for -> figExt2b


#===============================================================================
#bash
#n5 tCRE submit version
setwd(paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/n5_intersect/"))

system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.donor2bp.bed.gz -s | gzip > tCREsummit_donor.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.acceptor2bp.bed.gz -s | gzip > tCREsummit_acceptor.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.intron.bed.gz -s | gzip > tCREsummit_intron.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.exon.bed.gz -s | gzip > tCREsummit_exon.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.3UTR.bed.gz -s | gzip > tCREsummit_3UTR.bed.gz"))
system(paste0("bedtools intersect -c -wa -a ",SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz -b ",GENCODE_path,"gencode.v39.annotation.bed6.5UTR.bed.gz -s | gzip > tCREsummit_5UTR.bed.gz"))
#===========================================

summit=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.stranded_summit.bed.gz"), header=F, stringsAsFactors = F)
pathTSS=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/n5_intersect/")
files=list.files(path=pathTSS, pattern="tCREsummit_")
files.names=gsub(".bed.gz","",files)
for (i in 1:6){
  data=read.delim(paste0(pathTSS, files[i]),header=F, stringsAsFactors=F)
  data$label=paste0(data$V1,"_",data$V2,"_",data$V3,"_",data$V6)
  data$V7[which(data$V7>0)]=1
  summit=left_join(summit, data[,c(4,7)], by="V4", copy=F)}
colnames(summit)[c(7:12)]=files.names
write.table(summit, gzfile(paste0(pathTSS, "tCREsummit_all_info.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
gc()

pathTSS=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/n5_intersect/")
summit=read.delim(paste0(pathTSS, "tCREsummit_all_info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

summit$class2="Intergenic"
summit$class2[which(summit$tCREsummit_intron >0)]="Intron"
summit$class2[which(summit$tCREsummit_exon >0)]="Exon"
summit$class2[which(summit$tCREsummit_3UTR >0)]="3'UTR"
summit$class2[which(summit$tCREsummit_5UTR >0)]="5'UTR"
summit%>%group_by(class2)%>%dplyr::summarise(count=n())

n5bed%>%group_by(class2)%>%dplyr::summarise(percentRead=sum(V4)/sum(n5bed$V4))
n5bed%>%group_by(class2)%>%dplyr::summarise(percentCluster=length(unique(n5_TSScluster))/length(unique(n5bed$n5_TSScluster)))
n5bed%>%group_by(class2)%>%dplyr::summarise(percentCluster=length(unique(CREID))/length(unique(n5bed$CREID)))

n5bed1=n5bed[which(n5bed$ATAC == "noATAC" &  n5bed$promoter_type == "unclassed"),]
length(unique(n5bed1$n5_TSScluster)) #5324
length(unique(n5bed1$CREID)) #4139
length(unique(n5bed1$CREID[which(n5bed1$promoter_type2 != "unclassed")]))#255

n5bed1%>%group_by(class2)%>%dplyr::summarise(percentCluster=length(unique(n5_TSScluster))/length(unique(n5bed1$n5_TSScluster)))
n5bed1%>%group_by(class2)%>%dplyr::summarise(percentCluster=length(unique(CREID))/length(unique(n5bed1$CREID)))

n5bed1%>%group_by(class2)%>%dplyr::summarise(count=sum(V4),event=n(), cluster=length(unique(n5_TSScluster)))
n5bed2=n5bed[which(n5bed$ATAC == "noATAC" & n5bed$class2 == "3'UTR" & n5bed$promoter_type == "unclassed"),]
length(unique(n5bed2$n5_TSScluster)) #654
n5bed13=unique(n5bed1[,c(16,15)])%>%group_by(CREID)%>%dplyr::summarise(class2=paste(class2,collapse=";"))

#========================================
#bash
#to get CHROMHMM-promoter-type: directly link tCRE with the CHROMHMM result without thu aCRE
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
path_HMM=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/chromHMM/")
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",path_HMM, "markalone_iPS_k16.bed.gz > ontCAGE.Neuron_THP1.CRE.coord.markalone_k16ips.bed"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",path_HMM, "markalone_NSC_k16.bed.gz > ontCAGE.Neuron_THP1.CRE.coord.markalone_k16nsc.bed"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz  -b ",path_HMM, "markalone_Neuron_k16.bed.gz > ontCAGE.Neuron_THP1.CRE.coord.markalone_k16nrn.bed"))

#=========================================
#R
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
CRE.ips=read.delim("ontCAGE.Neuron_THP1.CRE.coord.markalone_k16ips.bed", header=F, stringsAsFactors = F, check.names = F)
CRE.nsc=read.delim("ontCAGE.Neuron_THP1.CRE.coord.markalone_k16nsc.bed", header=F, stringsAsFactors = F, check.names = F)
CRE.nrn=read.delim("ontCAGE.Neuron_THP1.CRE.coord.markalone_k16nrn.bed", header=F, stringsAsFactors = F, check.names = F)

CRE.ips1=CRE.ips%>%group_by(V4)%>%dplyr::summarise(type=paste(V16, collapse=";"))
CRE.nsc1=CRE.nsc%>%group_by(V4)%>%dplyr::summarise(type=paste(V16, collapse=";"))
CRE.nrn1=CRE.nrn%>%group_by(V4)%>%dplyr::summarise(type=paste(V16, collapse=";"))

CRE.ips1$type2="CTCF_alone"
CRE.ips1$type2[grep("enhancer", CRE.ips1$type)]="enhancer"
CRE.ips1$type2[grep("romoter", CRE.ips1$type)]="promoter"
CRE.ips1%>%group_by(type2)%>%dplyr::summarise(count=n())
CRE.nsc1$type2="CTCF_alone"
CRE.nsc1$type2[grep("enhancer", CRE.nsc1$type)]="enhancer"
CRE.nsc1$type2[grep("romoter", CRE.nsc1$type)]="promoter"
CRE.nsc1%>%group_by(type2)%>%dplyr::summarise(count=n())
CRE.nrn1$type2="CTCF_alone"
CRE.nrn1$type2[grep("enhancer", CRE.nrn1$type)]="enhancer"
CRE.nrn1$type2[grep("romoter", CRE.nrn1$type)]="promoter"
CRE.nrn1%>%group_by(type2)%>%dplyr::summarise(count=n())

CREa=full_join(CRE.ips1[,c(1,3)],CRE.nsc1[,c(1,3)],by=c("V4"="V4"),copy=F)
CREa=full_join(CREa,CRE.nrn1[,c(1,3)],by=c("V4"="V4"),copy=F)
colnames(CREa)[c(2:4)]=c("k16ips","k16nsc","k16nrn")
CREa$ChromHMM="CTCF_alone"
CREa$ChromHMM[which(CREa$k16ips == "enhancer" | CREa$k16nsc == "enhancer" | CREa$k16nrn == "enhancer")]="enhancer-like"
CREa$ChromHMM[which(CREa$k16ips == "promoter" | CREa$k16nsc == "promoter" | CREa$k16nrn == "promoter")]="promoter-like"
CREa%>%group_by(ChromHMM)%>%dplyr::summarise(count=n())

CRE=left_join(CRE,CREa[,c(1,5)],by=c("CREID"="V4"),copy=F)
CRE$ChromHMM[which(is.na(CRE$ChromHMM))]="unclassed"

CRE1=CRE[which(CRE$promoter_type == "unclassed" & CRE$Neuron_series == "Yes"),]
CRE1$cell_type="non-specific"
CRE1$cell_type[which(rowSums(CRE1[,c(45,46)]==0)==2)]="iPSC-specific"
CRE1$cell_type[which(rowSums(CRE1[,c(44,46)]==0)==2)]="NSC-specific"
CRE1$cell_type[which(rowSums(CRE1[,c(44,45)]==0)==2)]="Neuron-specific"
CRE1$cell_type[which(CRE1$typeStr=="gene_tss")]="GENCODE"

CRE1=left_join(CRE1, summit[,c(4,13)],by=c("CREID"="V4"),copy=F)
CRE1$ChromHMM[which(CRE1$ChromHMM != "unclassed")]="ChroHMM"
CRE1$ChromHMM[which(CRE1$ChromHMM == "unclassed")]="no_ChroHMM"
write.table(CRE1,gzfile(paste0(path_fig2_data,"tCRE_unclass_Neuron.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data
# for -> figExt2c&d


#================================
#CpG cover 200, num 40, TATA motif score >=3
cpg_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/CGI/")

CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
cpg_info=read.delim(paste0(cpg_path,"ontCAGE.Neuron_THP1.CRE.stranded_summit_5kb_extend.CpG.fakebed.tsv"), header=T, stringsAsFactors = F, check.names = F)
cpg_info0=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 4500),] #take the 1000bp around summit
cpg_info0$cover1=5500-cpg_info0$locS
cpg_info0$cover2=cpg_info0$locE-4500
cpg_info0$cover = apply(cpg_info0[,c(8,9)], 1, min)
cpg_info0$cover[which(cpg_info0$cover>1000)]=1000
cpg_info0=cpg_info0[which(cpg_info0$cpgNum >=40),]
cpg_info0=cpg_info0%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info0=cpg_info0[which(cpg_info0$cover >200),]
cpg_info2=cpg_info0%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)
CRE$any_CpG_island="No"
CRE$any_CpG_island[which(CRE$CREID %in% unique(cpg_info2$CREID))]="Yes"
cpg_info1=cpg_info[-which(cpg_info$locS>5000 | cpg_info$locE < 4500),] #only take the 500bp before summit
cpg_info1$cover1=5000-cpg_info1$locS
cpg_info1$cover2=cpg_info1$locE-4500
cpg_info1$cover = apply(cpg_info1[,c(8,9)], 1, min)
cpg_info1$cover[which(cpg_info1$cover>500)]=500
cpg_info1=cpg_info1[which(cpg_info1$cpgNum >=40),]
cpg_info1=cpg_info1%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1=cpg_info1[which(cpg_info1$cover >200),]
cpg_info2=cpg_info1%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)#this doesnt change here
CRE$upstream_CpG_island="No"
CRE$upstream_CpG_island[which(CRE$CREID %in% unique(cpg_info2$CREID))]="Yes"
cpg_info1a=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 5000),] #only take the 500bp after summit
cpg_info1a$cover1=5500-cpg_info1a$locS
cpg_info1a$cover2=cpg_info1a$locE-5000
cpg_info1a$cover = apply(cpg_info1a[,c(8,9)], 1, min)
cpg_info1a$cover[which(cpg_info1a$cover>500)]=500
cpg_info1a=cpg_info1a[which(cpg_info1a$cpgNum >=40),]
cpg_info1a=cpg_info1a%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1a=cpg_info1a[which(cpg_info1a$cover >200),]
cpg_info2a=cpg_info1a%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)
CRE$downstream_CpG_island="No"
CRE$downstream_CpG_island[which(CRE$CREID %in% cpg_info2a$CREID)]="Yes"

tata_info=read.delim(paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/TATAbox/ontCAGE.Neuron_THP1.CRE.stranded_summit_50bb_extend.TBP.fakebed.tsv"), header=T, stringsAsFactors = F)
tata_info1=tata_info[which(tata_info$motif_score >=3),]
tata_info1=tata_info1[which(tata_info1$locS>15 & tata_info1$locS < 24),] #only take the 15bp region
tata_info2=tata_info1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)

CRE$TATA_box="No"
CRE$TATA_box[which(CRE$CREID %in% tata_info2$CREID)]="Yes"
write.table(CRE, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#rename directionality and super enhancer

CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE$orientation=gsub("convergent","2D",CRE$orientation)
CRE$orientation=gsub("bidirectional","2D",CRE$orientation)
CRE$orientation=gsub("unidirectional","1D",CRE$orientation)
CRE$orientation[which(CRE$fwd_rev_count <5)]="Others"
CRE%>%group_by(promoter_type,orientation)%>%dplyr::summarise(count=n())
count.CRE=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count/output/count_matrix/ontCAGE.Neuron_THP1.count.txt"), header=T, stringsAsFactors =F, check.names = F)
count.CRE$iPS=rowSums(count.CRE[,c(22,23)])
count.CRE$NSC=rowSums(count.CRE[,c(2,3)])
count.CRE$Neuron=rowSums(count.CRE[,c(4,5)])
CRE=left_join(CRE, count.CRE[,c(1,24:26)],by="CREID", copy=F)
CRE$SE_all[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like")]="TE"
CRE$SE_all[which(!is.na(CRE$SE_source) & CRE$Neuron_series=="Yes")]="SE"
CRE$SE_all[which(CRE$Neuron_series=="No")]="Others" 
CRE%>%group_by(promoter_type,SE_all)%>%dplyr::summarise(count=n())
CRE$SE_iPSC="Others"
CRE$SE_iPSC[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like" & CRE$iPS >0)]="TE"
CRE$SE_iPSC[which(CRE$Neuron_series=="Yes" & CRE$SE.ips==1 & CRE$SE_iPSC=="TE")]="SE"
CRE$SE_NSC="Others"
CRE$SE_NSC[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like" & CRE$NSC >0)]="TE"
CRE$SE_NSC[which(CRE$Neuron_series=="Yes" & CRE$SE.nsc==1 & CRE$SE_NSC=="TE")]="SE"
CRE$SE_Neuron="Others"
CRE$SE_Neuron[which(CRE$Neuron_series=="Yes" & CRE$promoter_type == "enhancer-like" & CRE$Neuron >0)]="TE"
CRE$SE_Neuron[which(CRE$Neuron_series=="Yes" & CRE$SE.neuron==1 & CRE$SE_Neuron=="TE")]="SE"
write.table(CRE, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#define major strand
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE$chr=sapply(strsplit(CRE$CREID,"_"),"[",1)
CRE$start=as.numeric(sapply(strsplit(CRE$CREID,"_"),"[",2))
CRE$end=as.numeric(sapply(strsplit(CRE$CREID,"_"),"[",3))
CRE$strand=sapply(strsplit(CRE$CREID,"_"),"[",4)
library(GenomicRanges)
gr <- GRanges(seqnames = CRE$chr, ranges = IRanges(start = CRE$start, end = CRE$end))
reduced_gr <- GenomicRanges::reduce(gr, min.gapwidth=1)
overlap_hits <- findOverlaps(gr, reduced_gr)
CRE$region_ID = subjectHits(overlap_hits)
CREa=CRE[,c(1,4,70,71)]%>%group_by(region_ID)%>%dplyr::slice_max(score)
CREa=CREa%>%group_by(region_ID)%>%dplyr::slice_max(strand) #take the "+" > "." > "-"
CRE$representative="No"
CRE$representative[which(CRE$CREID %in% CREa$CREID)]="Yes"
write.table(CRE, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#===============================================================================
#for minor strand tCRE
#CpG cover 200, num 40, TATA motif score >=3
CREsummit4001=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.minorstrand.complete.summit5000.bed.gz"), header=F, stringsAsFactors = F)
CREanno=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREanno1=CREanno[which(CREanno$CREID %in% CREsummit4001$V4),]
write.table(CREanno1, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.minorstrand.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

cpg_path=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/CGI/")
cpg_info=read.delim(paste0(cpg_path,"ontCAGE.Neuron_THP1.CRE.minorstrand.stranded_summit_5kb_extend.CpG.fakebed.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpg_info0=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 4500),] #take the 1000bp surround summit
cpg_info0$cover1=5500-cpg_info0$locS
cpg_info0$cover2=cpg_info0$locE-4500
cpg_info0$cover = apply(cpg_info0[,c(8,9)], 1, min)
cpg_info0$cover[which(cpg_info0$cover>1000)]=1000
cpg_info0=cpg_info0[which(cpg_info0$cpgNum >=40),]
cpg_info0=cpg_info0%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info0=cpg_info0[which(cpg_info0$cover >200),]
cpg_info2=cpg_info0%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)#missed 16
CREanno1$any_CpG_island="No"
CREanno1$any_CpG_island[which(CREanno1$CREID %in% unique(cpg_info2$CREID))]="Yes"
cpg_info1=cpg_info[-which(cpg_info$locS>5000 | cpg_info$locE < 4500),] #only take the 500bp before summit
cpg_info1$cover1=5000-cpg_info1$locS
cpg_info1$cover2=cpg_info1$locE-4500
cpg_info1$cover = apply(cpg_info1[,c(8,9)], 1, min)
cpg_info1$cover[which(cpg_info1$cover>500)]=500
cpg_info1=cpg_info1[which(cpg_info1$cpgNum >=40),]
cpg_info1=cpg_info1%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1=cpg_info1[which(cpg_info1$cover >200),]
cpg_info2=cpg_info1%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)#this doesnt change here
CREanno1$upstream_CpG_island="No"
CREanno1$upstream_CpG_island[which(CREanno1$CREID %in% unique(cpg_info2$CREID))]="Yes"
cpg_info1a=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 5000),] #only take the 500bp after summit
cpg_info1a$cover1=5500-cpg_info1a$locS
cpg_info1a$cover2=cpg_info1a$locE-5000
cpg_info1a$cover = apply(cpg_info1a[,c(8,9)], 1, min)
cpg_info1a$cover[which(cpg_info1a$cover>500)]=500
cpg_info1a=cpg_info1a[which(cpg_info1a$cpgNum >=40),]
cpg_info1a=cpg_info1a%>%group_by(CREID)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1a=cpg_info1a[which(cpg_info1a$cover >200),]
cpg_info2a=cpg_info1a%>%group_by(CREID)%>%dplyr::slice_max(cpgNum)
CREanno1$downstream_CpG_island="No"
CREanno1$downstream_CpG_island[which(CREanno1$CREID %in% cpg_info2a$CREID)]="Yes"

tata_info=read.delim(paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/TATAbox/ontCAGE.Neuron_THP1.CRE.minorstrand_summit_50bb_extend.TBP.fakebed.tsv.gz"), header=T, stringsAsFactors = F)
tata_info1=tata_info[which(tata_info$motif_score >=3),]
tata_info1=tata_info1[which(tata_info1$locS>15 & tata_info1$locS < 24),] #only take the 15bp region
tata_info2=tata_info1%>%group_by(CREID)%>%dplyr::slice_max(motif_score)
CREanno1$TATA_box="No"
CREanno1$TATA_box[which(CREanno1$CREID %in% tata_info2$CREID)]="Yes"
write.table(CREanno1, gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.minorstrand.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#continue with all the tCREs
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE$CpGTATA="Null"
CRE$CpGTATA[which(CRE$any_CpG_island=="Yes" & CRE$TATA_box=="Yes")]="Mix"
CRE$CpGTATA[which(CRE$any_CpG_island=="Yes" & CRE$TATA_box=="No")]="CGI"
CRE$CpGTATA[which(CRE$any_CpG_island=="No" & CRE$TATA_box=="Yes")]="TATA"
CRE%>%group_by(CpGTATA)%>%dplyr::summarise(count=n())

CREminor=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.minorstrand.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CREminor$CpGTATA="Null"
CREminor$CpGTATA[which(CREminor$any_CpG_island=="Yes" & CREminor$TATA_box=="Yes")]="Mix"
CREminor$CpGTATA[which(CREminor$any_CpG_island=="Yes" & CREminor$TATA_box=="No")]="CGI"
CREminor$CpGTATA[which(CREminor$any_CpG_island=="No" & CREminor$TATA_box=="Yes")]="TATA"
CREminor1=CREminor[which(CREminor$strand == "."),]
CREminor2=CREminor[which(CREminor$strand != "."),]
CREminor2=CREminor2[-which(CREminor2$region_ID %in% CREminor1$region_ID),]
CREminor2=CREminor2%>%group_by(region_ID)%>%dplyr::slice_max(strand)
CREminor=rbind(CREminor1,CREminor2)

CRE1=CRE[which(CRE$representative == "Yes"),]
CRE1%>%group_by(CpGTATA)%>%dplyr::summarise(count=n())

CRE=CRE%>%group_by(region_ID)%>%dplyr::summarise(n5_string=paste(CREID, collapse=";"),
                                                 promoter_type=paste(unique(promoter_type), collapse=";"),
                                                 strand=paste(unique(strand),collapse=";"), count=n())

CRE1a=CRE[-grep(";",CRE$promoter_type),]
CRE1b=CRE1a[which(CRE1a$count<=2),]
CRE1b=left_join(CRE1b,CRE1[,c(72,86)],by="region_ID",copy=F)
CRE1b=left_join(CRE1b,CREminor[,c(72,86)],by="region_ID",copy=F, suffix=c("_main","_minor"))

CRE1b$antisense_tCRE="No"
CRE1b$antisense_tCRE[which(!is.na(CRE1b$CpGTATA_minor))]="Yes"
write.table(CRE1b,gzfile(paste0(path_fig2_data,"tCRE_minor_strand_summary.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig2/data

#===============================================================================
#which tCREs are more conserved
#compare with conservation
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE$length=as.numeric(sapply(strsplit(CRE$CREID,"_"),"[",3))-as.numeric(sapply(strsplit(CRE$CREID,"_"),"[",2))
#
CRE=CRE[which(CRE$region1001_rep == "Yes"),]

CREa=CRE[which(CRE$SE_all != "Others"),c("CREID","SE_all","phastCon17_mean_1001")]
CREi=CRE[which(CRE$SE_iPSC != "Others"),c("CREID","SE_iPSC","phastCon17_mean_1001")]
CREn=CRE[which(CRE$SE_NSC != "Others"),c("CREID","SE_NSC","phastCon17_mean_1001")]
CREne=CRE[which(CRE$SE_Neuron != "Others"),c("CREID","SE_Neuron","phastCon17_mean_1001")]
CREa$group4="Neuron-series"
CREi$group4="iPSC"
CREn$group4="NSC"
CREne$group4="Neuron"
colnames(CREa)=colnames(CREi)
colnames(CREn)=colnames(CREi)
colnames(CREne)=colnames(CREi)
data2=rbind(CREa,CREi,CREn,CREne)
#data2=melt(data2, id=c(1,2,4,5))

data3=CRE[which(CRE$orientation != "Others"), c("CREID","orientation","phastCon17_mean_1001","promoter_type")]
data2$group2="Super enhancer"
data3$group2="Directionality"
colnames(data3)=colnames(data2)
data2=rbind(data2,data3)
data2=data2[-which(data2$group4 %in% c("unclassed","CTCF-alone","Neuron-series")),]
colnames(data2)[2]="feature"

data2=left_join(data2, CRE[,c("CREID","CGI","dCGI","TATA","region1001_ID")], by="CREID", copy=F)
#these CGI, dCGI and TATA are mutually exclusive, as in Ex2f

write.table(data2, gzfile(paste0(path_fig2_data,"filtered_tCRE_for_conservation.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
#stored in [primary_folder]/fig2/data

#===============================================================================
# find distance between enhancer and promoter
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CRE.p=read.delim(paste0(SCAFE_path, "ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.CRE.coord.p.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
CRE.p=CRE.p[which(CRE.p$V4 %in% CRE$CREID[which(CRE$promoter_type != "promoter-like")]),]
CRE.p$V19=abs(CRE.p$V19)
CRE.p=unique(CRE.p[,c(4,19)])
colnames(CRE.p)[2]="Distance_PLScCRE"
CRE=left_join(CRE,CRE.p,by=c("CREID"="V4"),copy=F)

#===============================================================================
#define histone marks linked to each CRE
#bash
#Add H3K27me3
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
path_HMM=paste0(primary_folder,"code_n_data/Fig2_CRE_analysis/chromHMM/")
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"iPS-K27m3_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.iPS-K27m3.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"NSC-K27m3_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.NSC-K27m3.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"NRN-K27m3_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.NRN-K27m3.bed.gz"))

#Add H3k27ac
setwd(paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"iPS-K27ac_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.iPS-K27ac.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"NSC-K27ac_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.NSC-K27ac.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ontCAGE.Neuron_THP1.CRE.coord.bed.gz -b ",path_HMM,"NRN-K27ac_pval0.01.300K.bfilt.narrowPeak.hammock.bed.gz | gzip > ontCAGE.Neuron_THP1.CRE.coord.NRN-K27ac.bed.gz"))
#=======================================

pathk=paste0(SCAFE_path,"ontCAGE.Neuron_THP1/bed/")
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
ipsme3=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.iPS-K27m3.bed.gz"), header=F, stringsAsFactors = F)
nscme3=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.NSC-K27m3.bed.gz"), header=F, stringsAsFactors = F)
nrnme3=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.NRN-K27m3.bed.gz"), header=F, stringsAsFactors = F)
CRE$K27ME3_iPS="No"
CRE$K27ME3_iPS[which(CRE$CREID %in% unique(ipsme3$V4))]="Yes"
CRE$K27ME3_NSC="No"
CRE$K27ME3_NSC[which(CRE$CREID %in% unique(nscme3$V4))]="Yes"
CRE$K27ME3_Neuron="No"
CRE$K27ME3_Neuron[which(CRE$CREID %in% unique(nrnme3$V4))]="Yes"
CRE%>%group_by(K27ME3_iPS,K27ME3_NSC,K27ME3_Neuron)%>%dplyr::summarise(count=n())

#
ipsac=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.iPS-K27ac.bed.gz"), header=F, stringsAsFactors = F)
nscac=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.NSC-K27ac.bed.gz"), header=F, stringsAsFactors = F)
nrnac=read.delim(paste0(pathk,"ontCAGE.Neuron_THP1.CRE.coord.NRN-K27ac.bed.gz"), header=F, stringsAsFactors = F)
CRE$K27Ac_iPS="No"
CRE$K27Ac_iPS[which(CRE$CREID %in% unique(ipsac$V4))]="Yes"
CRE$K27Ac_NSC="No"
CRE$K27Ac_NSC[which(CRE$CREID %in% unique(nscac$V4))]="Yes"
CRE$K27Ac_Neuron="No"
CRE$K27Ac_Neuron[which(CRE$CREID %in% unique(nrnac$V4))]="Yes"
CRE%>%group_by(K27Ac_iPS,K27Ac_NSC,K27Ac_Neuron)%>%dplyr::summarise(count=n())
CRE%>%group_by(K27Ac_iPS,K27ME3_iPS)%>%dplyr::summarise(count=n())
write.table(CRE, paste0(SCAFE_path,gzfile("ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#poised, active, bivalent enhancer enrichment with CGI and TATA

#iPSC state alone
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
CRE$CGI="No"
CRE$CGI[which(CRE$any_CpG_island == "Yes" )]="Yes"
CRE$CGI[which(CRE$TATA_box == "Yes")]="No"
CRE$dCGI="No"
CRE$dCGI[which(CRE$downstream_CpG_island == "Yes" & CRE$upstream_CpG_island == "No" & CRE$TATA_box == "No")]="Yes"
CRE$TATA="No"
CRE$TATA[which(CRE$any_CpG_island == "No" & CRE$TATA_box == "Yes")]="Yes"
CRE$CGInap="No"
CRE$CGInap[which(CRE$CGI == "Yes" & CRE$Distance_PLScCRE >= 2000)]="Yes"
CRE$CGIap="No"
CRE$CGIap[which(CRE$CGI == "Yes" & CRE$Distance_PLScCRE < 2000)]="Yes"

write.table(CRE,gzfile(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#================================================
#put the data in [primary_folder]/fig2/data/"
RLE.CRE=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count/output/count_matrix/ontCAGE.Neuron_THP1.RLE.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
write.table(RLE.CRE,gzfile(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.RLE.tsv.gz")),col.names=T, row.names=T, sep="\t", quote=F)
write.table(CRE,gzfile(paste0(path_fig2_data,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)







