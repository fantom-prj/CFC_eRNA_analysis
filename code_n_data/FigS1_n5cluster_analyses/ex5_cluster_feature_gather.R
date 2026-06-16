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
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/")
CGI_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/CGI/")
TATA_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/TATAbox/")
phast_path=paste0(primary_folder,"code_n_data/FigS1_n5cluster_analyses/phastCon/")

#===============================================================================
#ex5 cluster version
#Add feature to ex5_cluster
#define CGI and TATA (CGI cover 200nt, num 40, TATA motif score >=3)

n5cluster=read.delim(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

cpg_info=read.delim(paste0(CGI_path,"Neuron_THP1.S3.end5.summit.table5.5000bp_extend.CpG.fakebed.tsv"), header=T, stringsAsFactors = F, check.names = F)
cpg_info0=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 4500),] #take the 1000bp surround summit
cpg_info0$cover1=5500-cpg_info0$locS
cpg_info0$cover2=cpg_info0$locE-4500
cpg_info0$cover = apply(cpg_info0[,c(9,10)], 1, min)
cpg_info0$cover[which(cpg_info0$cover>1000)]=1000
cpg_info0=cpg_info0[which(cpg_info0$cpgNum >=40),]
cpg_info0=cpg_info0%>%group_by(n5_string)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info0=cpg_info0[which(cpg_info0$cover >200),]
cpg_info2=cpg_info0%>%group_by(n5_string)%>%dplyr::slice_max(cpgNum)#missed 3
n5cluster$any_CpG_island="No"
n5cluster$any_CpG_island[which(n5cluster$n5_string %in% unique(cpg_info2$n5_string))]="Yes"
cpg_info1=cpg_info[-which(cpg_info$locS>5000 | cpg_info$locE < 4500),] #only take the 500bp before summit
cpg_info1$cover1=5000-cpg_info1$locS
cpg_info1$cover2=cpg_info1$locE-4500
cpg_info1$cover = apply(cpg_info1[,c(9,10)], 1, min)
cpg_info1$cover[which(cpg_info1$cover>500)]=500
cpg_info1=cpg_info1[which(cpg_info1$cpgNum >=40),]
cpg_info1=cpg_info1%>%group_by(n5_string)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1=cpg_info1[which(cpg_info1$cover >200),]
cpg_info2=cpg_info1%>%group_by(n5_string)%>%dplyr::slice_max(cpgNum)#missed 3
n5cluster$upstream_CpG_island="No"
n5cluster$upstream_CpG_island[which(n5cluster$n5_string %in% cpg_info2$n5_string)]="Yes"
cpg_info1a=cpg_info[-which(cpg_info$locS>5500 | cpg_info$locE < 5000),] #only take the 500bp after summit
cpg_info1a$cover1=5500-cpg_info1a$locS
cpg_info1a$cover2=cpg_info1a$locE-5000
cpg_info1a$cover = apply(cpg_info1a[,c(9,10)], 1, min)
cpg_info1a$cover[which(cpg_info1a$cover>500)]=500
cpg_info1a=cpg_info1a[which(cpg_info1a$cpgNum >=40),]
cpg_info1a=cpg_info1a%>%group_by(n5_string)%>%dplyr::mutate(cover=sum(cover),count=n())
cpg_info1a=cpg_info1a[which(cpg_info1a$cover >200),]
cpg_info2a=cpg_info1a%>%group_by(n5_string)%>%dplyr::slice_max(cpgNum)#missed 2
n5cluster$downstream_CpG_island="No"
n5cluster$downstream_CpG_island[which(n5cluster$n5_string %in% cpg_info2a$n5_string)]="Yes"
tata_info=read.delim(paste0(TATA_path,"Neuron_THP1.S3.end5.summit.table5.50bp_extend.TBP.fakebed.tsv"), header=T, stringsAsFactors = F)
tata_info1=tata_info[which(tata_info$motif_score >=3),]
tata_info1=tata_info1[which(tata_info1$locS>15 & tata_info1$locS < 24),] #only take the 15bp region
tata_info2=tata_info1%>%group_by(n5_string)%>%dplyr::slice_max(motif_score)
n5cluster$TATA_box="No"
n5cluster$TATA_box[which(n5cluster$n5_string %in% tata_info2$n5_string)]="Yes"

con17=fread(paste0(phast_path,"end5cluster_phastconALL_summary.tsv.gz"), header=T)
con17=con17[which(con17$label=="17way"),]
n5cluster=left_join(n5cluster, con17[,c(2,18,23)], by="n5_string",copy=F)

write.table(n5cluster, gzfile(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)


#===============================================================================
#Add data from matching CRE: directionality, promoter-typing, SE/TE
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

n5cluster=read.delim(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5cluster=right_join(n5cluster, unique(table5[,c("n5_string", "CREID")]),by=c("V4"="n5_string"),copy=F)
n5cluster=left_join(n5cluster, CRE[,c(1:53)], by="CREID",copy=F)

#===============================================================================
#Add data from matching CRE: histone marks inherit from CRE region
CRE=read.delim(paste0(SCAFE_path,"ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, check.names = F, stringsAsFactors = F)

n5cluster=left_join(n5cluster,CRE[,c("CREID","K27ME3_iPS","K27ME3_NSC","K27ME3_Neuron")], by="CREID",copy=F)
n5cluster[which(n5cluster$CpG_island == "No" & n5cluster$TATA_box == "No"),]%>%group_by(downstream_CpG_island,K27ME3_iPS)%>%dplyr::summarise(count=n())
n5cluster=left_join(n5cluster,CRE[,c("CREID","K27Ac_iPS","K27Ac_NSC","K27Ac_Neuron")], by="CREID",copy=F)
n5cluster[which(n5cluster$promoter_type == "enhancer-like" &n5cluster$CpG_island == "No" & n5cluster$TATA_box == "No"),]%>%group_by(K27Ac_iPS,K27ME3_iPS,downstream_CpG_island)%>%dplyr::summarise(count=n())
write.table(n5cluster, gzfile(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#set 1kb region

#n5cluster$length=n5cluster$V3-n5cluster$V2
n5cluster$region1001_start=n5cluster$V7-500
n5cluster$region1001_end=n5cluster$V8+500
#identify overlap 1001 region 
library(GenomicRanges)
gr <- GRanges(seqnames = n5cluster$V1, ranges = IRanges(start = n5cluster$region1001_start, end = n5cluster$region1001_end))
reduced_gr <- reduce(gr)
overlap_hits <- findOverlaps(gr, reduced_gr)
n5cluster$region1001_ID = subjectHits(overlap_hits)
#
library(GenomicRanges)
gr <- GRanges(seqnames = n5cluster$V1, ranges = IRanges(start = n5cluster$V2+1, end = n5cluster$V3))
reduced_gr <- reduce(gr)
overlap_hits <- findOverlaps(gr, reduced_gr)
n5cluster$clustermerge_ID = subjectHits(overlap_hits)

n5cluster1=n5cluster%>%group_by(region1001_ID)%>%dplyr::slice_max(V5)
n5cluster$region1001_rep="No"
n5cluster$region1001_rep[which(n5cluster$n5_string %in% n5cluster1$n5_string)]="Yes"

n5cluster1=n5cluster%>%group_by(clustermerge_ID)%>%dplyr::slice_max(V5)
n5cluster$representative="No"
n5cluster$representative[which(n5cluster$n5_string %in% n5cluster1$n5_string)]="Yes"
write.table(n5cluster, gzfile(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
n5cluster1=n5cluster[which(n5cluster$region1001_rep == "Yes"),]
n5clustera=n5cluster1[which(n5cluster1$SE_all != "Others"),c("n5_string",  "SE_all",  "mean_1001",  "mean_up500")]
n5clusteri=n5cluster1[which(n5cluster1$SE_iPSC != "Others"),c("n5_string",  "SE_iPSC",  "mean_1001",  "mean_up500")]
n5clustern=n5cluster1[which(n5cluster1$SE_NSC != "Others"),c("n5_string",  "SE_NSC",  "mean_1001",  "mean_up500")]
n5clusterne=n5cluster1[which(n5cluster1$SE_Neuron != "Others"),c("n5_string",  "SE_Neuron",  "mean_1001",  "mean_up500")]
n5clustera$group4="Neuron-series"
n5clusteri$group4="iPSC"
n5clustern$group4="NSC"
n5clusterne$group4="Neuron"
colnames(n5clustera)=colnames(n5clusteri)
colnames(n5clustern)=colnames(n5clusteri)
colnames(n5clusterne)=colnames(n5clusteri)
data2=rbind(n5clustera,n5clusteri,n5clustern,n5clusterne)
#data2=melt(data2, id=c(1,2,4,5))

data3=n5cluster1[which(n5cluster1$orientation != "Others"), c("n5_string", "orientation","mean_1001","mean_up500","promoter_type")]
data2$group2="Super enhancer"
data3$group2="Directionality"
colnames(data3)=colnames(data2)
data2=rbind(data2,data3)

data4=data2
data4=data4[-which(data4$group4 %in% c("unclassed","CTCF-alone","Neuron-series")),]
data4$group4=factor(data4$group4,levels=c("enhancer-like","promoter-like","iPSC","NSC","Neuron"))
data4$mean_up500 = as.numeric(data4$mean_up500)
data4$mean_1001 = as.numeric(data4$mean_1001)

data4=left_join(data4, n5cluster[,c("n5_string","CGI","dCGI","TATA")], by="n5_string", copy=F)
#keep only 1D 2D
data4=data4[which(data4$group2=="Directionality" & data4$group4 %in% c("promoter-like","enhancer-like")),]
data4$group=paste0(substr(data4$group4,1,1),"_",data4$SE_iPSC)
write.table(data4, gzfile(paste0(path_fig2_data,"filtered_ex5cluster_for_conservation.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#ex5_cluster base adjacent to promoter
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/"))

n5cluster=read.delim(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5_bed=n5cluster[which(n5cluster$promoter_type == "enhancer-like"),c(1:6)]
n5_bedp=n5cluster[which(n5cluster$promoter_type == "promoter-like"),c(1:6)]
write.table(n5_bed[order(n5_bed$V1,n5_bed$V2),],gzfile("enhancer_ex5_cluster.bed.gz"),col.names=F, row.names=F,sep="\t", quote=F)
write.table(n5_bedp[order(n5_bedp$V1,n5_bedp$V2),],gzfile("promoter_ex5_cluster.bed.gz"),col.names=F, row.names=F,sep="\t", quote=F)

#=====
# against promoter-like cluster
# not used for adjacent-promoter definition
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/"))
system("bedtools closest -a enhancer_ex5_cluster.bed.gz -b promoter_ex5_cluster.bed.gz -D a | gzip > enhancer_promoter_distance_ex5_cluster.bed.gz")
distance=read.delim("enhancer_promoter_distance_ex5_cluster.bed.gz",header=F, stringsAsFactors = F, check.names=F)
distance=unique(distance[,c(4,13)])
colnames(distance)[2]="distance_ex5_p"
n5cluster=left_join(n5cluster,distance, by=c("n5_string"="V4"),copy=F)

#=====
# against promoter cCRE from SCREEN
# used for adjacent-promoter definition
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/"))
system(paste0("bedtools closest -a enhancer_ex5_cluster.bed.gz -b ",primary_folder,"code_n_data/n5_regions/GRCh38-PLS.all.promoter.sort.bed.gz -D a | gzip > enhancer_promoter_distance_cCRE.bed.gz"))
pCREdistance=read.delim("enhancer_promoter_distance_cCRE.bed.gz",header=F, stringsAsFactors = F, check.names=F)
pCREdistance=pCREdistance%>%group_by(V4)%>%slice_max(V13)
pCREdistance=unique(pCREdistance[,c(4,13)])
colnames(pCREdistance)[2]="distance_cCRE_PLS"
n5cluster=left_join(n5cluster,pCREdistance, by=c("n5_string"="V4"),copy=F)

write.table(n5cluster, gzfile(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
#H3K27 state and regulatory elements
n5cluster=read.delim(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n5cluster%>%group_by(downstream_CpG_island, TATA_box)%>%dplyr::summarise(count=n())
n5cluster$group="Un-marked"
n5cluster$group[which(n5cluster$K27ME3_iPS=="Yes")]="Repressed"
n5cluster$group[which(n5cluster$K27Ac_iPS=="Yes")]="Active"
n5cluster$group[which(n5cluster$K27ME3_iPS=="Yes" & n5cluster$K27Ac_iPS=="Yes")]="Co-marked"
n5cluster$CGI="No"
n5cluster$CGI[which(n5cluster$any_CpG_island == "Yes")]="Yes"
n5cluster$CGI[which(n5cluster$TATA_box == "Yes")]="No"
n5cluster$dCGI="No"
n5cluster$dCGI[which(n5cluster$downstream_CpG_island == "Yes" & n5cluster$upstream_CpG_island == "No" & n5cluster$TATA_box == "No")]="Yes"
n5cluster$TATA="No"
n5cluster$TATA[which(n5cluster$downstream_CpG_island == "No" & n5cluster$upstream_CpG_island == "No" & n5cluster$TATA_box == "Yes")]="Yes"
n5cluster$CGInap="No"
n5cluster$CGInap[which(n5cluster$CGI == "Yes" & abs(n5cluster$distance_cCRE_PLS) >= 2000)]="Yes"
n5cluster$CGIap="No"
n5cluster$CGIap[which(n5cluster$CGI == "Yes" & abs(n5cluster$distance_cCRE_PLS) < 2000)]="Yes"

#===============================================================================
#expression of ex5_cluster from iPSC
n5_matrix=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/count_n5cluster/output/count_matrix/ontCAGE.Neuron_THP1.count.txt"), header=T, stringsAsFactors = F, check.names = F)

n5_matrix=n5_matrix[which(n5_matrix$CREID%in% n5cluster$n5_string),]
n5_matrix$ex5_iPSC=rowSums(n5_matrix[,c(22,23)])
n5_matrix$ex5_NSC=rowSums(n5_matrix[,c(2,3)])
n5_matrix$ex5_Neuron=rowSums(n5_matrix[,c(4,5)])
n5_matrix$ex5_THP1=rowSums(n5_matrix[,c(6:13)])
n5_matrix$ex5_dTHP1=rowSums(n5_matrix[,c(14:21)])
n5cluster=left_join(n5cluster,n5_matrix[,c(1,24:28)],by=c("n5_string"="CREID"),copy=F)

n5cluster$CpGTATA="Null"
n5cluster$CpGTATA[which(n5cluster$any_CpG_island == "Yes")]="CGI"
n5cluster$CpGTATA[which(n5cluster$CGIap == "Yes")]="CGIap"
n5cluster$CpGTATA[which(n5cluster$CGInap == "Yes")]="CGInap"
n5cluster$CpGTATA[which(n5cluster$TATA_box == "Yes")]="TATA"
n5cluster$CpGTATA[which(n5cluster$any_CpG_island == "Yes" & n5cluster$CpGTATA=="TATA")]="Others"

write.table(n5cluster, gzfile(paste0(SCAFE_path,"end5.cluster.info.p.e.se.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===============================================================================
write.table(n5cluster,gzfile(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#available also in [primary_folder]/fig2/data







# Analyses not included in the manuscript
#===============================================================================
#perform fisher's exact for CGI, Null, TATA, towards 2D, SE

n5cluster$CpGTATA="Null"
n5cluster$CpGTATA[which(n5cluster$any_CpG_island == "Yes")]="CGI"
n5cluster$CpGTATA[which(n5cluster$TATA_box == "Yes")]="TATA"
n5cluster$CpGTATA[which(n5cluster$any_CpG_island == "Yes" & n5cluster$CpGTATA=="TATA")]="Mix"
n5cluster1=n5cluster[which(n5cluster$promoter_type %in% c("promoter-like","enhancer-like")),]
#n5cluster1=n5cluster1[which(n5cluster1$orientation %in% c("2D","1D")),]
n5cluster1=n5cluster1[which(n5cluster1$CpGTATA != "Mix"),]
n5cluster1$value="Yes"
FE1=spread(n5cluster1[,c("n5_string","promoter_type","orientation","SE_all","CpGTATA","value")],key=5,value=6)
FE1[is.na(FE1)]="No"

FE2=reshape2::melt(FE1[,c("promoter_type","orientation","SE_all","CGI","Null","TATA")], id=c(1:3))
colnames(FE2)[c(2,3,4,5)]=c("2D","SE","variable2","value2")
FE3=reshape2::melt(FE2, id=c(1,4:5))
FE3=FE3[-which(FE3$variable=="SE" & FE3$value %in% c("No","Others")),]
FE3=FE3[-which(FE3$variable=="2D" & FE3$value %in% c("Others")),]
FE3$value[which(FE3$value == "2D")]="Yes"
FE3$value[which(FE3$value == "1D")]="No"
FE3$value[which(FE3$value == "SE")]="Yes"
FE3$value[which(FE3$value == "TE")]="No"

FE3$result=paste0(FE3$value,"_",FE3$value2)
FE4=FE3%>%group_by(promoter_type, variable, variable2, result)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

FE5=spread(FE4[,c(1:5)],key=4, value=5)
colnames(FE5)=c("promoter_type","group","variable","no_no","no_yes","yes_no","yes_yes")
for(i in 1:9){
  GSEATasting <- matrix(c(FE5$no_no[i], FE5$no_yes[i], FE5$yes_no[i], FE5$yes_yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  FE5$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  FE5$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
FE5$FE_logOR=log(FE5$OR)
FE5$sig_level="ns"
FE5$sig_level[which(FE5$p.val<0.05)]="*"
FE5$sig_level[which(FE5$p.val<0.01)]="**"
FE5$sig_level[which(FE5$p.val<0.001)]="***"
FE5$group=factor(FE5$group, levels=c("2D","SE"))
FE5$sig_level=factor(FE5$sig_level, levels=c("ns","*","**","***"))
write.table(FE5, "GCI_TATA_2D_SE.FE.result.tsv", col.names=T, row.names=F, sep="\t", quote=F)
FE5$variable=factor(FE5$variable, levels=c("TATA","Null","CGI"))

ggplot() + 
  labs(x=NULL, y = "n5_cluster group",title= "Enrichment") +
  facet_grid(rows=vars(promoter_type))+
  geom_point(data=FE5, mapping=aes(y=variable, x=group, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0, limits=c(-0.8,0.8))+ 
  scale_size_manual(values=c(0.2,1.8))+
  theme1
#Figure not included

#===============================================================================
#state across cell type

n5cluster$group="Others"
n5cluster$group[which(n5cluster$K27ME3_iPS=="Yes")]="Poised"
n5cluster$group[which(n5cluster$K27Ac_iPS=="Yes")]="Active"
n5cluster$group[which(n5cluster$K27ME3_iPS=="Yes" & n5cluster$K27Ac_iPS=="Yes")]="Bivalent"
n5cluster$group[which(n5cluster$iPS == 0)]="Undetected"
n5cluster$group_NSC="Others"
n5cluster$group_NSC[which(n5cluster$K27ME3_NSC=="Yes")]="Poised"
n5cluster$group_NSC[which(n5cluster$K27Ac_NSC=="Yes")]="Active"
n5cluster$group_NSC[which(n5cluster$K27ME3_NSC=="Yes" & n5cluster$K27Ac_NSC=="Yes")]="Bivalent"
n5cluster$group_NSC[which(n5cluster$NSC == 0)]="Undetected"
n5cluster$group_Neuron="Others"
n5cluster$group_Neuron[which(n5cluster$K27ME3_Neuron=="Yes")]="Poised"
n5cluster$group_Neuron[which(n5cluster$K27Ac_Neuron=="Yes")]="Active"
n5cluster$group_Neuron[which(n5cluster$K27ME3_Neuron=="Yes" & n5cluster$K27Ac_Neuron=="Yes")]="Bivalent"
n5cluster$group_Neuron[which(n5cluster$Neuron == 0)]="Undetected"
n5cluster6=n5cluster[which(n5cluster$promoter_type == "enhancer-like"),]
aa=n5cluster6%>%group_by(group)%>%dplyr::summarise(cell="iPSC",count=n())%>%dplyr::mutate(percent=count/sum(count))
ab=n5cluster6%>%group_by(group_NSC)%>%dplyr::summarise(cell="NSC",count=n())%>%dplyr::mutate(percent=count/sum(count))
ac=n5cluster6%>%group_by(group_Neuron)%>%dplyr::summarise(cell="Neuron",count=n())%>%dplyr::mutate(percent=count/sum(count))
colnames(ab)=colnames(aa)
colnames(ac)=colnames(aa)
aa=rbind(aa,ab,ac)
aa$cell=factor(aa$cell,levels=c("iPSC","NSC","Neuron"))
aa$group=factor(aa$group, levels=c("Active","Bivalent","Poised","Others","Undetected"))
ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enhancer chromatin state", fill=NULL) +
  scale_fill_npg()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1))+
  geom_histogram(data=aa, mapping=aes(y=percent, x=cell, fill=group), size=0.25, color = "black", stat="identity", alpha=1) + 
  theme1
#Figure not included

#===========================

