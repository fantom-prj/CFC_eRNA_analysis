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
spliceAI_input=paste0(primary_folder,"code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_input/")
spliceAI_output=paste0(primary_folder,"code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_output/")
path_fig6_data=paste0(primary_folder,"fig6/data/")

#===============================================================================
#generate vcf file for spliceAI
#GWAS data was downloaded from "http://www.mulinlab.org/causaldb/documentaiton.html" and liftover to hg38
#file not provided, please download it from the webpage
fine=fread("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/gwas_eqtl/chung_gwas/v2.0/credible_set.hg38.txt", header=T, stringsAsFactors = F)
fine=fine[,c(1,27,3,5:6)]
colnames(fine)=c("#CHROM",  "POS",     "ID",  "REF", "ALT") 
fine$POS=fine$POS+1 #1-base
fine$QUAL="."
fine$FILTER="."
fine$INFO="."
fine=unique(fine)
nrow(unique(fine[,c(1,2)]))#917107
write.table(fine[order(fine$`#CHROM`,fine$POS),], paste0(spliceAI_input,"credible_set.hg38.vcf"), col.names=T, row.names=F, sep="\t", quote=F)

#vcf file stored in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_input


#eQTL, pip3, only from brain or iPSC
#data from https://www.ebi.ac.uk/eqtl/
#downloadable files are not provided
data=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/eqtl_LD_pip3.tsv", header=T, stringsAsFactors =F, check.names = F)
data$V1=sapply(strsplit(data$variant,"_"),"[",1)
data$V2=sapply(strsplit(data$variant,"_"),"[",2)
data$V3=sapply(strsplit(data$variant,"_"),"[",3)
data$V4=sapply(strsplit(data$variant,"_"),"[",4)
data$V1=gsub("chr","",data$V1)
data$V2=as.numeric(data$V2)+1 #1-base
data=data[,c(21,22,5,23,24)]
data=unique(data)
colnames(data)=c("#CHROM",  "POS",     "ID",  "REF", "ALT") 
data$QUAL="."
data$FILTER="."
data$INFO="."
nrow(unique(data[,c(1,2)])) #92918
write.table(data[order(data$`#CHROM`,data$POS),], paste0(spliceAI_input,"eqtl.hg38.vcf"), col.names=T, row.names=F, sep="\t", quote=F)

#vcf file stored in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_input


#eQTL, pip3, all samples
#data from https://www.ebi.ac.uk/eqtl/
#downloadable files are not provided
data=read.delim("/analysisdata/fantom6/Interactome/single_cell_wallace/eqtl/allsamples_eqtl_LD_pip3.tsv", header=T, stringsAsFactors =F, check.names = F)
data$V1=sapply(strsplit(data$variant,"_"),"[",1)
data$V2=sapply(strsplit(data$variant,"_"),"[",2)
data$V3=sapply(strsplit(data$variant,"_"),"[",3)
data$V4=sapply(strsplit(data$variant,"_"),"[",4)
data$V1=gsub("chr","",data$V1)
data$V2=as.numeric(data$V2)+1 #1-base
data=data[,c(21,22,5,23,24)]
data=unique(data)
colnames(data)=c("#CHROM",  "POS",     "ID",  "REF", "ALT") 
data$QUAL="."
data$FILTER="."
data$INFO="."
nrow(unique(data[,c(1,2)])) #367989
write.table(data[order(data$`#CHROM`,data$POS),], paste0(spliceAI_input,"allsample_eqtl.hg38.vcf"), col.names=T, row.names=F, sep="\t", quote=F)

#vcf file stored in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_input

#===============================================================================
#these lines were added to the vcf files manually:
##fileformat=VCFv4.2
##reference=GRCh38/hg38
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>

#===============================================================================

#revise spliceAI gene annotation file to include novel transcript model
eg=read.delim(paste0(spliceAI_input,"grch38.txt"), header=T, stringsAsFactors = F, check.names = F)
bed12=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
gtf=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
gtf$transcriptID=sapply(strsplit(gtf$V9,";"),"[",2)
gtf$transcriptID=gsub(" transcript_id ","",gtf$transcriptID)
gtfe=gtf[which(gtf$V3 == "exon"),]
gtft=gtf[which(gtf$V3 == "transcript"),]
gtfe$V4=gtfe$V4-1
gtfe1=gtfe%>%group_by(transcriptID)%>%dplyr::summarise(EXON_START=paste(V4,collapse=","),EXON_END=paste(V5,collapse=","))
gtfe1$EXON_START=paste0(gtfe1$EXON_START,",")
gtfe1$EXON_END=paste0(gtfe1$EXON_END,",")
gtft$V4=gtft$V4-1
final=left_join(gtft[,c(10,1,7,4,5)],gtfe1, by="transcriptID", copy=F)
final$V1=gsub("chr","",final$V1)
colnames(final)=colnames(eg)
write.table(final,paste0(spliceAI_input,"table5pENST.txt"), col.names=T, row.names=F, sep="\t", quote=F)
final=read.delim(paste0(spliceAI_input,"table5pENST.txt"), header=T)

#the annotation file located in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_input
# -> this file "table5pENST.txt" was used to replace grch38.txt

#=========================================
#bash
#getfasta for spliceAI

setwd(spliceAI_input)
system("for file in *1.count.bed.gz; do bedtools getfasta -s -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed \"$file\" |cut -f4,7 > \"${file%.count.bed.gz}.tsv\" ; done")
system("gzip *.tsv")

#Running on separate GPUs the predictions for the VCFs (thus using CUDA_VISIBLE_DEVICES)
setwd(spliceAI_input)
system("CUDA_VISIBLE_DEVICES=0 nohup python3 __main__.py -I allsample_eqtl.hg38.vcf -O spliceai_out/out_allsample_eqtl.hg38.vcf -R genome.fa -A grch38 > run0.log 2>&1 &")
system("CUDA_VISIBLE_DEVICES=0 nohup python3 __main__.py -I eqtl.hg38.vcf -O spliceai_out/out_eqtl.hg38.vcf -R genome.fa -A grch38 > run0.log 2>&1 &")
system("CUDA_VISIBLE_DEVICES=1 nohup python3 __main__.py -I credible_set.hg38.vcf -O spliceai_out/out_credible_set.hg38.vcf -R genome.fa -A grch38 > run1.log 2>&1 &")

#output of spliceAI located in [primary_folder]/code_n_data/Fig6_GWAS_eQTL_SNP/spliceAI_snp_output

#===============================================================================
#parse spliceAI result
eqtlai=read.delim(paste0(spliceAI_output,"out_eqtl.hg38.vcf"), header=T, stringsAsFactors = F, check.names = F, skip=28)
eqtlai$INFO=gsub("SpliceAI=","",eqtlai$INFO)
eqtlai=separate_rows(eqtlai, INFO, sep=",")
eqtlai=eqtlai%>%separate_wider_delim(INFO, delim = "|", names = c("ALLELE", "transcriptID", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"), too_few = "align_start")
eqtlai1=eqtlai[which(eqtlai$DS_AG > 0.5 | eqtlai$DS_AL > 0.5 | eqtlai$DS_DG > 0.5 | eqtlai$DS_DL > 0.5),]
eqtlai1$group="eQTL"
eqtlaiall=read.delim(paste0(spliceAI_output,"out_allsample_eqtl.hg38.vcf"), header=T, stringsAsFactors = F, check.names = F, skip=28)
eqtlaiall$INFO=gsub("SpliceAI=","",eqtlaiall$INFO)
eqtlaiall=separate_rows(eqtlaiall, INFO, sep=",")
eqtlaiall=eqtlaiall%>%separate_wider_delim(INFO, delim = "|", names = c("ALLELE", "transcriptID", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"), too_few = "align_start")
eqtlaiall1=eqtlaiall[which(eqtlaiall$DS_AG > 0.5 | eqtlaiall$DS_AL > 0.5 | eqtlaiall$DS_DG > 0.5 | eqtlaiall$DS_DL > 0.5),]
eqtlaiall1$group="eQTL_all"
gwasai=read.delim(paste0(spliceAI_output,"out_credible_set.hg38.vcf"), header=T, stringsAsFactors = F, check.names = F, skip=28)
gwasai$INFO=gsub("SpliceAI=","",gwasai$INFO)
gwasai=separate_rows(gwasai, INFO, sep=",")
gwasai=gwasai%>%separate_wider_delim(INFO, delim = "|", names = c("ALLELE", "transcriptID", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"), too_few = "align_start")
gwasai1=gwasai[which(gwasai$DS_AG > 0.5 | gwasai$DS_AL > 0.5 | gwasai$DS_DG > 0.5 | gwasai$DS_DL > 0.5),]
gwasai1$group="GWAS"

#===============================================================================
#Limited to the SJs in our list, remove cryptic junction
both=rbind(eqtlai1,eqtlaiall1,gwasai1)

#===
#export for Table S15
both1=both[which(both$group %in% c("eQTL_all", "GWAS")),]
both1$group[which(both1$group=="eQTL_all")]="eQTL"
both1$region="Acceptor"
both1$region[which(both1$DS_DL > 0.5 | both1$DS_DG >0.5)]="Donor"
both1$group2="ENST"
both1$group2[grep("ONTT",both1$transcriptID)]="ONTT"
write.table(both1, paste0(primary_folder,"supplementary_table/TableS15.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

#================
both$DS_AL=as.numeric(both$DS_AL)*(-1)
both$DS_DL=as.numeric(both$DS_DL)*(-1)

bothm1=reshape2::melt(both[,c(1:13,18)],id=c(1:9,14))
bothm1$variable=gsub("DS_","",bothm1$variable)
bothm2=reshape2::melt(both[,c(1:9,14:18)],id=c(1:9,14))
bothm2$variable=gsub("DP_","",bothm2$variable)
bothm2$POS2=bothm2$POS+as.numeric(bothm2$value)
colnames(bothm2)[12]="distance_SNP"
bothm=left_join(bothm1,bothm2[,c("POS","ID","transcriptID", "variable","group","distance_SNP","POS2")],by=c("POS","ID","transcriptID", "variable","group"),copy=F)
bothma=bothm[which(abs(as.numeric(bothm$value)) > 0.5),]
bothma$region="Acceptor"
bothma$region[which(bothma$variable %in% c("DG","DL"))]="Donor"

final=read.delim(paste0(spliceAI_input,"table5pENST.txt"), header=T, stringsAsFactors =F, check.names = F)
final1=final[which(final$`#NAME` %in% unique(bothma$transcriptID)),]
final2=separate_rows(final1, EXON_START, EXON_END, sep=",")
final2=left_join(final2, unique(bothma[,c(9,3,14)]), by=c("#NAME"="transcriptID"), copy=F)
final2$start_diff=abs(final2$POS2-as.numeric(final2$EXON_START))
final2$end_diff=abs(final2$POS2-as.numeric(final2$EXON_END))
final3=final2[which(final2$start_diff<10 | final2$end_diff<10),]
final3$region="Acceptor"
final3$region[which(final3$STRAND=="+" & final3$end_diff <10)]="Donor"
final3$region[which(final3$STRAND=="-" & final3$start_diff <10)]="Donor"

bothma =left_join(bothma, unique(final3[,c(1,8,9,12,10,11,3)]), by=c("transcriptID"="#NAME", "ID","POS2","region"), copy=F)
bothma1=bothma[which(!is.na(bothma$STRAND)),]

bothma1$group2="ENST"
bothma1$group2[grep("ONTT",bothma1$transcriptID)]="ONTT"
bothma1$SNP_pos=paste0(bothma1$`#CHROM`,"_",bothma1$POS)

bothma1=left_join(bothma1,table5[,c("model_ID","transcript_novelty")], by=c("transcriptID"="model_ID"),copy=F)
bothma1$transcript_novelty[which(is.na(bothma1$transcript_novelty))]="ENST"
write.table(bothma1,gzfile(paste0(path_fig6_data,"spliceAI_eQTL_GWAS_acceprot_donor_gain_loss.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#file located in [primary_folder]/fig6/data/
#-> combine with SNP from exon data in GWAS_eQTL_SNP.R
#===============================================================================


