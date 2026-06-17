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
path_fig3_data=paste0(primary_folder,"fig3/data/")
SALA_path=paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/")
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/iPSchro/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.all/")

#===============================================================================
# The following folders were removed after incorporation:
# table0_gene, table4_gene, CPAT, Input 
# only provide upon request
#===============================================================================

#===============================================================================
setwd(paste0(SALA_path,"Input_iPS_chro"))
system("sh ./bam_to_bed/pool/00_pool_bed.sh")
system("sh ./CTES_clusters2/00_run_end3_cluster.sh")
system("sh ./CTES_clusters2/00_run.transcript_bed_to_end_bed_bigwig.sh")
system("sh ./CTSS_clusters2/00_cp_CTSS_clusters.sh")
system("perl ./junction_extractor2/batch_run_interactome_long_read.junction_extractor.pl")
system("perl ./junction_extractor2/pool/pool_junction_extractor_info.pl")
system("sh ./input_end3_end5_junct_bed/00_prepare_input_bed.sh")
system("sh ./other_junctions/00_get_junction_list.sh")
# and prepare ./GENCODE_info manually

#===============================================================================
#run SALA assembler
system(paste0("sh ",SALA_path,"transcript.sh"))

#===============================================================================
path1=paste0(SALA_path,"transcript/zenbu/")
path2=paste0(SALA_path,"transcript/log/")
path3=paste0(SALA_path,"transcript/bed/")
path4=paste0(SALA_path,"transcript/CPAT/")
gene0_path=paste0(SALA_path,"table0_gene/iPSchro.table4ref.disable_ref_chain_bound_gene_anno_10percent/")
gene4_path=paste0(SALA_path,"table4_gene/T4_10percent/")

setwd(path2)

######log table
setwd(path2)
transcript_info=read.delim("iPSchro.table4ref.model.info.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

#group for non-detectable, standard and permissive
read_info=fread("iPSchro.table4ref.trnscpt.info.tsv.gz", header=T, stringsAsFactors = F, select=c(1,4,7))
read_info$trnscpt_ID=gsub("WTC-11-NGN2-hiPSC_chromatin_1_rep", "Ic", read_info$trnscpt_ID)
read_info$rep = substr(read_info$trnscpt_ID, start = 1, stop = 3)

read_info2=read_info[which(read_info$set_ID %in% unique(transcript_info$full_set_ID)),c(2,4)]%>%group_by(rep,set_ID)%>%dplyr::summarise(count=n())
read_info2=read_info2[-grep("EN",read_info2$rep),]
read_info2=read_info2[-grep("ONT",read_info2$rep),]
read_info3=spread(read_info2, key=1, value=3)
read_info3=left_join(transcript_info[,c(1,4)],read_info3, by=c("full_set_ID"="set_ID"),copy=F)
read_info3[is.na(read_info3)]=0
colnames(read_info3)=c("model_ID","full_length_set","iPSC_chromatin_1","iPSC_chromatin_2")
write.table(read_info3, gzfile("full_length_support_matrix.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

read_infoa=read_info[which(read_info$set_ID %in% unique(transcript_info$full_set_ID)),]
read_infoa=read_infoa[-grep("EN",read_infoa$trnscpt_ID),]
read_infoa=read_infoa[-grep("ONT",read_infoa$trnscpt_ID),]
write.table(read_infoa[,c(1:3)], gzfile("full_length_support_readID_modelID_pair.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
rm(read_info2)
rm(read_info3)
rm(read_infoa)

read_info1a=read_info[,c(4,3)]%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=n())
rm(read_info)
read_info1a=read_info1a[-grep("EN",read_info1a$rep),]
read_info1a=read_info1a[-grep("ONT",read_info1a$rep),]
read_info1a=separate_rows(read_info1a, model_ID_str, sep = ";", convert = FALSE)
read_info1a=read_info1a%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=sum(count))
read_info2a=spread(read_info1a, key=1, value=3)
colnames(read_info2a)=c("model_ID","iPSC_chromatin_1","iPSC_chromatin_2")

read_info2a=left_join(transcript_info[,c(1,4)],read_info2a, by="model_ID", copy=F)
read_info2a=read_info2a[,-2]
read_info2a[is.na(read_info2a)]=0
write.table(read_info2a, gzfile("partial_length_support_matrix.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
#######################

read_info3=fread("full_length_support_matrix.tsv.gz", header=T, )
transcript_info=left_join(transcript_info,read_info3[,c(1,3:4)], by="model_ID",copy=F)

transcript_info$source="permissive"
transcript_info$source[which(rowSums(transcript_info[,c(19,20)] >=5)>=2)]="standard"

transcript_info$source[grep("ENST",transcript_info$model_ID)]="fulllength_gencode"
transcript_info$source[grep("ONT",transcript_info$model_ID)]="fulllength_ONT"
transcript_info$source[intersect(grep("ENST",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count>0))]="partial_gencode"
transcript_info$source[intersect(grep("ENST",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count==0))]="non_detectable_gencode"
transcript_info$source[intersect(grep("ONT",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count>0))]="partial_ONT"
transcript_info$source[intersect(grep("ONT",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count==0))]="non_detectable_ONT"
transcript_info%>%group_by(source)%>%dplyr::summarise(count=n())

#========================================
#run gene annotator
system(paste0("sh ",SALA_path,"gene0.sh"))

#========================================
#add temparory geneID
gene.info=read.delim(paste0(gene0_path,"log/iPSchro.table4ref.disable_ref_chain_bound_gene_anno_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene.info%>%group_by(g_assign)%>%dplyr::summarise(count=n())
length(unique(gene.info$gene_ID))

colnames(gene.info)[c(7,8)]=c("IN1_gene_ID","IN1_gene_name")
transcript_info=left_join(transcript_info,gene.info[,c(1,7,8)],by="model_ID",copy=F)
transcript_info$gene_novelty="novel"
transcript_info$gene_novelty[grep("ENSG",transcript_info$IN1_gene_ID)]="GENCODE"
transcript_info$gene_novelty[grep("ONTG",transcript_info$IN1_gene_ID)]="ONT"
transcript_info$transcript_novelty="novel"
transcript_info$transcript_novelty[grep("ENST",transcript_info$model_ID)]="GENCODE"
transcript_info$transcript_novelty[grep("ONTT",transcript_info$model_ID)]="ONT"
transcript_info%>%group_by(gene_novelty,transcript_novelty)%>%dplyr::summarise(count=n())

##add transcript ratio##
transcript_info$iPSC_chromatin=transcript_info$iPSC_chromatin_1+transcript_info$iPSC_chromatin_2
transcript_info$iPSC_chromatin[is.na(transcript_info$iPSC_chromatin)]=0
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_iPSC_chromatin=iPSC_chromatin/sum(iPSC_chromatin))
transcript_info$T_ratio_iPSC_chromatin[which(transcript_info$T_ratio_iPSC_chromatin == "NaN")]=0

write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info=read.delim("table0.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

#=======================================================
#generate transcript model 3n and 5n bed
table0_model=read.delim(paste0(path3,"iPSchro.table4ref.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)

options(scipen=999)
table0_model3n=table0_model
table0_model5n=table0_model
table0_model5n$V3[which(table0_model5n$V6=="+")]=table0_model5n$V2[which(table0_model5n$V6=="+")]+1
table0_model5n$V2[which(table0_model5n$V6=="-")]=table0_model5n$V3[which(table0_model5n$V6=="-")]-1
write.table(table0_model5n[order(table0_model5n$V1,table0_model5n$V2),c(1:6)],gzfile(paste0(path3,"iPSchro.table4ref.model.5n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
table0_model3n$V3[which(table0_model3n$V6=="-")]=table0_model3n$V2[which(table0_model3n$V6=="-")]+1
table0_model3n$V2[which(table0_model3n$V6=="+")]=table0_model3n$V3[which(table0_model3n$V6=="+")]-1
write.table(table0_model3n[order(table0_model3n$V1,table0_model3n$V2),c(1:6)],gzfile(paste0(path3,"iPSchro.table4ref.model.3n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(table0_model[order(table0_model$V1,table0_model$V2),],gzfile(paste0(path4,"iPSchro.table0.bed12.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#=========================================================
#bash
#bed12tobed6
setwd(path1)
system("bed12ToBed6 -i iPSchro.table0.bed12.bed.gz | gzip > iPSchro.table0.bed6.bed.gz")
#=========================================================
bed6=read.delim(paste0(path4,"iPSchro.table0.bed6.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
bed6$length=bed6$V3-bed6$V2
bed6a=bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(),transcript_length=sum(length))
transcript_info=left_join(transcript_info, bed6a, by=c("model_ID"="V4"), copy=F)

#========================================================
# Internal priming prediction
chrom_table=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/chrom.tsv"), header=T, stringsAsFactors = F, check.names = F)
options(scipen=999)
t0bed12 <- read.delim(paste0(path3,"iPSchro.table4ref.model.3n.bed.gz"), header=F, stringsAsFactors=F)
t0bed12$V4 <- paste0(t0bed12$V1,"_",t0bed12$V2,"_",t0bed12$V3,"_",t0bed12$V6)
t0bed12 <- t0bed12%>%group_by(V4)%>%dplyr::mutate(V5=n())
t0bed12a3 <- unique(t0bed12[,c(1:6)])
write.table(t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),],gzfile(paste0(path3,"iPSchro.table4ref.model.n3.bed6.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

t0bed12a3$V2 <- t0bed12a3$V2-50
t0bed12a3$V3 <- t0bed12a3$V3+50
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V2>=0),]
t0bed12a3 <- left_join(t0bed12a3,chrom_table,by=c("V1"="Chromosome"),copy=F)
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V3 < t0bed12a3$Length),]
t0bed12a3 <- t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),]
write.table(t0bed12a3[,c(1:6)], gzfile(paste0(path3,"iPSchro.table4ref.model.n3.n101.bed6.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#==========================
#bash
#getfasta
setwd(path3)
system("bedtools getfasta -s -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed iPSchro.table4ref.model.n3.n101.bed6.bed.gz | gzip > iPSchro.table4ref.model.n3.n101.FASTA.gz")
#please download the fasta file separately:
#wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

#======================================
setwd(path3)
extracted_df=read.delim("iPSchro.table4ref.model.n3.n101.FASTA.gz", header=F, stringsAsFactors = F, check.names = F)

for (i in 1:101){extracted_df[,i+7]=substr(extracted_df$V7, start = i, stop = i)}
extracted_df$fracA08=(rowSums(extracted_df[,c(59:66)] == "A"))/8
extracted_df$fracA16=(rowSums(extracted_df[,c(59:74)] == "A"))/16
extracted_df$fracA20=(rowSums(extracted_df[,c(59:78)] == "A"))/20
extracted_df$internal_prime="no"
extracted_df$internal_prime[which(extracted_df$fracA08 > 0.75)]="fracA08"
extracted_df$internal_prime[which(extracted_df$fracA16 > 0.5)]="fracA16"
with_IP=100-(sum(extracted_df$V5[which(extracted_df$internal_prime == "no")])/sum(extracted_df$V5))*100
print(paste0("% of transcript model with potential internal priming: ", signif(with_IP,3),"%"))
extracted_df=extracted_df[,c(1:7,109:112)]
write.table(extracted_df,gzfile(paste0(path2,"potential_internal_prime_TES.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
rm(extracted_df)
rm(t0bed12, t0bed12a3)

#=================================
#add n3 string
data=separate_rows(transcript_info[,c(1,18)],full_set_bound_str,sep="_")
data=data[grep("T",data$full_set_bound_str),]
colnames(data)[2]="n3_string"
transcript_info=left_join(transcript_info,data,by="model_ID",copy=F)

#add 3' end feature: internal priming, ivano PAS, poly(A) prediction
internal_prime_sample=read.delim(paste0(path2,"potential_internal_prime_TES.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n3_ref <- read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n3.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
need=unique(paste0(n3_ref$V1,"_",n3_ref$V2,"_",n3_ref$V3,"_",n3_ref$V6))
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime=="no")] <- "No"
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime!="No")] <- "Yes"
internal_prime_sample$GENCODE <- "No"
internal_prime_sample$GENCODE[which(internal_prime_sample$V4%in% need)] <- "Yes"
internal_prime_sample <- internal_prime_sample[,c(4,12,11)]
colnames(internal_prime_sample) <- c("label","GENCODE","internal_priming")

table0_model3n=read.delim(paste0(path3,"iPSchro.table4ref.model.3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table0_model3n$label=paste0(table0_model3n$V1,"_",table0_model3n$V2,"_",table0_model3n$V3,"_",table0_model3n$V6)
table0_model3n=left_join(table0_model3n, internal_prime_sample, by ="label", copy=F)

length(which(table0_model3n$internal_priming == "Yes")) #66105 model
length(unique(table0_model3n$label[which(table0_model3n$internal_priming == "Yes")])) #51471 CES affect 63984 model
table0_model3n$internal_priming[which(table0_model3n$GENCODE == "Yes" & table0_model3n$internal_priming == "Yes")]="Yes;GENCODE"
table0_model3n%>%group_by(internal_priming)%>%dplyr::summarise(count=n())
transcript_info=left_join(transcript_info, table0_model3n[c(4,9)], by=c("model_ID" = "V4"), copy=F)

#======================================

#get n5 and n3 from table5+ENST
bed12=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.bed12.bed.bgz"), header=F, stringsAsFactors = F)
bed12$V2ex=bed12$V2+1
bed12$V3ex=bed12$V3-1
bed12$V4=substr(bed12$V4, start = 1, stop = 4)
#bed12a=bed12[grep("ENS",bed12$V4),]
bed12_n5a=unique(bed12[which(bed12$V6=="+"),c(1,2,13,6,4)])
bed12_n5b=unique(bed12[which(bed12$V6=="-"),c(1,14,3,6,4)])
colnames(bed12_n5b)=colnames(bed12_n5a)
bed12_n5=rbind(bed12_n5a,bed12_n5b)
bed12_n3a=unique(bed12[which(bed12$V6=="-"),c(1,2,13,6,4)])
bed12_n3b=unique(bed12[which(bed12$V6=="+"),c(1,14,3,6,4)])
colnames(bed12_n3b)=colnames(bed12_n3a)
bed12_n3=rbind(bed12_n3a,bed12_n3b)
remove(bed12_n5a,bed12_n5b,bed12_n3a,bed12_n3b)
bed12_n5=bed12_n5%>%group_by(V1,V2,V2ex,V6)%>%dplyr::slice_min(V4)
bed12_n3=bed12_n3%>%group_by(V1,V2,V2ex,V6)%>%dplyr::slice_min(V4)
bed12_n5$V5=0
bed12_n5=bed12_n5[order(bed12_n5$V1, bed12_n5$V2),c(1,2,3,5,6,4)]
bed12_n5$V4=paste0(bed12_n5$V4,"_CTSS_",1:nrow(bed12_n5))
bed12_n3$V5=0
bed12_n3=bed12_n3[order(bed12_n3$V1, bed12_n3$V2),c(1,2,3,5,6,4)]
bed12_n3$V4=paste0(bed12_n3$V4,"_CTES_",1:nrow(bed12_n3))
write.table(bed12_n5,gzfile(paste0(path3,"table5pENST.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(bed12_n3,gzfile(paste0(path3,"table5pENST.n3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#add GENCODE 3'end to ontCAGE n3 clusters
setwd(path3)
options(scipen=999)
end3=read.delim("iPSchro.table4ref.end3.bed.bgz", header=F, stringsAsFactors = F, check.names = F)
write.table(end3[order(end3$V1,end3$V2),],gzfile(paste0(path1,"iPSchro.table4ref.end3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(end3[order(end3$V1,end3$V2),c(1:6)],gzfile("iPSchro.table4ref.end3.cluster.region.bed6.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#=================================
#bash
#bedtools 
setwd(path3)
system("bedtools intersect -s -wa -wb -a iPSchro.table4ref.end3.cluster.region.bed6.bed.gz -b table5pENST.n3.bed.gz | gzip > iPSchro.table4ref.end3.cluster.region.gencode_t53n.bed.gz")

#================================
setwd(path3)
table0_model2=read.delim("iPSchro.table4ref.end3.cluster.region.gencode_t53n.bed.gz", header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n3_GENCODE_ONT="No"
transcript_info$n3_GENCODE_ONT[which(transcript_info$n3_string %in% unique(table0_model2$V4))]="Yes"
transcript_info%>%group_by(n3_GENCODE_ONT)%>%summarise(count=n())

transcript_info$n3_support="cluster"
transcript_info$n3_support[grep("XT",transcript_info$n3_string)]="non_cluster"
transcript_info%>%group_by(internal_priming,n3_support)%>%dplyr::summarise(count=n())
transcript_info%>%group_by(n3_GENCODE_ONT,n3_support)%>%dplyr::summarise(count=n())

data0=transcript_info[,c(1,32:33)]
data0$n3_support[which(data0$n3_support=="cluster")]="n3_cluster"
data0$n3_support[which(data0$n3_support=="non_cluster")]=NA
data0$n3_GENCODE_ONT[which(data0$n3_GENCODE_ONT == "Yes")]="GENCODE_ONT"
data0$n3_GENCODE_ONT[which(data0$n3_GENCODE_ONT == "No")]=NA
data0=melt(data0, id=1)
data0=data0[which(!is.na(data0$value)),]
data01=data0%>%group_by(model_ID)%>%dplyr::summarise(n3_support=paste(value,collapse=" & "))
transcript_info=transcript_info[,-which(colnames(transcript_info) == "n3_support")]

transcript_info=left_join(transcript_info, data01, by="model_ID", copy=F)
transcript_info$n3_support[which(is.na(transcript_info$n3_support))]="no_support"
transcript_info$n3_support[which(transcript_info$internal_priming == "Yes")]="internal_priming"
transcript_info%>%group_by(n3_support)%>%dplyr::summarise(count=n())

#===================================================================================================================
###assign tCRE to SCAFE supported transcript
#manage n5 cluster bed
setwd(path3)
options(scipen=999)
end5=read.delim("iPSchro.table4ref.end5.bed.bgz", header=F, stringsAsFactors = F, check.names = F)
write.table(end5[order(end5$V1,end5$V2),],gzfile(paste0(path1,"iPSchro.table4ref.end5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(end5[order(end5$V1,end5$V2),c(1:6)],gzfile("iPSchro.table4ref.end5.cluster.region.bed6.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

##bedtools intersect
#========================================
setwd(path3)
system(paste0("bedtools intersect -s -wa -wb -a iPSchro.table4ref.end5.cluster.region.bed6.bed.gz -b ",SCAFE_path, "bed/ontCAGE.all.cluster.coord.bed.gz | gzip > iPSchro.table4ref.end5.cluster.region.cluster.bed.gz"))
system("bedtools intersect -s -wa -wb -a iPSchro.table4ref.end5.cluster.region.bed6.bed.gz -b table5pENST.n5.bed.gz | gzip > iPSchro.table4ref.end5.cluster.region.gencode_t55n.bed.gz")

#========================================
#use cluster to cluster intersect and link to tCRE
setwd(path3)
transcript_info$n5_string=sapply(strsplit(transcript_info$full_set_bound_str,"_"),"[",1)
table0_model4=read.delim("iPSchro.table4ref.end5.cluster.region.cluster.bed.gz", header=F, stringsAsFactors = F, check.names = F) 
cluster=read.delim(paste0(SCAFE_path,"log/ontCAGE.all.cluster.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table0_model4=left_join(table0_model4,cluster[,c(1,16)],by=c("V10"="clusterID"),copy=F)
table0_model4a=unique(table0_model4[,c(4,10,19)])%>%group_by(V4)%>%dplyr::summarise(TSScluster=paste(unique(V10), collapse=";"),
                                                                                    CREID=paste(unique(CREID), collapse=";"))
grep(";",table0_model4a$CREID)

CRE=read.delim(paste0(primary_folder,"code_n_data/SCAFE/iPSchro/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.all.CRE.info.p.e.se.tsv"), header=T, stringsAsFactors = F, check.names = F)
table0_model4a=left_join(table0_model4a,CRE[,c(1,31,32,34,39)], by="CREID",copy=F)
transcript_info=left_join(transcript_info, table0_model4a, by=c("n5_string"="V4"),copy=F)

table0_model5=read.delim("iPSchro.table4ref.end5.cluster.region.gencode_t55n.bed.gz", header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n5_GENCODE_ONT="no"
transcript_info$n5_GENCODE_ONT[which(transcript_info$n5_string %in% unique(table0_model5$V4))]="yes"
transcript_info%>%group_by(n5_GENCODE_ONT,promoter_type)%>%summarise(count=n())

data1=transcript_info[,c(1,37,41)]
data1$promoter_type[!is.na(data1$promoter_type)]="SCAFE"
data1$n5_GENCODE_ONT[which(data1$n5_GENCODE_ONT == "yes")]="GENCODE_ONT"
data1$n5_GENCODE_ONT[which(data1$n5_GENCODE_ONT == "no")]=NA
data1=melt(data1, id=1)
data1=data1[which(!is.na(data1$value)),]
data2=data1%>%group_by(model_ID)%>%dplyr::summarise(n5_support=paste(value,collapse=" & "))
transcript_info=left_join(transcript_info, data2, by="model_ID", copy=F)
transcript_info$n5_support[which(is.na(transcript_info$n5_support))]="no_support"
transcript_info%>%group_by(n5_support)%>%dplyr::summarise(count=n())

#check the 5n and 3n
transcript_infok=transcript_info[grep("ENST",transcript_info$model_ID),]
transcript_infok%>%group_by(n5_support)%>%dplyr::summarise(count=n())
transcript_infok%>%group_by(n3_support)%>%dplyr::summarise(count=n())

#add CPAT coding potential
#=================================================
##bed12 getfasta
setwd(path1)
system("bedtools getfasta -s -nameOnly -split -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed iPSchro.table0.bed12.bed.gz > iPSchro.table0.fasta")

##run CPAT
setwd(path4)
system(paste0("cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData --top-orf=5 -g ",path1,"iPSchro.table0.fasta -o output1"))
system("rm output1.ORF_seqs.fa")
system("gzip *.tsv")

setwd(path1)
system("gzip iPSchro.table0.fasta")

#========================================

#add CPAT coding potential
setwd(path2)

CPAT=read.delim(paste0(pth4,"output1.ORF_prob.best.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CPAT$seq_ID=sapply(strsplit(CPAT$seq_ID,"\\("),"[",1)
CPAT$seq_ID=gsub("ONTTC","ONTTc",CPAT$seq_ID)
transcript_info=left_join(transcript_info,CPAT[,c(1,8,11)], by=c("model_ID"="seq_ID"),copy=F)
transcript_info$Coding_prob[is.na(transcript_info$Coding_prob)]=0
transcript_info$CPAT_class="coding"
transcript_info$CPAT_class[which(transcript_info$Coding_prob < 0.364)]="non-coding"
transcript_info%>%group_by(CPAT_class)%>%dplyr::summarise(count=n())
write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#assign transcript class
#class from gencode
gtf=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.gtf.gz"), header=F, stringsAsFactors = F, check.names = F)
gtf=gtf[which(gtf$V3 == "transcript"),]
gtf$gene=sapply(strsplit(gtf$V9, "; "),"[",1)
gtf$transcript=sapply(strsplit(gtf$V9, "; "),"[",2)
gtf$gene_type=sapply(strsplit(gtf$V9, "; "),"[",3)
gtf$transcript_type=sapply(strsplit(gtf$V9, "; "),"[",5)
gtf$gene=gsub("gene_id ","",gtf$gene)
gtf$transcript=gsub("transcript_id ","",gtf$transcript)
gtf$gene_type=gsub("gene_type ","",gtf$gene_type)
gtf$transcript_type=gsub("transcript_type ","",gtf$transcript_type)
gtf=gtf[,c(10:13)]
write.table(gtf, gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.transcript.gene.class.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info=left_join(transcript_info, unique(gtf[,c(1,3)]), by=c("IN1_gene_ID"="gene"),copy=F)
transcript_info$Gencode_ONT_geneClass2="NA"
transcript_info$Gencode_ONT_geneClass2[which(!is.na(transcript_info$gene_type))]="others"
transcript_info$Gencode_ONT_geneClass2[which(transcript_info$gene_type == "lncRNA")]="lncRNA"
transcript_info$Gencode_ONT_geneClass2[which(transcript_info$gene_type == "protein_coding")]="protein_coding"

transcript_info=left_join(transcript_info, gtf[,c(2,4)], by=c("model_ID"="transcript"),copy=F)
transcript_info$Gencode_ONT_transcriptClass2="NA"
transcript_info$Gencode_ONT_transcriptClass2[which(!is.na(transcript_info$transcript_type))]="others"
transcript_info$Gencode_ONT_transcriptClass2[which(transcript_info$transcript_type == "lncRNA")]="lncRNA"
transcript_info$Gencode_ONT_transcriptClass2[which(transcript_info$transcript_type == "protein_coding")]="protein_coding"
transcript_info%>%group_by(Gencode_ONT_geneClass2)%>%dplyr::summarise(count=n())
transcript_info%>%group_by(Gencode_ONT_transcriptClass2)%>%dplyr::summarise(count=n())

transcript_info$Novel_transcriptClass = "NA"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel")] = "others"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel" & transcript_info$CPAT_class == "non-coding")] = "ncRNA"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel" & transcript_info$CPAT_class == "non-coding"  & transcript_info$transcript_length > 200)] = "lncRNA"
transcript_info%>%group_by(Gencode_ONT_geneClass2,Novel_transcriptClass)%>%dplyr::summarise(count=n())
write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

##see which ENST and ENSG we have re-adjusted 5'end, report every ENST changed, filter them afterwards. For ENSG, need to re-do after table 4
transcript_info1aa=transcript_info[which(transcript_info$source == "fulllength_gencode" | transcript_info$source == "fulllength_ONT"),]
transcript_model=read.delim(paste0(path3,"iPSchro.table4ref.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
transcript_model_genecode=transcript_model[which(transcript_model$V4 %in% transcript_info1aa$model_ID),]

gencode_tran=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.bed12.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
gencode_tran1=gencode_tran[which(gencode_tran$V4 %in% transcript_info1aa$model_ID),]
gencode_tran1=left_join(gencode_tran1[,c(1:6)], transcript_model_genecode[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new"))
gencode_tran1=left_join(gencode_tran1,transcript_info1aa[,c(1,31:33,42)],by=c("V4"="model_ID"), copy=F)
gencode_tran1$n5_adjust="5n_adjust"
gencode_tran1$n3_adjust="3n_adjust"
gencode_tran1$n5_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="+"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="-"))]=NA
gencode_tran1$n3_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="-"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="+"))]=NA
gencode_tran2=melt(gencode_tran1[,c(4,16,17)], id=1)
gencode_tran2=gencode_tran2[which(!is.na(gencode_tran2$value)),]
gencode_tran2=gencode_tran2%>%group_by(V4)%>%dplyr::summarise(ENST_adjust=paste(value,collapse=" & "))
gencode_tran1=left_join(gencode_tran1,gencode_tran2, by="V4", copy=F)
gencode_tran1$ENST_adjust[is.na(gencode_tran1$ENST_adjust)]="No"

gencode_tran1n5=gencode_tran1[c(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="+"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="-")),]
gencode_tran1n5%>%group_by(n5_support)%>%dplyr::summarise(count=n())
gencode_tran1n3=gencode_tran1[c(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="-"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="+")),]
gencode_tran1n3%>%group_by(n3_support)%>%dplyr::summarise(count=n())

#===============================================================================
#We dont change the n5 and n3 according to the chromatin bound data

#===============================================================================
#gene level novel geneclass, temporary for table 1
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_overall=full_qry_count/sum(full_qry_count))
transcript_infoa=transcript_info[which(transcript_info$Novel_transcriptClass %in% c("ncRNA","lncRNA") & transcript_info$gene_novelty == "novel"),]
transcript_infoa=transcript_infoa%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_overall2=full_qry_count/sum(full_qry_count))
transcript_infob=transcript_infoa%>%group_by(IN1_gene_ID)%>%dplyr::summarise(ncRNA_rate=sum(T_ratio_overall), weighted_average_length=sum(T_ratio_overall2*transcript_length))
transcript_infob$Novel_geneClass="others"
transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5)]="ncRNA"
transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5 & transcript_infob$weighted_average_length > 200)]="lncRNA"
transcript_infob%>%group_by(Novel_geneClass)%>%dplyr::summarise(count=n())
transcript_info=left_join(transcript_info,transcript_infob[,c(1,4)], by="IN1_gene_ID", copy=F)
transcript_info$Novel_geneClass[which(is.na(transcript_info$Novel_geneClass) & transcript_info$gene_novelty == "novel")]="others"
transcript_info%>%group_by(Novel_geneClass, Novel_transcriptClass)%>%dplyr::summarise(count=n())

write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info=read.delim("table0.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

transcript_info%>%group_by(promoter_type)%>%dplyr::summarise(total=sum(full_qry_count))

#ncRNA class for table0 -> all transcript first
#table5a=transcript_info[which(transcript_info$Novel_transcriptClass %in% c("lncRNA","ncRNA")),]
#table5b=transcript_info[which(transcript_info$Gencode_transcriptClass2 %in% c("lncRNA")),]
#table5=rbind(table5a, table5b)

#novel isoform consider as "genomic" or not
#genomic->does not overlap with any GENCODE acceptor nor donor, overlap with both exon and intron of any GENCODE transcripts
table6=transcript_info[which(transcript_info$transcript_novelty == "novel"),]
options(scipen=999)
table0.bed12=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/iPSchro.table4ref/zenbu/iPSchro.table0.bed12.bed.gz", header=F, stringsAsFactors = F, check.names = F)
table6.bed12=table0.bed12[which(table0.bed12$V4 %in% table6$model_ID),]
write.table(table6.bed12[order(table6.bed12$V1, table6.bed12$V2),], gzfile("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/iPSchro.table4ref/zenbu/table6.onlyNovelTranscript.bed12.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)
table0.end5=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/iPSchro.table4ref/bed/iPSchro.table4ref.model.5n.bed.gz", header=F, stringsAsFactors = F, check.names = F)
table6.end5=table0.end5[which(table0.end5$V4 %in% table6$model_ID),]
write.table(table6.end5[order(table6.end5$V1, table6.end5$V2),], gzfile("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/iPSchro.table4ref/zenbu/table6.onlyNovelTranscript.5n.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
setwd(path2)
write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#pre-filter: internal priming gone
transcript_info1a=transcript_info[which(transcript_info$source != "non_detectable_gencode" & transcript_info$source != "non_detectable_ONT"),]
transcript_info1c=transcript_info1a[which(transcript_info1a$internal_priming != "Yes" ),]
transcript_info1d=transcript_info1c[which(transcript_info1c$source != "partial_gencode" & transcript_info1c$source != "partial_ONT"),]

write.table(transcript_info1c, gzfile("table1.remove_undetect_gencode_and_internal_prime.573k.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info1b=transcript_info[which(transcript_info$source == "non_detectable_gencode" | transcript_info$source == "non_detectable_ONT"),]
write.table(transcript_info1b, gzfile("table1.only_undetect_gencode_ONT.277K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info1b1=transcript_info[which(transcript_info$source == "partial_gencode" | transcript_info$source == "partial_ONT"),]
write.table(transcript_info1b1, gzfile("table1.only_partial_gencode_ONT.55K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info1b2=transcript_info1a[which(transcript_info1a$internal_priming == "Yes"),]
write.table(transcript_info1b2, gzfile("table1.only_internal_priming.52K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)


#table4 for downstream analysis, all include SCAFE supported 5' end
transcript_info4a=transcript_info1c[which(transcript_info1c$transcript_novelty != "novel" & transcript_info1c$gene_novelty != "novel"),] #GENCODE with SCAFE support
transcript_info4a=transcript_info4a[grep("SCAFE", transcript_info4a$n5_support),]
transcript_info4b=transcript_info1c[intersect(grep("SCAFE", transcript_info1c$n5_support),which(transcript_info1c$transcript_novelty == "novel" & transcript_info1c$gene_novelty != "novel" & transcript_info1c$source=="standard")),] 
transcript_info4b=transcript_info4b[which(transcript_info4b$T_ratio_iPSC_chromatin >= 0.1 ),] #novel isoform
transcript_info4c=transcript_info1c[intersect(grep("SCAFE", transcript_info1c$n5_support),which(transcript_info1c$transcript_novelty == "novel" & transcript_info1c$gene_novelty == "novel")),] #novel gene
transcript_info4=rbind(transcript_info4a,transcript_info4b,transcript_info4c)
write.table(transcript_info4, gzfile("table4.chimeric.76K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info4%>%group_by(gene_novelty,transcript_novelty,CPAT_class)%>%dplyr::summarise(count=n())

#bed12 for each table
table0.bed12=read.delim(paste0(path1,"iPSchro.table0.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table1.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1c$model_ID),]
table2.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info2$model_ID),]
table4.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info4$model_ID),]
write.table(table1.bed12, gzfile(paste0(path1,"iPSchro.table1.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table2.bed12, gzfile(paste0(path1,"iPSchro.table2.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table4.bed12, gzfile(paste0(path1,"iPSchro.table4.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare files for table4 gene annotation
setwd(paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/table4_gene/"))
table0.bed12=read.delim(paste0(path3,"iPSchro.table4ref.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
table0.bed12a=table0.bed12[union(union(grep("ENST",table0.bed12$V4),grep("ONT",table0.bed12$V4)),which(table0.bed12$V4 %in% transcript_info4$model_ID)),]
write.table(table0.bed12a[order(table0.bed12a$V1,table0.bed12a$V2),],gzfile("table4wENST.bed12.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#========================================================================
setwd(paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/table4_gene/"))
system("zcat table4wENST.bed12.bed.gz | bgzip > table4wENST.bed12.bed.bgz")
system("tabix -p bed table4wENST.bed12.bed.bgz")

#========================================================================
setwd(paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/table4_gene/"))
table0=read.delim(paste0(path2,"iPSchro.table4ref.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table0a=table0[which(table0$model_ID %in% table0.bed12a$V4),]
write.table(table0a, gzfile("table4wENST.info.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#=====================
#run gene annotator
system(paste0("sh ",SALA_path,"gene4.sh"))


#=======================================================================================================================
#Re-run all gene related properties
transcript_info4=read.delim(paste0(path2,"table4.chimeric.76K.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t4gene=read.delim(paste0(gene4_path,"log/T4_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t4gene%>%group_by(g_assign)%>%dplyr::summarise(count=n())
colnames(t4gene)[c(7,8)]=c("T4_gene_ID","T4_gene_name")
length(unique(t4gene$T4_gene_ID))
transcript_info4=left_join(transcript_info4,t4gene[,c(1,7,8)], by="model_ID", copy=F)
transcript_info4$T4_gene_novelty="novel"
transcript_info4$T4_gene_novelty[grep("ENSG",transcript_info4$T4_gene_ID)]="GENCODE"
transcript_info4$T4_gene_novelty[grep("ONTG",transcript_info4$T4_gene_ID)]="ONT"

#GENCODE geneClass
anno=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5pENST.transcript.gene.class.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(anno)[3]="T4_Gencode_ONT_geneCalss"
transcript_info4=left_join(transcript_info4, unique(anno[,c(1,3)]), by=c("T4_gene_ID"="gene"),copy=F)
transcript_info4$T4_Gencode_ONT_geneCalss2="NA"
transcript_info4$T4_Gencode_ONT_geneCalss2[which(!is.na(transcript_info4$T4_Gencode_ONT_geneCalss))]="others"
transcript_info4$T4_Gencode_ONT_geneCalss2[which(transcript_info4$T4_Gencode_ONT_geneCalss == "lncRNA")]="lncRNA"
transcript_info4$T4_Gencode_ONT_geneCalss2[which(transcript_info4$T4_Gencode_ONT_geneCalss == "protein_coding")]="protein_coding"

#===============================================================================
#gene level novel geneclass
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_overall=full_qry_count/sum(full_qry_count))
transcript_infoa=transcript_info4[which(transcript_info4$Novel_transcriptClass %in% c("ncRNA","lncRNA") & transcript_info4$T4_gene_novelty == "novel"),]
transcript_infoa=transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_overall2=full_qry_count/sum(full_qry_count))
transcript_infob=transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::summarise(ncRNA_rate=sum(T_ratio_overall), weighted_average_length=sum(T4_T_ratio_overall2*transcript_length))
transcript_infob$T4_Novel_geneClass="others"
transcript_infob$T4_Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5)]="ncRNA"
transcript_infob$T4_Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5 & transcript_infob$weighted_average_length > 200)]="lncRNA"
transcript_infob%>%group_by(T4_Novel_geneClass)%>%dplyr::summarise(count=n())
transcript_info4=left_join(transcript_info4,transcript_infob[,c(1,4)], by="T4_gene_ID", copy=F)
transcript_info4$T4_Novel_geneClass[which(is.na(transcript_info4$T4_Novel_geneClass) & transcript_info4$T4_gene_novelty == "novel")]="others"
transcript_info4%>%group_by(T4_Novel_geneClass, Novel_transcriptClass)%>%dplyr::summarise(count=n())

write.table(transcript_info4, gzfile(paste0(path2,"table4.chimeric.76K.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=========================
#re-run transcript ratio to exclude novel transcript model newly assigned to ENSG

transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_iPS_chromatin=iPSC_chromatin/sum(iPSC_chromatin))
transcript_info4$T4_T_ratio_iPS_chromatin[which(transcript_info4$T4_T_ratio_iPS_chromatin == "NaN")]=0

#remove newly assigned ENSG with 10% & 5 reads cutoff
table4ENSG=transcript_info4[which(transcript_info4$T4_gene_novelty != "novel" & transcript_info4$transcript_novelty == "novel"),]
table4ENSG1=table4ENSG[which(table4ENSG$T4_T_ratio_iPS_chromatin >= 0.1 ),]
table4ENSG1=table4ENSG1[which(table4ENSG1$source != "permissive"),] 
table4off=setdiff(table4ENSG$model_ID,table4ENSG1$model_ID)

table4a=transcript_info4[which(transcript_info4$model_ID %in% table4off),]
table5=transcript_info4[-which(transcript_info4$model_ID %in% table4off),]

#gene-base promoter_type
table5=table5%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_gene_promoter_type=paste(unique(promoter_type),collapse=";"))
#promoter>enhancer>unclass
table5$T4_gene_promoter_type[grep("promoter",table5$T4_gene_promoter_type)]="promoter-like"
table5$T4_gene_promoter_type[grep("enhancer",table5$T4_gene_promoter_type)]="enhancer-like"
unique(table5[,c(53,61)])%>%group_by(T4_gene_promoter_type)%>%dplyr::summarise(count=n())

write.table(table5,gzfile(paste0(path2,"table5.chimeric.76K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

table0.bed12=read.delim(paste0(path1,"iPSchro.table0.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.bed12=table0.bed12[which(table0.bed12$V4 %in% table5$model_ID),]
write.table(table5.bed12, gzfile(paste0(path1,"iPSchro.table5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)


#make the pie chart#
table5$transcript_novelty[which(table5$T4_gene_novelty == "GENCODE" & table5$transcript_novelty == "novel")]="Novel isoform"
table5t=table5%>%group_by(T4_gene_novelty,transcript_novelty,promoter_type)%>%dplyr::summarise(count=n())
table5g=unique(table5[,c(85,87, 98)])
table5g1=table5g%>%group_by(T4_gene_novelty,T4_gene_promoter_type)%>%dplyr::summarise(count=n())

#==================================================================
#build gtf in another Rscript

#==================================================================
#subclass of ncRNA
table5=read.delim(paste0(path2,"table5.chimeric.76K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5%>%group_by(T4_gene_novelty,transcript_novelty)%>%dplyr::summarise(count=n())

table5a=table5[which(table5$T4_Gencode_ONT_geneClass2=="lncRNA" & table5$Gencode_ONT_transcriptClass2 =="lncRNA"),]
table5b=table5[which(table5$T4_Gencode_ONT_geneClass2=="lncRNA" & table5$Novel_transcriptClass %in% c("lncRNA")),]
table5c=table5[which(table5$T4_gene_novelty == "novel" & table5$T4_Novel_geneClass %in% c("lncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","ncRNA")),]
table5d=table5[which(table5$T4_gene_novelty == "novel" & table5$T4_Novel_geneClass %in% c("ncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","ncRNA")),]
table5a$T4_ncRNA_source="ENST_ONT_lncRNA"
table5b$T4_ncRNA_source="novel_lncRNA_isoform"
table5c$T4_ncRNA_source="novel_lncRNA"
table5d$T4_ncRNA_source="novel_ncRNA"
table7=rbind(table5a,table5b,table5c,table5d)

options(scipen=999)
table0.bed12=read.delim(paste0(path1,"iPSchro.table0.bed12.bed.gz"),header=F, stringsAsFactors = F, check.names = F)
table7.bed12=table0.bed12[which(table0.bed12$V4 %in% table7$model_ID),]
write.table(table7.bed12[order(table7.bed12$V1,table7.bed12$V2),],gzfile(paste0(path1,"table7_ncRNA.bed12.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
#prepare n5
table7.bed12$V3[which(table7.bed12$V6 == "+")]=table7.bed12$V2[which(table7.bed12$V6 == "+")]+1
table7.bed12$V2[which(table7.bed12$V6 == "-")]=table7.bed12$V3[which(table7.bed12$V6 == "-")]-1
write.table(table7.bed12[order(table7.bed12$V1, table7.bed12$V2),c(1:6)], gzfile(paste0(path1,"table7_ncRNA.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#use #ENST with updated range previously prepared in SALA_Neuron_THP1.R

#copy to current folder
setwd(path1)
system("bed12ToBed6 -i table7_ncRNA.bed12.bed.gz | gzip > table7_ncRNA.bed6.bed.gz")

#bedtools #start site overlap with mRNA exon or not
system("bedtools intersect -wa -s -c -a table7_ncRNA.n5.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | gzip > table7_ncRNA.n5.inter.mRNApseudo.exon.bed.gz")

#bedtools #start site overlap with mRNA transcript
system("bedtools intersect -wa -s -c -a table7_ncRNA.n5.bed.gz -b ENST.T4updated.mRNA_pseudo.bed12.bed.gz | gzip > table7_ncRNA.n5.inter.mRNApseudo.transcript.bed.gz")

#bedtools #whole transcript overlap with mRNA transcript #sense and anti-sense #no hit mean intergenic
system("bedtools intersect -wa -s -c -a table7_ncRNA.bed12.bed.gz -b ENST.T4updated.mRNA_pseudo.bed12.bed.gz | cut -f-6,13- | gzip > table7_ncRNA.bed12.inter.mRNApseudo.transcript.bed.gz")
system("bedtools intersect -wa -S -c -a table7_ncRNA.bed12.bed.gz -b ENST.T4updated.mRNA_pseudo.bed12.bed.gz | cut -f-6,13- | gzip > table7_ncRNA.bed12.antisense.inter.mRNApseudo.transcript.bed.gz")

#bedtools #exon overlap with mRNA exon #hit mean sense overlap
system("bedtools intersect -wa -wb -s -a table7_ncRNA.bed6.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | bedtools overlap -i stdin -cols 2,3,8,9 | gzip > table7_ncRNA.bed6.inter.mRNApseudo.exon.bed.gz")
system("bedtools intersect -wa -wb -S -a table7_ncRNA.bed6.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | bedtools overlap -i stdin -cols 2,3,8,9 | gzip > table7_ncRNA.bed6.antisense.inter.mRNApseudo.exon.bed.gz")

#bedtools #whole transcript overlap with mRNA exon #no hit mean sense intronic/intergenic
system("bedtools intersect -wa -s -c -a table7_ncRNA.bed12.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | gzip > table7_ncRNA.bed12.inter.mRNApseudo.exon.bed.gz")
system("bedtools intersect -wa -S -c -a table7_ncRNA.bed12.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | gzip > table7_ncRNA.bed12.antisense.inter.mRNApseudo.exon.bed.gz")

system("bedtools closest -a table7_ncRNA.n5.bed.gz  -b ENST.T4updated.mRNA_pseudo.n5.bed.gz -D a -S | gzip > table7_ncRNA.n5.divergent.bed.gz")
#===================================================
setwd(path1)
div=read.delim("table7_ncRNA.n5.divergent.bed.gz", header=F, stringsAsFactors = F, check.names = F)
div=div[which(div$V11 != "."),]
div=div[which(abs(div$V13) <= 2000),]

startExon=read.delim("table7_ncRNA.n5.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
startExon=startExon[which(startExon$V7>0),]
startTranscript=read.delim("table7_ncRNA.n5.inter.mRNApseudo.transcript.bed.gz", header=F, stringsAsFactors = F, check.names = F)
startTranscript=startTranscript[which(startTranscript$V7>0),]
fullTranscript=read.delim("table7_ncRNA.bed12.inter.mRNApseudo.transcript.bed.gz", header=F, stringsAsFactors = F, check.names = F)
fullTranscript=fullTranscript[which(fullTranscript$V7>0),]
fullTranscriptAs=read.delim("table7_ncRNA.bed12.antisense.inter.mRNApseudo.transcript.bed.gz", header=F, stringsAsFactors = F, check.names = F)
fullTranscriptAs=fullTranscriptAs[which(fullTranscriptAs$V7>0),]
fullExon=read.delim("table7_ncRNA.bed12.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
fullExon=fullExon[which(fullExon$V13>0),]
fullExonAs=read.delim("table7_ncRNA.bed12.antisense.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
fullExonAs=fullExonAs[which(fullExonAs$V13>0),]
exonExon=read.delim("table7_ncRNA.bed6.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
exonExon=exonExon%>%group_by(V4,V10)%>%dplyr::summarise(overlap=sum(V13))
ref=read.delim("ENST.T4updated.mRNA_pseudo.bed6.bed.gz", header=F, stringsAsFactors = F)
ref=ref%>%group_by(V4)

exonExon=left_join(exonExon,table7[,c(1,29)],by=c("V4"="model_ID"),copy=F)
exonExon$overlap_rate=exonExon$overlap/exonExon$transcript_length
exonExonAs=read.delim("table7_ncRNA.bed6.antisense.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
exonExonAs=exonExonAs%>%group_by(V4,V10)%>%dplyr::summarise(overlap=sum(V13))
exonExonAs=left_join(exonExonAs,table7[,c(1,29)],by=c("V4"="model_ID"),copy=F)
exonExonAs$overlap_rate=exonExonAs$overlap/exonExonAs$transcript_length

table7$divergent="No"
table7$divergent[which(table7$model_ID %in% unique(div$V4))] ="Yes"
table7$startGexon="No"
table7$startGexon[which(table7$model_ID %in% unique(startExon$V4))]="Yes"
table7$startGtranscript="No"
table7$startGtranscript[which(table7$model_ID %in% unique(startTranscript$V4))]="Yes"
table7$fullGexon="No"
table7$fullGexon[which(table7$model_ID %in% unique(fullExon$V4))]="Yes"
table7$fullGexonAs="No"
table7$fullGexonAs[which(table7$model_ID %in% unique(fullExonAs$V4))]="Yes"
table7$exonGexon01="No"
table7$exonGexon01[which(table7$model_ID %in% unique(exonExon$V4[which(exonExon$overlap_rate>0.1)]))]="Yes"
table7$exonGexonAs01="No"
table7$exonGexonAs01[which(table7$model_ID %in% unique(exonExonAs$V4[which(exonExonAs$overlap_rate>0.1)]))]="Yes"
table7$fullGtranscript="No"
table7$fullGtranscript[which(table7$model_ID %in% unique(fullTranscript$V4))]="Yes"
table7$fullGtranscriptAs="No"
table7$fullGtranscriptAs[which(table7$model_ID %in% unique(fullTranscriptAs$V4))]="Yes"

table7$T4_ncRNA_subclass="Unknown"
table7$T4_ncRNA_subclass[which(table7$startGexon=="No" & table7$startGtranscript=="No" & table7$fullGexon=="No" & table7$fullGtranscript=="No" & table7$fullGtranscriptAs=="No")]="Intergenic"
table7$T4_ncRNA_subclass[which(table7$startGexon=="No" & table7$startGtranscript=="No" & table7$fullGexon=="No" & table7$fullGtranscript=="No" & table7$fullGtranscriptAs=="Yes" & table7$exonGexonAs01=="Yes")]="Antisense"
table7$T4_ncRNA_subclass[which(table7$startGexon=="No" & table7$startGtranscript=="No" & table7$fullGexon=="No" & table7$fullGtranscript=="No" & table7$fullGtranscriptAs=="Yes" & table7$exonGexonAs01=="No")]="Antisense_other"
table7$T4_ncRNA_subclass[which(table7$startGexon=="No" & table7$startGtranscript=="No" & table7$fullGexon=="No" & table7$fullGtranscript=="No" & table7$fullGtranscriptAs=="Yes" & table7$fullGexonAs=="No")]="Antisense_intronic"
table7$T4_ncRNA_subclass[which(table7$fullGtranscript=="Yes" & table7$exonGexon01=="Yes" )]="Sense_overlap_RNA"
table7$T4_ncRNA_subclass[which(table7$fullGtranscript=="Yes" & table7$exonGexon01=="No" )]="Sense_RNA_other"
table7$T4_ncRNA_subclass[which(table7$startGexon=="No" & table7$startGtranscript=="Yes" & table7$fullGexon=="No" & table7$fullGtranscript=="Yes")]="Sense_intronic"
table7$T4_ncRNA_subclass[which(table7$divergent=="Yes"  & table7$exonGexonAs01=="No" )]="Divergent"

table7%>%group_by(T4_ncRNA_subclass)%>%dplyr::summarise(count=n())
maindata_table7=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/log/table7.ncRNA.99K.subclass.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
table7a=table7[which(table7$gene_novelty != "novel"),c(1,29,72)]
table7a=left_join(table7a,maindata_table7[,c(1,57,110)], by="model_ID", copy=F)
table7b=table7a[which(table7a$T4_ncRNA_subclass.x != table7a$T4_ncRNA_subclass.y),] #302
table7=left_join(table7,maindata_table7[,c(1,110)], by="model_ID", copy=F)
table7$T4_ncRNA_subclass.x[which(table7$gene_novelty != "novel" & !is.na(table7$T4_ncRNA_subclass.y))]=table7$T4_ncRNA_subclass.y[which(table7$gene_novelty != "novel" & !is.na(table7$T4_ncRNA_subclass.y))]
table7=table7[,c(1:72)]
colnames(table7)[72]="T4_ncRNA_subclass"

write.table(table7,gzfile(paste0(path2,"table7.ncRNA.32K.subclass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
table5=left_join(table5, table7[,c(1,62,72)], by="model_ID", copy=F)
write.table(table5,gzfile(paste0(path2,"table5.chimeric.76K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

table7s=table7[,c(1,72,62:71)]
table7s1=table7s[which(table7s$T4_ncRNA_subclass == "Sense_overlap_RNA"),]



#3n PAS=====================================================
#bedtools closest -a read.collapsed.CESup25.filtered.neuron.all.sort.bed -b /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/full_length_read/3nbed/PAS/parsed.GENCODE_polyA_signal.score_3.plus.filtered.re_name.3n.bed.gz -s -D a | cut -f-7,12-15 | gzip > read.collapsed.CESup25.filtered.neuron.all.sort.PAS.tsv.gz

prediction.neuron=read.csv("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/n3_Ivano/transcriptome_predictions_TES_all_chr.csv", header=T)
prediction.thp1=read.csv("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/n3_Ivano/transcriptome_predictions_TES_all_chr_THP1.csv", header=T)
commonp=intersect(prediction.neuron$cluster_id,prediction.thp1$cluster_id)
prediction=left_join(prediction.neuron[which(prediction.neuron$cluster_id %in% commonp),],prediction.thp1[which(prediction.thp1$cluster_id %in% commonp),], by="cluster_id",copy=F)
prediction1=prediction[which(prediction$tes_predicted.x != prediction$tes_predicted.y),]

#intron count for splicing efficiency
intron=read.delim(paste0(path3,"iPSchro.table4ref.junct.bed.bgz"), header=F, stringsAsFactors = F)
intron1=intron[which(intron$V5 !=0),]
write.table(intron1[order(intron1$V1,intron1$V2),],gzfile(paste0(path1,"iPSchro.table4ref.junct.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)




