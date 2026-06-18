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
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/")
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")
GENCODE_path=paste0(primary_folder,"code_n_data/GENCODEv39/")

#===============================================================================
setwd(paste0(SALA_path,"Input_Neuron_THP1"))
system("perl ./bam_to_bed/bam_to_bed.pl")
system("sh ./bam_to_bed/pool/00_pool_bed.sh")
system("sh ./CTES_clusters/00_run_end3_cluster.sh")
system("sh ./CTES_clusters/00_run.transcript_bed_to_end_bed_bigwig.sh")
system("sh ./CTSS_clusters/00_cp_CTSS_clusters.sh")
system("perl ./junction_extractor/batch_run_interactome_long_read.junction_extractor.pl")
system("perl ./junction_extractor/pool/pool_junction_extractor_info.pl")
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
gene0_path=paste0(SALA_path,"table0_gene/iPSC_NSC_Neuron.S3.disable_ref_chain_bound_gene_anno_10percent/")
gene4_path=paste0(SALA_path,"table4_gene/Neuron_THP1_T4_10percent/")

setwd(path2)

#gather files for bedtools#
######

#=======================================================
#generate transcript model 3n and 5n position bed
setwd(path2)
table0_model=read.delim(paste0(path3,"Neuron_THP1.S3.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
options(scipen=999)
table0_model3n=table0_model
table0_model5n=table0_model
table0_model5n$V3[which(table0_model5n$V6=="+")]=table0_model5n$V2[which(table0_model5n$V6=="+")]+1
table0_model5n$V2[which(table0_model5n$V6=="-")]=table0_model5n$V3[which(table0_model5n$V6=="-")]-1
write.table(table0_model5n[order(table0_model5n$V1,table0_model5n$V2),c(1:6)],gzfile(paste0(path3,"Neuron_THP1.S3.model.5n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
table0_model3n$V3[which(table0_model3n$V6=="-")]=table0_model3n$V2[which(table0_model3n$V6=="-")]+1
table0_model3n$V2[which(table0_model3n$V6=="+")]=table0_model3n$V3[which(table0_model3n$V6=="+")]-1
write.table(table0_model3n[order(table0_model3n$V1,table0_model3n$V2),c(1:6)],gzfile(paste0(path3,"Neuron_THP1.S3.model.3n.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
write.table(table0_model[order(table0_model$V1,table0_model$V2),],gzfile(paste0(path1,"Neuron_THP1.table0.bed12.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

#n3_cluster
options(scipen=999)
end3=read.delim(paste0(path3,"Neuron_THP1.S3.end3.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
write.table(end3[order(end3$V1,end3$V2),],gzfile(paste0(path1,"Neuron_THP1.S3.end3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(end3[order(end3$V1,end3$V2),c(1:6)],gzfile(paste0(path3,"Neuron_THP1.S3.end3.cluster.region.bed6.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#n5_cluster
end5=read.delim(paste0(path3,"Neuron_THP1.S3.end5.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
write.table(end5[order(end5$V1,end5$V2),],gzfile(paste0(path1,"Neuron_THP1.S3.end5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(end5[order(end5$V1,end5$V2),c(1:6)],gzfile(paste0(path3,"Neuron_THP1.S3.end5.cluster.region.bed6.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
#bedtools
#bed12tobed6
setwd(path1)
system("bed12ToBed6 -i Neuron_THP1.table0.bed12.bed.gz | gzip > Neuron_THP1.table0.bed6.bed.gz")

#detect GENCODE end from n3 cluster 
setwd(path3)
system("bedtools intersect -s -wa -wb -a Neuron_THP1.S3.end3.cluster.region.bed6.bed.gz -b gencode.3n.bed.gz | gzip > Neuron_THP1.S3.end3.cluster.region.gencode3n.bed.gz")

#detect GENCODE 5n and SCAFE cluster from n5 cluster
setwd(path3)
system(paste0("bedtools intersect -s -wa -wb -a Neuron_THP1.S3.end5.cluster.region.bed6.bed.gz -b ",SCAFE_path,"ontCAGE.Neuron_THP1.cluster.coord.bed.gz | gzip > Neuron_THP1.S3.end5.cluster.region.cluster.bed.gz"))
system("bedtools intersect -s -wa -wb -a Neuron_THP1.S3.end5.cluster.region.bed6.bed.gz -b gencode.5n.bed.gz | gzip > Neuron_THP1.S3.end5.cluster.region.gencode5n.bed.gz")

#intersect transcript n5 position with SCAFE TSSCluster
setwd(path3)
system(paste0("bedtools intersect -s -wa -wb -a Neuron_THP1.S3.model.5n.bed.gz -b ",SCAFE_path,"ontCAGE.Neuron_THP1.cluster.coord.bed.gz | gzip > Neuron_THP1.S3.model.5n.SCAFE_TSScluster.bed.gz"))

#Coding potential - bed12 getfasta & CPAT
setwd(path1)
system("bedtools getfasta -s -nameOnly -split -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed Neuron_THP1.table0.bed12.bed.gz > Neuron_THP1.table0.fasta")

setwd(path4)
system(paste0("cpat.py -x Human_Hexamer.tsv -d Human_logitModel.RData --top-orf=5 -g ",path1,"Neuron_THP1.table0.fasta -o output2"))
system("rm output2.ORF_seqs.fa")
system("gzip *.tsv")

setwd(path1)
system("gzip Neuron_THP1.table0.fasta")

#===============================================================================
# Internal priming prediction
chrom_table=read.delim(paste0(SALA_path,"Input_Neuron_THP1/GENCODE_info/chrom.tsv"), header=T, stringsAsFactors = F, check.names = F)
options(scipen=999)
t0bed12 <- read.delim(paste0(path3,"Neuron_THP1.S3.model.bed.bgz"), header=F, stringsAsFactors=F)
t0bed12$V4 <- paste0(t0bed12$V1,"_",t0bed12$V2,"_",t0bed12$V3,"_",t0bed12$V6)
t0bed12 <- t0bed12%>%group_by(V4)%>%dplyr::mutate(V5=n())
t0bed12a3 <- unique(t0bed12[,c(1:6)])
write.table(t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),],gzfile(paste0(path3,"Neuron_THP1.S3.model.n3.bed6.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)

t0bed12a3$V2 <- t0bed12a3$V2-50
t0bed12a3$V3 <- t0bed12a3$V3+50
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V2>=0),]
t0bed12a3 <- left_join(t0bed12a3,chrom_table,by=c("V1"="Chromosome"),copy=F)
t0bed12a3 <- t0bed12a3[which(t0bed12a3$V3 < t0bed12a3$Length),]
t0bed12a3 <- t0bed12a3[order(t0bed12a3$V1,t0bed12a3$V2),]
write.table(t0bed12a3[,c(1:6)], gzfile(paste0(path3,"Neuron_THP1.S3.model.n3.n101.bed6.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#==========================
#bash
#getfasta
setwd(path3)
system("bedtools getfasta -s -bedOut -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed Neuron_THP1.S3.model.n3.n101.bed6.bed.gz | gzip > Neuron_THP1.S3.model.n3.n101.FASTA.gz")
#please download the fasta file separately:
#wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz

#======================================
setwd(path3)
extracted_df=read.delim("Neuron_THP1.S3.model.n3.n101.FASTA.gz", header=F, stringsAsFactors = F, check.names = F)

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

#===============================================================================
#log table
setwd(path2)
transcript_info=read.delim("Neuron_THP1.S3.model.info.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

#group for non-detectable, standard and permissive
read_info=fread("Neuron_THP1.S3.trnscpt.info.tsv.gz", header=T, stringsAsFactors = F, select=c(1,4,7))
read_info$rep = substr(read_info$trnscpt_ID, start = 1, stop = 2)
read_info$rep[grep("P",read_info$trnscpt_ID)] = substr(read_info$trnscpt_ID[grep("P",read_info$trnscpt_ID)], start = 1, stop = 4)
read_info$rep[grep("D",read_info$trnscpt_ID)] = substr(read_info$trnscpt_ID[grep("D",read_info$trnscpt_ID)], start = 1, stop = 4)

read_info2=read_info[which(read_info$set_ID %in% unique(transcript_info$full_set_ID)),c(2,4)]%>%group_by(rep,set_ID)%>%dplyr::summarise(count=n())
read_info2=read_info2[-grep("EN",read_info2$rep),]
read_info3=spread(read_info2, key=1, value=3)
read_info3=left_join(transcript_info[,c(1,4)],read_info3, by=c("full_set_ID"="set_ID"),copy=F)
read_info3[is.na(read_info3)]=0
read_info3=read_info3[,c(1,2,11,12,23,24,13,14,3:10,15:22)]
colnames(read_info3)=c("model_ID","full_length_set","iPSC_1","iPSC_2","NSC_1","NSC_2","Neuron_1","Neuron_2",
                       "DMSO24.noPAP_1","DMSO24.noPAP_2","DMSO24.PAP_1","DMSO24.PAP_2","DMSO96.noPAP_1","DMSO96.noPAP_2","DMSO96.PAP_1","DMSO96.PAP_2",
                       "PMA24.noPAP_1","PMA24.noPAP_2","PMA24.PAP_1","PMA24.PAP_2","PMA96.noPAP_1","PMA96.noPAP_2","PMA96.PAP_1","PMA96.PAP_2")
write.table(read_info3, gzfile("full_length_support_matrix.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

read_infoa=read_info[which(read_info$set_ID %in% unique(transcript_info$full_set_ID)),]
read_infoa=read_infoa[-grep("EN",read_infoa$trnscpt_ID),]
write.table(read_infoa[,c(1:3)], gzfile("full_length_support_readID_modelID_pair.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
rm(read_info2)
rm(read_info3)
rm(read_infoa)

read_info1a=read_info[,c(4,3)]%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=n())
rm(read_info)
read_info1a=read_info1a[-grep("EN",read_info1a$rep),]
read_info1a=separate_rows(read_info1a, model_ID_str, sep = ";", convert = FALSE)
read_info1a=read_info1a%>%group_by(rep,model_ID_str)%>%dplyr::summarise(count=sum(count))
read_info2a=spread(read_info1a, key=1, value=3)
read_info2a=read_info2a[,c(1,10,11,22,23,12,13,2:9,14:21)]
colnames(read_info2a)=c("model_ID","iPSC_1","iPSC_2","NSC_1","NSC_2","Neuron_1","Neuron_2",
                       "DMSO24.noPAP_1","DMSO24.noPAP_2","DMSO24.PAP_1","DMSO24.PAP_2","DMSO96.noPAP_1","DMSO96.noPAP_2","DMSO96.PAP_1","DMSO96.PAP_2",
                       "PMA24.noPAP_1","PMA24.noPAP_2","PMA24.PAP_1","PMA24.PAP_2","PMA96.noPAP_1","PMA96.noPAP_2","PMA96.PAP_1","PMA96.PAP_2")

read_info2a=left_join(transcript_info[,c(1,4)],read_info2a, by="model_ID", copy=F)
read_info2a=read_info2a[,-2]
read_info2a[is.na(read_info2a)]=0
write.table(read_info2a, gzfile("partial_length_support_matrix.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
#######################

read_info3=fread("full_length_support_matrix.tsv.gz", header=T, )
transcript_info=left_join(transcript_info,read_info3[,c(1,3:24)], by="model_ID",copy=F)
colnames(transcript_info)[c(25:40)]=paste0("THP1.",colnames(transcript_info)[c(25:40)])

transcript_info$source="permissive"
transcript_info$source[which(rowSums(transcript_info[,c(19,20)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(21,22)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(23,24)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(25,26)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(27,28)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(29,30)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(31,32)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(33,34)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(35,36)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(37,38)] >=5)>=2)]="standard"
transcript_info$source[which(rowSums(transcript_info[,c(39,40)] >=5)>=2)]="standard"

transcript_info$source[grep("ENST",transcript_info$model_ID)]="fulllength_gencode"
transcript_info$source[intersect(grep("ENST",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count>0))]="partial_gencode"
transcript_info$source[intersect(grep("ENST",transcript_info$model_ID),which(transcript_info$full_qry_count==0 & transcript_info$partial_qry_count==0))]="non_detectable_gencode"
transcript_info%>%group_by(source)%>%dplyr::summarise(count=n())


#========================================
#run gene annotator
system(paste0("sh ",SALA_path,"gene0.sh"))

#========================================
#add temparory geneID

gene.info=read.delim(paste0(gene0_path,"log/iPSC_NSC_Neuron.full.disable_ref_chain_bound_gene_anno_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene.info%>%group_by(g_assign)%>%dplyr::summarise(count=n())
length(unique(gene.info$gene_ID))

colnames(gene.info)[c(7,8)]=c("IN1_gene_ID","IN1_gene_name")
transcript_info=left_join(transcript_info,gene.info[,c(1,7,8)],by="model_ID",copy=F)
transcript_info$gene_novelty="novel"
transcript_info$gene_novelty[grep("ENSG",transcript_info$IN1_gene_ID)]="GENCODE"
transcript_info$transcript_novelty="novel"
transcript_info$transcript_novelty[grep("ENST",transcript_info$model_ID)]="GENCODE"
transcript_info%>%group_by(gene_novelty,transcript_novelty)%>%dplyr::summarise(count=n())

##add transcript ratio##
transcript_info$iPS=transcript_info$iPSC_1+transcript_info$iPSC_2
transcript_info$iPS[is.na(transcript_info$iPS)]=0
transcript_info$NSC=transcript_info$NSC_1+transcript_info$NSC_2
transcript_info$NSC[is.na(transcript_info$NSC)]=0
transcript_info$Neuron=transcript_info$Neuron_1+transcript_info$Neuron_2
transcript_info$Neuron[is.na(transcript_info$Neuron)]=0
transcript_info$THP1_DMSO=rowSums(transcript_info[,c(25:32)])
transcript_info$THP1_DMSO[is.na(transcript_info$THP1_DMSO)]=0
transcript_info$THP1_PMA=rowSums(transcript_info[,c(33:40)])
transcript_info$THP1_PMA[is.na(transcript_info$THP1_PMA)]=0

transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_iPS=iPS/sum(iPS))
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_NSC=NSC/sum(NSC))
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_Neuron=Neuron/sum(Neuron))
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_THP1_DMSO=THP1_DMSO/sum(THP1_DMSO))
transcript_info=transcript_info%>%group_by(IN1_gene_ID)%>%dplyr::mutate(T_ratio_THP1_PMA=THP1_PMA/sum(THP1_PMA))

transcript_info$T_ratio_iPS[which(transcript_info$T_ratio_iPS == "NaN")]=0
transcript_info$T_ratio_NSC[which(transcript_info$T_ratio_NSC == "NaN")]=0
transcript_info$T_ratio_Neuron[which(transcript_info$T_ratio_Neuron == "NaN")]=0
transcript_info$T_ratio_THP1_DMSO[which(transcript_info$T_ratio_THP1_DMSO == "NaN")]=0
transcript_info$T_ratio_THP1_PMA[which(transcript_info$T_ratio_THP1_PMA == "NaN")]=0

write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info=read.delim("table0.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

#==============================================
#add transcript length and exon number
bed6=read.delim(paste0(path1,"Neuron_THP1.table0.bed6.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
bed6$length=bed6$V3-bed6$V2
bed6a=bed6%>%group_by(V4)%>%dplyr::summarise(n_exon=n(),transcript_length=sum(length))
transcript_info=left_join(transcript_info, bed6a, by=c("model_ID"="V4"), copy=F)

#extract n3 string
data=separate_rows(transcript_info[,c(1,18)],full_set_bound_str,sep="_")
data=data[grep("T",data$full_set_bound_str),]
colnames(data)[2]="n3_string"
transcript_info=left_join(transcript_info,data,by="model_ID",copy=F)

#add 3' end feature: internal priming & GENCODE 3' end
internal_prime_sample=read.delim(paste0(path2,"potential_internal_prime_TES.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
n3_ref <- read.delim(paste0(SALA_path,"Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n3.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
need=unique(paste0(n3_ref$V1,"_",n3_ref$V2,"_",n3_ref$V3,"_",n3_ref$V6))
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime=="no")] <- "No"
internal_prime_sample$internal_prime[which(internal_prime_sample$internal_prime!="No")] <- "Yes"
internal_prime_sample$GENCODE <- "No"
internal_prime_sample$GENCODE[which(internal_prime_sample$V4%in% need)] <- "Yes"
internal_prime_sample <- internal_prime_sample[,c(4,12,11)]
colnames(internal_prime_sample) <- c("label","GENCODE","internal_priming")

table0_model3n=read.delim(paste0(path3,"Neuron_THP1.S3.model.3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table0_model3n$label=paste0(table0_model3n$V1,"_",table0_model3n$V2,"_",table0_model3n$V3,"_",table0_model3n$V6)
table0_model3n=left_join(table0_model3n, internal_prime_sample, by ="label", copy=F)

length(which(table0_model3n$internal_priming == "Yes")) #449138 model
length(unique(table0_model3n$label[which(table0_model3n$internal_priming == "Yes")])) #180817 CES affect 449138 model
table0_model3n$internal_priming[which(table0_model3n$GENCODE == "Yes" & table0_model3n$internal_priming == "Yes")]="Yes;GENCODE"
table0_model3n%>%group_by(internal_priming)%>%dplyr::summarise(count=n())
transcript_info=left_join(transcript_info, table0_model3n[c(4,9)], by=c("model_ID" = "V4"), copy=F)

#GENCODE n3 cluster
#================================
table0_model2=read.delim(paste0(path3,"Neuron_THP1.S3.end3.cluster.region.gencode3n.bed.gz"), header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n3_GENCODE="No"
transcript_info$n3_GENCODE[which(transcript_info$n3_string %in% unique(table0_model2$V4))]="Yes"
transcript_info%>%group_by(n3_GENCODE)%>%summarise(count=n())

transcript_info$n3_support="cluster"
transcript_info$n3_support[grep("XT",transcript_info$n3_string)]="non_cluster"
transcript_info%>%group_by(internal_priming,n3_support)%>%dplyr::summarise(count=n())

data0=transcript_info[,c(1,60:61)]
data0$n3_support[which(data0$n3_support=="cluster")]="n3_cluster"
data0$n3_support[which(data0$n3_support=="non_cluster")]=NA
data0$n3_GENCODE[which(data0$n3_GENCODE == "Yes")]="GENCODE"
data0$n3_GENCODE[which(data0$n3_GENCODE == "No")]=NA
data0=melt(data0, id=1)
data0=data0[which(!is.na(data0$value)),]
data01=data0%>%group_by(model_ID)%>%dplyr::summarise(n3_support=paste(value,collapse=" & "))
transcript_info=transcript_info[,-which(colnames(transcript_info) == "n3_support")]

transcript_info=left_join(transcript_info, data01, by="model_ID", copy=F)
transcript_info$n3_support[which(is.na(transcript_info$n3_support))]="no_support"
transcript_info$n3_support[which(transcript_info$internal_priming == "Yes")]="internal_priming"
transcript_info%>%group_by(n3_support)%>%dplyr::summarise(count=n())
#===================================================================================================================

#GENCODE and SCAFE n5 cluster
#finally, use cluster to cluster intersect and link to tCRE
transcript_info$n5_string=sapply(strsplit(transcript_info$full_set_bound_str,"_"),"[",1)
table0_model4=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/Neuron_THP1.S3.end5.cluster.region.cluster.bed.gz", header=F, stringsAsFactors = F, check.names = F) 
cluster=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/log/ontCAGE.Neuron_THP1.cluster.info.tsv.gz",header=T, stringsAsFactors = F, check.names = F)
table0_model4=left_join(table0_model4,cluster[,c(1,16)],by=c("V10"="clusterID"),copy=F)
table0_model4a=unique(table0_model4[,c(4,10,19)])%>%group_by(V4)%>%dplyr::summarise(TSScluster=paste(unique(V10), collapse=";"),
                                                                                    CREID=paste(unique(CREID), collapse=";"))
grep(";",table0_model4a$CREID)
CRE=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv", header=T, stringsAsFactors = F, check.names = F)
table0_model4a=left_join(table0_model4a,CRE[,c(1,31,32,34,39)], by="CREID",copy=F)
transcript_info=left_join(transcript_info, table0_model4a, by=c("n5_string"="V4"),copy=F)

table0_model5=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/Neuron_THP1.S3.end5.cluster.region.gencode5n.bed.gz", header=F, stringsAsFactors = F, check.names = F) 
transcript_info$n5_GENCODE="no"
transcript_info$n5_GENCODE[which(transcript_info$n5_string %in% unique(table0_model5$V4))]="yes"
transcript_info%>%group_by(n5_GENCODE,promoter_type)%>%summarise(count=n())

data1=transcript_info[,c(1,65,69)]
data1$promoter_type[!is.na(data1$promoter_type)]="SCAFE"
data1$n5_GENCODE[which(data1$n5_GENCODE == "yes")]="GENCODE"
data1$n5_GENCODE[which(data1$n5_GENCODE == "no")]=NA
data1=melt(data1, id=1)
data1=data1[which(!is.na(data1$value)),]
data2=data1%>%group_by(model_ID)%>%dplyr::summarise(n5_support=paste(value,collapse=" & "))
transcript_info=left_join(transcript_info, data2, by="model_ID", copy=F)
transcript_info$n5_support[which(is.na(transcript_info$n5_support))]="no_support"
transcript_info%>%group_by(n5_support)%>%dplyr::summarise(count=n())
write.table(transcript_info, gzfile("table0.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#check the 5n and 3n
transcript_infok=transcript_info[grep("ENST",transcript_info$model_ID),]
transcript_infok%>%group_by(n5_support)%>%dplyr::summarise(count=n())
transcript_infok%>%group_by(n3_support)%>%dplyr::summarise(count=n())

#===============================================================================
#add CPAT coding potential
transcript_info=read.delim("table0.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
CPAT=read.delim(paste0(path4,"output2.ORF_prob.best.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
CPAT$seq_ID=sapply(strsplit(CPAT$seq_ID,"\\("),"[",1)
transcript_info=left_join(transcript_info,CPAT[,c(1,8,11)], by=c("model_ID"="seq_ID"),copy=F)
transcript_info$Coding_prob[is.na(transcript_info$Coding_prob)]=0
transcript_info$CPAT_class="coding"
transcript_info$CPAT_class[which(transcript_info$Coding_prob < 0.364)]="non-coding"

#===============================================================================
#assign transcript class
#class from gencode
gencode.gtf=fread(paste0(GENCODE_path,"gencode.v39.annotation.gtf.gz"), header=F, stringsAsFactors = F)
gtf=gtf[which(gtf$V3 == "transcript"),]
gtf$geneID=sapply(strsplit(gtf$V9,";"),"[",1)
gtf$transcriptID=sapply(strsplit(gtf$V9,";"),"[",2)
gtf$Gencode_geneClass=sapply(strsplit(gtf$V9,";"),"[",3)
gtf$Gencode_transcriptClass=sapply(strsplit(gtf$V9,";"),"[",5)
gtf$geneID=gsub("gene_id ","",gtf$geneID)
gtf$transcriptID=gsub(" transcript_id ","",gtf$transcriptID)
gtf$Gencode_geneClass=gsub(" gene_type ","",gtf$Gencode_geneClass)
gtf$Gencode_transcriptClass=gsub(" transcript_type ","",gtf$Gencode_transcriptClass)
print(gtf%>%group_by(Gencode_geneClass)%>%dplyr::summarise(count=n()),n=44)
write.table(gtf[c(10:13,1,4,5,7)],gzfile(paste0(path2,"GENCODE.gene.n.transcript.class.1base.tsv.gz")), col.names=T, row.names=F, sep="\t", quote= F)

transcript_info=left_join(transcript_info, unique(gtf[,c(10,12)]), by=c("IN1_gene_ID"="geneID"),copy=F)
transcript_info$Gencode_geneClass2="NA"
transcript_info$Gencode_geneClass2[which(!is.na(transcript_info$Gencode_geneClass))]="others"
transcript_info$Gencode_geneClass2[which(transcript_info$Gencode_geneClass == "lncRNA")]="lncRNA"
transcript_info$Gencode_geneClass2[which(transcript_info$Gencode_geneClass == "protein_coding")]="protein_coding"

transcript_info=left_join(transcript_info, gtf[,c(11,13)], by=c("model_ID"="transcriptID"),copy=F)
transcript_info$Gencode_transcriptClass2="NA"
transcript_info$Gencode_transcriptClass2[which(!is.na(transcript_info$Gencode_transcriptClass))]="others"
transcript_info$Gencode_transcriptClass2[which(transcript_info$Gencode_transcriptClass == "lncRNA")]="lncRNA"
transcript_info$Gencode_transcriptClass2[which(transcript_info$Gencode_transcriptClass == "protein_coding")]="protein_coding"
transcript_info%>%group_by(Gencode_geneClass2)%>%dplyr::summarise(count=n())
transcript_info%>%group_by(Gencode_transcriptClass2)%>%dplyr::summarise(count=n())

transcript_info$Novel_transcriptClass = "NA"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel")] = "others"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel" & transcript_info$CPAT_class == "non-coding")] = "ncRNA"
transcript_info$Novel_transcriptClass[which(transcript_info$transcript_novelty == "novel" & transcript_info$CPAT_class == "non-coding"  & transcript_info$transcript_length > 200)] = "lncRNA"
transcript_info%>%group_by(Gencode_geneClass2,Novel_transcriptClass)%>%dplyr::summarise(count=n())

#===============================================================================
#which ENST and ENSG we have re-adjusted 5'end, report every ENST changed, filter them afterwards. For ENSG, need to re-do after filtering
transcript_info1aa=transcript_info[which(transcript_info$source == "fulllength_gencode"),]
transcript_model=read.delim(paste0(path3,"Neuron_THP1.S3.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
transcript_model_genecode=transcript_model[which(transcript_model$V4 %in% transcript_info1aa$model_ID),]
gencode_tran=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.transcript.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
gencode_tran1=gencode_tran[which(gencode_tran$V4 %in% transcript_info1aa$model_ID),]
gencode_tran1=left_join(gencode_tran1[,c(1:6)], transcript_model_genecode[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new"))
gencode_tran1=left_join(gencode_tran1,transcript_info1aa[,c(1,59:61,70)],by=c("V4"="model_ID"), copy=F)
gencode_tran1$n5_adjust="5n_adjust"
gencode_tran1$n3_adjust="3n_adjust"
gencode_tran1$n5_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="+"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="-"))]=NA
gencode_tran1$n3_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="-"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="+"))]=NA
gencode_tran2=reshape2::melt(gencode_tran1[,c(4,16,17)], id=1)
gencode_tran2=gencode_tran2[which(!is.na(gencode_tran2$value)),]
gencode_tran2=gencode_tran2%>%group_by(V4)%>%dplyr::summarise(ENST_adjust=paste(value,collapse=" & "))
gencode_tran1=left_join(gencode_tran1,gencode_tran2, by="V4", copy=F)
gencode_tran1$ENST_adjust[is.na(gencode_tran1$ENST_adjust)]="No"
gencode_tran1a=gencode_tran1[which(gencode_tran1$internal_priming != "Yes"),] #do not adjust those with potential internal priming
transcript_info=left_join(transcript_info,gencode_tran1a[,c(4,18)], by=c("model_ID" = "V4"), copy=F)
#stat depends on which table starting from, here contain 49527 ENST, full-length support without internal priming is 44000
length(which(is.na(gencode_tran1a$n5_adjust)))/nrow(gencode_tran1a)
#0.3044482 have same 5'end
length(which(is.na(gencode_tran1a$n3_adjust)))/nrow(gencode_tran1a)
#0.3826495 have same 3'end
length(which(is.na(gencode_tran1a$n5_adjust) & is.na(gencode_tran1a$n3_adjust)))/nrow(gencode_tran1a)
#0.1526797 exactly the same
#NOW include internal primed 3' end in the adjusted 3' end, n=5749, these should be removed (table 0 -> table 1 (raw))
#
gene_model=read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/table0_gene/iPSC_NSC_Neuron.S3.disable_ref_chain_bound_gene_anno_10percent/bed/iPSC_NSC_Neuron.S3.disable_ref_chain_bound_gene_anno_10percent.gene.bed.bgz", header=F, stringsAsFactors = F, check.names = F)
gene_model$V4=sapply(strsplit(gene_model$V4,"\\|"),"[",1)
gencode_gene=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.gene.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
gencode_gene=left_join(gencode_gene[,c(1:6)], gene_model[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new")) #61533
gencode_gene=gencode_gene[which(!is.na(gencode_gene$V1_new)),] #61422, 111 ENSG missing
gencode_gene$n5_adjust="5n_adjust"
gencode_gene$n3_adjust="3n_adjust"
gencode_gene$n5_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="+"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="-"))]=NA
gencode_gene$n3_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="-"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="+"))]=NA
gencode_gene2=melt(gencode_gene[,c(4,12,13)], id=1)
gencode_gene2=gencode_gene2[which(!is.na(gencode_gene2$value)),]
gencode_gene2=gencode_gene2%>%group_by(V4)%>%dplyr::summarise(ENSG_adjust=paste(value,collapse=" & "))
gencode_gene=left_join(gencode_gene,gencode_gene2, by="V4", copy=F)
gencode_gene$ENSG_adjust[is.na(gencode_gene$ENSG_adjust)]="No"
transcript_info=left_join(transcript_info,gencode_gene[,c(4,14)], by=c("IN1_gene_ID" = "V4"), copy=F)

#===============================================================================
#gene level novel geneclass, temporary for raw table
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

#===============================================================================
#apply different filters
#pre-filter: internal priming gone
transcript_info1a=transcript_info[which(transcript_info$source != "non_detectable_gencode"),]
transcript_info1c=transcript_info1a[which(transcript_info1a$internal_priming != "Yes" ),]
transcript_info1d=transcript_info1c[which(transcript_info1c$source != "partial_gencode"),]

write.table(transcript_info1c, gzfile("table1.remove_undetect_gencode_and_internal_prime.2.51M.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info1b=transcript_info[which(transcript_info$source == "non_detectable_gencode"),]
write.table(transcript_info1b, gzfile("table1.only_undetect_gencode.151K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

transcript_info1b1=transcript_info[which(transcript_info$source == "partial_gencode"),]
write.table(transcript_info1b1, gzfile("table1.only_partial_gencode.34K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info1b1=read.delim("table1.only_partial_gencode.34K.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

transcript_info1b2=transcript_info1a[which(transcript_info1a$internal_priming == "Yes"),]
write.table(transcript_info1b2, gzfile("table1.only_internal_priming.430K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info1b2=read.delim("table1.only_internal_priming.430K.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

#table2, standard similar to TALON, 5 reads support
transcript_info2=transcript_info1c[which(transcript_info1c$source != "permissive"),]
write.table(transcript_info2, gzfile("table2.standard.153K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info2=read.delim("table2.standard.153K.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
transcript_info2%>%group_by(transcript_novelty,n5_support)%>%dplyr::summarise(count=n())
transcript_info2%>%group_by(gene_novelty,transcript_novelty)%>%dplyr::summarise(count=n())

#table4 for downstream analysis, all include SCAFE supported 5' end
transcript_info4a=transcript_info1d[which(transcript_info1d$transcript_novelty == "GENCODE" & transcript_info1d$gene_novelty == "GENCODE"),] #GENCODE with SCAFE support
transcript_info4a=transcript_info4a[grep("SCAFE", transcript_info4a$n5_support),]
transcript_info4b=transcript_info1d[intersect(grep("SCAFE", transcript_info1d$n5_support),which(transcript_info1d$transcript_novelty == "novel" & transcript_info1d$gene_novelty == "GENCODE" & transcript_info1d$source=="standard")),] 
transcript_info4b=transcript_info4b[which(transcript_info4b$T_ratio_iPS >= 0.1 | transcript_info4b$T_ratio_NSC >= 0.1 | transcript_info4b$T_ratio_Neuron >=0.1 | transcript_info4b$T_ratio_THP1_DMSO >= 0.1 | transcript_info4b$T_ratio_THP1_PMA >= 0.1),] #novel isoform
transcript_info4c=transcript_info1d[intersect(grep("SCAFE", transcript_info1d$n5_support),which(transcript_info1d$transcript_novelty == "novel" & transcript_info1d$gene_novelty == "novel")),] #novel gene
transcript_info4=rbind(transcript_info4a,transcript_info4b,transcript_info4c)
write.table(transcript_info4, gzfile("table4.chimeric.199K.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info4%>%group_by(gene_novelty,transcript_novelty,transcript_group,CPAT_class)%>%dplyr::summarise(count=n())

#bed12 for each table
table0.bed12=read.delim(paste0(path1,"Neuron_THP1.table0.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table1.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1c$model_ID),]
table2.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info2$model_ID),]
table4.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info4$model_ID),]
write.table(table1.bed12, gzfile(paste0(path1,"Neuron_THP1.table1.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table2.bed12, gzfile(paste0(path1,"Neuron_THP1.table2.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(table4.bed12, gzfile(paste0(path1,"Neuron_THP1.table4.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare input file for gene annotator using table 4 info
path6=paste0(SALA_path,"table4_gene/")
setwd(path6)
table0.bed12=read.delim(paste0(path3, "Neuron_THP1.S3.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
table0.bed12a=table0.bed12[union(grep("ENST",table0.bed12$V4),which(table0.bed12$V4 %in% transcript_info4$model_ID)),]
write.table(table0.bed12a[order(table0.bed12a$V1,table0.bed12a$V2),],gzfile("table4wENST.bed12.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#bash
setwd(path6)
system("zcat table4wENST.bed12.bed.gz | bgzip > table4wENST.bed12.bed.bgz")
system("tabix -p bed table4wENST.bed12.bed.bgz")
system("rm table4wENST.bed12.bed.gz")

#===============================================================================
table0=read.delim(paste0(path2,"Neuron_THP1.S3.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table0a=table0[which(table0$model_ID %in% table0.bed12a$V4),]
write.table(table0a, gzfile("table4wENST.info.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#=====================
#run gene annotator
system(paste0("sh ",SALA_path,"gene4.sh"))

#===============================================================================
#Re-run all gene related properties
transcript_info4=read.delim(paste0(path2,"table4.chimeric.199K.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t4gene=read.delim(paste0(gene4_path,"log/Neuron_THP1_T4_10percent.model.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
t4gene%>%group_by(g_assign)%>%dplyr::summarise(count=n())
colnames(t4gene)[c(7,8)]=c("T4_gene_ID","T4_gene_name")
length(unique(t4gene$T4_gene_ID))
transcript_info4=left_join(transcript_info4,t4gene[,c(1,7,8)], by="model_ID", copy=F)
transcript_info4$T4_gene_novelty="novel"
transcript_info4$T4_gene_novelty[grep("ENSG",transcript_info4$T4_gene_ID)]="GENCODE"

#GENCODE geneClass
gencode.anno=read.delim(paste0(path2,"GENCODE.gene.n.transcript.class.1base.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(gencode.anno)[3]="T4_Gencode_geneCalss"
transcript_info4=left_join(transcript_info4, unique(gencode.anno[,c(1,3)]), by=c("T4_gene_ID"="geneID"),copy=F)
transcript_info4$T4_Gencode_geneCalss2="NA"
transcript_info4$T4_Gencode_geneCalss2[which(!is.na(transcript_info4$T4_Gencode_geneCalss))]="others"
transcript_info4$T4_Gencode_geneCalss2[which(transcript_info4$T4_Gencode_geneCalss == "lncRNA")]="lncRNA"
transcript_info4$T4_Gencode_geneCalss2[which(transcript_info4$T4_Gencode_geneCalss == "protein_coding")]="protein_coding"

##ENSG with re-adjusted 3/5'end
gene_model=read.delim(paste0(gene4_path,"bed/Neuron_THP1_T4_10percent.gene.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
gene_model$V4=sapply(strsplit(gene_model$V4,"\\|"),"[",1)
gencode_gene=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.gene.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
gencode_gene=left_join(gencode_gene[,c(1:6)], gene_model[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new")) #61533
gencode_gene=gencode_gene[which(!is.na(gencode_gene$V1_new)),] #61420, 113 ENSG missing
gencode_gene$n5_adjust="5n_adjust"
gencode_gene$n3_adjust="3n_adjust"
gencode_gene$n5_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="+"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="-"))]=NA
gencode_gene$n3_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="-"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="+"))]=NA
gencode_gene2=melt(gencode_gene[,c(4,12,13)], id=1)
gencode_gene2=gencode_gene2[which(!is.na(gencode_gene2$value)),]
gencode_gene2=gencode_gene2%>%group_by(V4)%>%dplyr::summarise(T4_ENSG_adjust=paste(value,collapse=" & "))
gencode_gene=left_join(gencode_gene,gencode_gene2, by="V4", copy=F)
gencode_gene$T4_ENSG_adjust[is.na(gencode_gene$T4_ENSG_adjust)]="No"
transcript_info4=left_join(transcript_info4,gencode_gene[,c(4,14)], by=c("T4_gene_ID" = "V4"), copy=F)

#===============================================================================
#gene level novel geneclass
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_overall=full_qry_count/sum(full_qry_count))
transcript_infoa=transcript_info4[which(transcript_info4$Novel_transcriptClass %in% c("ncRNA","lncRNA") & transcript_info4$T4_gene_novelty == "novel"),]
transcript_infoa=transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_overall2=full_qry_count/sum(full_qry_count))
transcript_infob=transcript_infoa%>%group_by(T4_gene_ID)%>%dplyr::summarise(ncRNA_rate=sum(T4_T_ratio_overall), weighted_average_length=sum(T4_T_ratio_overall2*transcript_length))
transcript_infob$T4_Novel_geneClass="others"
transcript_infob$T4_Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5)]="ncRNA"
transcript_infob$T4_Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5 & transcript_infob$weighted_average_length > 200)]="lncRNA"
transcript_infob%>%group_by(T4_Novel_geneClass)%>%dplyr::summarise(count=n())
transcript_info4=left_join(transcript_info4,transcript_infob[,c(1,4)], by="T4_gene_ID", copy=F)
transcript_info4$T4_Novel_geneClass[which(is.na(transcript_info4$T4_Novel_geneClass) & transcript_info4$T4_gene_novelty == "novel")]="others"
transcript_info4%>%group_by(T4_Novel_geneClass, Novel_transcriptClass)%>%dplyr::summarise(count=n())

write.table(transcript_info4, gzfile(paste0(path2,"table4.chimeric.199K.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
transcript_info4=read.delim(paste0(path2,"table4.chimeric.199K.tsv.gz"), header=T, stringsAsFactor=F, check.names=F)
transcript_info4%>%group_by(T4_gene_novelty,transcript_novelty,transcript_group,CPAT_class)%>%dplyr::summarise(count=n())
transcript_info4%>%group_by(T4_gene_novelty,transcript_group)%>%dplyr::summarise(count=n())
transcript_info4%>%group_by(T4_gene_novelty,ncRNA_subclass)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

#===============================================================================
#re-run transcript ratio to exclude novel transcript model newly assigned to ENSG
length(which(transcript_info4$gene_novelty=="novel" & transcript_info4$T4_gene_novelty == "GENCODE"))
#4626
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_iPS=iPS/sum(iPS))
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_NSC=NSC/sum(NSC))
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_Neuron=Neuron/sum(Neuron))
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_THP1_DMSO=THP1_DMSO/sum(THP1_DMSO))
transcript_info4=transcript_info4%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_T_ratio_THP1_PMA=THP1_PMA/sum(THP1_PMA))
transcript_info4$T4_T_ratio_iPS[which(transcript_info4$T4_T_ratio_iPS == "NaN")]=0
transcript_info4$T4_T_ratio_NSC[which(transcript_info4$T4_T_ratio_NSC == "NaN")]=0
transcript_info4$T4_T_ratio_Neuron[which(transcript_info4$T4_T_ratio_Neuron == "NaN")]=0
transcript_info4$T4_T_ratio_THP1_DMSO[which(transcript_info4$T4_T_ratio_THP1_DMSO == "NaN")]=0
transcript_info4$T4_T_ratio_THP1_PMA[which(transcript_info4$T4_T_ratio_THP1_PMA == "NaN")]=0

length(which(transcript_info4$T4_gene_novelty == "GENCODE" & transcript_info4$transcript_novelty == "novel" & transcript_info4$T4_T_ratio_iPS <0.1 & transcript_info4$T4_T_ratio_NSC <0.1 & transcript_info4$T4_T_ratio_Neuron <0.1 & transcript_info4$T4_T_ratio_THP1_DMSO <0.1 & transcript_info4$T4_T_ratio_THP1_PMA <0.1))
#4183
length(which(transcript_info4$gene_novelty=="novel" & transcript_info4$T4_gene_novelty == "GENCODE" & transcript_info4$transcript_novelty == "novel" & transcript_info4$T4_T_ratio_iPS <0.1 & transcript_info4$T4_T_ratio_NSC <0.1 & transcript_info4$T4_T_ratio_Neuron <0.1 & transcript_info4$T4_T_ratio_THP1_DMSO <0.1 & transcript_info4$T4_T_ratio_THP1_PMA <0.1))
#4087
length(which(transcript_info4$gene_novelty=="novel" & transcript_info4$T4_gene_novelty == "GENCODE" & transcript_info4$T4_T_ratio_iPS <0.1 & transcript_info4$T4_T_ratio_NSC <0.1 & transcript_info4$T4_T_ratio_Neuron <0.1 & transcript_info4$T4_T_ratio_THP1_DMSO <0.1 & transcript_info4$T4_T_ratio_THP1_PMA <0.1 & transcript_info4$CPAT_class == "coding"))
#212
#remove newly assigned ENSG with 10% & 5 reads cutoff
table4ENSG=transcript_info4[which(transcript_info4$T4_gene_novelty == "GENCODE" & transcript_info4$transcript_novelty == "novel"),]
table4ENSG1=table4ENSG[which(table4ENSG$T4_T_ratio_iPS >= 0.1 | table4ENSG$T4_T_ratio_NSC >= 0.1 | table4ENSG$T4_T_ratio_Neuron >= 0.1 | table4ENSG$T4_T_ratio_THP1_DMSO >=0.1 | table4ENSG$T4_T_ratio_THP1_PMA >=0.1),]
table4ENSG1=table4ENSG1[which(table4ENSG1$source != "permissive"),] 
table4off=setdiff(table4ENSG$model_ID,table4ENSG1$model_ID)

table4a=transcript_info4[which(transcript_info4$model_ID %in% table4off),]
table5=transcript_info4[-which(transcript_info4$model_ID %in% table4off),]

#remove isoforms without confident ex3_cluster
table4b=table5[which(table5$T4_gene_novelty == "GENCODE" & table5$transcript_novelty == "Novel isoform" & table5$n3_support == "no_support"),]
table5=table5[-which(table5$model_ID %in% table4b$model_ID),]

table4a=rbind(table4a, table4b)
write.table(table4a, gzfile(paste0(path2,"table4.removed.4734.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#gene-base promoter_type
table5$promoter_type[which(table5$promoter_type=="unclassed" & table5$ATAC == "noATAC")]="excluded"
table5$promoter_type[which(table5$promoter_type=="CTCF-alone" & table5$ATAC == "noATAC")]="excluded"
table5=table5%>%group_by(T4_gene_ID)%>%dplyr::mutate(T4_gene_promoter_type=paste(unique(promoter_type),collapse=";"))
#promoter>enhancer>unclass
table5$T4_gene_promoter_type[grep("promoter",table5$T4_gene_promoter_type)]="promoter-like"
table5$T4_gene_promoter_type[grep("enhancer",table5$T4_gene_promoter_type)]="enhancer-like"
table5$T4_gene_promoter_type[grep("CTCF",table5$T4_gene_promoter_type)]="CTCF-alone"
table5$T4_gene_promoter_type[grep("unclassed",table5$T4_gene_promoter_type)]="unclassed"
unique(table5[,c(86,99)])%>%group_by(T4_gene_promoter_type)%>%dplyr::summarise(count=n())

write.table(table5,gzfile(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

table0.bed12=read.delim(paste0(path1,"Neuron_THP1.table0.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.bed12=table0.bed12[which(table0.bed12$V4 %in% table5$model_ID),]
write.table(table5.bed12, gzfile(paste0(path1,"Neuron_THP1.table5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#============================
#gene
colnames(table5)[c(89,90)]=c("T4_Gencode_geneClass","T4_Gencode_geneClass2")

table5gene=unique(table5[,c("T4_gene_ID","T4_gene_name","T4_gene_novelty","T4_Gencode_geneClass","T4_Gencode_geneClass2","T4_ENSG_adjust","T4_Novel_geneClass","T4_gene_promoter_type")])
table5gene$geneClass=paste0(table5gene$T4_gene_novelty,":",table5gene$T4_Novel_geneClass)
table5gene$geneClass[which(table5gene$T4_gene_novelty == "GENCODE")]=paste0(table5gene$T4_gene_novelty[which(table5gene$T4_gene_novelty == "GENCODE")],":",table5gene$T4_Gencode_geneClass2[which(table5gene$T4_gene_novelty == "GENCODE")])
table5gene%>%group_by(geneClass, T4_gene_promoter_type)%>%dplyr::summarise(count=n())
write.table(table5gene, paste0(path2,"table5.gene.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#add 3'end frequency
table5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=unique(table5[,c(1:105)])
table5.bed12=read.delim(paste0(path1,"Neuron_THP1.table5.bed12.bed.gz"), header=F,  stringsAsFactors = F)
table5.bed12$V2[which(table5.bed12$V6=="+")]=table5.bed12$V3[which(table5.bed12$V6=="+")]-1
table5.bed12$V3[which(table5.bed12$V6=="-")]=table5.bed12$V2[which(table5.bed12$V6=="-")]+1
table5.bed12$TESID=paste0(table5.bed12$V1,"_",table5.bed12$V2,"_",table5.bed12$V3,"_",table5.bed12$V6)
data3=read.delim(paste0(path2,"table0.TESbymodelID.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data3=data3%>%group_by(TES)%>%dplyr::summarise(TEScount=sum(TEScount))
table5.bed12=left_join(table5.bed12[,c(4,13)],data3, by=c("TESID"="TES"),copy=F)
table5.bed12$TES_recur="No"
table5.bed12$TES_recur[which(table5.bed12$TEScount >1)]="Yes"
transcript_model=read.delim(paste0(path2,"table0.TESbymodelID_modelonly.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
transcript_model$polyA="poly(A)"
transcript_model$polyA[which(transcript_model$polyA_prediction == "No" & transcript_model$PAS2 == "No" & transcript_model$GENCODE_polyA == "No")]="non-poly(A)"

table5.bed12=left_join(table5.bed12, transcript_model[,c(1,2,4,11,12,15,16)], by=c("TESID"="TES", "V4"="V4"), copy=F)
table5=left_join(table5, table5.bed12, by=c("model_ID"="V4"), copy=F)
table5%>%group_by(polyA,TES_recur)%>%dplyr::summarise(count=n())
rm(data3)
rm(transcript_model)
write.table(table5,paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
# RNA class for SALA Final
table5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5$T4_Novel_geneClass[which(table5$T4_Novel_geneClass == "ncRNA")]="short_ncRNA"
table5$Novel_transcriptClass[which(table5$Novel_transcriptClass == "ncRNA")]="short_ncRNA"

#simple RNA class both gene and transcript (ncRNA / others)
table5$T4_Novel_geneClass=paste0("Novel:",table5$T4_Novel_geneClass)
table5$T4_Gencode_geneClass2=paste0("GENCODE:",table5$T4_Gencode_geneClass2)
table5$T4_Novel_geneClass[which(table5$T4_Novel_geneClass == "Novel:NA")]=table5$T4_Gencode_geneClass2[which(table5$T4_Novel_geneClass == "Novel:NA")]
sum1=unique(table5[c(86,93)])%>%group_by(T4_Novel_geneClass)%>%dplyr::summarise(G_count=n())
#table5$Novel_transcriptClass[grep("ncRNA",table5$Novel_transcriptClass)]="non-coding"
table5$Novel_transcriptClass=paste0("Novel:",table5$Novel_transcriptClass)
table5$Gencode_transcriptClass2=paste0("GENCODE:",table5$Gencode_transcriptClass2)
table5$Novel_transcriptClass[which(table5$Novel_transcriptClass == "Novel:NA")]=table5$Gencode_transcriptClass2[which(table5$Novel_transcriptClass == "Novel:NA")]
sum2=table5%>%group_by(Novel_transcriptClass)%>%dplyr::summarise(T_count=n())
sum1$percent=sum1$G_count/sum(sum1$G_count)
sum2$percent=sum2$T_count/sum(sum2$T_count)
sum1$group="Gene"
sum2$group="Transcript"
colnames(sum1)[c(1,2)]=c("class","count")
colnames(sum2)[c(1,2)]=c("class","count")
sum1=rbind(sum1, sum2)
sum1[sum1=="Novel:others"]="Novel:potential_coding"
write.table(sum1, gzfile(paste0(path_fig3_data,"RNA_class_summary.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)


sum3=table5%>%group_by(T4_Novel_geneClass,Novel_transcriptClass)%>%dplyr::summarise(T_count=n())
sum3[sum3=="Novel:others"]="Novel:potential_coding"
sum3=left_join(sum3, sum1[c(1:6),c(1:2)], by=c("T4_Novel_geneClass"="class"), copy=F)
sum3a=sum3[which(sum3$T4_Novel_geneClass %in% c("GENCODE:lncRNA","Novel:lncRNA","Novel:short_ncRNA")),]
sum3a$include=c("Yes","No","Yes","No","No","Yes","No","Yes","Yes","No","Yes")
write.table(sum3a,gzfile(paste0(path_fig3_data,"ncRNA.class.plot.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#sub-classification: identify "ncRNA subclass" for SALA Final 
setwd(path2)
table5=read.delim("table5.chimeric.194K.remove.permissive.isoform.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

table5a=table5[which(table5$T4_Gencode_geneClass2=="lncRNA" & table5$Gencode_transcriptClass2 =="lncRNA"),]
table5b=table5[which(table5$T4_Gencode_geneClass2=="lncRNA" & table5$Novel_transcriptClass %in% c("lncRNA")),]
table5c=table5[which(table5$T4_gene_novelty == "novel" & table5$T4_Novel_geneClass %in% c("lncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","short_ncRNA")),]
table5d=table5[which(table5$T4_gene_novelty == "novel" & table5$T4_Novel_geneClass %in% c("short_ncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","short_ncRNA")),]
table5x=table5[which(table5$T4_Gencode_geneClass2=="protein_coding" & table5$Novel_transcriptClass =="lncRNA"),c(1,71:112)]
table5a$T4_ncRNA_source="ENST_lncRNA"
table5b$T4_ncRNA_source="novel_lncRNA_isoform"
table5c$T4_ncRNA_source="novel_lncRNA"
table5d$T4_ncRNA_source="novel_ncRNA"
table7=rbind(table5a,table5b,table5c,table5d)

options(scipen=999)
table0.bed12=read.delim(paste0(path1,"Neuron_THP1.table0.bed12.bed.gz"),header=F, stringsAsFactors = F, check.names = F)
table7.bed12=table0.bed12[which(table0.bed12$V4 %in% table7$model_ID),]
write.table(table7.bed12[order(table7.bed12$V1,table7.bed12$V2),],gzfile(paste0(path1,"table7_ncRNA.bed12.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
#prepare n5
table7.bed12$V3[which(table7.bed12$V6 == "+")]=table7.bed12$V2[which(table7.bed12$V6 == "+")]+1
table7.bed12$V2[which(table7.bed12$V6 == "-")]=table7.bed12$V3[which(table7.bed12$V6 == "-")]-1
write.table(table7.bed12[order(table7.bed12$V1, table7.bed12$V2),c(1:6)], gzfile(paste0(path1,"table7_ncRNA.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#prepare #ENST with updated range
options(scipen=999)
ENSTbed12=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.transcript.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.bed12=read.delim(paste0(path1,"Neuron_THP1.table5.bed12.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.bed12a=table5.bed12[grep("ENST",table5.bed12$V4),]
ENSTbed12a=ENSTbed12[-which(ENSTbed12$V4 %in% table5.bed12a$V4),]
ENSTbed12=rbind(ENSTbed12a,table5.bed12a)
write.table(ENSTbed12[order(ENSTbed12$V1, ENSTbed12$V2),], gzfile(paste0(path1,"ENST.T4updated.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#prepare mRNA & pseudo alone
gencode.anno=read.delim(paste0(path2,"GENCODE.gene.n.transcript.class.1base.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gencode.anno1=gencode.anno[which(gencode.anno$Gencode_transcriptClass == "protein_coding"),]
gencode.anno2=gencode.anno[grep("pseudogene", gencode.anno$Gencode_transcriptClass),]
ENSTbed12a=ENSTbed12[which(ENSTbed12$V4 %in% c(gencode.anno1$transcriptID, gencode.anno2$transcriptID)),]
write.table(ENSTbed12a[order(ENSTbed12a$V1, ENSTbed12a$V2),], gzfile(paste0(path1,"ENST.T4updated.mRNA_pseudo.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#prepare n5
ENSTbed12a$V3[which(ENSTbed12a$V6 == "+")]=ENSTbed12a$V2[which(ENSTbed12a$V6 == "+")]+1
ENSTbed12a$V2[which(ENSTbed12a$V6 == "-")]=ENSTbed12a$V3[which(ENSTbed12a$V6 == "-")]-1
write.table(ENSTbed12a[order(ENSTbed12a$V1, ENSTbed12a$V2),c(1:6)], gzfile(paste0(path1,"ENST.T4updated.mRNA_pseudo.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===================
#bedtools
setwd(path1)
system("bed12ToBed6 -i ENST.T4updated.mRNA_pseudo.bed12.bed.gz | gzip > ENST.T4updated.mRNA_pseudo.bed6.bed.gz")
system("zcat ENST.T4updated.mRNA_pseudo.bed12.bed.gz | bedparse introns | bed12ToBed6 -i stdin | gzip > ENST.T4updated.mRNA_pseudo.intron.bed6.bed.gz")
system("bedtools subtract -s -a ENST.T4updated.mRNA_pseudo.intron.bed6.bed.gz -b ENST.T4updated.mRNA_pseudo.bed6.bed.gz | gzip > ENST.T4updated.mRNA_pseudo.pure.intron.bed6.bed.gz")
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

exonExon=left_join(exonExon,table7[,c(1,57)],by=c("V4"="model_ID"),copy=F)
exonExon$overlap_rate=exonExon$overlap/exonExon$transcript_length
exonExonAs=read.delim("table7_ncRNA.bed6.antisense.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
exonExonAs=exonExonAs%>%group_by(V4,V10)%>%dplyr::summarise(overlap=sum(V13))
exonExonAs=left_join(exonExonAs,table7[,c(1,57)],by=c("V4"="model_ID"),copy=F)
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
write.table(table7[,c(1,100:110)],gzfile(paste0(path2,"table7.ncRNA.121K.subclass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
table5=left_join(table5, table7[,c(1,100,110)], by="model_ID", copy=F)
write.table(table5,gzfile(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
table7s=table7[,c(1,110,86:109)]
table7s1=table7s[which(table7s$T4_ncRNA_subclass == "Sense_overlap_RNA"),]

#gene level ncRNA subclass
#order -> divergent lncRNAs, sense overlap lncRNAs, sense intronic lncRNAs, antisense lncRNAs and intergenic lncRNAs 
table7=read.delim(paste0(path2,"table7.ncRNA.99K.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
unique(table7$T4_ncRNA_subclass)
table7$order=0
table7$order[which(table7$T4_ncRNA_subclass == "Intergenic")]=1
table7$order[which(table7$T4_ncRNA_subclass == "Antisense_other")]=2
table7$order[which(table7$T4_ncRNA_subclass == "Antisense_intronic")]=3
table7$order[which(table7$T4_ncRNA_subclass == "Antisense")]=4
table7$order[which(table7$T4_ncRNA_subclass == "Sense_RNA_other")]=5
table7$order[which(table7$T4_ncRNA_subclass == "Sense_intronic")]=6
table7$order[which(table7$T4_ncRNA_subclass == "Sense_overlap_RNA")]=7
table7$order[which(table7$T4_ncRNA_subclass == "Divergent")]=8
gene_data=unique(table7[,c(86,110,111)])%>%group_by(T4_gene_ID)%>%slice_max(order)
colnames(gene_data)[2]="T4_gene_ncRNA_subclass"
table7=left_join(table7, gene_data[,c(1,2)], by="T4_gene_ID",copy=F)
table5=left_join(table5, gene_data[,c(1,2)], by="T4_gene_ID", copy=F)
write.table(table5,gzfile(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(colnames(table5), "readme.tsv", row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#prepare transcript base 5' and 3' bed file
t1.bed12=read.delim(paste0(path1,"Neuron_THP1.table1.bed12.bed.gz"), header=F, stringsAsFactors = F)
#t2.bed12=read.delim(paste0(path1,"Neuron_THP1.table2.bed12.bed.gz"), header=F, stringsAsFactors = F)
t5.bed12=read.delim(paste0(path1,"Neuron_THP1.table5.bed12.bed.gz"), header=F, stringsAsFactors = F)

#count n5, n3 number
t1.bed12$n5=t1.bed12$V2
t1.bed12$n5[which(t1.bed12$V6=="-")]=t1.bed12$V3[which(t1.bed12$V6=="-")]
t1.bed12$n3=t1.bed12$V3
t1.bed12$n3[which(t1.bed12$V6=="-")]=t1.bed12$V2[which(t1.bed12$V6=="-")]
t1.bed12$n5_1=t1.bed12$n5+1
t1.bed12$n3_1=t1.bed12$n3+1
t1.bed12=left_join(t1.bed12, SALAt1[,c(1,86)], by=c("V4"="model_ID"),copy=F)
nrow(unique(t1.bed12[,c(1,13)])) #1269925
nrow(unique(t1.bed12[,c(1,14)])) #1087801
#t2.bed12$n5=t2.bed12$V2
#t2.bed12$n5[which(t2.bed12$V6=="-")]=t2.bed12$V3[which(t2.bed12$V6=="-")]
#t2.bed12$n3=t2.bed12$V3
#t2.bed12$n3[which(t2.bed12$V6=="-")]=t2.bed12$V2[which(t2.bed12$V6=="-")]
#t2.bed12$n5_1=t2.bed12$n5+1
#t2.bed12$n3_1=t2.bed12$n3+1
#t2.bed12=left_join(t2.bed12, SALAt1[,c(1,86)], by=c("V4"="model_ID"),copy=F)
#nrow(unique(t2.bed12[,c(1,13)])) #84519
#nrow(unique(t2.bed12[,c(1,14)])) #80147
t5.bed12$n5=t5.bed12$V2
t5.bed12$n5[which(t5.bed12$V6=="-")]=t5.bed12$V3[which(t5.bed12$V6=="-")]
t5.bed12$n3=t5.bed12$V3
t5.bed12$n3[which(t5.bed12$V6=="-")]=t5.bed12$V2[which(t5.bed12$V6=="-")]
t5.bed12$n5_1=t5.bed12$n5+1
t5.bed12$n3_1=t5.bed12$n3+1
t5.bed12=left_join(t5.bed12, SALAt5[,c(1,86)], by=c("V4"="model_ID"),copy=F)
nrow(unique(t5.bed12[,c(1,13)])) #150334
nrow(unique(t5.bed12[,c(1,14)])) #147529

write.table(t1.bed12[order(t1.bed12$V1,t1.bed12$n5),c(1,13,15,4,5,6,17)],gzfile(paste0(path1,"SALA_Raw.table1.N5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#write.table(t2.bed12[order(t2.bed12$V1,t2.bed12$n5),c(1,13,15,4,5,6,17)],gzfile(paste0(path1,"SALA_Read_filtered.table2.N5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(t5.bed12[order(t5.bed12$V1,t5.bed12$n5),c(1,13,15,4,5,6,17)],gzfile(paste0(path1,"SALA_Final.table5.N5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(t1.bed12[order(t1.bed12$V1,t1.bed12$n3),c(1,14,16,4,5,6,17)],gzfile(paste0(path1,"SALA_Raw.table1.N3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#write.table(t2.bed12[order(t2.bed12$V1,t2.bed12$n3),c(1,14,16,4,5,6,17)],gzfile(paste0(path1,"SALA_Read_filtered.table2.N3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(t5.bed12[order(t5.bed12$V1,t5.bed12$n3),c(1,14,16,4,5,6,17)],gzfile(paste0(path1,"SALA_Final.table5.N3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===================================
#intersect with external and data-driven 5' features
n5_path=paste0(primary_folder,"code_n_data/n5_regions/")
cluster_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")
setwd(path1)
cmd <- paste0(
  "for file in *.N5.bed.gz; do ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-ELS.all.enhancer.sort.bed.gz | gzip > \"${file%.bed.gz}.Ecount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-PLS.all.promoter.sort.bed.gz | gzip > \"${file%.bed.gz}.Pcount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"GRCh38-CTCF.sort.bed.gz | gzip > \"${file%.bed.gz}.Ccount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"F6_CAT.promoter.bed.gz | gzip > \"${file%.bed.gz}.CAGEcount.bed.gz\"; ",
  "bedtools intersect -c -a \"$file\" -b ",n5_path,"peaks.merged.all.main_chr.bed.gz | gzip > \"${file%.bed.gz}.ATACcount.bed.gz\" ;",
  "bedtools intersect -c -s -a \"$file\" -b ",cluster_path,"ontCAGE.Neuron_THP1.cluster.coord.bed.gz | gzip > \"${file%.bed.gz}.SCAFEcount.bed.gz\" ;",
  "done")
system(cmd)

#===============================================================================
#generate supp table
path_supp=paste0(primary_folder,"supplementary_table/")
SCAFE_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")

table5=read.delim(paste0(path2,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
CREanno=read.delim(paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1.CRE.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=left_join(table5[,-67], CREanno[,c(1,35)], by="CREID",copy=F)
table5=table5[,c(1:66,126,67:125)]
all_pas1=read.delim(paste0(path_supp,"TableS1_polyA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=left_join(table5[,-c(111,112,113)], all_pas1[,c(7,4,6,8:10,11)], by="TESID", copy=F)

table5$polyA="No"
table5$polyA[which(table5$polyA_predict_class=="significant poly(A)" | table5$PAS_motif_class == "Yes")]="Yes"
table5%>%group_by(polyA_predict_class,PAS_motif_class,polyA)%>%dplyr::summarise(count=n())
table5a=table5[,c(1:18,41:43,45:50,56:73,76:79,86:102,106:130)]
colnames(table5a)[c(26,27,53,54,61,62,66)]=c("THP1","dTHP1","T4_Gencode_geneClass","T4_Gencode_geneClass2","T_ratio_THP1","T_ratio_dTHP1","transcript_group")
colnames(table5a)=gsub("T4_","",colnames(table5a))
write.table(table5a,paste0(path_supp,"TableS4_finalTx.tsv"), col.names=T, row.names=F, sep="\t",quote=F)

table1=read.delim(paste0(path2,"table1.remove_undetect_gencode_and_internal_prime.2.51M.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
table1=left_join(table1[,-67], CREanno[,c(1,35)], by="CREID",copy=F)
table1=table1[,c(1:66,88,67:87)]
table1a=table1[,c(1:18,41:43,45:50,56:73,76:79)]
colnames(table1a)[c(26,27)]=c("THP1","dTHP1")
write.table(table1a,gzfile(paste0(path_supp,"TableS5_rawTx.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)






