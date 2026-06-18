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
path_fig3_data=paste0(primary_folder,"fig3/data/")
short_read_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/kallisto_short_t5_partialYes.ENST/")
long_read_bambu=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/bambu_long_t5_partialYes.ENST/")
SQANTI_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SQANTI/")

#===============================================================================
#quantification from short and long read
shortkallisto=read.delim(paste0(short_read_path,"Neuron.series.transcript.tpm.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
bambu_tcpm_new=read.delim(paste0(long_read_bambu,"CPM_transcript.txt.gz"), header=T, stringsAsFactors = F, check.names = F)

bambu_tcpm_new=bambu_tcpm_new[,-2]
bambu_tcpm_new[is.na(bambu_tcpm_new)]=0
bambu_tcpm_new$iPSC=rowMeans(bambu_tcpm_new[,c(2:3)])
bambu_tcpm_new$NSC=rowMeans(bambu_tcpm_new[,c(6:7)])
bambu_tcpm_new$Neuron=rowMeans(bambu_tcpm_new[,c(4:5)])
bambu_tcpm_new1=reshape2::melt(bambu_tcpm_new[,c(1,24:26)], id=1)
bambu_tcpm_new1=bambu_tcpm_new1[which(bambu_tcpm_new1$TXNAME %in% shortkallisto$target_id),]

shortkallisto$iPSC=rowMeans(shortkallisto[,c(2:3)])
shortkallisto$NSC=rowMeans(shortkallisto[,c(6:7)])
shortkallisto$Neuron=rowMeans(shortkallisto[,c(4:5)])
shortkallisto1=reshape2::melt(shortkallisto[,c(1,8:10)], id=1)
together=left_join(shortkallisto1, bambu_tcpm_new1, by=c("target_id"="TXNAME", "variable"="variable"),copy=F, suffix=c("_short","_long"))
colnames(together)[c(2,3,4)]=c("Cell","Kallisto","Bambu")
t_to_g=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.transcript_gene_link.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
together=left_join(together, t_to_g, by=c("target_id"="transcriptID"))
gene=together%>%group_by(geneID,Cell)%>%dplyr::summarise(Kallisto=sum(Kallisto), Bambu=sum(Bambu))

#add length and bin
gtf=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.gtf.gz"), header=F, stringsAsFactors = F)
gtf$length=gtf$V5-gtf$V4+1
gtfe=gtf[which(gtf$V3=="exon"),]
gtfe$transcript=sapply(strsplit(gtfe$V9,"; "),"[",2)
gtfe$transcript=gsub("transcript_id ","",gtfe$transcript)
gtfe1=gtfe%>%group_by(transcript)%>%dplyr::summarise(length=sum(length))
together=left_join(together, gtfe1, by=c("target_id"="transcript"),copy=F)
together$length_bin="<=1kb"
together$length_bin[which(together$length>1000)]="1-2kb"
together$length_bin[which(together$length>2000)]="2-3kb"
together$length_bin[which(together$length>3000)]="3-4kb"
together$length_bin[which(together$length>4000)]=">4kb"

write.table(together,gzfile(paste0(path_fig3_data,"short_long_transcript_TPM.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(gene,gzfile(paste0(path_fig3_data,"short_long_gene_TPM.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5g & h, 3e

#===============================================================================
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
quant_count1=read.delim(paste0(short_read_path, "Neuron.series.gene.count.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
quant_count=read.delim(paste0(short_read_path, "Neuron.series.transcript.count.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
quant_count1$group="GENCODE"
quant_count1$group[grep("ONT",quant_count1$geneID)]="Novel"
quant_count1$detection="Yes"
quant_count1$detection[which(rowSums(quant_count1[,c(2:7)])==0)]="No"
quant_count1=left_join(quant_count1, unique(table5[,c(86,99,110)]), by=c("geneID"="T4_gene_ID"),copy=F)
write.table(quant_count1,gzfile(paste0(path_fig3_data,"Kallisto_Neuron.series.gene.count.matrix.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
aa1=quant_count1[which(!is.na(quant_count1$T4_gene_promoter_type)),]%>%group_by(group,detection)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
bb1=quant_count1[which(!is.na(quant_count1$T4_gene_promoter_type)),]%>%group_by(group,detection,T4_gene_ncRNA_subclass)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aa1$group2="Gene"
quant_count$group="GENCODE"
quant_count$group[grep("ONT",quant_count$target_id)]="Novel"
quant_count$detection="Yes"
quant_count$detection[which(rowSums(quant_count[,c(2:7)])==0)]="No"
quant_count=left_join(quant_count, unique(table5[,c(1,65,101)]), by=c("target_id"="model_ID"),copy=F)
write.table(quant_count,gzfile(paste0(path_fig3_data,"Kallisto_Neuron.series.transcript.count.matrix.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5e

#===============================================================================
# content of novel isoform
table0=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table0.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
table0a=separate_rows(table0[,c(1,18)],full_set_bound_str, sep="_")
table0a2=table0a[grep("F",table5a$full_set_bound_str),]
table0a3=table0a[grep("T",table5a$full_set_bound_str),]
table0b2=table0a2[grep("ENST",table0a2$model_ID),]
table0b3=table0a3[grep("ENST",table0a3$model_ID),]

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
table5a=separate_rows(table5[,c(1,18)],full_set_bound_str, sep="_")
table5a2=table5a[grep("F",table5a$full_set_bound_str),]
table5a3=table5a[grep("T",table5a$full_set_bound_str),]
table5c2=table5a2[-grep("ENST",table5a2$model_ID),]
table5c3=table5a3[-grep("ENST",table5a3$model_ID),]

table5c2$novel_ex5="Yes"
table5c2$novel_ex5[which(table5c2$full_set_bound_str %in% unique(table0b2$full_set_bound_str))]="No"
table5c3$novel_ex3="Yes"
table5c3$novel_ex3[which(table5c3$full_set_bound_str %in% unique(table0b3$full_set_bound_str))]="No"


table5d=left_join(table5c2,table5c3, by="model_ID",copy=F)
#table5d1=table5c%>%group_by(model_ID)%>%dplyr::summarise(n_novel_sj=length(which(novel_sj=="Yes")),novel_sj=paste(full_set_bound_str[which(novel_sj=="Yes")],collapse="_"))
#table5d=left_join(table5d,table5d1, by="model_ID",copy=F)
colnames(table5d)[c(2,4)]=c("ex5_cluster","ex3_cluster")

#incorporate SQANTI result of SALA_finalized transcriptome
#result generated from [primary_folder]/code_n_data/transcript_model_analyses_Fig3/compare_diff_assemblers.R
SQANTI3.table5=read.delim(paste0(SQANTI_path,"SALA_finalized/SALA_table5_IP_vs_gencode/SALA_table5_IP_vs_gencode_classification.txt.gz"), header=T, stringsAsFactors = F, check.names = F)

table5d=left_join(table5d,table5[,c("model_ID","n_exon","transcript_length","T4_Gencode_geneCalss2","Coding_prob","CPAT_class","polyA","n3_support","n5_support")], by="model_ID", copy=F)
table5d%>%group_by(T4_Gencode_geneCalss2)%>%dplyr::summarise(count=n())
table5d=left_join(table5d, SQANTI3.table5[,c(1,6,8)], by=c("model_ID"="isoform"), copy=F)
table5d[which(table5d$T4_Gencode_geneCalss2 == "protein_coding"),]%>%group_by(structural_category)%>%dplyr::summarise(count=n())

#get the n5 n3 string for all the ENST, no match is new
options(scipen=999)
gencode12=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.bed.bgz"), header=F, stringsAsFactors = F)
gencode12_5=gencode12[,c(1:6)]
gencode12_3=gencode12[,c(1:6)]
gencode12_5$V3[which(gencode12_5$V6 == "+")]=gencode12_5$V2[which(gencode12_5$V6 == "+")]+1
gencode12_5$V2[which(gencode12_5$V6 == "-")]=gencode12_5$V3[which(gencode12_5$V6 == "-")]-1
gencode12_3$V3[which(gencode12_3$V6 == "-")]=gencode12_3$V2[which(gencode12_3$V6 == "-")]+1
gencode12_3$V2[which(gencode12_3$V6 == "+")]=gencode12_3$V3[which(gencode12_3$V6 == "+")]-1
write.table(gencode12_5[order(gencode12_5$V1,gencode12_5$V2),],gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n5.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)
write.table(gencode12_3[order(gencode12_3$V1,gencode12_3$V2),],gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n3.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#=======================================
#bash
#intersect gencode start and end sites with ex5_cluster and ex3_cluster
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
system(paste0("bedtools intersect -wa -wb -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n5.bed.gz -b Neuron_THP1.S3.end5.bed.bgz -s | gzip > gencode.v39.n5_string.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n3.bed.gz -b Neuron_THP1.S3.end3.bed.bgz -s | gzip > gencode.v39.n3_string.bed.gz"))
#=======================================

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
gencode_n5_string=read.delim("gencode.v39.n5_string.bed.gz", header=F, stringsAsFactors = F, check.names = F)
gencode_n3_string=read.delim("gencode.v39.n3_string.bed.gz", header=F, stringsAsFactors = F, check.names = F)

#add polyA prediction
gc_pas=read.delim(paste0(primary_folder,"code_n_data/Fig1_read_analyses/TES/all_read_n3_ont_gencode_prediction20240407.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gencode_n3_string$TESID=paste0(gencode_n3_string$V1,"_",gencode_n3_string$V2,"_",gencode_n3_string$V3,"_",gencode_n3_string$V6)
gencode_n3_string=left_join(gencode_n3_string[,c(1:6,10,19)],gc_pas[,c(8,12)], by=c("TESID"="label"),copy=F)
length(which(gencode_n3_string$GENCODE_polyA == "No"))

table5d=left_join(table5d, table0[,c(1, 56,57)], by=c("associated_transcript"="model_ID"),copy=F, suffix=c("","_ENST"))
table5d=left_join(table5d, gencode_n5_string[,c(4,10)], by=c("associated_transcript"="V4"),copy=F)
table5d=left_join(table5d, gencode_n3_string[,c(4,7,9)], by=c("associated_transcript"="V4"),copy=F, suffix=c("_n5","_n3"))
colnames(table5d)[c(13,18,19,20)]=c("polyA","n5_string_ENST","n3_string_ENST","polyA_ENST")

#===============================================================================
#combine structural_category from SQANTI with alternative TSS/TES
table5d$alt_TSS="No"
table5d$alt_TSS[which(table5d$ex5_cluster != table5d$n5_string_ENST)]="Yes"
table5d$alt_TES="No"
table5d$alt_TES[which(table5d$ex3_cluster != table5d$n3_string_ENST)]="Yes"
table5d$class=table5d$structural_category
table5d$class[which(table5d$structural_category == "full-splice_match" & table5d$alt_TSS=="Yes")] <- "alternative_TSS"
table5d$class[which(table5d$structural_category == "full-splice_match" & table5d$alt_TES=="Yes")] <- "alternative_TES"
table5d$class[which(table5d$structural_category == "full-splice_match" & table5d$alt_TSS=="Yes" & table5d$alt_TES=="Yes")] <- "alternative_TSS&TES"
table5d$class[which(table5d$structural_category == "incomplete-splice_match" & table5d$alt_TSS=="Yes")] <- "ISM_alternative_TSS"
table5d$class[which(table5d$structural_category == "incomplete-splice_match" & table5d$alt_TES=="Yes")] <- "ISM_alternative_TES"
table5d$class[which(table5d$structural_category == "incomplete-splice_match" & table5d$alt_TSS=="Yes" & table5d$alt_TES=="Yes")] <- "ISM_alternative_TSS&TES"
table5d$class[which(table5d$structural_category %in% c("antisense","fusion","genic","intergenic" ))] <- "novel_not_in_catalog"
table5d$structural_category[which(table5d$structural_category %in% c("antisense","fusion","genic","intergenic" ))] <- "novel_not_in_catalog"
table5d$end_class="n.a."
table5d$end_class[which(table5d$alt_TSS=="Yes")] <- "alt_TSS"
table5d$end_class[which(table5d$alt_TES=="Yes")] <- "alt_TES"
table5d$end_class[which(table5d$alt_TSS=="Yes" & table5d$alt_TES=="Yes")] <- "alt_TSS&TES"
write.table(table5d,gzfile(paste0(path_fig3_data,"novel.transcript.SQANTI3.table5.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
table5d=read.delim(paste0(path_fig3_data,"novel.transcript.SQANTI3.table5.tsv.gz"), header=T, stringsAsFactors = F)

#any of alt TES stop at donor site?
options(scipen=999)
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info"))
gencode_intron=read.delim("gencode.v39.annotation.bed6.intron.bed.gz", header=F, stringsAsFactors = F)
gencode_intron=gencode_intron%>%group_by(V1,V2,V3,V6)%>%summarise(V5=n())
gencode_intron$V4=paste0(gencode_intron$V1,"_",gencode_intron$V2,"_",gencode_intron$V3,"_",gencode_intron$V6)
gencode_intron=gencode_intron[,c(1:3,6,5,4)]
gencode_donor=gencode_intron
gencode_acceptor=gencode_intron
gencode_donor$V3[which(gencode_donor$V6 == "+")]=gencode_donor$V2[which(gencode_donor$V6 == "+")]
gencode_donor$V2[which(gencode_donor$V6 == "+")]=gencode_donor$V2[which(gencode_donor$V6 == "+")]-1
gencode_donor$V2[which(gencode_donor$V6 == "-")]=gencode_donor$V3[which(gencode_donor$V6 == "-")]
gencode_donor$V3[which(gencode_donor$V6 == "-")]=gencode_donor$V3[which(gencode_donor$V6 == "-")]+1
gencode_acceptor$V3[which(gencode_acceptor$V6 == "-")]=gencode_acceptor$V2[which(gencode_acceptor$V6 == "-")]
gencode_acceptor$V2[which(gencode_acceptor$V6 == "-")]=gencode_acceptor$V2[which(gencode_acceptor$V6 == "-")]-1
gencode_acceptor$V2[which(gencode_acceptor$V6 == "+")]=gencode_acceptor$V3[which(gencode_acceptor$V6 == "+")]
gencode_acceptor$V3[which(gencode_acceptor$V6 == "+")]=gencode_acceptor$V3[which(gencode_acceptor$V6 == "+")]+1
write.table(gencode_donor[order(gencode_donor$V1,gencode_donor$V2),],gzfile("gencode.v39.donor_minus1.bed.gz"),col.names=F, row.names=F, sep="\t", quote=F)
write.table(gencode_acceptor[order(gencode_acceptor$V1,gencode_acceptor$V2),],gzfile("gencode.v39.acceptor_plus1.bed.gz"),col.names=F, row.names=F, sep="\t", quote=F)

#=======================================
#bash
#intersect TES with donor site
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info"))
system(paste0("bedtools closest -a GENCODEv39.transcript.n3.bed.gz -b gencode.v39.donor_minus1.bed.gz -s -D a | gzip > ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed/gencode.v39.n3_closest.donor_minus1.bed.gz"))

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
system(paste0("bedtools closest -a Neuron_THP1.S3.TES.table5.1bp.bed -b ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/gencode.v39.donor_minus1.bed.gz -s -D a | gzip > table5.n3_closest.donor_minus1.bed.gz"))
#=======================================

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
gencode_n3_donor=read.delim("gencode.v39.n3_closest.donor_minus1.bed.gz", header=F, stringsAsFactors = F, check.names = F)
table5_n3_donor=read.delim("table5.n3_closest.donor_minus1.bed.gz", header=F, stringsAsFactors = F, check.names = F)
gencode_n3_donor=unique(gencode_n3_donor[which(abs(gencode_n3_donor$V13) == 0),c(4,13)])
table5_n3_donor=unique(table5_n3_donor[which(abs(table5_n3_donor$V13) == 0),c(4,13)])
table5_n3_donor=left_join(table5_n3_donor,table5[,c(1,106)],by=c("V4"="TESID"),copy=F)

table5d$TES_donor_minus1="No"
table5d$TES_donor_minus1[which(table5d$model_ID %in% table5_n3_donor$model_ID)]="Yes"

#=====
# Add ORF info
SQANTI3.table5=read.delim(paste0(SQANTI_path,"SALA_final/SALA_final_classification.txt.gz"), header=T, stringsAsFactors = F, check.names = F)
table5d=left_join(table5d,SQANTI3.table5[,c(1,31:36,47)], by=c("model_ID"="isoform"), copy=F)
write.table(table5d,gzfile(paste0(path_fig3_data,"novel.transcript.SQANTI3.table5.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

table5d$ENST_FL_detected=NA
table5d$ENST_FL_detected[which(table5d$associated_transcript!= "novel")]="No"
table5d$ENST_FL_detected[which(table5d$associated_transcript %in% table5$model_ID)]="Yes"

table5e=table5d[which(table5d$T4_Gencode_geneCalss2 == "protein_coding"),]
table5e=table5e[which(table5e$n3_support != "no_support"),]
write.table(table5e,gzfile(paste0(path_fig3_data,"novel.isoform.protein_coding.SQANTI3.table5.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5i & j



#IsoformSwitchAnalyzeR
#================================================================================================================
#bambu -> this one with gtf version with all partially detectable ENST add back
#====================================================================================
#Open Targets: downloaded from https://platform.opentargets.org/ (Koscielny et al. 2017)
#Refseq gtf (release from 2024_08) was downloaded from NCBI
opendata1=read.delim(paste0(long_read_bambu,"combined_disease_association_02.tsv.gz"),header=T, stringsAsFactors=F, check.names = F)
refseq=read.delim(paste0(primary_folder,"code_n_data/GENCODEv39/gencode.v39.metadata.RefSeq.gz"), header=F, stringsAsFactors = F)
refseq=refseq%>%group_by(V1)%>%dplyr::summarise(refseq_RNA=paste(unique(V2),collapse=";",recycle0 =F),refseq_protein=paste(unique(V3),collapse=";",recycle0 =F))

library(IsoformSwitchAnalyzeR)

bambu_t=fread(paste0(long_read_bambu,"counts_transcript.txt"), header=T)

study=data.frame(colnames(bambu_t))
colnames(study)[1]="sampleID"
study$sampleID=gsub("Set18-4_","",study$sampleID)
study$sampleID=gsub("Set19-2_","",study$sampleID)
study$sampleID=gsub("_concatenated.sorted","",study$sampleID)
study$sampleID=gsub("_v4.labeled_filtered.sorted","",study$sampleID)
study$sampleID=gsub("iPSC_iPSC","iPSC",study$sampleID)
study$sampleID=gsub("NSC_NSC","NSC",study$sampleID)
study$sampleID=gsub("Neuron_Neuron","Neuron",study$sampleID)
study$sampleID=gsub("-","_",study$sampleID)

study$condition=sapply(strsplit(study$sampleID,"_"),"[",1)
study$condition[c(13:16)]="PMA24"
study$condition[c(21:24)]="PMA96"
study=study[-c(1:2),]
colnames(bambu_t)[c(3:24)]=study$sampleID

#==========================================================
bambu_t1=bambu_t[,-2]
write.table(bambu_t1,gzfile(paste0(long_read_bambu,"transcript.count.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

bambu_t1=bambu_t1[,c(1:7)]
colnames(bambu_t1)[1]="isoform_id"
study=study[c(1:6),]

#============================================
#prepare transcript fasta
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/"))
system("bedtools getfasta -s -split -nameOnly -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed table5.final.partial_yes_detected.alone.bed12.bed.gz > table5.final.partial_yes_detected.alone.fasta")

fastaa=fread(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.fasta"), header=F)
fastaa$V1=gsub("\\(\\+\\)","",fastaa$V1)
fastaa$V1=gsub("\\(\\-\\)","",fastaa$V1)
write.table(fastaa, gzfile(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.fasta.gz")), col.names =F, row.names=F, sep="\t", quote=F)

#isoform switch - use the original files containing rRNA
aSwitchList <- importRdata(
  isoformCountMatrix   = bambu_t1,
  designMatrix         = study,
  isoformExonAnnoation = paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.gtf.gz"),
  isoformNtFasta     = paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.fasta.gz"),
  fixStringTieAnnotationProblem = TRUE,
  showProgress = FALSE)

#=======
switchResult=isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  alpha = 0.05,
  dIFcutoff = 0.1)

extractSwitchSummary(switchResult)
switchResult_df=data.frame(switchResult$isoformSwitchAnalysis)
ips_nsc=switchResult_df[which(switchResult_df$condition_1=="iPSC" & switchResult_df$condition_2 == "NSC"),]
ips_nrn=switchResult_df[which(switchResult_df$condition_1=="iPSC" & switchResult_df$condition_2 == "Neuron"),]
nrn_nsc=switchResult_df[which(switchResult_df$condition_1=="Neuron" & switchResult_df$condition_2 == "NSC"),]
colnames(ips_nsc)[c(6:8)]=paste0("iPS_NSC_",colnames(ips_nsc)[c(6:8)])
colnames(ips_nsc)[c(9,10)]=c("iPSC_IF","NSC_IF")
colnames(ips_nrn)[c(6:8)]=paste0("iPS_NRN_",colnames(ips_nrn)[c(6:8)])
colnames(nrn_nsc)[c(6:8)]=paste0("NSC_NRN_",colnames(nrn_nsc)[c(6:8)])
colnames(nrn_nsc)[c(9)]=c("NRN_IF")
nrn_nsc$NSC_NRN_dIF=nrn_nsc$NSC_NRN_dIF*(-1)
combine=left_join(ips_nsc[,c(3,6:10)], nrn_nsc[,c(3,9,6:8)], by="isoform_id", copy=F)
combine=left_join(combine, ips_nrn[,c(3,6:8)], by="isoform_id", copy=F)
combine=combine[,c(1,5:7,2:4,8:13)]
write.table(combine, gzfile(paste0(long_read_bambu,"isoform.switch.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#add CPM from bambu
bambu_tcpm=fread(paste0(long_read_bambu,"CPM_transcript.txt.gz"), header=T)
bambu_tcpm=bambu_tcpm[,c(1,3:8)]
bambu_tcpm$iPSC_CPM=bambu_tcpm$`iPSC_Set18-4_iPSC-PAP_rep1_concatenated.sorted`+bambu_tcpm$`iPSC_Set18-4_iPSC-PAP_rep2_concatenated.sorted`
bambu_tcpm$NSC_CPM=bambu_tcpm$`NSC_Set18-4_NSC-PAP_rep1_concatenated.sorted`+bambu_tcpm$`NSC_Set18-4_NSC-PAP_rep2_concatenated.sorted`
bambu_tcpm$Neuron_CPM=bambu_tcpm$`Neuron_Set18-4_Neuron-PAP_rep1_concatenated.sorted`+bambu_tcpm$`Neuron_Set18-4_Neuron-PAP_rep2_concatenated.sorted`
bambu_tcpm=bambu_tcpm[,c(1,8:10)]

combine=combine[,c(1:13)]
combine=left_join(combine, table5[,c(1,77,78,79,80,86,71,72,73,89,113)], by=c("isoform_id"="model_ID"),copy=F)
combine=left_join(combine,opendata1[,c(10,1,11,9)], by=c("T4_gene_ID"="geneID"),copy=F)
combine=left_join(combine,refseq, by=c("isoform_id"="V1"),copy=F)
combine=left_join(combine,bambu_tcpm,by=c("isoform_id"="TXNAME"), copy=F)
combine=combine[,c(29:31,1:28)]
write.table(combine, gzfile(paste0(long_read_bambu,"isoform.switch.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
combine=read.delim(paste0(long_read_bambu,"isoform.switch.tsv.gz"),header=T, stringsAsFactors=F, check.names = F)
combine=combine[which(!is.na(combine$T4_gene_ID)),]
combine$iPS_NSC_max_CPM=apply(combine[,c(1:2)], 1, max)
combine$NSC_NRN_max_CPM=apply(combine[,c(2:3)], 1, max)
combine$iPS_NRN_max_CPM=apply(combine[,c(1,3)], 1, max)

combine=combine[,-which(colnames(combine) %in% c("ENST_adjust","ENSG_adjust"))]

needk=union(union(which(abs(combine$iPS_NSC_dIF)>0.5 & combine$iPS_NSC_padj<0.05 & combine$iPS_NSC_max_CPM>2),
                  which(abs(combine$NSC_NRN_dIF)>0.5 & combine$NSC_NRN_padj<0.05 & combine$NSC_NRN_max_CPM>2)),
            which(abs(combine$iPS_NRN_dIF)>0.5 & combine$iPS_NRN_padj<0.05 & combine$iPS_NRN_max_CPM>2))
combine2=combine[needk,]
combine2=left_join(combine2, table5[,c("model_ID","n_exon","transcript_length","n5_string","n3_string","n5_support","n3_support")], by=c("isoform_id"="model_ID"), copy=F)
combine2=combine2[grep("ENS",combine2$T4_gene_ID),]
combine2$ONTT_involve="No"
combine2$ONTT_involve[grep("ONTT",combine2$isoform_id )]="Yes"
combine2%>%group_by(ONTT_involve)%>%dplyr::summarise(count=n())

write.table(combine2[order(combine2$maxScore, decreasing = T),], gzfile(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#stored in [primary_folder]/fig3/data
# for -> Table S9

#=================================
#selected those change with both sig up and down in a gene
acom=combine%>%group_by(T4_gene_ID)%>%dplyr::summarise(count=n(), count_up=length(which(iPS_NSC_dIF>0.5 & iPS_NSC_padj<0.05 & iPS_NSC_max_CPM>2)), count_down=length(which(iPS_NSC_dIF< (-0.5) & iPS_NSC_padj<0.05 & iPS_NSC_max_CPM>2)))
bcom=combine%>%group_by(T4_gene_ID)%>%dplyr::summarise(count=n(), count_up=length(which(NSC_NRN_dIF>0.5 & NSC_NRN_padj<0.05 & NSC_NRN_max_CPM>2)), count_down=length(which(NSC_NRN_dIF< (-0.5) & NSC_NRN_padj<0.05 & NSC_NRN_max_CPM>2)))
ccom=combine%>%group_by(T4_gene_ID)%>%dplyr::summarise(count=n(), count_up=length(which(iPS_NRN_dIF>0.5 & iPS_NRN_padj<0.05 & iPS_NRN_max_CPM>2)), count_down=length(which(iPS_NRN_dIF< (-0.5) & iPS_NRN_padj<0.05 & iPS_NRN_max_CPM>2)))

acom=acom[which(acom$count_up>0 & acom$count_down >0),]
bcom=bcom[which(bcom$count_up>0 & bcom$count_down >0),]
ccom=ccom[which(ccom$count_up>0 & ccom$count_down >0),]

acom1=combine[grep("ONT",combine$isoform_id),]%>%group_by(T4_gene_ID)%>%dplyr::summarise(count_up_ONT=length(which(iPS_NSC_dIF>0.5 & iPS_NSC_padj<0.05 & iPS_NSC_max_CPM>2)), count_down_ONT=length(which(iPS_NSC_dIF< (-0.5) & iPS_NSC_padj<0.05 & iPS_NSC_max_CPM>2)))
bcom1=combine[grep("ONT",combine$isoform_id),]%>%group_by(T4_gene_ID)%>%dplyr::summarise(count_up_ONT=length(which(NSC_NRN_dIF>0.5 & NSC_NRN_padj<0.05 & NSC_NRN_max_CPM>2)), count_down_ONT=length(which(NSC_NRN_dIF< (-0.5) & NSC_NRN_padj<0.05 & NSC_NRN_max_CPM>2)))
ccom1=combine[grep("ONT",combine$isoform_id),]%>%group_by(T4_gene_ID)%>%dplyr::summarise(count_up_ONT=length(which(iPS_NRN_dIF>0.5 & iPS_NRN_padj<0.05 & iPS_NRN_max_CPM>2)), count_down_ONT=length(which(iPS_NRN_dIF< (-0.5) & iPS_NRN_padj<0.05 & iPS_NRN_max_CPM>2)))

acom1=acom1[intersect(which(acom1$count_up_ONT>0 | acom1$count_down_ONT>0),which(acom1$T4_gene_ID %in% acom$T4_gene_ID)),]
bcom1=bcom1[intersect(which(bcom1$count_up_ONT>0 | bcom1$count_down_ONT>0),which(bcom1$T4_gene_ID %in% bcom$T4_gene_ID)),]
ccom1=ccom1[intersect(which(ccom1$count_up_ONT>0 | ccom1$count_down_ONT>0),which(ccom1$T4_gene_ID %in% ccom$T4_gene_ID)),]

needd=union(union(acom$T4_gene_ID, bcom$T4_gene_ID), ccom$T4_gene_ID)
needd1=union(union(acom1$T4_gene_ID, bcom1$T4_gene_ID), ccom1$T4_gene_ID)
combine4=combine[which(combine$T4_gene_ID %in% needd),]
combine4$iPS_NSC_hit="No"
combine4$iPS_NSC_hit[which(combine4$T4_gene_ID %in% acom$T4_gene_ID)]="Yes"
combine4$NSC_NRN_hit="No"
combine4$NSC_NRN_hit[which(combine4$T4_gene_ID %in% bcom$T4_gene_ID)]="Yes"
combine4$iPS_NRN_hit="No"
combine4$iPS_NRN_hit[which(combine4$T4_gene_ID %in% ccom$T4_gene_ID)]="Yes"
combine4$ONTT_involve="No"
combine4$ONTT_involve[combine4$T4_gene_ID %in% needd1]="Yes"

neede=union(union(which(abs(combine4$iPS_NSC_dIF)>0.5 & combine4$iPS_NSC_padj<0.05 & combine4$iPS_NSC_max_CPM>2),
                  which(abs(combine4$NSC_NRN_dIF)>0.5 & combine4$NSC_NRN_padj<0.05 & combine4$NSC_NRN_max_CPM>2)),
            which(abs(combine4$iPS_NRN_dIF)>0.5 & combine4$iPS_NRN_padj<0.05 & combine4$iPS_NRN_max_CPM>2))
combine5=combine4[neede,]
combine5=left_join(combine5, table5[,c("model_ID","n_exon","transcript_length","n5_string","n3_string","n5_support","n3_support")], by=c("isoform_id"="model_ID"), copy=F)
write.table(combine5[order(combine5$maxScore, decreasing = T),], gzfile(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5k & l, 3f


#=================================
#if novel start site initiate from simple repeat
#hipstr reference was downloaded from https://github.com/HipSTR-Tool/HipSTR-references
hipstr=read.delim("hg38.hipstr_reference.bed.gz", header=F, stringsAsFactors = F)
#=================================
#bash
#intersect this hipstr with table5 transcript model 5' end and gencode 5' end
setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
system(paste0("bedtools intersect -wa -wb -a ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/Input_Neuron_THP1/GENCODE_info/GENCODEv39.transcript.n5.bed.gz -b ",long_read_bambu,"hg38.hipstr_reference.bed.gz | gzip > gencode.v39.n5_hipstr.bed.gz"))
system(paste0("bedtools intersect -wa -wb -a Neuron_THP1.S3.model.5n.bed.gz -b ",long_read_bambu,"hg38.hipstr_reference.bed.gz | gzip > Neuron_THP1.S3.model.5n_hipstr.bed.gz"))
#================================

setwd(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/bed"))
gencode_hip=read.delim("gencode.v39.n5_hipstr.bed.gz", header=F, stringsAsFactors = F, check.names = F)
table0_hip=read.delim("Neuron_THP1.S3.model.5n_hipstr.bed.gz", header=F, stringsAsFactors = F, check.names = F)
length(intersect(table0_hip$V4,gencode_hip$V4))
table0_hip=left_join(table0_hip, table5[,c(1,62,64,65,86,87,90,110)], by=c("V4"="model_ID"),copy=F)
table5_hip=table0_hip[which(!is.na(table0_hip$T4_gene_ID)),]
k0=table0_hip%>%group_by(V13,V6)%>%dplyr::summarise(count=n())
k5=table5_hip%>%group_by(V13,V6)%>%dplyr::summarise(count=n())

table5e=read.delim(paste0(path_fig3_data,"novel.isoform.protein_coding.SQANTI3.table5.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5e=left_join(table5e, table5_hip[,c(4,13)], by=c("model_ID"="V4"),copy=F)
table5e=left_join(table5e, table5_hip[,c(4,13)], by=c("associated_transcript"="V4"),copy=F)
colnames(table5e)[c(34,35)]=c("TSS_hipSTR","TSS_hipSTR_ENST")
write.table(table5e, paste0(path_fig3_data,"novel.isoform.protein_coding.SQANTI3.table5.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)

combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
combine2=left_join(combine2, table5e[,c(1,23:25,34)], by=c("isoform_id"="model_ID"),copy=F)
combine2=left_join(combine2 ,SQANTI3.table5[,c(1,31:36,47)], by=c("isoform_id"="isoform"), copy=F)
combine2$class[grep("ENS",combine2$isoform_id)]="GENCODE"
combine2$class[which(is.na(combine2$class))]="novel_in_catalog"
colnames(combine2)[42]="structural_category"
combine2$end_class[which(is.na(combine2$end_class))]="n.a."

combine5=read.delim(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
combine5=left_join(combine5, table5e[,c(1,23:25,34)], by=c("isoform_id"="model_ID"),copy=F)
combine5=left_join(combine5 ,SQANTI3.table5[,c(1,31:36,47)], by=c("isoform_id"="isoform"), copy=F)
combine5$class[grep("ENS",combine5$isoform_id)]="GENCODE"
combine5$class[which(is.na(combine5$class))]="novel_in_catalog"
colnames(combine5)[46]="structural_category"
combine5$end_class[which(is.na(combine5$end_class))]="n.a."

write.table(combine2, gzfile(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(combine5, gzfile(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#add ATAC count into ex5_cluster
#scATAC raw data is provided in DDBJ
library(Signac)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(GenomicRanges)
library(future)
library(chromVAR)

peak_ex5 <- getPeaks(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/zenbu/Neuron_THP1.S3.end5.bed.gz"), sort_peaks = TRUE)
IDList <- c("iPS","NSC","NRN")
md.List <- list()
for(i in IDList){
  tmp <- read.table(
    file = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/singlecell.csv"), 
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1)[-1, ]
  md.List[[i]] <- tmp}

for(i in IDList){md.List[[i]] <- md.List[[i]][md.List[[i]]$passed_filters > 500, ]}

## create fragment objects
frags.List <- list()
for(i in IDList){
  tmp <- CreateFragmentObject(
    path = paste0("/analysisdata/fantom6/Interactome/scATACcellranger_Kouno/",i,"/outs/fragments.tsv.gz"),
    cells = rownames(md.List[[i]]))
  frags.List[[i]] <- tmp}
plan("multicore", workers = 6)
options(future.globals.maxSize = 50 * 1024^3) 
counts.List <- list()
for(i in IDList){
  tmp <-  FeatureMatrix(
    fragments = frags.List[[i]],
    features = peak_ex5,
    cells = rownames(md.List[[i]]))
  counts.List[[i]] <- tmp}

c1=as.data.frame(rowSums(as.matrix(counts.List[["iPS"]])))
c2=as.data.frame(rowSums(as.matrix(counts.List[["NSC"]])))
c3=as.data.frame(rowSums(as.matrix(counts.List[["NRN"]])))
c1$peakID=rownames(c1)
c2$peakID=rownames(c2)
c3$peakID=rownames(c3)
c0=full_join(c1,c2, by="peakID",copy=F)
c0=full_join(c0,c3, by="peakID",copy=F)
c0$peakID=rownames(as.matrix(counts.List[["iPS"]]))
colnames(c0)[c(1,3,4)]=c("ATACcount_iPSC","ATACcount_NSC","ATACcount_Neuron")
c0$chr=sapply(strsplit(c0$peakID,"-"),"[",1)
c0$start=as.numeric(sapply(strsplit(c0$peakID,"-"),"[",2))-1
c0$end=as.numeric(sapply(strsplit(c0$peakID,"-"),"[",3))
c0=left_join(c0,kk1,by=c("chr"="V1","start"="V2","end"="V3"),copy=F)
c0a=separate_rows(c0, V4, sep="\\|")
write.table(c0[,c(5:8,2,1,3,4)],gzfile(paste0(path_fig3_data,"ex5cluster_ATACcount_celltype.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
combine5=read.delim(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)                   
ex5_ATACcount=read.delim(paste0(path_fig3_data,"ex5cluster_ATACcount_celltype.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(ex5_ATACcount)[c(6:8)]=paste0("ex5_",colnames(ex5_ATACcount)[c(6:8)])
ex5_ATACcount$ex5_ATAC_CPM_iPSC=ex5_ATACcount$ex5_ATACcount_iPSC/sum(ex5_ATACcount$ex5_ATACcount_iPSC)*1000000
ex5_ATACcount$ex5_ATAC_CPM_NSC=ex5_ATACcount$ex5_ATACcount_NSC/sum(ex5_ATACcount$ex5_ATACcount_NSC)*1000000
ex5_ATACcount$ex5_ATAC_CPM_Neuron=ex5_ATACcount$ex5_ATACcount_Neuron/sum(ex5_ATACcount$ex5_ATACcount_Neuron)*1000000

combine2=left_join(combine2, ex5_ATACcount[,c(4,9:11)], by=c("n5_string"="V4"),copy=F)
combine5=left_join(combine5, ex5_ATACcount[,c(4,9:11)], by=c("n5_string"="V4"),copy=F)

#===============================================================================
#add back quantification from short read
combine2=left_join(combine2, shortkallisto[,c(1,8:10)] ,by=c("isoform_id"="target_id"), copy=F)
combine2=combine2[,c(1:3,51:53,4:50)]
colnames(combine2)[c(1:6)]=c("iPSC_CPM_longread","NSC_CPM_longread","Neuron_CPM_longread","iPSC_TPM_shortread","NSC_TPM_shortread","Neuron_TPM_shortread")
combine5=left_join(combine5, shortkallisto[,c(1,8:10)] ,by=c("isoform_id"="target_id"), copy=F)
combine5=combine5[,c(1:3,54:56,4:53)]
colnames(combine5)[c(1:6)]=c("iPSC_CPM_longread","NSC_CPM_longread","Neuron_CPM_longread","iPSC_TPM_shortread","NSC_TPM_shortread","Neuron_TPM_shortread")

write.table(combine2, gzfile(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(combine5, gzfile(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(combine2, paste0(primary_folder,"/supplementary_table/TableS9_isoform.switch.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5k & l, 3f, Table S9

#===============================================================================
#prepare bed12 for zenbu
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/")
bed12=read.delim(paste0(SALA_path,"all_gtf_file/table5.final.partial_yes_detected.alone.bed12.bed.gz"), header=F, stringsAsFactors = F)
bed12=left_join(bed12, bambu_tcpm, by=c("V4"="TXNAME"),copy=F)
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
bed12=left_join(bed12, table5[,c(1,94:96)], by=c("V4"="model_ID"),copy=F)
bed12i=bed12
bed12i$V4=paste0(bed12i$V4," | CPM:",bed12i$iPSC_CPM)
write.table(bed12i[order(bed12i$V1,bed12i$V2),c(1:4,16,6:12)],gzfile(paste0(SALA_path,"zenbu/detectable_table5_partial_yes/table5.final.partial_yes_detected.alone.iPSCV5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
bed12s=bed12
bed12s$V4=paste0(bed12s$V4," | CPM:",bed12s$NSC_CPM)
write.table(bed12s[order(bed12s$V1,bed12s$V2),c(1:4,17,6:12)],gzfile(paste0(SALA_path,"zenbu/detectable_table5_partial_yes/table5.final.partial_yes_detected.alone.NSCV5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
bed12n=bed12
bed12n$V4=paste0(bed12n$V4," | CPM:",bed12n$Neuron_CPM)
write.table(bed12n[order(bed12n$V1,bed12n$V2),c(1:4,18,6:12)],gzfile(paste0(SALA_path,"zenbu/detectable_table5_partial_yes/table5.final.partial_yes_detected.alone.NRNV5.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)



