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
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_default/")
path1=paste0(SALA_path,"default/transcript/zenbu/")
path2=paste0(SALA_path,"default/transcript/log/")
path3=paste0(SALA_path,"default/transcript/bed/")
path4=paste0(SALA_path,"default/transcript/cpat/")
gene0_path=paste0(SALA_path,"default/table0_gene/")
gene4_path=paste0(SALA_path,"default/table4_gene/")
GENCODE_path=paste0(primary_folder,"code_n_data/GENCODEv39/")

#===============================================================================
#run SALA default parameter [https://github.com/fantom-prj/SALA]
# input files derived from SALA Final (paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/")

setwd(SALA_path)
system("sh transcript.sh")
system("sh gene0.sh")
system("sh count.sh")
system("sh filter.sh")
system("sh gene4.sh")

#===============================================================================
# The following folders were removed after incorporation:
# table0_gene, table4_gene, CPAT, Input 
# only provide upon request
#===============================================================================


#===============================================================================
# perform analyses not included in SALA
#sub-classification: identify "ncRNA subclass" for SALA Final 
setwd(path2)
table5=read.delim("Neuron_THP1.S3.table4_filtered.noIP.All_Ref_updated.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
unique(table5$ref_source)
unique(table5$Ref_geneClass)

table5a=table5[which(table5$Ref_geneClass2=="lncRNA" & table5$Ref_transcriptClass2 =="lncRNA" & table5$ref_source== "fulllength_ref"),]
table5b=table5[which(table5$Ref_geneClass2=="lncRNA" & table5$Novel_transcriptClass %in% c("lncRNA")),]
table5c=table5[which(table5$T4_gene_novelty == "Novel" & table5$Novel_geneClass %in% c("lncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","ncRNA")),]
table5d=table5[which(table5$T4_gene_novelty == "Novel" & table5$Novel_geneClass %in% c("ncRNA") & table5$Novel_transcriptClass %in% c("lncRNA","ncRNA")),]

table5a$T4_ncRNA_source="ENST_lncRNA"
table5b$T4_ncRNA_source="novel_lncRNA_isoform"
table5c$T4_ncRNA_source="novel_lncRNA"
table5d$T4_ncRNA_source="novel_ncRNA"
table7=rbind(table5a,table5b,table5c,table5d)

options(scipen=999)
table5.bed12=read.delim(paste0(path3,"Neuron_THP1.S3.table4.bed12.bed.gz"), header=F,  stringsAsFactors = F)
table7.bed12=table5.bed12[which(table5.bed12$V4 %in% table7$model_ID),]
write.table(table7.bed12[order(table7.bed12$V1,table7.bed12$V2),],gzfile(paste0(path3,"table7_ncRNA.bed12.bed.gz")),col.names=F, row.names=F, sep="\t", quote=F)
#prepare n5
table7.bed12$V3[which(table7.bed12$V6 == "+")]=table7.bed12$V2[which(table7.bed12$V6 == "+")]+1
table7.bed12$V2[which(table7.bed12$V6 == "-")]=table7.bed12$V3[which(table7.bed12$V6 == "-")]-1
write.table(table7.bed12[order(table7.bed12$V1, table7.bed12$V2),c(1:6)], gzfile(paste0(path3,"table7_ncRNA.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#prepare #ENST with updated range
options(scipen=999)
ENSTbed12=read.delim(paste0(GENCODE_path,"gencode.v39.annotation.transcript.bed.gz"), header=F, stringsAsFactors = F, check.names = F)
table5.bed12a=table5.bed12[grep("ENST",table5.bed12$V4),]
ENSTbed12a=ENSTbed12[-which(ENSTbed12$V4 %in% table5.bed12a$V4),]
ENSTbed12=rbind(ENSTbed12a,table5.bed12a)
#write.table(ENSTbed12[order(ENSTbed12$V1, ENSTbed12$V2),], gzfile(paste0(path3,"ENST.T4updated.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#prepare mRNA & pseudo alone
gencode.anno=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/GENCODE.gene.n.transcript.class.1base.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gencode.anno1=gencode.anno[which(gencode.anno$Gencode_transcriptClass == "protein_coding"),]
gencode.anno2=gencode.anno[grep("pseudogene", gencode.anno$Gencode_transcriptClass),]
ENSTbed12a=ENSTbed12[which(ENSTbed12$V4 %in% c(gencode.anno1$transcriptID, gencode.anno2$transcriptID)),]
write.table(ENSTbed12a[order(ENSTbed12a$V1, ENSTbed12a$V2),], gzfile(paste0(path3,"ENST.T4updated.mRNA_pseudo.bed12.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#prepare n5
ENSTbed12a$V3[which(ENSTbed12a$V6 == "+")]=ENSTbed12a$V2[which(ENSTbed12a$V6 == "+")]+1
ENSTbed12a$V2[which(ENSTbed12a$V6 == "-")]=ENSTbed12a$V3[which(ENSTbed12a$V6 == "-")]-1
write.table(ENSTbed12a[order(ENSTbed12a$V1, ENSTbed12a$V2),c(1:6)], gzfile(paste0(path3,"ENST.T4updated.mRNA_pseudo.n5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#===================
#bedtools
setwd(path3)
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
setwd(path3)
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

exonExon=left_join(exonExon,table7[,c("model_ID","transcript_length")],by=c("V4"="model_ID"),copy=F)
exonExon$overlap_rate=exonExon$overlap/exonExon$transcript_length
exonExonAs=read.delim("table7_ncRNA.bed6.antisense.inter.mRNApseudo.exon.bed.gz", header=F, stringsAsFactors = F, check.names = F)
exonExonAs=exonExonAs%>%group_by(V4,V10)%>%dplyr::summarise(overlap=sum(V13))
exonExonAs=left_join(exonExonAs,table7[,c("model_ID","transcript_length")],by=c("V4"="model_ID"),copy=F)
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
write.table(table7,gzfile(paste0(path2,"table7.ncRNA.41K.subclass.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
table5=left_join(table5, table7[,c("model_ID","T4_ncRNA_source","T4_ncRNA_subclass")], by="model_ID", copy=F)
write.table(table5,gzfile(paste0(path2,"Neuron_THP1.S3.table4_filtered.noIP.All_Ref_updated.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#gene level ncRNA subclass
#order -> divergent lncRNAs, sense overlap lncRNAs, sense intronic lncRNAs, antisense lncRNAs and intergenic lncRNAs 
table7=read.delim(paste0(path2,"table7.ncRNA.41K.subclass.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
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
gene_data=unique(table7[,c("T4_gene_ID","T4_ncRNA_subclass","order")])%>%group_by(T4_gene_ID)%>%slice_max(order)
colnames(gene_data)[2]="T4_gene_ncRNA_subclass"
table7=left_join(table7, gene_data[,c(1,2)], by="T4_gene_ID",copy=F)
table5=left_join(table5, gene_data[,c(1,2)], by="T4_gene_ID", copy=F)

# add transcript_group
table5$transcript_group="ENST"
table5$transcript_group[which(table5$transcript_novelty != "Ref" & table5$gene_novelty != "Ref")]="Transcript_from_novel_gene"
table5$transcript_group[which(table5$transcript_novelty != "Ref" & table5$gene_novelty == "Ref")]="Novel_isoform"
table5%>%group_by(transcript_group)%>%dplyr::summarise(n())
write.table(table5,gzfile(paste0(path2,"Neuron_THP1.S3.table4_filtered.noIP.All_Ref_updated.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
#prepare transcript base 5' and 3' bed file
t5.bed12=read.delim(paste0(path3,"Neuron_THP1.S3.table4.bed12.bed.gz"), header=F,  stringsAsFactors = F)
t5.bed12=t5.bed12[which(t5.bed12$V4 %in% table5$model_ID[which(table5$ref_source %in% c("fulllength_ref","novel_transcript"))]),]


#count n5, n3 number
t5.bed12$n5=t5.bed12$V2
t5.bed12$n5[which(t5.bed12$V6=="-")]=t5.bed12$V3[which(t5.bed12$V6=="-")]
t5.bed12$n3=t5.bed12$V3
t5.bed12$n3[which(t5.bed12$V6=="-")]=t5.bed12$V2[which(t5.bed12$V6=="-")]
t5.bed12$n5_1=t5.bed12$n5+1
t5.bed12$n3_1=t5.bed12$n3+1
t5.bed12=left_join(t5.bed12, table5[,c("model_ID","transcript_group")], by=c("V4"="model_ID"),copy=F)
nrow(unique(t5.bed12[,c(1,13)])) #52249
nrow(unique(t5.bed12[,c(1,14)])) #69078

write.table(t5.bed12[order(t5.bed12$V1,t5.bed12$n5),c(1,13,15,4,5,6,17)],gzfile(paste0(path3,"SALA_Default.table5.N5.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)
write.table(t5.bed12[order(t5.bed12$V1,t5.bed12$n3),c(1,14,16,4,5,6,17)],gzfile(paste0(path3,"SALA_Default.table5.N3.bed.gz")), col.names=F, row.names=F, sep="\t", quote=F)

#===================================
#intersect with external and data-driven 5' features
n5_path=paste0(primary_folder,"code_n_data/n5_regions/")
cluster_path=paste0(primary_folder,"code_n_data/SCAFE/CFC_Neuron_THP1/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/")
setwd(path3)
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


#============================================================
### include information in raw table
# build gtf from raw table
#=====
# for SALA default
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_default/")
path1=paste0(SALA_path,"default/transcript/zenbu/")
path2=paste0(SALA_path,"default/transcript/log/")
path3=paste0(SALA_path,"default/transcript/bed/")

setwd(path3)
gene0_path=paste0(SALA_path,"default/table0_gene/")

transcript_info_final=read.delim(paste0(path2,"Neuron_THP1.S3.table0_raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
transcript_info_final=transcript_info_final[which(transcript_info_final$internal_priming != "Yes"),]
transcript_info_final$Novel_transcriptClass <- "NA"
transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" )] <- "others"
transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" & transcript_info_final$CPAT_class == "non-coding")] <- "ncRNA"
transcript_info_final$Novel_transcriptClass[which(transcript_info_final$transcript_novelty == "Novel" & transcript_info_final$CPAT_class == "non-coding"  & transcript_info_final$transcript_length > 200)] <- "lncRNA"
transcript_info_final%>%group_by(Novel_transcriptClass)%>%dplyr::summarise(count=n())

transcript_info_final <- transcript_info_final%>%group_by(IN1_gene_ID)%>%dplyr::mutate(overall_T_ratio=full_qry_count/sum(full_qry_count))
transcript_infoa <- transcript_info_final[which(transcript_info_final$Novel_transcriptClass %in% c("ncRNA","lncRNA") & transcript_info_final$gene_novelty == "Novel"),]
transcript_infoa <- transcript_infoa%>%group_by(IN1_gene_ID)%>%dplyr::mutate(overall_T_ratio2=full_qry_count/sum(full_qry_count))
transcript_infob <- transcript_infoa%>%group_by(IN1_gene_ID)%>%dplyr::summarise(ncRNA_rate=sum(overall_T_ratio), weighted_average_length=sum(overall_T_ratio2*transcript_length))
transcript_infob$Novel_geneClass <- "others"
transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5)] <- "ncRNA"
transcript_infob$Novel_geneClass[which(transcript_infob$ncRNA_rate > 0.5 & transcript_infob$weighted_average_length > 200)] <- "lncRNA"
transcript_info_final <- left_join(transcript_info_final,transcript_infob[,c(1,4)], by="IN1_gene_ID", copy=F)
transcript_info_final$Novel_geneClass[which(is.na(transcript_info_final$Novel_geneClass) & transcript_info_final$gene_novelty == "Novel")] <- "others"

write.table(transcript_info_final, gzfile(paste0(path2,"Neuron_THP1.S3.table0_raw.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#=====
##Ref gene with re-adjusted 3/5'end
gene_model <- read.delim(paste0(gene0_path,"iPSC_NSC_Neuron.S3.t0.gene/bed/iPSC_NSC_Neuron.S3.t0.gene.gene.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
gene_model$V4 <- sapply(strsplit(gene_model$V4,"\\|"),"[",1)
gencode_gene <- read.delim("/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/resources/GENCODE_V39/gene.bed.gz", header=F, stringsAsFactors = F, check.names = F)
gencode_gene <- left_join(gencode_gene[,c(1:6)], gene_model[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new")) #61533
gencode_gene <- gencode_gene[which(!is.na(gencode_gene$V1_new)),]
gencode_gene$n5_adjust <- "5n_adjust"
gencode_gene$n3_adjust <- "3n_adjust"
gencode_gene$n5_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="+"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="-"))] <- NA
gencode_gene$n3_adjust[union(which(gencode_gene$V2_ori==gencode_gene$V2_new & gencode_gene$V6_ori=="-"), which(gencode_gene$V3_ori==gencode_gene$V3_new & gencode_gene$V6_ori=="+"))] <- NA
gencode_gene2 <- reshape2::melt(gencode_gene[,c(4,12,13)], id=1)
gencode_gene2 <- gencode_gene2[which(!is.na(gencode_gene2$value)),]
gencode_gene2 <- gencode_gene2%>%group_by(V4)%>%dplyr::summarise(Ref_gene_adjust=paste(value,collapse=" & "))
gencode_gene <- left_join(gencode_gene,gencode_gene2, by="V4", copy=F)
gencode_gene$Ref_gene_adjust[is.na(gencode_gene$Ref_gene_adjust)] <- "No"
transcript_info_final <- left_join(transcript_info_final,gencode_gene[,c(4,14)], by=c("IN1_gene_ID" = "V4"), copy=F)

#=====
##Ref transcript with re-adjusted 3/5'end
transcript_info1aa <- transcript_info_final[which(transcript_info_final$ref_source == "fulllength_ref"),]
transcript_model <- read.delim("Neuron_THP1.S3.model.bed.bgz", header=F, stringsAsFactors = F, check.names = F)
transcript_model_genecode <- transcript_model[which(transcript_model$V4 %in% transcript_info1aa$model_ID),]
gencode_tran <- read.delim("/osc-fs_home/yip/CFC_seq_paper_fig_data/code_n_data/SALA/SALA_from_github/resources/GENCODE_V39/transcript.bed.bgz", header=F, stringsAsFactors = F, check.names = F)
gencode_tran1 <- gencode_tran[which(gencode_tran$V4 %in% transcript_info1aa$model_ID),]
gencode_tran1 <- left_join(gencode_tran1[,c(1:6)], transcript_model_genecode[,c(1:6)], by="V4", copy=F, suffix=c("_ori","_new"))
gencode_tran1$n5_adjust <- "5n_adjust"
gencode_tran1$n3_adjust <- "3n_adjust"
gencode_tran1$n5_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="+"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="-"))] <- NA
gencode_tran1$n3_adjust[union(which(gencode_tran1$V2_ori==gencode_tran1$V2_new & gencode_tran1$V6_ori=="-"), which(gencode_tran1$V3_ori==gencode_tran1$V3_new & gencode_tran1$V6_ori=="+"))] <- NA
gencode_tran2 <- reshape2::melt(gencode_tran1[,c(4,12,13)], id=1)
gencode_tran2 <- gencode_tran2[which(!is.na(gencode_tran2$value)),]
gencode_tran2 <- gencode_tran2%>%group_by(V4)%>%dplyr::summarise(Ref_transcript_adjust=paste(value,collapse=" & "))
gencode_tran1 <- left_join(gencode_tran1,gencode_tran2, by="V4", copy=F)
gencode_tran1$Ref_transcript_adjust[is.na(gencode_tran1$Ref_transcript_adjust)] <- "No"
transcript_info_final <- left_join(transcript_info_final,gencode_tran1[,c(4,14)], by=c("model_ID" = "V4"), copy=F)
write.table(transcript_info_final, gzfile(paste0(path2,"Neuron_THP1.S3.table0_raw.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
rm(gene_model,gencode_gene,transcript_model,gencode_tran)

#build gtf
#=============
path2=paste0(SALA_path,"default/transcript/log/")
path3=paste0(SALA_path,"default/transcript/bed/")

setwd(path3)
gene0_path=paste0(SALA_path,"default/table0_gene/")

gencode.gtf=fread(paste0(GENCODE_path,"gencode.v39.annotation.gtf.gz"), header=F, stringsAsFactors = F)

table0.bed12=read.delim("Neuron_THP1.S3.model.bed.bgz",header=F, stringsAsFactors = F, check.names = F)
table0.bed6=read.delim("Neuron_THP1.S3.model.bed6.bed.gz",header=F, stringsAsFactors = F, check.names = F)
table0.genebed=read.delim(paste0(gene0_path,"iPSC_NSC_Neuron.S3.t0.gene/bed/iPSC_NSC_Neuron.S3.t0.gene.gene.bed.bgz"),header=F, stringsAsFactors = F, check.names = F)
transcript_info1a=read.delim(paste0(path2,"Neuron_THP1.S3.table0_raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
transcript_info1a=transcript_info1a[which(transcript_info1a$internal_priming != "Yes"),]

transcript_info1aa=transcript_info1a[which(transcript_info1a$Ref_transcript_adjust != "No" & !is.na(transcript_info1a$Ref_transcript_adjust)),]
transcript_info1bb=transcript_info1a[which(transcript_info1a$Ref_gene_adjust != "No" & !is.na(transcript_info1a$Ref_gene_adjust)),]

gencode.gtf.g=gencode.gtf[which(gencode.gtf$V3 == "gene"),]
gencode.gtf.t=gencode.gtf[which(gencode.gtf$V3 == "transcript"),]
gencode.gtf.e=gencode.gtf[which(gencode.gtf$V3 == "exon"),]

library(tidyverse)
gencode.gtf.g$gene_ID=sapply(strsplit(gencode.gtf.g$V9," "),"[",2)
gencode.gtf.g$gene_ID=gsub(";","",gencode.gtf.g$gene_ID)
gencode.gtf.g$gene_type=sapply(strsplit(gencode.gtf.g$V9," "),"[",4)
gencode.gtf.g$gene_type=gsub(";","",gencode.gtf.g$gene_type)
gencode.gtf.g$gene_name=sapply(strsplit(gencode.gtf.g$V9," "),"[",6)
gencode.gtf.g$gene_name=gsub(";","",gencode.gtf.g$gene_name)
gencode.gtf.g=gencode.gtf.g %>% mutate(across(c(gene_ID,gene_type,gene_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ga=gencode.gtf.g[-which(gencode.gtf.g$gene_ID %in% unique(transcript_info1bb$IN1_gene_ID)),] 
gencode.gtf.ga=gencode.gtf.ga[,c(1:8,10:12)] #component1
gencode.gtf.ga$V2="GENCODE"
gencode.gtf.ga$gene_novelty="GENCODEv39"

table0.genebed$gene_ID=sapply(strsplit(table0.genebed$V4,"\\|"),"[",1)
table1.genebed=table0.genebed[which(table0.genebed$gene_ID %in% unique(transcript_info1a$IN1_gene_ID)),]
table1.genebed$label="ONTCAGE"
table1.genebed$type="gene"
table1.genebed=table1.genebed[,c(1,14,15,2,3,5,6,5,13)]
table1.genebed$V2=table1.genebed$V2+1
table1.genebed$V5="."
table1.genebed$V5.1="."
table1a.genebed=table1.genebed[-grep("ENSG",table1.genebed$gene_ID),]
table1a.genebed=left_join(table1a.genebed,unique(transcript_info1a[,c("IN1_gene_ID","Novel_geneClass")]), by=c("gene_ID"="IN1_gene_ID"), copy=F)
table1a.genebed$gene_name=sapply(strsplit(table1a.genebed$gene_ID,"\\."),"[",1)
table1a.genebed$gene_novelty="novel"
table1b.genebed=table1.genebed[which(table1.genebed$gene_ID %in% unique(transcript_info1bb$IN1_gene_ID)),]
table1b.genebed=left_join(table1b.genebed,gencode.gtf.g[,c(10:12)],by="gene_ID",copy=F)
table1b.genebed$novelty="GENCODE_updated"
table1b.genebed$label="GENCODE"
colnames(table1a.genebed)=colnames(gencode.gtf.ga)
colnames(table1b.genebed)=colnames(gencode.gtf.ga)
table1.genebed=rbind(table1a.genebed,table1b.genebed,gencode.gtf.ga)
table1.genebed$V9=paste0("gene_id \"",table1.genebed$gene_ID,"\"; gene_type \"",table1.genebed$gene_type,"\"; gene_name \"",table1.genebed$gene_name,"\"; gene_novelty \"",table1.genebed$gene_novelty,"\";")
#write.table(table1.genebed,gzfile(paste0(path3,"table1pENST.gene.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)
#

gencode.gtf.t$transcript_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",4)
gencode.gtf.t$transcript_ID=gsub(";","",gencode.gtf.t$transcript_ID)
gencode.gtf.t$gene_ID=sapply(strsplit(gencode.gtf.t$V9," "),"[",2)
gencode.gtf.t$gene_ID=gsub(";","",gencode.gtf.t$gene_ID)
gencode.gtf.t$transcript_type=sapply(strsplit(gencode.gtf.t$V9," "),"[",10)
gencode.gtf.t$transcript_type=gsub(";","",gencode.gtf.t$transcript_type)
gencode.gtf.t$transcript_name=sapply(strsplit(gencode.gtf.t$V9," "),"[",12)
gencode.gtf.t$transcript_name=gsub(";","",gencode.gtf.t$transcript_name)
gencode.gtf.t=gencode.gtf.t %>% mutate(across(c(gene_ID,transcript_ID,transcript_type,transcript_name),~ map_chr(.x, ~ gsub("\"", "", .x))))
gencode.gtf.ta=gencode.gtf.t[-which(gencode.gtf.t$transcript_ID %in% transcript_info1aa$model_ID),] 
gencode.gtf.ta=gencode.gtf.ta[,c(1:8,10:13)] #component1
gencode.gtf.ta$V2="GENCODE"
gencode.gtf.ta$transcript_novelty="GENCODEv39"

table1.bed12=table0.bed12[which(table0.bed12$V4 %in% transcript_info1a$model_ID),]
table1.bed12$label="ONTCAGE"
table1.bed12$type="transcript"
table1.bed12=table1.bed12[,c(1,13,14,2,3,5,6,5,4)]
table1.bed12$V2=table1.bed12$V2+1
table1.bed12$V5="."
table1.bed12$V5.1="."
table1a.bed12=table1.bed12[-grep("ENST",table1.bed12$V4),]
table1a.bed12=left_join(table1a.bed12, transcript_info1a[,c("model_ID","IN1_gene_ID","Novel_transcriptClass")], by=c("V4"="model_ID"),copy=F)
table1a.bed12$transcript_name=sapply(strsplit(table1a.bed12$V4,"\\."),"[",1)
table1a.bed12$transcript_novelty="novel"
table1b.bed12=table1.bed12[which(table1.bed12$V4 %in% transcript_info1aa$model_ID),]
table1b.bed12=left_join(table1b.bed12,gencode.gtf.t[,c(10:13)],by=c("V4"="transcript_ID"),copy=F) #take geneID, transcript_type, transcript_name
table1b.bed12$transcript_novelty="GENCODE_updated"
table1b.bed12$label="GENCODE"
colnames(table1a.bed12)=colnames(gencode.gtf.ta)
colnames(table1b.bed12)=colnames(gencode.gtf.ta)
table1.bed12=rbind(table1a.bed12,table1b.bed12,gencode.gtf.ta)
table1.bed12=left_join(table1.bed12,table1.genebed[,c(9:12)], by="gene_ID", copy=F) #get back gene info from last table
table1.bed12$V9=paste0("gene_id \"",table1.bed12$gene_ID,"\"; transcript_id \"",table1.bed12$transcript_ID,"\"; gene_type \"",table1.bed12$gene_type,"\"; gene_name \"",table1.bed12$gene_name,"\"; transcript_type \"",table1.bed12$transcript_type,"\"; transcript_name \"",table1.bed12$transcript_name,"\"; gene_novelty \"",table1.bed12$gene_novelty,"\"; transcript_novelty \"",table1.bed12$transcript_novelty,"\";")
#write.table(table1.bed12,gzfile(paste0(path3,"table1pENST.transcript.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

gencode.gtf.e$transcript_ID=sapply(strsplit(gencode.gtf.e$V9," "),"[",4)
gencode.gtf.e$transcript_ID=gsub(";","",gencode.gtf.e$transcript_ID)
gencode.gtf.e$exon_number=sapply(strsplit(gencode.gtf.e$V9," "),"[",14)
gencode.gtf.e$exon_number=as.numeric(gsub(";","",gencode.gtf.e$exon_number))
gencode.gtf.e$exon_id=sapply(strsplit(gencode.gtf.e$V9," "),"[",16)
gencode.gtf.e$exon_id=gsub(";","",gencode.gtf.e$exon_id)
gencode.gtf.e=gencode.gtf.e %>% mutate(across(c(transcript_ID,exon_id),~ map_chr(.x, ~ gsub("\"", "", .x))))
exon.list=unique(gencode.gtf.e[,c(10:12)])
gencode.gtf.ea=gencode.gtf.e[-which(gencode.gtf.e$transcript_ID %in% transcript_info1aa$model_ID),]  
gencode.gtf.ea=gencode.gtf.ea[,c(1:8,10:11)] #component1

table1.bed6=table0.bed6[which(table0.bed6$V4 %in% transcript_info1a$model_ID),]
table1.bed6$label="ONTCAGE"
table1.bed6$type="exon"
table1.bed6=table1.bed6[,c(1,7,8,2,3,5,6,5,4)]
table1.bed6$V2=table1.bed6$V2+1
table1.bed6$V5="."
table1.bed6$V5.1="."
table1.bed6p=table1.bed6[which(table1.bed6$V6 == "+"),]%>%group_by(V4)%>%dplyr::arrange(V2)%>%dplyr::mutate(exon_number = 1: n())
table1.bed6n=table1.bed6[which(table1.bed6$V6 == "-"),]%>%group_by(V4)%>%dplyr::arrange(desc(V2))%>%dplyr::mutate(exon_number = 1: n())
table1.bed6=rbind(table1.bed6p,table1.bed6n)

table1a.bed6=table1.bed6[-grep("ENST",table1.bed6$V4),]
table1b.bed6=table1.bed6[which(table1.bed6$V4 %in% transcript_info1aa$model_ID),]
colnames(table1a.bed6)=colnames(gencode.gtf.ea)
colnames(table1b.bed6)=colnames(gencode.gtf.ea)
table1.bed6=rbind(table1a.bed6,table1b.bed6,gencode.gtf.ea)

table1.bed6=left_join(table1.bed6, table1.bed12[,c(9:12,14:15)], by="transcript_ID", copy=F)
table1.bed6$V9=paste0("gene_id \"",table1.bed6$gene_ID,"\"; transcript_id \"",table1.bed6$transcript_ID,"\"; gene_type \"",table1.bed6$gene_type,"\"; gene_name \"",table1.bed6$gene_name,"\"; transcript_type \"",table1.bed6$transcript_type,"\"; transcript_name \"",table1.bed6$transcript_name,"\"; exon_number \"",table1.bed6$exon_number,"\";")

options(scipen=999)
final=rbind(table1.genebed[,c(1:8,13)],table1.bed12[,c(1:8,17)], table1.bed6[,c(1:8,16)])
final=final[order(final$V1,final$V4),]
write.table(final,gzfile(paste0(path2,"table0.All_Ref.gtf.gz")), col.names=F, row.names=F, sep="\t", quote=F)

setwd(path2)
system("zcat table0.All_Ref.gtf.gz | bedparse gtf2bed | gzip > table0.All_Ref.bed12.bed.gz")


#===============================================================================
path_supp=paste0(primary_folder,"supplementary_table/")

table5=read.delim(paste0(path2,"Neuron_THP1.S3.table4_filtered.noIP.All_Ref_updated.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
colnames(table5)[c(48,49)]=c("THP1","dTHP1")
table5=table5[which(table5$ref_source %in% c("fulllength_ref","novel_transcript")),]
write.table(table5,paste0(path_supp,"TableS6_SALA_default_Tx.tsv"), col.names=T, row.names=F, sep="\t",quote=F)
# -> corresponding gtf: table4.Detected_Fulllength_Ref.gtf.gz

#===============================================================================
# all tmp folders were removed,
# all read base info and bed file and fasta files were removed

