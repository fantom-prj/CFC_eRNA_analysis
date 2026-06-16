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

################
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig3_data=paste0(primary_folder,"fig3/data/")
CAT_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_FANTOMCAT/")
lncbook_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_lncbook/")
gencode_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_gencodev47/")
refseq_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/Compare_refseq2024/")

output_path=paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/SALA_compare_exisiting_databases/compare_4Database/")

#===============================================================================
#number of transcripts
table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
cfc_anno=table5[grep("ONT",table5$model_ID),c(1,86,89,93,102)]
cfc_anno$group[grep("ncRNA",cfc_anno$group)]="ncRNA"
cfc_anno$group[which(cfc_anno$group != "ncRNA")]="others"
cfc_anno1=cfc_anno%>%group_by(group)%>%dplyr::summarise(count=n())

FCAT_anno=read.delim(paste0(CAT_path,"transcript_gene_link.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
FCAT_anno%>%group_by(geneClass)%>%dplyr::summarise(count=n())
FCAT_anno$geneClass[which(FCAT_anno$geneClass %in% c("sense_overlap_RNA", "short_ncRNA", "lncRNA_antisense", "lncRNA_divergent", "lncRNA_intergenic", "lncRNA_sense_intronic"))]="ncRNA"
FCAT_anno$geneClass[which(FCAT_anno$geneClass != "ncRNA")]="others"
FCAT_anno1=FCAT_anno%>%group_by(geneClass)%>%dplyr::summarise(count=n())

refseq2024_anno=read.delim(paste0(refseq_path,"gene_transcript_link.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
refseq2024_anno$gbkey[which(refseq2024_anno$gbkey != "ncRNA")]="others"
refseq2024_anno1=refseq2024_anno%>%group_by(gbkey)%>%dplyr::summarise(count=n())

lncbook_anno=read.delim(paste0(lncbook_path,"gene_transcript.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
lncbook_anno1=lncbook_anno%>%group_by(transcript_type)%>%dplyr::summarise(count=n())

gencodev47_anno=read.delim(paste0(gencode_path,"gene_transcript_link.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
aa=gencodev47_anno%>%group_by(transcriptType)%>%dplyr::summarise(count=n())
gencodev47_anno$transcriptType[which(gencodev47_anno$transcriptType == "lncRNA")]="ncRNA"
gencodev47_anno$transcriptType[which(gencodev47_anno$transcriptType != "ncRNA")]="others"
gencodev47_anno1=gencodev47_anno%>%group_by(transcriptType)%>%dplyr::summarise(count=n())

Tnumber=data.frame(cbind(database=c("CFC_novel","FANTOM_CAT","LncBook","Refseq","GENCODEv47"),
                         transcript=c(145951,705472,323945,181361,385659),
                         ncRNA=c(109562,205315,323945,34304,189177)))
Tnumber$others=as.numeric(Tnumber$transcript)-as.numeric(Tnumber$ncRNA)
Tnumber1=reshape2::melt(Tnumber, id=c(1,2))
write.table(Tnumber1,gzfile(paste0(path_fig3_data,"transcript_number_5dataset.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#put together FCAT & lncbook & refseq & GENCODEv47
setwd(output_path)
CATT=read.delim("ONTT_FCATT_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
CATG=read.delim("ENSG_ONTG_FCATT_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

lncbookT=read.delim("ONTT_lncbookT_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
lncbookG=read.delim("ENSG_ONTG_lncbookT_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

refseq2024T=read.delim("ONTT_refseq2024T_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
colnames(refseq2024T)[2]="refseq2024T_transcriptID"
refseq2024G=read.delim("ENSG_ONTG_refseq2024T_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

gencodev47T=read.delim("ONTT_gencodev47T_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)
colnames(gencodev47T)[2]="gencodev47T_transcriptID"
gencodev47G=read.delim("ENSG_ONTG_gencodev47T_match_transcript.tsv.gz", header=T, stringsAsFactors = F, check.names = F)

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)

transcript=table5[grep("ONTT",table5$model_ID),c(1,102)]
transcript$group[grep("ncRNA",transcript$group)]="ncRNA"
transcript$group[which(transcript$group != "ncRNA")]="others"
transcript%>%group_by(group)%>%dplyr::summarise(count=n())
transcript=left_join(transcript,CATT,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,lncbookT,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,refseq2024T,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,gencodev47T,by=c("model_ID"="model_ID_str"),copy=F)
transcript[is.na(transcript)]=0
indices <- transcript[, 3:6] != 0
transcript[, c(3:6)][indices] = 1
colnames(transcript)[c(3:6)]=c("FANTOM_CAT","LncBook","Refseq","GENCODEv47")
transcript$CFC_novel=1
transcript$same_model=""
write.table(transcript, gzfile(paste0(path_fig3_data,"transcript.upsetter.tsv.gz")), col.names=T, row.names=F, quote=F, sep="\t")
#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5a

gene=unique(table5[grep("ONTG",table5$T4_gene_ID),c(86,93,99)])
gene=left_join(gene,CATG,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,lncbookG,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,refseq2024G,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,gencodev47G,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene[is.na(gene)]=0
indices <- gene[, 4:7] != 0
gene[, c(4:7)][indices] = 1
colnames(gene)[c(4:7)]=c("FANTOM_CAT","LncBook","Refseq","GENCODEv47")
gene$CFC_novel=1
gene$same_model=""
write.table(gene, gzfile(paste0(path_fig3_data,"gene.upsetter.tsv.gz")), col.names=T, row.names=F, quote=F, sep="\t")
#stored in [primary_folder]/fig3/data
# for -> Fig. Ext5b


#===============================================================================
#linked annotated transcript/gene ID to novel transcript/gene
FCAT_anno=read.delim(paste0(CAT_path,"transcript_gene_link.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
FCAT_anno$geneClass[which(FCAT_anno$geneClass %in% c("sense_overlap_RNA", "short_ncRNA", "lncRNA_antisense", "lncRNA_divergent", "lncRNA_intergenic", "lncRNA_sense_intronic"))]="ncRNA"
FCAT_anno$geneClass[which(FCAT_anno$geneClass != "ncRNA")]="others"
lncbook_anno=read.delim(paste0(lncbook_path,"gene_transcript.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
lncbook_anno$transcript_type="ncRNA"
refseq2024_anno=read.delim(paste0(refseq_path,"gene_transcript_link.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
refseq2024_anno$gbkey[which(refseq2024_anno$gbkey != "ncRNA")]="others"
refseq2024_anno1=refseq2024_anno%>%group_by(gbkey)%>%dplyr::summarise(count=n())
gencodev47_anno=read.delim(paste0(gencode_path,"gene_transcript_link.tsv.gz"), header=T, check.names = F,stringsAsFactors = F)
gencodev47_anno$transcriptType[which(gencodev47_anno$transcriptType == "lncRNA")]="ncRNA"
gencodev47_anno$transcriptType[which(gencodev47_anno$transcriptType != "ncRNA")]="others"
colnames(FCAT_anno)[3]="transcriptType"
colnames(lncbook_anno)[6]="transcriptType"
colnames(refseq2024_anno)[4]="transcriptType"
all_anno=rbind(FCAT_anno[,c(1,3)],lncbook_anno[,c(2,6)],refseq2024_anno[,c(2,4)],gencodev47_anno[,c(2,3)])

transcript=table5[grep("ONTT",table5$model_ID),c(1,102)]
transcript$group[grep("ncRNA",transcript$group)]="ncRNA"
transcript$group[which(transcript$group != "ncRNA")]="others"
transcript=left_join(transcript,CATT,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,lncbookT,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,refseq2024T,by=c("model_ID"="model_ID_str"),copy=F)
transcript=left_join(transcript,gencodev47T,by=c("model_ID"="model_ID_str"),copy=F)
transcript1=reshape2::melt(transcript[,c(1,3:6)], id=1)
transcript1=transcript1[which(!is.na(transcript1$value)),]
transcript1$variable=gsub("_transcriptID","",as.character(transcript1$variable))
transcript1=separate_rows(transcript1, value, sep=";")
transcript1=left_join(transcript1, all_anno, by=c("value"="transcriptID"),copy=F)
transcript2=transcript1%>%group_by(model_ID)%>%dplyr::summarise(database=paste(variable,collapse=";"), transcriptID=paste(value,collapse=";"), transcriptType=paste(transcriptType,collapse=";"))
transcript=left_join(transcript, transcript2, by="model_ID",copy=F)
transcript$database[which(is.na(transcript$database))]="Novel"
#transcript3=transcript%>%group_by(database)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

write.table(transcript1, gzfile(paste0(output_path,"all.match.transcript_long_format.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
transcript1=read.delim(paste0(output_path,"all.match.transcript_long_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
transcript5=transcript1%>%group_by(model_ID, variable)%>%dplyr::summarise(value=paste(value,collapse=";"))
transcript6=spread(transcript5, key=2, value=3)
transcript6=transcript6[,c(1,2,4,5,3)]
colnames(transcript6)[c(2:5)]=c("FANTOM_CAT","LncBook","Refseq","GENCODEv47")
write.table(transcript6, gzfile(paste0(output_path,"all.match.transcript_short_format.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig3/data
# for -> supplementary table S7


#do the same for gene
gene=unique(table5[grep("ONTG",table5$T4_gene_ID),c(86,93)])
gene=left_join(gene,CATG,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,lncbookG,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,refseq2024G,by=c("T4_gene_ID"="gene_ID"),copy=F)
gene=left_join(gene,gencodev47G,by=c("T4_gene_ID"="gene_ID"),copy=F)
# dont use gene annotation
gene1=reshape2::melt(gene[,c(1,3:6)], id=1)
gene1=gene1[which(!is.na(gene1$value)),]
gene1$variable=gsub("_transcriptID","",as.character(gene1$variable))
gene2=spread(gene1, key=2, value=3)
gene2=gene2[,c(1,2,4,5,3)]
colnames(gene2)[c(2:5)]=c("FANTOM_CAT","LncBook","Refseq","GENCODEv47")
write.table(gene2, gzfile("all.match.gene_short_format.tsv.gz"), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig3/data
# for -> supplementary table S7

#========
#also make one for gene to gene
gene2=read.delim(paste0(output_path,"all.match.gene_short_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene3=reshape2::melt(gene2, id=1)
gene3=separate_rows(gene3,value,sep=";")
gene3a=left_join(gene3[which(gene3$variable=="FANTOM_CAT" & !is.na(gene3$value)),], FCAT_anno[,c(1:2)], by=c("value"="transcriptID"),copy=F)
gene3b=left_join(gene3[which(gene3$variable=="LncBook" & !is.na(gene3$value)),], lncbook_anno[,c(1:2)], by=c("value"="transcriptID"),copy=F)
gene3c=left_join(gene3[which(gene3$variable=="Refseq" & !is.na(gene3$value)),], refseq2024_anno[,c(1:2)], by=c("value"="transcriptID"),copy=F)
gene3d=left_join(gene3[which(gene3$variable=="GENCODEv47" & !is.na(gene3$value)),], gencodev47_anno[,c(1:2)], by=c("value"="transcriptID"),copy=F)

gene3a=gene3a%>%group_by(T4_gene_ID)%>%dplyr::summarise(FANTOM_CAT=paste(unique(geneID),collapse=";"))
gene3b=gene3b%>%group_by(T4_gene_ID)%>%dplyr::summarise(LncBook=paste(unique(geneID),collapse=";"))
gene3c=gene3c%>%group_by(T4_gene_ID)%>%dplyr::summarise(Refseq=paste(unique(geneID),collapse=";"))
gene3d=gene3d%>%group_by(T4_gene_ID)%>%dplyr::summarise(GENCODEv47=paste(unique(geneID),collapse=";"))

gene4=full_join(gene3a,gene3b, by="T4_gene_ID", copy=F)
gene4=full_join(gene4,gene3c, by="T4_gene_ID", copy=F)
gene4=full_join(gene4,gene3d, by="T4_gene_ID", copy=F)
write.table(gene4[order(gene4$T4_gene_ID),],gzfile(paste0(output_path,"all.match.genetogene_short_format.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#=======================================================================
#get expression matrix from bambu
bambu_tcpm_new=read.delim(paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/bambu_long_t5_partialYes.ENST/CPM_transcript.txt.gz"), header=T, stringsAsFactors = F, check.names = F)

bambu_tcpm_new=bambu_tcpm_new[,-2]
bambu_tcpm_new[is.na(bambu_tcpm_new)]=0
bambu_tcpm_new$iPSC=rowMeans(bambu_tcpm_new[,c(2:3)])
bambu_tcpm_new$NSC=rowMeans(bambu_tcpm_new[,c(6:7)])
bambu_tcpm_new$Neuron=rowMeans(bambu_tcpm_new[,c(4:5)])
bambu_tcpm_new$THP1=rowMeans(bambu_tcpm_new[,c(8:11,16:19)])
bambu_tcpm_new$dTHP1=rowMeans(bambu_tcpm_new[,c(12:15, 20:23)])
t_to_g=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.transcript_gene_link.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
bambu_tcpm_new=left_join(bambu_tcpm_new, t_to_g, by=c("TXNAME"="transcriptID"))
bambu_tcpm_new1=bambu_tcpm_new[,c(29,24:28)]%>%group_by(geneID)%>%summarise(across(everything(), sum, na.rm = TRUE))

#calculate gini index
library(DescTools)

bambu_tcpm_new1$gini=apply(bambu_tcpm_new1[,c(2:6)], 1, Gini)
bambu_tcpm_new1$maxExp=apply(bambu_tcpm_new1[,c(2:6)], 1, max)
bambu_tcpm_new1=bambu_tcpm_new1[which(bambu_tcpm_new1$geneID %in% gene$T4_gene_ID),]
bambu_tcpm_new1$group="Un-annotated"
bambu_tcpm_new1$group[which(bambu_tcpm_new1$geneID %in% gene2$T4_gene_ID)]="Annotated"
bambu_tcpm_new1%>%group_by(group)%>%dplyr::summarise(median_gini=median(gini,na.rm=T), median_maxExp=median(maxExp))
write.table(bambu_tcpm_new1, gzfile(paste0(path_fig3_data,"bambu_expression_annot_unannot.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)


#===============================================================================
path_supp=paste0(primary_folder,"supplementary_table/")
transcript6=read.delim(paste0(output_path,"all.match.transcript_short_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene2=read.delim(paste0(output_path,"all.match.gene_short_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene4=read.delim(paste0(output_path,"all.match.genetogene_short_format.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
colnames(transcript6)=c("CFCseq_transcriptID","FANTOM_CAT_transcriptID","LncBook_transcriptID","Refseq_transcriptID","GENCODEv47_transcriptID")
colnames(gene2)=c("CFCseq_geneID","FANTOM_CAT_transcriptID","LncBook_transcriptID","Refseq_transcriptID","GENCODEv47_transcriptID")
colnames(gene4)=c("CFCseq_geneID","FANTOM_CAT_geneID","LncBook_geneID","Refseq_geneID","GENCODEv47_geneID")
write.table(transcript6,gzfile(paste0(path_supp,"TableS7a_novelT.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
write.table(gene2,gzfile(paste0(path_supp,"TableS7b_novelG.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
write.table(gene4,gzfile(paste0(path_supp,"TableS7c_novelG_G.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

#stored in [primary_folder]/fig3/data
# for -> supplementary table S7




