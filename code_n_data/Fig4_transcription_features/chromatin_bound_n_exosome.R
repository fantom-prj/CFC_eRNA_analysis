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

#===============================================================================
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig4_data=paste0(primary_folder,"fig4/data/")
exo_path=paste0(primary_folder,"code_n_data/Fig4_transcription_features/exosome_sensitivity/")
chromatin_path=paste0(primary_folder,"code_n_data/Fig4_transcription_features/quantification_chr_total/")
SALA_path=paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/transcript/")
RBP_path=paste0(exo_path,"RBP_result/")
SALA_foler=paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/transcript/")

#quantification of chromatin-bound and total RNA from iPSC
#===============================================================================
library(bambu)
#run bambu
#bam files please download from DDBJ
results_dir=paste0(chromatin_path,"bambu")
fa.file <- "/analysisdata/fantom6/Interactome/CFC_THP1new/SALA/resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
bam <- read.delim("/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/exosome_sensitivity/iPSC_total_chr_TC_bam.path.txt", header = F, stringsAsFactors = F, check.names = F)
test.bam <- bam$V2
annotations <- prepareAnnotations(paste0(SALA_path,"gtf/table5_partial_yes_detected.alone_allNeuron_THP1t5.gtf.gz"))
ID_link=unique(data.frame(cbind(mcols(annotations)$GENEID, bambu:::assignGeneIds(annotations, GRangesList())$GENEID)))
ID_link=ID_link%>%group_by(X2)%>%dplyr::summarise(geneID=paste(X1,collapse=";"))
mcols(annotations)$GENEID <- bambu:::assignGeneIds(annotations, GRangesList())$GENEID
dir.create(paste0(results_dir,"/rcOut"), recursive=TRUE)
print("starting running bambu...")
se <- bambu(reads =  test.bam, 
            annotations = annotations, 
            genome = fa.file, 
            discovery = FALSE, 
            opt.discovery = list(min.exonDistance = 0), 
            rcOutDir = paste0(results_dir,"/rcOut"),returnDistTable=TRUE)
saveRDS(se, paste0(results_dir,"/se.rds"))
writeBambuOutput(se, results_dir)
print(paste0("finish running bambu. results located in ",results_dir))

gene_count <- read.delim(paste0(results_dir,"/counts_gene.txt"), header=T)
transcript_count <- read.delim(paste0(results_dir,"/counts_transcript.txt"), header=T)
transcript_CMP <- read.delim(paste0(results_dir,"/CPM_transcript.txt"), header=T)
gene_count <- left_join(gene_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_count <- left_join(transcript_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_CMP <- left_join(transcript_CMP, ID_link, by=c("GENEID"="X2"), copy=F)
write.table(gene_count, paste0(results_dir,"/counts_gene.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(transcript_count, paste0(results_dir,"/counts_transcript.txt"), col.names=T, row.names=F, sep="\t", quote=F)
write.table(transcript_CMP, paste0(results_dir,"/CPM_transcript.txt"), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/code_n_data/Fig4_transcription_features/quantification_chr_total/bambu

#===============================================================================
#edgeR
transcript_count <- read.delim(paste0(results_dir,"/counts_transcript.txt"), header=T)
transcript_count$iPSC_total_r1=rowSums(transcript_count[,c(3:5)])
transcript_count$iPSC_total_r2=rowSums(transcript_count[,c(6:8)])
transcript_count$iPSC_chr_r1=transcript_count$WTC.11.NGN2.hiPSC_chromatin_1_rep1_labeled.sorted
transcript_count$iPSC_chr_r2=transcript_count$WTC.11.NGN2.hiPSC_chromatin_1_rep2_labeled.sorted

link=read.delim(paste0(SALA_path,"gtf/table5_partial_yes_detected.alone_allNeuron_THP1t5.transcript_gene_link.tsv"),header=T, stringsAsFactors = F, check.names = F)
transcript_count=left_join(transcript_count[,c(1,12:15)], link, by=c("TXNAME"="transcriptID"), copy=F)
gene_count=transcript_count%>%group_by(geneID)%>%dplyr::summarise(across(2:5, sum))
gene_count=gene_count[which(rowSums(gene_count[,c(2:5)] >0)>1),]
gene_count=data.frame(gene_count)
gene_count=gene_count[which(!is.na(gene_count$geneID)),]
rownames(gene_count)=gene_count$geneID
gene_count=gene_count[,c(2:5)]

transcript_count=transcript_count[which(rowSums(transcript_count[,c(2:5)] >0)>1),]
rownames(transcript_count)=transcript_count$TXNAME
transcript_count=transcript_count[,c(2:5)]

TS=c(rep("total",2),rep("chr",2))
TS <- factor(TS, levels=c("total","chr"))
counttable1=gene_count
my_data1= DGEList(counts=counttable1, group=TS)
my_data1 <- calcNormFactors(my_data1, method="RLE")
my_data1 <- estimateDisp(my_data1)
my_data1$common.dispersion
plotBCV(my_data1)
levels(my_data1$samples$group)
lrt = exactTest(my_data1, pair=c("total","chr"))
lrt$table$FDR=p.adjust(lrt$table$PValue, method="BH")
exp.df = lrt$table
exp.df$geneID=row.names(exp.df)
gene_count$geneID=rownames(gene_count)
gene_count=left_join(gene_count, exp.df, by="geneID", copy=F)
write.table(gene_count,paste0(chromatin_path,"edgeR_bambu_genebase_total_chr_ipsc.tsv"),col.names=T, row.names=F, sep="\t", quote=F)


counttable1=transcript_count
my_data1= DGEList(counts=counttable1, group=TS)
my_data1 <- calcNormFactors(my_data1, method="RLE")
my_data1 <- estimateDisp(my_data1)
my_data1$common.dispersion
plotBCV(my_data1)
levels(my_data1$samples$group)
lrt = exactTest(my_data1, pair=c("total","chr"))
lrt$table$FDR=p.adjust(lrt$table$PValue, method="BH")
exp.df = lrt$table
exp.df$transcriptID=row.names(exp.df)
transcript_count$transcriptID=rownames(transcript_count)
transcript_count=left_join(transcript_count, exp.df, by="transcriptID", copy=F)
write.table(transcript_count,paste0(chromatin_path,"edgeR_bambu_Txbase_total_chr_ipsc.tsv"),col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/code_n_data/Fig4_transcription_features/quantification_chr_total


#===================================
#for transcript group
chr_table5=read.delim(paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/transcript/log/table5.chimeric.76K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
unique(chr_table5[,c(53,59,61)])%>%group_by(T4_gene_promoter_type_CHR,T4_Novel_geneClass_CHR)%>%dplyr::summarise(count=n())
unique(chr_table5[,c(53,57,61)])%>%group_by(T4_gene_promoter_type_CHR,T4_Gencode_ONT_geneClass2)%>%dplyr::summarise(count=n())

chr_table5%>%group_by(source)%>%dplyr::summarise(count=n())
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
table5a=table5[which(table5$T4_gene_ID %in% chr_table5$T4_gene_ID),]

chr_table5=left_join(chr_table5,unique(table5a[,c(86,99,90,93)]), by="T4_gene_ID",copy=F, suffix=c("_CHR","_NEURON_SERIES"))
chr_table5$gene_group="others"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "promoter-like"),grep("ncRNA",chr_table5$T4_Novel_geneClass_NEURON_SERIES))]="ONTG.p_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "promoter-like"),grep("ncRNA",chr_table5$T4_Gencode_geneCalss2))]="ENSG.p_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "enhancer-like"),grep("ncRNA",chr_table5$T4_Novel_geneClass_NEURON_SERIES))]="ONTG.e_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "enhancer-like"),grep("ncRNA",chr_table5$T4_Gencode_geneCalss2))]="ENSG.e_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "unclassed"),grep("ncRNA",chr_table5$T4_Novel_geneClass_NEURON_SERIES))]="ONTG.other_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_NEURON_SERIES == "unclassed"),grep("ncRNA",chr_table5$T4_Gencode_geneCalss2))]="ENSG.other_ncRNA"
#chr_table5$novel_gene_group="others"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_CHR == "promoter-like"),grep("ncRNA",chr_table5$T4_Novel_geneClass_CHR))]="CHRG.p_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_CHR == "enhancer-like"),grep("ncRNA",chr_table5$T4_Novel_geneClass_CHR))]="CHRG.e_ncRNA"
chr_table5$gene_group[intersect(which(chr_table5$T4_gene_promoter_type_CHR == "unclassed"),grep("ncRNA",chr_table5$T4_Novel_geneClass_CHR))]="CHRG.other_ncRNA"

chr_table5$transcript_group="others"
chr_table5$transcript_group[intersect(which(chr_table5$promoter_type== "promoter-like"),grep("ncRNA",chr_table5$T4_ncRNA_source))]="p_ncRNA"
chr_table5$transcript_group[intersect(which(chr_table5$promoter_type== "enhancer-like"),grep("ncRNA",chr_table5$T4_ncRNA_source))]="e_ncRNA"
chr_table5$transcript_group[intersect(which(chr_table5$promoter_type== "unclassed"),grep("ncRNA",chr_table5$T4_ncRNA_source))]="other_ncRNA"
chr_table5$transcript_group[which(chr_table5$Gencode_ONT_transcriptClass2== "protein_coding")]="mRNA"
chr_table5%>%group_by(transcript_group)%>%dplyr::summarise(count=n())



chr_table5=left_join(chr_table5,table5[,c(1,102)], by="model_ID",copy=F)
chr_table5$group[grep("CHR",chr_table5$model_ID)]=paste0(chr_table5$promoter_type[grep("CHR",chr_table5$model_ID)],"_",chr_table5$T4_ncRNA_source[grep("CHR",chr_table5$model_ID)])
chr_table5$group[grep("_NA",chr_table5$group)]="others"
chr_table5$group=gsub("_isoform","",chr_table5$group)
chr_table5$group=gsub("nhancer-like_novel_l","_",chr_table5$group)
chr_table5$group=gsub("romoter-like_novel_l","_",chr_table5$group)
chr_table5$group=gsub("nhancer-like_novel_","_",chr_table5$group)
chr_table5$group=gsub("romoter-like_novel_","_",chr_table5$group)
chr_table5$group=gsub("unclassed_novel_l","other_",chr_table5$group)
chr_table5$group=gsub("unclassed_novel_","other_",chr_table5$group)
#undetected from total dataset
chr_table5$Gencode_ONT_transcriptClass2 [which(is.na(chr_table5$group))]

chr_table5a=unique(chr_table5[,c(53,67)])%>%group_by(gene_group)%>%dplyr::summarise(count=n())
chr_table5a$group=sapply(strsplit(chr_table5a$gene_group,"\\."),"[",1)
chr_table5a$group2=sapply(strsplit(chr_table5a$gene_group,"\\."),"[",2)
chr_table5a$group3="Detected"
chr_table5a$group3[which(chr_table5a$group == "CHRG")]="Novel"
chr_table5a=chr_table5a[which(chr_table5a$group != "others"),]
chr_table5a$group=factor(chr_table5a$group,levels=c("ENSG","ONTG","CHRG"))
chr_table5a$group2=factor(chr_table5a$group2,levels=c("p_ncRNA","e_ncRNA","other_ncRNA"))
chr_table5a$label=chr_table5a$count
chr_table5a$label[which(chr_table5a$group3=="Detected")]=NA
write.table(chr_table5a,gzfile(paste0(path_fig4_data,"chr_bound_iPSC_gene_count.plot.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)
#for ex6k
#=====================================


#calculate exosome sensitivity 
#(TPMExosome-suppressed - TPMControl)/TPMExosome-suppressed
#=====================================
#bash
#kallisto
#bam file of EXOSC3 KD in iPSC can be found in DDBJ
setwd(exo_path)
system("perl kallisto_quantify_exosome.pl")
#=====================================

path1=paste0(exo_path,"kallisto_result_SALA_table5_partialYes.ENST_chro/")
files=list.files(pattern="abundance", path=path1, recursive =T)
files.names1=unique(sapply(strsplit(files, "\\/"),"[",1))

gene_transcript=read.delim(paste0(SALA_path,"gtf/table5_partial_yes_detected.alone_allNeuron_THP1t5.transcript_gene_link.tsv"), header=T, stringsAsFactors = F, check.names = F)
quant=read.delim(paste0(path1,files[1]), header=T, stringsAsFactors = F)
quant$target_id=sapply(strsplit(quant$target_id,"\\("),"[",1)
quant_tpm=quant[,c(1,5)]
quant_count=quant[,c(1,4)]
colnames(quant_tpm)[1+1]=files.names1[1]
colnames(quant_count)[1+1]=files.names1[1]
for (i in 2:length(files)){
  quant1=read.delim(paste0(path1,files[i]), header=T, stringsAsFactors = F)
  quant1$target_id=sapply(strsplit(quant1$target_id,"\\("),"[",1)
  quant_tpm=left_join(quant_tpm,quant1[,c(1,5)],by="target_id")
  quant_count=left_join(quant_count,quant1[,c(1,4)],by="target_id")
  colnames(quant_tpm)[1+i]=files.names1[i]
  colnames(quant_count)[1+i]=files.names1[i]}

setwd(exo_path)
write.table(quant_tpm,gzfile("ExoKD.transcript.tpm.matrix.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)
write.table(quant_count,gzfile("ExoKD.transcript.count.matrix.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)
quant_tpm=left_join(quant_tpm,gene_transcript, by=c("target_id"="transcriptID"),copy=F)
quant_count=left_join(quant_count,gene_transcript, by=c("target_id"="transcriptID"),copy=F)
quant_tpm1=quant_tpm[,c(2:24)]%>%group_by(geneID)%>%dplyr::summarise(across(everything(), sum))
quant_count1=quant_count[,c(2:24)]%>%group_by(geneID)%>%dplyr::summarise(across(everything(), sum))
write.table(quant_tpm1,gzfile("ExoKD.gene.tpm.matrix.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)
write.table(quant_count1,gzfile("ExoKD.gene.count.matrix.tsv.gz"), row.names=F, col.names=T, sep="\t", quote=F)

quant_count.tran=read.delim("ExoKD.transcript.count.matrix.tsv.gz", row.names=1, header=T, stringsAsFactors = F, check.names = F)
quant_count.gene=read.delim("ExoKD.gene.count.matrix.tsv.gz", row.names=1, header=T, stringsAsFactors = F, check.names = F)
quant_count.tran=quant_count.tran[,c(7,8,11,12)]
quant_count.gene=quant_count.gene[,c(7,8,11,12)]
d <- DGEList(counts=quant_count.tran)
RLE <- calcNormFactors(d, method="RLE")
quant_RLE.tran=data.frame(cpm(RLE, normalized.lib.sizes=TRUE))
d <- DGEList(counts=quant_count.gene)
RLE <- calcNormFactors(d, method="RLE")
quant_RLE.gene=data.frame(cpm(RLE, normalized.lib.sizes=TRUE))

quant_RLE.tran$exo_supressed=(quant_RLE.tran$iPS_EX3si2_rep1+quant_RLE.tran$iPS_EX3si2_rep2)/2
quant_RLE.tran$control=(quant_RLE.tran$iPS_NC1_rep3+quant_RLE.tran$iPS_NC1_rep4)/2
quant_RLE.tran$exo_sensitivity=(quant_RLE.tran$exo_supressed-quant_RLE.tran$control)/quant_RLE.tran$exo_supressed
quant_RLE.tran$exo_sensitivity[which(rowSums(quant_RLE.tran[,c(1:4)])==0)]=NA
quant_RLE.tran$exo_sensitivity[which(quant_RLE.tran$exo_sensitivity<=0)]=0
length(which(quant_RLE.tran$exo_sensitivity >0.5 )) #68009
length(which(!is.na(quant_RLE.tran$exo_sensitivity))) #178444

quant_RLE.tran$group="GENCODE"
quant_RLE.tran$group[grep("CHR",rownames(quant_RLE.tran))]="CHR_iPSC"
quant_RLE.tran$group[grep("ONT",rownames(quant_RLE.tran))]="ONTG"
quant_RLE.tran$detection="No"
quant_RLE.tran$detection[which(rowSums(quant_RLE.tran[,c(1:4)])>=0.1)]="Yes"
write.table(quant_RLE.tran,gzfile("ExoKD.transcript.RLE.matrix.tsv.gz"), col.names=T, row.names=T, sep="\t", quote=F)

quant_RLE.gene$exo_supressed=(quant_RLE.gene$iPS_EX3si2_rep1+quant_RLE.gene$iPS_EX3si2_rep2)/2
quant_RLE.gene$control=(quant_RLE.gene$iPS_NC1_rep3+quant_RLE.gene$iPS_NC1_rep4)/2
quant_RLE.gene$exo_sensitivity=(quant_RLE.gene$exo_supressed-quant_RLE.gene$control)/quant_RLE.gene$exo_supressed
quant_RLE.gene$exo_sensitivity[which(rowSums(quant_RLE.gene[,c(1:4)])==0)]=NA
quant_RLE.gene$exo_sensitivity[which(quant_RLE.gene$exo_sensitivity<=0)]=0
length(which(quant_RLE.gene$exo_sensitivity >0.5)) #24753
length(which(!is.na(quant_RLE.gene$exo_sensitivity))) #61037

quant_RLE.gene$group="GENCODE"
quant_RLE.gene$group[grep("CHR",rownames(quant_RLE.gene))]="CHR_iPSC"
quant_RLE.gene$group[grep("ONT",rownames(quant_RLE.gene))]="ONTG"
quant_RLE.gene$detection="No"
quant_RLE.gene$detection[which(rowSums(quant_RLE.gene[,c(1:4)])>=0.1)]="Yes"
write.table(quant_RLE.gene,gzfile("ExoKD.gene.RLE.matrix.tsv.gz"), col.names=T, row.names=T, sep="\t", quote=F)

#exosome sensitivity at gene level and transcript level using RLE normalized matrix are stored in [primary_folder]/code_n_data/Fig4_transcription_features/exosome_sensitivity
# for -> Fig4 k-m Ext6o


#===============================================================================
#add back DEG total to chr - gene base
# gene based was not used in the manuscript
exp_gene=read.delim(paste0(chromatin_path,"edgeR_bambu_genebase_total_chr_ipsc.tsv"), header=T, stringsAsFactors = F, check.names = F)
quant_RLE.gene1=quant_RLE.gene[which(quant_RLE.gene$detection=="Yes"),]
quant_RLE.gene1$geneID=rownames(quant_RLE.gene1)

exp.df1=inner_join(exp_gene, quant_RLE.gene1, by="geneID", copy=F)
exp.df1=left_join(exp.df1,unique(chr_table5[,c("T4_gene_ID","gene_group")]), by=c("geneID"="T4_gene_ID"), copy=F)
exp.df1$bin=round(exp.df1$logFC)
kk=exp.df1[grep("e_",exp.df1$gene_group),]%>%group_by(bin)%>%dplyr::summarise(count=n(), med=median(exo_sensitivity))
exp.df1%>%group_by(gene_group)%>%dplyr::summarise(count=n())

exp.df1$bin[which(exp.df1$bin>5)]=5
exp.df1$bin[which(exp.df1$bin<(-3))]=(-3)
exp.df1$group2="Others"
exp.df1$group2[grep("e_",exp.df1$gene_group)]="e_ncRNA"
exp.df1$group2[grep("p_",exp.df1$gene_group)]="p_ncRNA"
exp.df8=exp.df1[which(exp.df1$group2!="Others" & exp.df1$logCPM>(-2)),]%>%group_by(group2, bin)%>%dplyr::summarise(count=n(), median=median(exo_sensitivity))

ggplot()+
  geom_boxplot(data=exp.df1[which(exp.df1$group2!="Others" & exp.df1$logCPM > (-2)),], mapping=aes(x=as.factor(bin), y=exo_sensitivity), size=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA)+
  facet_grid(cols=vars(group2))+
  scale_x_discrete(labels=c("-3 or less","-2","-1","0","1","2","3","4","5 or greater"))+
  labs(y="Exosome sensitivity", x="Log2FC (chromatin-bound/total)", title="Exosome sensitivty of iPSC ncRNAs")+
  geom_smooth(data = exp.df8, aes(x = as.numeric(as.factor(bin)), y = median), method = "loess", color = "blue", size=0.25) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1, lineheight = 0.7))
#gene-base not used in the manuscript

#===============================================================================
#add back DEG total to chr - transcript base
exp_tx=read.delim(paste0(chromatin_path,"edgeR_bambu_Txbase_total_chr_ipsc.tsv"), header=T, stringsAsFactors = F, check.names = F)
quant_RLE.tran1=quant_RLE.tran[which(quant_RLE.tran$detection=="Yes"),]
quant_RLE.tran1$group=factor(quant_RLE.tran1$group, levels=c("GENCODE","ONTG","CHR_iPSC"))
quant_RLE.tran1$exo_sensitivity_bin=round(as.numeric(quant_RLE.tran1$exo_sensitivity)*100)/100
quant_RLE.tran1$transcriptID=rownames(quant_RLE.tran1)

exp.df2=inner_join(exp_tx, quant_RLE.tran1, by="transcriptID", copy=F)
exp.df2=left_join(exp.df2,unique(chr_table5[,c("model_ID","transcript_group")]), by=c("transcriptID"="model_ID"), copy=F)
exp.df2=left_join(exp.df2,unique(table5[,c("model_ID","group")]), by=c("transcriptID"="model_ID"), copy=F)
exp.df2$transcript_group[which(is.na(exp.df2$transcript_group))]=exp.df2$group.y[which(is.na(exp.df2$transcript_group))]
exp.df2$bin=round(exp.df2$logFC)
kk=exp.df2[grep("e_",exp.df2$transcript_group),]%>%group_by(bin)%>%dplyr::summarise(count=n(), med=median(exo_sensitivity))
exp.df2%>%group_by(transcript_group)%>%dplyr::summarise(count=n())

exp.df2$bin[which(exp.df2$bin>6)]=6
exp.df2$bin[which(exp.df2$bin<(-3))]=(-3)
exp.df2$group2="Others"
exp.df2$group2[grep("e_",exp.df2$transcript_group)]="e_ncRNA"
exp.df2$group2[grep("p_",exp.df2$transcript_group)]="p_ncRNA"
write.table(exp.df2,gzfile(paste0(path_fig4_data,"Transcript_chr_total_DET_exosome_iPSC.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig4/data
# for -> Fig. Ext6m


#fisher's exact for RBP binding
#===============================================================================
#RBP from transcript exon postar
#please download the RBP data from POSTAR3 (http://111.198.139.65/RBP.html)
#==================
#bash
#bedtools intersect between POSTAR and SALA final transcriptome exon
setwd(RBP_path)
system(paste0("bedtools bed12tobed6 -i ",primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.bed12.bed.gz | gzip > table5.final.partial_yes_detected.exon.bed6.bed.gz"))
system("bedtools intersect -wa -wb -a table5.final.partial_yes_detected.exon.bed6.bed.gz -b /analysisdata/fantom6/Interactome/resources/POSTAR3.human.bed.gz | gzip > table5.final.partial_yes_detected.exon.POSTAR.bed.gz")
#==================
setwd(RBP_path)
RBP=read.delim("table5.final.partial_yes_detected.exon.POSTAR.bed.gz", header=F, stringsAsFactors = F)
RBP=unique(RBP[,c(4,10)])

table5=read.delim(paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"),header=T, stringsAsFactors = F, check.names=F)
table5g=unique(table5[,c("T4_gene_ID","promoter_type","CPAT_class","gene_group")])
RBP=left_join(RBP, table5[,c("model_ID","T4_gene_ID","promoter_type","CPAT_class","gene_group")], by=c("V4"="model_ID"),copy=F)
RBP=RBP[which(RBP$gene_group %in% c("p_ncRNA","e_ncRNA","mRNA")),]
RBP=unique(RBP[,c(2:6)])

RBP_cnc=RBP%>%group_by(V10)%>%dplyr::summarise(coding_bind=sum(gene_group=="mRNA"),non_coding_bind=sum(gene_group!="mRNA"))
RBP_cnc$coding_nobind=length(which(table5g$gene_group =="mRNA"))-RBP_cnc$coding_bind
RBP_cnc$non_coding_nobind=length(which(table5g$gene_group %in% c("p_ncRNA","e_ncRNA")))-RBP_cnc$non_coding_bind
write.table(RBP_cnc,"coding_non_coding_genebase.rbp_binding.txt", col.names=T, row.names=F, sep="\t", quote=F)

colnames(RBP_cnc)[c(2:5)]=c("Yes_Yes","Yes_No","No_Yes","No_No")

RBP_enp=RBP%>%group_by(V10)%>%dplyr::summarise(pro_bind=sum(promoter_type=="promoter-like"),enh_bind=sum(promoter_type=="enhancer-like"))
RBP_enp$pro_nobind=length(which(table5$CPAT_class =="coding"))-RBP_enp$pro_bind
RBP_enp$enh_nobind=length(which(table5$CPAT_class =="non-coding"))-RBP_enp$enh_bind


RBP_cnc=read.delim("coding_non_coding_trnscrpts.rbp_binding.txt", header=T, stringsAsFactors = F)
RBP_cnc$Yes_Yes=RBP_cnc$exon_prom
RBP_cnc$Yes_No=RBP_cnc$exon_prom/RBP_cnc$fraction_prom-RBP_cnc$exon_prom
RBP_cnc$No_Yes=RBP_cnc$exon_enh
RBP_cnc$No_No=RBP_cnc$exon_enh/RBP_cnc$fraction_enh-RBP_cnc$exon_enh
RBP_cnc$total=rowSums(RBP_cnc[,c(8:11)])
for(i in 1:nrow(RBP_cnc)){
  GSEATasting <- matrix(c(RBP_cnc$No_No[i], RBP_cnc$No_Yes[i], RBP_cnc$Yes_No[i], RBP_cnc$Yes_Yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  RBP_cnc$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  RBP_cnc$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
RBP_cnc$FE_logOR=log10(RBP_cnc$OR)
RBP_cnc$sig_level="ns"
RBP_cnc$sig_level[which(RBP_cnc$p.val<0.05)]="*"
RBP_cnc$sig_level[which(RBP_cnc$p.val<0.01)]="**"
RBP_cnc$sig_level[which(RBP_cnc$p.val<0.001)]="***"
RBP_cnc$sig_level=factor(RBP_cnc$sig_level, levels=c("ns","*","**","***"))
RBP_cnc$group="coding_poential"
RBP_cnc$FDR=p.adjust(RBP_cnc$p.val, method="BH")

RBP_enp=read.delim("prom_enh_like_trnscrpts.rbp_binding_greater_in_enh.txt", header=T, stringsAsFactors = F)
RBP_enp$Yes_Yes=RBP_enp$exon_prom
RBP_enp$Yes_No=RBP_enp$exon_prom/RBP_enp$fraction_prom-RBP_enp$exon_prom
RBP_enp$No_Yes=RBP_enp$exon_enh
RBP_enp$No_No=RBP_enp$exon_enh/RBP_enp$fraction_enh-RBP_enp$exon_enh
RBP_enp$total=rowSums(RBP_enp[,c(8:11)])
for(i in 1:nrow(RBP_enp)){
  GSEATasting <- matrix(c(RBP_enp$No_No[i], RBP_enp$No_Yes[i], RBP_enp$Yes_No[i], RBP_enp$Yes_Yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  RBP_enp$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  RBP_enp$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
RBP_enp$FE_logOR=log10(RBP_enp$OR)
RBP_enp$sig_level="ns"
RBP_enp$sig_level[which(RBP_enp$p.val<0.05)]="*"
RBP_enp$sig_level[which(RBP_enp$p.val<0.01)]="**"
RBP_enp$sig_level[which(RBP_enp$p.val<0.001)]="***"
RBP_enp$sig_level=factor(RBP_enp$sig_level, levels=c("ns","*","**","***"))
RBP_enp$group="promoter_type"
RBP_enp$FDR=p.adjust(RBP_enp$p.val, method="BH")

RBP_trans_exon=rbind(RBP_cnc,RBP_enp)
write.table(RBP_trans_exon,gzfile(paste0(path_fig4_data,"RBP.nc.enhancer.FE.result.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
#stored in [primary_folder]/fig4/data
# for -> Fig. Ext6n

#===============================================================================
# transcript length of chromatin bound iPSC transcript
table5=read.delim(paste0(SALA_foler,"log/table5.chimeric.76K.remove.permissive.isoform.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_infoa=fread(paste0(SALA_foler,"log/full_length_support_readID_modelID_pair.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
read_infoa1=read_infoa[which(read_infoa$model_ID_str %in% table5$model_ID),]
path1=paste0(SALA_foler,"tmp/")
files=list.files(path=path1, pattern="trnscpt.qry.bed.bgz", recursive=T)
files2=list.files(path=path1, pattern="trnscpt.qry.exon.bed6.bed.gz", recursive=T)
files.names=sapply(strsplit(files,".trnscpt.qry.bed.bgz"),"[",1)

for ( i in 1: length(files)){
  bed12=read.delim(paste0(path1,files[i]), header=F, stringsAsFactors = F, check.names = F)
  bed12$V4=gsub("WTC-11-NGN2-hiPSC_chromatin_1_rep","Ic",bed12$V4)
  bed12$V4=sapply(strsplit(bed12$V4,"\\|"),"[",1)
  bed12=bed12[which(bed12$V4 %in% read_infoa1$trnscpt_ID),c(4,1,2,3,6)]
  bed6=read.delim(paste0(path1,files2[i]), header=F, stringsAsFactors = F, check.names = F)
  bed6$V4=gsub("WTC-11-NGN2-hiPSC_chromatin_1_rep","Ic",bed6$V4)
  bed6$V4=sapply(strsplit(bed6$V4,"\\|"),"[",1)
  bed6$length=bed6$V3-bed6$V2
  bed6a=bed6[which(bed6$V4 %in% read_infoa1$trnscpt_ID),]%>%group_by(V4)%>%dplyr::summarise(n_exon=n(),length=sum(length))
  info=left_join(bed12,bed6a,by="V4",copy=F)
  info=left_join(info,read_infoa1,by=c("V4"="trnscpt_ID"),copy=F)
  write.table(info, gzfile(paste0(path1,files.names[i],".read.info.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)}
rm(bed6)
rm(bed12)
gc()

files=list.files(path=path1, pattern="read.info.tsv.gz", recursive=T)
files=files[-23] #remove chrM
i=1
data=fread(paste0(path1,files[i]), stringsAsFactors = F, select=c(1:7,9))
data=left_join(data,table5[,c(1,34,36,53)],by=c("model_ID_str"="model_ID"),copy=F)
data$genomic_range=data$V3-data$V2
data1=data%>%group_by(model_ID_str)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
data1a=data%>%group_by(CREID,V6)%>%dplyr::summarise(n5_var=var(V2), n3_var=var(V3),read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
data1b=data%>%group_by(T4_gene_ID)%>%dplyr::summarise(n5_var=var(V3), n3_var=var(V2),read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
data1c=data%>%group_by(n5_string)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())

for (i in 2:length(files)){
  data=fread(paste0(path1,files[i]), stringsAsFactors = F, select=c(1:7,9))
  data=left_join(data,table5[,c(1,34,36,53)],by=c("model_ID_str"="model_ID"),copy=F)
  data$genomic_range=data$V3-data$V2
  data2=data%>%group_by(model_ID_str)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
  data2a=data%>%group_by(CREID,V6)%>%dplyr::summarise(n5_var=var(V2), n3_var=var(V3),read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
  data2b=data%>%group_by(T4_gene_ID)%>%dplyr::summarise(n5_var=var(V3), n3_var=var(V2),read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
  data2c=data%>%group_by(n5_string)%>%dplyr::summarise(read_median_length=median(length), MAD=mad(length), read_mean_length=mean(length), SD=sd(length), median_exon=median(n_exon), median_range=median(genomic_range), n_read=n())
  data1=rbind(data1,data2)
  data1a=rbind(data1a,data2a)
  data1b=rbind(data1b,data2b)
  data1c=rbind(data1c,data2c)}

data1a=left_join(data1a, table5cre, by=c("CREID"="CREID", "V6"="strand"),copy=F)
k1=data1a%>%group_by(gene_group)%>%dplyr::summarise(median=median(read_median_length),exon=median(median_exon),count=n())
k2=data1a%>%group_by(promoter_type)%>%dplyr::summarise(median=median(read_median_length),exon=median(median_exon),count=n())
data1b=left_join(data1b, table5g, by="T4_gene_ID", copy=F)
k3=data1b%>%group_by(gene_group)%>%dplyr::summarise(median=median(read_median_length),exon=median(median_exon),count=n())

write.table(data1,gzfile(paste0(SALA_foler,"length/table8.length.exon.info.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
write.table(data1a,gzfile(paste0(SALA_foler,"length/table8.length.exon.info.CREversion.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
write.table(data1b,gzfile(paste0(SALA_foler,"length/table8.length.exon.info.GENEversion.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)
write.table(data1c,gzfile(paste0(SALA_foler,"length/table8.length.exon.info.ex5_cluster.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

data1=read.delim(paste0(SALA_foler,"length/table8.length.exon.info.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)
data1=left_join(data1,table5[,c(1,113)],by=c("model_ID_str"="model_ID"), copy=F)
data1=data1[which(!is.na(data1$length_diff)),]
write.table(data1, gzfile(paste0(path_fig4_data,"chr_bound_iPSC_length_diff.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#stored in [primary_folder]/fig4/data
# for -> Fig. Ext6l



