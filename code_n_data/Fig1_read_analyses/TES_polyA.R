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
path_fig1_data=paste0(primary_folder,"fig1/data/")
read_path=paste0(primary_folder,"code_n_data/Fig1_read_analyses/")
TES_path=paste0(read_path,"TES/")
GENCODE_path=paste0(primary_folder,"code_n_data/GENCODEv39/")
SALA_path=paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/")

#===============================================================================
#random forest polyA prediction result
ont_pas=read.delim(paste0(TES_path,"fantom_random_forest20240407/output/annotated.3pclusters.gtf.gz"), skip=3, header=F, stringsAsFactors = F)
gc_pas=read.delim(paste0(TES_path,"fantom_random_forest20240407/output/annotated.3pclusters.gencode.gtf.gz"), skip=3, header=F, stringsAsFactors = F)
all_pas=rbind(ont_pas,gc_pas)
#this include all TES (from raw Neuron-THP1 dataset, & raw chromatin-bound iPSC dataset, & all GENCODE v39)

all_pas[c('cluster_id', 'polyAsignalClass', 'polyA.Signal', 'probability_TRUE', 'inRef', "class")] <- str_split_fixed(all_pas$V9, ';', 6)
all_pas$cluster_id=gsub("cluster_id ","",all_pas$cluster_id)
all_pas$polyAsignalClass=gsub(" polyAsignalClass ","",all_pas$polyAsignalClass)
all_pas$polyA.Signal=gsub(" polyA.Signal ","",all_pas$polyA.Signal)
all_pas$probability_TRUE=gsub(" probability_TRUE ","",all_pas$probability_TRUE)
all_pas$inRef=gsub(" inRef ","",all_pas$inRef)
all_pas$class=gsub(" class ","",all_pas$class)
all_pas=all_pas[,c(1,4,5,7,13:15)]

all_pas$label=paste0(all_pas$V1,"_",all_pas$V5-1,"_",all_pas$V5,"_",all_pas$V7)

#===============================================================================
#TES position bed, at 25 nt upstream the end site 
options(scipen=999)
all3n=all_pas[,c(1,2,3,12,8,4)]
all3n$V4=all3n$V4-1
all3n$V4[which(all3n$V7 == "+")]=all3n$V4[which(all3n$V7 == "+")]-25
all3n$V5[which(all3n$V7 == "+")]=all3n$V5[which(all3n$V7 == "+")]-25
all3n$V4[which(all3n$V7 == "-")]=all3n$V4[which(all3n$V7 == "-")]+25
all3n$V5[which(all3n$V7 == "-")]=all3n$V5[which(all3n$V7 == "-")]+25
write.table(all3n[order(all3n$V1,all3n$V4),],gzfile("all3n_up25.bed.gz"), col.names=F, row.names=F, sep="\t", quote=F)

#===============================================================================
#run genome wide PAS identification [scanMotifGenomeWide.pl PAS.motif hg38 -bed -keepAll -p 10 > output.bed]
setwd(TES_path)
system("perl 20240317_run_HOMER_PAS_motifs.pl")

#========================
#process the output
#homer bed output is 1-based!
gPAS=fread(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.bed.gz"), header=F)
gPAS=gPAS[which(nchar(gPAS$V1)<=5),] #main chromaosome only
gPAS$V2=gPAS$V2-1
gPAS$V4=paste0("PAS",1:nrow(gPAS))
write.table(gPAS, gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.bed.gz")), row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#getfasta
setwd(paste0(TES_path,"version20240317/"))
system("bedtools getfasta -bedOut -s -fi /analysisdata/fantom6/Interactome/resources/minimap2_mapping/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -bed GENCODE_polyA_signal.motif_rng.main.bed.gz | gzip > GENCODE_polyA_signal.motif_rng.main.FASTA.bed.gz")
# fasta file removed afterwards
#========================
#obtain the motif seq
gPAS=fread(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.bed.gz"), header=F)
gPAS1=gPAS[which(nchar(gPAS$V7)>6),]
gPAS=gPAS[which(nchar(gPAS$V7)==6),]

unique(gPAS$V7[which(gPAS$V5>6)]) # "AATAAA" "ATTAAA" 
unique(gPAS$V7[which(gPAS$V5>4 & gPAS$V5<6)]) # "AGTAAA" "TATAAA" "AATATA" 
unique(gPAS$V7[which(gPAS$V5>3 & gPAS$V5<4)]) # "ATTATA" "TTTAAA" "AAAAAA" "AATTAA" "AATAAG" "AATACA" "AAGAAA" "ACTAAA" "AATGAA" "CATAAA" "AATAAT" "GATAAA" "AATAGA" 
options(scipen=999)
write.table(gPAS[which(gPAS$V5>3), c(1:3,7,5,6)], gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.zenbu.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)
gPAS$V2[which(gPAS$V=="+")]=gPAS$V3[which(gPAS$V=="+")]-1
gPAS$V3[which(gPAS$V=="-")]=gPAS$V2[which(gPAS$V=="-")]+1
write.table(gPAS[which(gPAS$V5>6), c(1:3,7,5,6)], gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.6.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)
write.table(gPAS[which(gPAS$V5>4 & gPAS$V5<6), c(1:3,7,5,6)], gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.4_6.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)
write.table(gPAS[which(gPAS$V5>3 & gPAS$V5<4), c(1:3,7,5,6)], gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.3_4.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)
write.table(gPAS[which(gPAS$V5>1 & gPAS$V5<3), c(1:3,7,5,6)], gzfile(paste0(TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.1_3.bed.gz")),row.names=F, col.names=F, sep="\t", quote=F)

#===============================================================================
#remove intermediate file
setwd(paste0(TES_path,"version20240317/"))
system("rm GENCODE_polyA_signal.motif_rng.main.FASTA.bed.gz")
system("rm GENCODE_polyA_signal.motif_rng.bed.gz")
system("rm GENCODE_polyA_signal.motif_rng.main.bed.gz")

#bedtools closest link our TES to genome wide PAS with different score
setwd(TES_path)
system(paste0("bedtools closest -a all3n_up25.sort.bed.gz -b ",TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.6.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS6.bed.gz"))
system(paste0("bedtools closest -a all3n_up25.sort.bed.gz -b ",TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.4_6.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS46.bed.gz"))
system(paste0("bedtools closest -a all3n_up25.sort.bed.gz -b ",TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.3_4.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS34.bed.gz"))
system(paste0("bedtools closest -a all3n_up25.sort.bed.gz -b ",TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.FASTA.1_3.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS13.bed.gz"))
system(paste0("bedtools closest -a all3n_up25.sort.bed.gz -b ",TES_path,"version20240317/GENCODE_polyA_signal.motif_rng.main.zenbu.bed.gz -s -D a | gzip > all3n_up25.PASstandard.bed.gz"))

#remove intermediate file
setwd(paste0(TES_path,"version20240317/"))
system("rm GENCODE_polyA_signal.motif_rng.main.*.bed.gz")

#===============================================================================
#parse the files to include higher motif score first, followed by lower score
setwd(TES_path)
files=list.files(path=TES_path, pattern="all3n_up25.PAS")
files=files[-c(1,6)]
files
"all3n_up25.PAS13.bed.gz" 
"all3n_up25.PAS34.bed.gz" 
"all3n_up25.PAS46.bed.gz" 
"all3n_up25.PAS6.bed.gz" 

all3n=fread(files[4], header=F, stringsAsFactors = F)
all3n=all3n[which(all3n$V10>= (-20) & all3n$V10 <= 20),]
all3n=all3n[which(all3n$V7 != "."),]
all3n$V10=all3n$V10-25
colnames(all3n)[c(7,8,10)]=c("PAS_motif","motif_score","PAS_distance")
all3n=all3n[which(all3n$PAS_distance>= (-35)),]
all3n=all3n%>%group_by(V4)%>%dplyr::slice(which.max(PAS_distance))
all_pas=left_join(all_pas, all3n[,c(4,7,8,10)], by=c("label"="V4"),copy=F)
all_pas0=all_pas[which(!is.na(all_pas$PAS_motif)),]
all_pas1=all_pas[which(is.na(all_pas$PAS_motif)),c(1:12)]

for (i in c(3,2,1)){
  all3n=fread(files[i], header=F, stringsAsFactors = F)
  all3n=all3n[which(all3n$V10>= (-20) & all3n$V10 <= 20),]
  all3n=all3n[which(all3n$V7 != "."),]
  all3n$V10=all3n$V10-25
  colnames(all3n)[c(7,8,10)]=c("PAS_motif","motif_score","PAS_distance")
  all3n=all3n[which(all3n$PAS_distance>= (-35)),]
  all3n=all3n%>%group_by(V4)%>%dplyr::slice(which.max(PAS_distance))
  all_pas1=left_join(all_pas1, all3n[,c(4,7,8,10)], by=c("label"="V4"),copy=F)
  all_pas0=rbind(all_pas0,all_pas1)
  all_pas1=all_pas0[which(is.na(all_pas0$PAS_motif)),c(1:12)]
  all_pas0=all_pas0[which(!is.na(all_pas0$PAS_motif)),]}

unique(all_pas0$motif_score)
length(which(all_pas0$PAS == "No" & all_pas0$motif_score > 4))
length(which(all_pas1$PAS == "Yes"))
all_pas1=left_join(all_pas1, all3n[,c(4,7,8,10)], by=c("label"="V4"),copy=F)
length(which(is.na(all_pas1$PAS_motif)))
all_pas=rbind(all_pas0,all_pas1)

gencode.polyA=read.delim(paste0(GENCODE_path,"gencode.v39.polyAs.site.bed.gz"), header=F, stringsAsFactors = F)
gencode.polyA$TES=paste0(gencode.polyA$V1,"_",gencode.polyA$V2,"_",gencode.polyA$V3,"_",gencode.polyA$V6)
all_pas$GENCODE_polyA="No"
all_pas$GENCODE_polyA[which(all_pas$label %in% gencode.polyA$TES)]="Yes"

write.table(all_pas, gzfile(paste0(TES_path,"all_read_n3_ont_gencode_prediction20240407.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#prepare for plotting roc curve
#THP-1 read base here
all.thp1=fread(paste0(TES_path,"all_thp1.tsv.gz"), header=T, stringsAsFactors =F, check.names=F)
all.thp1=all.thp1[,c(1:8,11,12)]
all.thp1$cell="THP-1"
all.thp1$cell[grep("PMA",all.thp1$V4)]="dTHP-1"
all.thp1$cell[which(all.thp1$PAtailing == "woPAP")]=paste0(all.thp1$cell[which(all.thp1$PAtailing == "woPAP")],".woPAP")
all.thp2=all.thp1%>%group_by(cell,label,PAtailing,internal_prime2)%>%dplyr::summarise(TES_count=n())

all.thp2=left_join(all.thp2,all_pas[,c(8:20)], by="label", copy=F)
write.table(all.thp2, gzfile(paste0(TES_path,"all_thp1.PAS.collapsed.tsv.gz")), col.names=T, row.names =F, quote=F, sep="\t")

#===============================================================================
#define presence of polyA from THP-1 data

all.thp2=all.thp1[-grep("D9P2",all.thp1$readname),]
all.thp1_ROCprep=all.thp2%>%group_by(label,PAtailing)%>%dplyr::summarise(count=n())
all.thp1_ROCprep=spread(all.thp1_ROCprep, key=2, value=3)
all.thp1_ROCprep[is.na(all.thp1_ROCprep)]=0
colSums(all.thp1_ROCprep[,c(2,3)])
all.thp1_ROCprep$PAP_TPM=all.thp1_ROCprep$PAP/sum(all.thp1_ROCprep$PAP)*1000000
all.thp1_ROCprep$woPAP_TPM=all.thp1_ROCprep$woPAP/sum(all.thp1_ROCprep$woPAP)*1000000
all.thp1_ROCprep1=all.thp1_ROCprep[which(all.thp1_ROCprep$PAP>=3 | all.thp1_ROCprep$woPAP >=3),]
all.thp1_ROCprep1$log2FC=log2(all.thp1_ROCprep1$PAP_TPM/all.thp1_ROCprep1$woPAP_TPM)

length(which(all.thp1_ROCprep1$log2FC == Inf & all.thp1_ROCprep1$PAP>=3)) #22540
length(which(all.thp1_ROCprep1$log2FC <= 0 & all.thp1_ROCprep1$woPAP>=6)) #74918

all.thp1_ROCprep2=rbind(all.thp1_ROCprep1[which(all.thp1_ROCprep1$log2FC == Inf & all.thp1_ROCprep1$PAP>=3),],
                        all.thp1_ROCprep1[which(all.thp1_ROCprep1$log2FC <= 0 & all.thp1_ROCprep1$woPAP>=6),])

all.thp1_ROCprep2=left_join(all.thp1_ROCprep2, unique(all.thp2[,c(8,12,22:25)]), by ="label", copy=F)
all.thp1_ROCprep2$exp_polyA="Yes"
all.thp1_ROCprep2$exp_polyA[which(all.thp1_ROCprep2$log2FC == Inf)]="No"
write.table(all.thp1_ROCprep2, gzfile(paste0(path_fig1_data,"prediction_THP1_roc_curve_prep.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#==================================
library(pROC)
#roc curve using FLAM-seq data as ground true
all_pas1=all_pas[which(all_pas$GENCODE_polyA == "No"),]
roc_data1 <- roc(all_pas1$inRef, all_pas1$probability_TRUE, levels=c("TRUE","FALSE"))
auc1=auc(roc_data1) #0.9412
coords(roc_data1, 0.175, best.method="youden") #95.03%
coords(roc_data1, 0.655, best.method="youden") #95.06%
coords(roc_data1, "best", best.method="youden") #0.3191026

#roc curve using THP1 long read result as ground true
all.thp1_ROC=read.delim(paste0(path_fig1_data,"prediction_THP1_roc_curve_prep.tsv.gz"), header=T, check.names=F, stringsAsFactors = F)
roc_data_2 <- roc(all.thp1_ROC$exp_polyA[which(all.thp1_ROC$internal_prime2 =="no")], as.numeric(all.thp1_ROC$probability_TRUE[which(all.thp1_ROC$internal_prime2 =="no")]), levels=c("Yes","No"))
auc2=auc(roc_data_2) #0.932
coords(roc_data_2, 0.085, best.method="youden") ##-> ~95%
coords(roc_data_2, 0.6, best.method="youden") ##-> ~95%
coords(roc_data_2, "best", best.method="youden") ##-> 0.2357589

#prepare roc curve plot table
roc_df1 <- data.frame(
  specificity = roc_data1$specificities,
  sensitivity = roc_data1$sensitivities,
  threshold = roc_data1$thresholds,
  group="FLAM-seq")
roc_df1a=roc_df1[sample(nrow(roc_df1),80000),]

roc_df2 <- data.frame(
  specificity = roc_data_2$specificities,
  sensitivity = roc_data_2$sensitivities,
  threshold = roc_data_2$thresholds,
  group="THP-1")

roc_df1=rbind(roc_df1a,roc_df2)
saveRDS(roc_df1,paste0(path_fig1_data,"roc_curve.RDS"))
#===============================================================================
#define classifier and motif result

all_pas$polyA_prediction3="No"
all_pas$polyA_prediction3[which(all_pas$class == "significant poly(A)")] ="Yes"

all_pas$GENCODE_polyA[which(is.na(all_pas$GENCODE_polyA))]="No"
all_pas$PAS2="No"
all_pas$PAS2[which(all_pas$motif_score>3)]="Yes"

write.table(all_pas, gzfile(paste0(TES_path,"all_read_n3_ont_gencode_prediction20240407.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#transfer the data to read base all TES
all_pas=read.delim(paste0(TES_path,"all_read_n3_ont_gencode_prediction20240407.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data3=read.delim(paste0(SALA_path,"transcript/log/table0.TESbymodelID.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data3=data3[,c(1:7)]
data3=left_join(data3,all_pas[,c("probability_TRUE", "inRef","class", "label", "PAS_motif",  "motif_score", "PAS_distance", "GENCODE_polyA", "polyA_prediction3" ,"PAS2")],by=c("TES"="label"),copy=F)
write.table(data3, gzfile(paste0(SALA_path,"transcript/log/table0.TESbymodelID.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

#===
#annotate TES with transcript model annotation, from raw output and filtered transcriptome
transcript_model=read.delim(paste0(SALA_path,"transcript/bed/Neuron_THP1.S3.model.bed.bgz"), header=F, stringsAsFactors = F, check.names = F)
transcript_model$start3=transcript_model$V2
transcript_model$end3=transcript_model$V2+1
transcript_model$start3[which(transcript_model$V6=="+")]=transcript_model$V3[which(transcript_model$V6=="+")]-1
transcript_model$end3[which(transcript_model$V6=="+")]=transcript_model$V3[which(transcript_model$V6=="+")]
transcript_model$TES=paste0(transcript_model$V1,"_",transcript_model$start3,"_",transcript_model$end3,"_",transcript_model$V6)
transcript_model=left_join(transcript_model[,c(4,15)],unique(all_pas[,c("probability_TRUE", "inRef","class", "label", "PAS_motif",  "motif_score", "PAS_distance", "GENCODE_polyA", "polyA_prediction3" ,"PAS2")]),by=c("TES"="label"),copy=F)
write.table(transcript_model, gzfile(paste0(SALA_path,"transcript/log/table0.TESbymodelID_modelonly.tsv.gz")), row.names=F, col.names=T, sep="\t", quote=F)

transcript_model%>%group_by(polyA_prediction3, PAS2)%>%dplyr::summarise(count=n())
transcript_model.t5=transcript_model[which(transcript_model$V4 %in% table5$model_ID),]
transcript_model.t5%>%group_by(polyA_prediction3, PAS2)%>%dplyr::summarise(count=n())
transcript_model.t5$polyA="No"
transcript_model.t5$polyA[which(transcript_model.t5$polyA_prediction3 == "Yes" | transcript_model.t5$PAS2 =="Yes")]="Yes"

table5=left_join(table5, transcript_model.t5[,c(1,5,11:12)], by=c("model_ID"="V4"),copy=F, suffix=c("","_20240407"))
write.table(table5, gzfile(paste0(SALA_path,"transcript/log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz")),col.names=T, row.names=F, sep="\t", quote=F)

#===============================================================================
#Prepare table S1
all_pas=read.delim(paste0(TES_path,"all_read_n3_ont_gencode_prediction20240407.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all_pas1=all_pas[,c("V1","V4","V7","probability_TRUE","inRef","class","label","PAS_motif","motif_score","PAS_distance","PAS2","GENCODE_polyA")]
colnames(all_pas1)[c(1,2,3,6,7,11)]=c("chr","position","strand","polyA_predict_class","TESID","PAS_motif_class")
write.table(all_pas1,gzfile(paste0(path_fig1_data,"TableS1_polyA.tsv.gz")), col.names=T, row.names=F, sep="\t",quote=F)

table5=left_join(table5,all_pas1[,c("TESID","polyA_predict_class")],by="TESID", copy=F)
table5%>%group_by(TES_recur,class,PAS2)%>%dplyr::summarise(count=n())
table5=table5[,c(1:110,118,112:117)]
colnames(table5)[c(111,112)]=c("polyA_predict_class","PAS_motif_class")
write.table(table5, gzfile("log/table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), col.names=T, row.names=F, sep="\t",quote=F)

#prepare bed file for zenbu
all_pas2=all_pas1[,c("chr","TESID","strand","polyA_predict_class","PAS_motif_class")]
all_pas2$V2=as.numeric(sapply(strsplit(all_pas2$TESID,"_"),"[",2))
all_pas2$V3=as.numeric(sapply(strsplit(all_pas2$TESID,"_"),"[",3))
all_pas2$score=0
all_pas2$score[which(all_pas2$polyA_predict_class == "significant poly(A)" & all_pas2$PAS_motif_class == "No")]=1
all_pas2$score[which(all_pas2$polyA_predict_class != "significant poly(A)" & all_pas2$PAS_motif_class == "Yes")]=2
all_pas2$score[which(all_pas2$polyA_predict_class == "significant poly(A)" & all_pas2$PAS_motif_class == "Yes")]=3
all_pas2=all_pas2[,c("chr","V2","V3","TESID","score","strand")]
transcript_model=read.delim(paste0(SALA_path,"transcript/log/table0.TESbymodelID_modelonly.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all_pas3=all_pas2[which(all_pas2$TESID %in% unique(transcript_model$TES)),]
write.table(all_pas3[order(all_pas3$chr,all_pas3$V2),], gzfile(paste0(path_fig1_data,"zenbu_rawTranscript_polyA.bed.gz")),col.names=F, row.names=F, sep="\t",quote=F)

#for chromatin bound iPSC
iPSChroTES=read.delim(paste0(primary_folder,"code_n_data/SALA/iPSC_chromatin_bound/transcript/log/potential_internal_prime_TES.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all_pas4=all_pas2[which(all_pas2$TESID %in% unique(iPSChroTES$V4)),]
write.table(all_pas4[order(all_pas4$chr,all_pas4$V2),], gzfile(paste0(path_fig1_data,"zenbu_rawTranscript_iPSChro_polyA.bed.gz")),col.names=F, row.names=F, sep="\t",quote=F)







