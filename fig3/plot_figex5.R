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
path_fig3=paste0(primary_folder,"fig3/out/")
path_fig3_data=paste0(primary_folder,"fig3/data/")

#=====================
theme1=theme(  panel.background = element_rect(fill = "white", color = "white", linewidth = 0.25, linetype = "solid"),
               panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid', color = "grey90"), 
               axis.text.x = element_text(color ="black"),
               axis.text.y = element_text(color ="black"),
               plot.title = element_text(hjust = 0.5,  color ="black", margin = margin(0.1,0.1,0.1,0.1, "cm")),
               plot.subtitle = element_text(color="black", hjust=0.5),
               text = element_text(size=6),
               strip.text = element_text(size=6),
               legend.background = element_rect(fill="white", color="black", linewidth=0.25),
               legend.title = element_text(hjust=0.5),
               legend.text = element_text(lineheight = 0.6, margin = margin(l = 1, unit = "pt")),
               legend.key.size = unit(0.2, 'cm'),
               legend.margin = margin(0.02,0.02,0.02,0.02, "cm"),
               axis.line = element_line(linewidth = 0.25, colour = "black"),
               axis.ticks = element_line(linewidth = 0.25,colour = "black"),
               strip.background =element_rect(fill="grey90"),
               strip.text.x = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               strip.text.y = element_text(margin = margin(0.02,0.02,0.02,0.02, "cm")),
               legend.key = element_blank(),
               legend.position = "right",
               legend.box.margin = margin(t=-10, r=0, b=0, l=-10),
               plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

#======================
#Needed
#ex5a
transcript=read.delim(paste0(path_fig3_data,"transcript.upsetter.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)

transcripta=transcript[which(rowSums(transcript[,3:6]==1)==0),]
transcriptb=transcript[which(rowSums(transcript[,3:6]==1)>0),]
transcripta$same_model=colnames(transcripta)[7]
for (i in 1:nrow(transcriptb)){
  number=which(transcriptb[i,c(3:7)] == 1)
  transcriptb$same_model[i]=list(colnames(transcriptb)[number+2])}
transcript=rbind(transcripta,transcriptb)
transcript_n=transcript%>%group_by(same_model)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
transcript_n$label=paste0(signif(transcript_n$percent,2)*100,"%")
transcript_n$label[which(transcript_n$percent<0.005)]=NA
options(scipen=999)
transcript$group=factor(transcript$group, levels=c("others","ncRNA"))
library(ggupset)
library(ComplexUpset)
ex5a1=ggplot(transcript, aes(x = same_model)) +
  geom_bar(fill="black", linewidth=0.25, alpha=1,  width=0.8,  stat="count") + 
  scale_y_continuous(limits=c(0,150000))+
  geom_text(data=transcript_n,aes(y=count+500, label=label), vjust=0, hjust=0, angle=25, size=2 )+
  scale_x_upset()+
  labs(x="", y = "Counts", title="Same transcript model") +
  theme1+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.line.size = 0.25,
                   combmatrix.label.text = element_text(color ="black", size=6),
                   combmatrix.label.extra_spacing = 0.3,
                   combmatrix.label.make_space = FALSE)
pdf(paste0(path_fig3,"ex5a.upseter.same.transcript.pdf"), width = 2.9, height = 1.8)
print(ex5a1)
dev.off()

transcript_n2=reshape2::melt(transcript[,c(1:7)], id=c(1,2))
transcript_n3=transcript_n2%>%group_by(group,variable)%>%dplyr::summarise(count=sum(as.numeric(value)))
transcript_n4=transcript_n2%>%group_by(variable)%>%dplyr::summarise(count=sum(as.numeric(value)))%>%dplyr::mutate(percent=count/146154)
transcript_n4$label=paste0(signif(transcript_n4$percent,2)*100,"%")
transcript_n3$variable=factor(transcript_n3$variable, levels=c("Refseq","LncBook","GENCODEv47","FANTOM_CAT"))

ex5a2=ggplot(transcript_n3[which(transcript_n3$variable!= "CFC_novel"),], aes(x = variable)) +
  geom_bar(aes( y=count), fill="black", linewidth=0.1, alpha=1,  width=0.8,  stat="identity") + 
  scale_fill_npg()+
  coord_flip()+
  scale_y_reverse(limits=c(7500,0))+
  scale_x_discrete(position="top")+
  geom_text(data=transcript_n4[which(transcript_n4$variable!= "CFC_novel"),],aes(y=count-150, label=label), hjust=(1), size=2.2 )+
  labs(x=NULL, y =NULL, title=NULL, fill="Class") +
  theme1+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
pdf(paste0(path_fig3,"ex5a.upseter.same.transcript.side.pdf"), width = 1.9, height = 0.6)
print(ex5a2)
dev.off()

#======================
#Needed
#ex5b
gene=read.delim(paste0(path_fig3_data,"gene.upsetter.tsv.gz"), header=T, stringsAsFactors =  F, check.names=F)

genea=gene[which(rowSums(gene[,4:7]==1)==0),]
geneb=gene[which(rowSums(gene[,4:7]==1)>0),]
genea$same_model=colnames(genea)[8]
for (i in 1:nrow(geneb)){
  number=which(geneb[i,c(4:8)] == 1)
  geneb$same_model[i]=list(colnames(geneb)[number+3])}
gene=rbind(genea,geneb)
gene_n=gene%>%group_by(same_model)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
gene_n$label=paste0(signif(gene_n$percent,2)*100,"%")
gene_n$label[which(gene_n$percent<0.01)]=NA
options(scipen=999)
library(ggupset)
library(ComplexUpset)

ex5b1=ggplot(gene, aes(x = same_model)) +
  geom_bar(fill="black", linewidth=0.5, alpha=0.9,  width=0.8,  stat="count") + 
  labs(x="", y = "Counts", title="Same gene model", fill="Class") +
  scale_fill_npg()+
  scale_y_continuous(limits=c(0,25000))+
  geom_text(data=gene_n,aes(y=count+50, label=label), vjust=0, hjust=0, angle=25, size=2 )+
  scale_x_upset()+
  theme1+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.panel.point.size = 1.5,
                   combmatrix.panel.line.size = 0.25,
                   combmatrix.label.text = element_text(color ="black", size=6),
                   combmatrix.label.extra_spacing = 0.3,
                   combmatrix.label.make_space = FALSE)
pdf(paste0(path_fig3,"ex5b.upseter.same.gene.pdf"), width = 2.9, height = 1.8)
print(ex5b1)
dev.off()

gene_n2=reshape2::melt(gene[,c(1,4:8)], id=c(1))
gene_n3=gene_n2%>%group_by(variable)%>%dplyr::summarise(count=sum(as.numeric(value)))%>%dplyr::mutate(percent=count/39425)
gene_n3$label=paste0(signif(gene_n3$percent,2)*100,"%")
gene_n3$variable=factor(gene_n3$variable, levels=c("Refseq","GENCODEv47","LncBook","FANTOM_CAT"))

ex5b2=ggplot(gene_n3[which(gene_n3$variable!= "CFC_novel"),], aes(x = variable)) +
  geom_bar(aes(y=count), fill="black", linewidth=0.1, alpha=1,  width=0.8,  stat="identity") + 
  scale_fill_npg()+
  coord_flip()+
  scale_y_reverse(limits=c(17000,0))+
  scale_x_discrete(position="top")+
  geom_text(data=gene_n3[which(gene_n3$variable!= "CFC_novel"),],aes(y=count-300, label=label), hjust=1, size=2.2 )+
  labs(x=NULL, y =NULL, title=NULL, fill="Class") +
  theme1+
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
pdf(paste0(path_fig3,"ex5b2.upseter.same.gene.side.pdf"), width = 1.9, height = 0.6)
print(ex5b2)
dev.off()

#======================
#Needed
#ex5c
Tnumber1=read.delim(paste0(path_fig3_data,"transcript_number_5dataset.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

Tnumber1$database=factor(Tnumber1$database, levels=c("GENCODEv47","Refseq","LncBook","FANTOM_CAT","CFC_novel"))
Tnumber1$variable=factor(Tnumber1$variable, levels=c("others","ncRNA"))

ex5c=ggplot(Tnumber1, aes(x = database, y=as.numeric(value))) +
  geom_bar(aes(fill=variable), alpha=1,  width=0.8, color="black", linewidth=0.25, stat="identity") + 
  scale_y_continuous(breaks=c(0, 200000, 400000, 600000), labels=c(0, "200k", "400k", "600k"))+
  scale_fill_manual(values=c("white","grey"))+
  coord_flip()+
  labs(x=NULL, y = "Counts", title="Transcript content from different databases", fill=NULL) +
  theme1+
  theme(legend.position = "bottom")
pdf(paste0(path_fig3,"ex5c.transcript_number_5dataset.pdf"), width = 2.6, height = 1.5)
print(ex5c)
dev.off() 

#======================
#Needed
#ex5d
expression1=read.delim(paste0(path_fig3_data,"bambu_expression_annot_unannot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

ex5d1=ggplot(expression1, aes(x =group, y = gini, fill=group)) +
  geom_boxplot(linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.6, outlier.shape = NA, alpha=1) + 
  scale_y_continuous(limits=c(0, 1.25), breaks=c(0,0.5,1))+
  labs(x=NULL, y = "Gini index", title="Gene\nspecificity") +
  geom_signif(comparisons = list(c("Un-annotated", "Annotated")), y_position =1.1, map_signif_level=TRUE, na.rm = TRUE, test = "wilcox.test", tip_length = 0, size = 0.25, color = "black", textsize = 2.5) +
  theme1+ theme(legend.position = "none", axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))+
  scale_fill_manual(values = c("#8491B4FF","#DC0000FF"))

ex5d2=ggplot(expression1, aes(x =group, y = maxExp, fill=group)) +
  geom_boxplot(linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.6, outlier.shape = NA, alpha=1) + 
  labs(x=NULL, y = "Max. Expression (TPM)", title="Gene\nexpression") +
  coord_cartesian(ylim=c(0, 2))+
  geom_signif(y_position = 1.5, xmin = 1, xmax = 2, annotation = "***", tip_length = 0, size = 0.25, color = "black", textsize = 2.5) +
  theme1+ theme(legend.position = "none", axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))+
  scale_fill_manual(values = c("#8491B4FF","#DC0000FF"))

pdf(paste0(path_fig3,"ex5d.gini_maxExp.pdf"), width = 1.5, height = 1.9)
grid.arrange(ex5d1,ex5d2, ncol=2, widths=c(1,1))
dev.off()

#======================
#Needed
#ex5ef
quant_count1=read.delim(paste0(path_fig3_data,"Kallisto_Neuron.series.gene.count.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
quant_count=read.delim(paste0(path_fig3_data,"Kallisto_Neuron.series.transcript.count.matrix.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

aa1=quant_count1[which(!is.na(quant_count1$T4_gene_promoter_type)),]%>%group_by(group,detection)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
bb1=quant_count1[which(!is.na(quant_count1$T4_gene_promoter_type)),]%>%group_by(group,detection,T4_gene_ncRNA_subclass)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aa1$group2="Gene"

aa2=quant_count[which(!is.na(quant_count$promoter_type)),]%>%group_by(group,detection)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
bb2=quant_count[which(!is.na(quant_count$promoter_type)),]%>%group_by(group,detection,T4_ncRNA_subclass)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
aa2$group2="Transcript"
aa=rbind(aa1, aa2)
aa$label=paste0(signif(aa$percent*100,2),"%")
aa$label[which(aa$detection == "No")]=NA

ex5e=ggplot(aa, aes(x=group, y=percent, fill=detection)) + 
  scale_fill_manual(values=c("white","#E64B35FF"), guide=NULL)+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  labs(x=NULL , y ="% of models" , title ="Model detection by\nshort-read RNA-seq",  fill=NULL)+
  facet_grid(cols=vars(group2),  scales="free_x", space="free_x")+
  geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
  geom_text(data=aa, mapping=aes(x=group, y=percent, group=detection, label=label), size=1.8, color="black", angle=90 ,position = position_stack(vjust = 0.4, reverse=F))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex5e.model_detection_short_read_RNAseq.pdf"), width = 1.1, height = 1.9)
print(ex5e)
dev.off() 

bb1$group2="Gene"
bb2$group2="Transcript"
colnames(bb1)=colnames(bb2)
bb=rbind(bb1, bb2)
bb=bb[which(bb$group == "Novel"),]
bb$T4_ncRNA_subclass[which(bb$T4_ncRNA_subclass == "excluded")]="Potential_coding"
bb$T4_ncRNA_subclass[which(is.na(bb$T4_ncRNA_subclass))]="Potential_coding"
bb$detection=factor(bb$detection, levels=c("Yes","No"))

ex5f=ggplot(bb, aes(x=detection, y=percent, fill=T4_ncRNA_subclass)) + 
  scale_fill_npg()+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  labs(y ="% of models", x="Detected by short-read RNA-seq", title ="RNA subclass of\nCFC-seq novel models",  fill=NULL)+
  facet_grid(cols=vars(group2),  scales="free_x", space="free_x")+
  geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
  theme1
pdf(paste0(path_fig3,"ex5f.ncRNA_subclass_detection_short_read_RNAseq.pdf"), width = 1.9, height = 1.7)
print(ex5f)
dev.off() 


#======================
#Needed
#ex5g
together=read.delim(paste0(path_fig3_data,"short_long_transcript_TPM.tsv.gz"), header=T, stringsAsFactors = F, sep="\t", check.names = F)

together=together[which(rowSums(together[,c(3,4)])>0),]
together$Cell=factor(together$Cell, levels=c("iPSC","NSC","Neuron"))
together2=together%>%group_by(Cell)%>%summarise(r=cor.test(log10(Kallisto+0.01),log10(Bambu+0.01), method="pearson")$estimate)

ex5g=ggplot()+
  geom_point_rast(data = together, aes(x = log10(Kallisto + 0.01), y = log10(Bambu + 0.01), color = Cell), alpha = 0.1, shape = 19, size = 0.1) +
  facet_grid(cols=vars(Cell))+
  scale_color_npg(guide=NULL)+
  geom_text(data=together2, mapping=aes(x=2.5, y=4.2, label=paste0("r = ",signif(r,3))),size=2, color="black")+
  labs(title="Transcript quantification between short-read and long-read", y="Long-read (log10TPM)", x="Short-read (log10TPM)")+
  theme1
pdf(paste0(path_fig3,"ex5g.Transcript_quantification_short_long.pdf"), width = 3.4, height = 1.6)
print(ex5g)
dev.off() 

#=================
#Needed
#ex5h
together=read.delim(paste0(path_fig3_data,"short_long_transcript_TPM.tsv.gz"), header=T, stringsAsFactors = F, sep="\t", check.names = F)

together=together[which(rowSums(together[,c(3,4)])>0),]
colnames(together)[c(3,4)]=c("short-read","long-read")
together1=reshape2::melt(together,id=c(1,2,5,6,7))
together1$length_bin=factor(together1$length_bin, levels=c("<=1kb","1-2kb","2-3kb","3-4kb",">4kb"))
together1$Cell=factor(together1$Cell,levels=c("iPSC","NSC","Neuron"))

ex5h=ggplot() + 
  labs(x=NULL, y = "Expression (log10TPM)",title= "Transcript quantification", fill="Library") +
  scale_fill_npg()+
  facet_wrap(vars(Cell),nrow=2, ncol=2)+
  coord_cartesian(ylim=c(-2,3))+
  geom_boxplot(data=together1, mapping=aes(y=log10(value+0.01), x=length_bin, fill=variable), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.7) + 
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex5h.expression_by_length_short.long.pdf"), width = 2.5, height = 1.9)
print(ex5h)
dev.off()

#======================
#Needed
#ex5i
table5e=read.delim(paste0(path_fig3_data,"novel.isoform.protein_coding.SQANTI3.table5.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

table5f=table5e%>%group_by(structural_category,end_class)%>%dplyr::summarise(count=n())
table5f$end_class=factor(table5f$end_class, levels=c("alt_TSS","alt_TSS&TES","alt_TES","n.a."))

ex5i1=ggplot(table5f, aes(x=structural_category, y=count, fill=end_class)) + 
  scale_fill_manual(values=c("#E64B35FF","#3C5488FF","#4DBBD5FF","grey"))+
  labs(y ="Number of models", x=NULL, title ="Class of novel isoforms\nfrom coding genes",  fill=NULL)+
  geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

table5e%>%group_by(CPAT_class)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
table5g=table5e%>%group_by(structural_category,CPAT_class)%>%dplyr::summarise(count=n())
table5g$CPAT_class=factor(table5g$CPAT_class, levels=c("non-coding","coding"))

ex5i2=ggplot(table5g, aes(x=structural_category, y=count, fill=CPAT_class)) + 
  scale_fill_manual(values=c("white","grey"))+
  labs(y ="Number of models", x=NULL, title ="Class of novel isoforms\nfrom coding genes",  fill=NULL)+
  geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

pdf(paste0(path_fig3,"ex5i.novel.isoform.protein_coding.SQANTI3.table5.pdf"), width = 3.65, height = 1.7)
grid.arrange(ex5i1,ex5i2, ncol=2, widths=c(1.85,1.8))
dev.off() 


#======================
#Needed
#ex5j
table5e=read.delim(paste0(path_fig3_data,"novel.isoform.protein_coding.SQANTI3.table5.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

table5h=table5e[which(table5e$alt_TES == "Yes"),]%>%group_by(structural_category,polyA,polyA_ENST,TES_donor_minus1)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
table5h$group=paste0(table5h$polyA_ENST,"_to_",table5h$polyA)
table5h$group=factor(table5h$group, levels=c("No_to_No","Yes_to_Yes","No_to_Yes","Yes_to_No"))
table5h$TES_donor_minus1=factor(table5h$TES_donor_minus1, levels=c("Yes","No"))

ex5j=ggplot(table5h, aes(x=group, y=count, alpha=TES_donor_minus1)) + 
  facet_grid(cols=vars(structural_category))+
  scale_fill_npg()+
  scale_alpha_manual(values=c(0.2,0.8))+
  labs(y ="Number of models", x="Poly(A) positive", title ="Novel isoforms with alternative TES",  alpha="end at\ndonor site")+
  geom_bar(linewidth=0.25, stat="identity", color="black", fill="#4DBBD5FF")+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig3,"ex5j.novel.isoform.protein_coding.SQANTI3.table5.polyA.pdf"), width = 2.5, height = 1.7)
print(ex5j)
dev.off() 


#======================
#Needed
#ex5kl
combine5=read.delim(paste0(path_fig3_data,"isoform.switch.534.hit.isoform.alone.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

combine6=combine5%>%group_by(T4_gene_ID)%>%dplyr::summarise(n_TSS=length(unique(n5_string)),n_TES=length(unique(n3_string)),
                                                            T4_Gencode_geneCalss=unique(T4_Gencode_geneCalss), iPS_NSC_hit=unique(iPS_NSC_hit),NSC_NRN_hit=unique(NSC_NRN_hit),iPS_NRN_hit=unique(iPS_NRN_hit),
                                                            Novel_transcriptClass=paste(unique(Novel_transcriptClass),collapse=";"),
                                                            Gencode_transcriptClass2=paste(unique(Gencode_transcriptClass2),collapse=";"),
                                                            n_exon=paste(unique(n_exon),collapse=";"))

n12=length(which(combine6$iPS_NSC_hit=="Yes" & combine6$NSC_NRN_hit=="Yes"))
n23=length(which(combine6$NSC_NRN_hit=="Yes" & combine6$iPS_NRN_hit=="Yes"))
n13=length(which(combine6$iPS_NSC_hit=="Yes" & combine6$iPS_NRN_hit=="Yes"))
n123=length(which(combine6$iPS_NSC_hit=="Yes" & combine6$NSC_NRN_hit=="Yes" & combine6$iPS_NRN_hit=="Yes"))
library(VennDiagram)
grid.newpage()
ex5k=draw.triple.venn(length(which(combine6$iPS_NSC_hit=="Yes")), length(which(combine6$NSC_NRN_hit=="Yes")), length(which(combine6$iPS_NRN_hit=="Yes")), n12, n23, n13, n123, category=c("iPS vs NSC","NSC vs NRN","iPS vs NRN"),
                    cex = 0.7, cat.cex = 0.75, cat.col = "black", lwd=c(0.25,0.25,0.25))
pdf(paste0(path_fig3,"ex5k.isoform_switch.venn.pdf"), width = 1.1, height = 1.1)
grid.draw(ex5k)
dev.off()

#=======
n1=which(combine6$n_TSS>1)
n2=which(combine6$n_TES>1)
n3=grep(";", combine6$n_exon)
n12=intersect(n1,n2)
n23=intersect(n2,n3)
n13=intersect(n1,n3)
n123=intersect(n1,n23)
grid.newpage()
ex5l=draw.triple.venn(length(n1), length(n2), length(n3), length(n12), length(n23), length(n13), length(n123), category=c("TSS","TES","Exon"),
                    cex = 0.7, cat.cex = 0.75, cat.col = "black", lwd=c(0.25,0.25,0.25))
pdf(paste0(path_fig3,"ex5l.isoform_switch.venn2.pdf"), width = 1.1, height = 1.1)
grid.draw(ex5l)
dev.off()


