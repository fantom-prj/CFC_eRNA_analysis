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
path_fig1=paste0(primary_folder,"fig1/out/")
path_fig1_data=paste0(primary_folder,"fig1/data/")

#=======
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

#============================================================================================
#Needed
#ex1c
total_read=read.delim(paste0(path_fig1_data,"read_number_plot.tsv.gz"), header=T, stringsAsFactors = F, check.names=F, row.names=1)

colnames(total_read)[c(1:7)]=c("iPSC_PAT","NSC_PAT","Neuron_PAT","THP-1_PAT","THP-1_noPAT","dTHP-1_PAT","dTHP-1_noPAT")
total_read=data.frame(t(total_read))
total_read$raw.reads=total_read$raw.reads-total_read$aligned.reads
total_read$lib=row.names(total_read)
total_read1=reshape::melt(total_read[,c(4,1,3)],id=1)
total_read1$lib=factor(total_read1$lib, levels=c("iPSC_PAT","NSC_PAT","Neuron_PAT","THP-1_PAT","THP-1_noPAT","dTHP-1_PAT","dTHP-1_noPAT"))
total_read1$variable=as.character(total_read1$variable)
total_read1$variable[c(1:7)]="incomplete /\nnot mappable"
total_read1$variable[c(8:14)]="mappable"
total_read1$variable=factor(total_read1$variable, levels=c("incomplete /\nnot mappable", "mappable"))

ex1c=ggplot(total_read1, aes(x=lib, y=value, fill=variable)) + 
  scale_fill_manual(values=c("white", "grey"))+
  scale_y_continuous(breaks=c(0,25000000,50000000,75000000,100000000), labels=c("0", "25M", "50M", "75M", "100M"))+
  labs(fill=NULL, x=NULL, y ="Number of reads", title="Raw to mappable reads")+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.7)+
  theme1+theme(axis.text.x = element_text(color ="black", angle = 25, hjust=1, vjust=1),legend.position=c(0.8,0.8))
pdf(paste0(path_fig1,"ex1c.read.number.pdf"), width = 1.5, height = 1.6)
print(ex1c)
dev.off() 


#============================================================================================
#Needed
#ex1d
CRE.hit.rate=read.delim(paste0(path_fig1_data,"CRE.hit.rate.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

CRE.hit.rate$label=sapply(strsplit(CRE.hit.rate$name,"_EN"),"[",1)
CRE.hit.rate$label[grep("iPS",CRE.hit.rate$label)]="CFC-seq"
CRE.hit.rate$label=gsub("TOS","TSO",CRE.hit.rate$label)
CRE.hit.rate$label=gsub("Ont","ONT",CRE.hit.rate$label)
CRE.hit.rate$label=gsub("Pacbio","PB",CRE.hit.rate$label)
CRE.hit.rate$label=gsub("CAP_trap","CAP-trap",CRE.hit.rate$label)
CRE.hit.rate$label=gsub("uncap_deplete","uncap-deplete",CRE.hit.rate$label)

CRE.hit.rate1=CRE.hit.rate%>%group_by(group,label)%>%dplyr::summarise(rate=mean(rate), rate2=mean(rate2))
CRE.hit.rate1$group[which(CRE.hit.rate1$group == "cCRE")]="SCREEN cCRE"
CRE.hit.rate1$group[which(CRE.hit.rate1$group == "F5_CAGE")]="F5 CAGE"
CRE.hit.rate1$group=factor(CRE.hit.rate1$group,levels=c("SCREEN cCRE","ATAC","F5 CAGE"))
CRE.hit.rate1$label=factor(CRE.hit.rate1$label, levels=c("CFC-seq","CAP-trap_PB","uncap-deplete_PB","TSO_ONT","R2C2_ONT","dRNA_ONT"))

ex1d=ggplot(CRE.hit.rate1[which(CRE.hit.rate1$group %in% c("SCREEN cCRE","ATAC","F5 CAGE")),], aes(x=reorder(label,-rate), y=rate, fill=label)) + 
  #scale_fill_manual(values=c("other"="grey","current"="black"),guide="none")+
  scale_fill_npg()+
  labs(x="Different long read protocol on iPSC (WTC11)", y ="Percentage of reads", title="mRNA reads start from annotated regions", fill="Protocol")+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  geom_bar(stat="identity")+
  facet_grid(cols = vars(group))+
  theme1+theme(axis.text.x = element_blank())
pdf(paste0(path_fig1,"ex1d.5end_coverage.pdf"), width = 3, height = 1.6)
print(ex1d)
dev.off() 

#================
#Needed
#ex1e
all6=read.delim(paste0(path_fig1_data,"non_polyA_gene_group.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
all6$group=factor(all6$group, levels=c("Neuron.series_PAT", "THP-1.series_PAT", "THP-1.series_noPAT"))

ex1e=ggplot(all6[which(all6$class != "others"),], aes(x=group, y=percent, fill=class)) + 
  scale_fill_npg()+
  labs(x=NULL, y ="Percentage of read", title="Known non-\npoly(A) gene", fill=NULL)+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.7)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig1,"ex1e.Neuron_THP1_nonployAclass.pdf"), width = 1.6, height = 1.6)
print(ex1e)
dev.off() 


#=====================
#Needed
#ex1f
IPplot=read.delim(paste0(path_fig1_data,"Neuron_THP1_IP.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

IPplot$PAtailing=factor(IPplot$PAtailing, levels=c("Neuron.series_PAT", "THP-1.series_PAT","THP-1.series_noPAT"))
IPplot$internal_prime2=factor(IPplot$internal_prime2, levels=c("yes","no"))

ex1f=ggplot(IPplot, aes(x=PAtailing, y=percent, fill=internal_prime2)) + 
  scale_fill_manual(values=c("white","grey"), guide="none")+
  labs(x=NULL, y ="%Read without internal priming", title="Internal priming")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth = 0.25)+
  theme1 + theme(axis.text.x = element_text(color ="black", angle = 25, hjust=1, vjust=1))
pdf(paste0(path_fig1,"ex1f.Neuron_THP1_IP.pdf"), width = 1.3, height = 1.6)
print(ex1f)
dev.off() 


#=====================
#Needed
#ex1g
plotTT9=read.delim(paste0(path_fig1_data,"Neuron_THP1_PAS_geneClass.tsv.gz"), header=T, check.names =F, stringsAsFactors = F)
#excluded potential internal primed

plotTT9$geneClass_T1=gsub("lncRNA","ncRNA",plotTT9$geneClass_T1)
plotTT9$geneClass_T1=factor(plotTT9$geneClass_T1, levels=c("mRNA","p_ncRNA", "e_ncRNA"))
plotTT9$group=factor(plotTT9$group, levels=c("Neuron.series_PAT", "THP-1.series_PAT", "THP-1.series_noPAT"))
plotTT9a=plotTT9%>%group_by(group,geneClass_T1,PAS_distance2)%>%dplyr::summarise(percent=sum(percent))
plotTT10=plotTT9[which(plotTT9$PAS_distance2 >= -1.544068 & plotTT9$PAS_distance2 <=(-0.69897)),]%>%group_by(group, geneClass_T1)%>%dplyr::summarise(count=sum(percent))

ex1g=ggplot(plotTT9a, aes(x=PAS_distance2, y=percent, fill=geneClass_T1)) + 
  scale_fill_npg(guide="none")+
  scale_x_continuous(breaks=c(-3,-1,0,1,3), labels=c("-1000","-10","","10","1000"), limits=c(-3,3))+
  labs(x="distance from read 3' end (nt) log10 scale", y ="Percentage of reads", title="Location of polyadenylation signal", fill="RNA class")+
  facet_grid(cols = vars(group), rows = vars(geneClass_T1), scales = "fixed", space = "fixed")+
  scale_y_continuous(labels = scales::percent, breaks=c(0,0.2,0.4))+
  geom_bar(stat="identity", position = "identity", width=0.15)+
  geom_label(data=plotTT10, aes(x=(-0.69897), y=0.35, label=paste0(signif(count*100, 3),"%")), size=2.2, alpha=0.4, hjust=0, fill="white", label.size=0)+
  geom_vline(xintercept=c(-1.544068, -0.69897), linetype="dashed", color = "black", linewidth=0.2)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig1,"ex1g.Neuron_THP1_PAS_geneClass.pdf"), width = 2.8, height = 1.9)
print(ex1g)
dev.off() 


#=====================
#Needed
#ex1h
allsummary=read.delim(paste0(path_fig1_data,"neuron_THP1_TES_genomic_location.summary.tsv.gz"), header=T, stringsAsFactors =F, check.names=F)
allsummary$label[which(allsummary$class2=="Intergenic/others")]=NA

ex1h=ggplot(allsummary[which(allsummary$geneClass_T1 == "mRNA"),], aes(y=percent2, x=PAtailing, fill=class2)) + 
  scale_fill_npg()+
  labs(x=NULL, y ="Percentage of all read", title="mRNA without PAS motif", fill=NULL)+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, alpha=0.7)+
  geom_text(aes(label=label), position = position_stack(vjust = 0.5), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))

pdf(paste0(path_fig1,"ex1h.mRNA_without_PAS.pdf"), width = 1.7, height = 1.9)
print(ex1h)
dev.off()

#=====================
#Needed
#ex1i
all7c=read.delim(paste0(path_fig1_data,"noPAT_PATspecific_TES_location.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

all7d=all7c%>%group_by(PAtailing, group2, geneClass_T1,class2)%>%dplyr::summarise(count=sum(TES_count))%>%dplyr::mutate(percent=count/sum(count))
all7d=left_join(all7d, allsummary3, by=c("PAtailing","geneClass_T1"),copy=F)
all7d$percent2=all7d$count/all7d$total
all7d$group2=factor(all7d$group2, levels=c("no PAT","PAT-specific: Partial end", "PAT-specific: Others"))
all7d$class2=factor(all7d$class2, levels=c("GENCODE_3n","Intron","Donor_site","Exon","3'UTR","Intergenic/others"))
all7d$scale=all7d$percent/all7d$percent2

all7e=all7c%>%group_by(PAtailing, geneClass_T1,class2)%>%dplyr::summarise(count=sum(TES_count))%>%dplyr::mutate(percent=count/sum(count))
all7e=left_join(all7e, allsummary3, by=c("PAtailing","geneClass_T1"),copy=F)
all7e$percent2=all7e$count/all7e$total
all7e$PAtailing=factor(all7e$PAtailing, levels=c("woPAP","PAP"))
all7e$class2=factor(all7e$class2, levels=c("GENCODE_3n","Intron","Donor_site","Exon","3'UTR","Intergenic/others"))
all7e$scale=all7e$percent/all7e$percent2

#fac1=unique(all7d$scale[which(all7d$group2== "no PAT" & all7d$geneClass_T1 == "mRNA")])
ex1i1=ggplot(all7d[which(all7d$group2== "no PAT" & all7d$geneClass_T1 == "mRNA"),], aes(y=percent2, x=1, fill=class2)) + 
  scale_fill_npg()+
  labs(x=NULL, y ="Percentage of all read", title="no PAT" ,fill=NULL)+
  scale_y_continuous(labels = scales::percent)+
  #scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~ . / fac1, name = "Percentage of all read", labels = scales::percent))+
  geom_bar(stat="identity", color="black", linewidth=0.25, alpha=0.7)+
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#fac4=unique(all7e$scale[which(all7e$PAtailing == "PAP" & all7e$geneClass_T1 == "mRNA")])
ex1i2=ggplot(all7e[which(all7e$PAtailing == "PAP" & all7e$geneClass_T1 == "mRNA"),], aes(y=percent2, x=1, fill=class2)) + 
  scale_fill_npg()+
  labs(x=NULL, y ="Percentage of all read", title="PAT-specific", fill=NULL)+
  scale_y_continuous(labels = scales::percent)+
  #scale_y_continuous(labels = scales::percent, sec.axis = sec_axis(~ . / fac4, name = "Percentage of all read", labels = scales::percent))+
  geom_bar(stat="identity", color="black", linewidth=0.25, alpha=0.7)+
  theme1+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

pdf(paste0(path_fig1,"ex1i.noPAT_PATspecific_TES_location.pdf"), width = 1.2, height = 1.9)
grid.arrange(arrangeGrob(ex1i1, ex1i2, ncol=1, nrow = 2, heights = c(1,1)))
dev.off() 

result4=read.delim(paste0(path_fig1_data,"THP1_w.wo.PAT_ATGC_sumnmary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

result4$group1=sapply(strsplit(result4$group,"_"),"[",1)
result4$group1[which(result4$group1 == "noPAT.noPAS")]="no PAT"
result4$group1[which(result4$group1 == "PAT.noPAS")]="PAT-specific"
result4$group2=sapply(strsplit(result4$group,"_"),"[",2)
result4=result4[which(result4$group2 != "exclude GENCODE.singular"),]

colnames(result4)[c(3:6)]=c("A","U","G","C")
mresult4=reshape2::melt(result4[,c(2:8)], id=c(6,7,1))
ex1i3=ggplot() + 
  labs(x="Stranded distance from the 3'end (nt)",  y="Sequence composition (%)", color=NULL)+
  geom_line(data=mresult4, mapping=aes(x=nt, color=variable, y=value), linewidth=0.3)+
  scale_color_npg()+
  scale_y_continuous(labels = scales::percent, limits=c(0.15,0.4), breaks=c(0.2,0.3,0.4))+
  scale_x_continuous(limits=c(-30,30))+
  facet_wrap(vars(group1,group2), nrow=2, ncol=2, scales = "free")+
  theme1
pdf(paste0(path_fig1,"ex1i3.THP1_noPAS_ATGC.pdf"), width = 2.2, height = 1.9)
print(ex1i3)
dev.off() 

#=====================
#Needed
#ex1j
all7c=read.delim(paste0(path_fig1_data,"noPAT_PATspecific_TES_location.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

all7q=all7c[which(all7c$geneClass_T1 == "mRNA" & all7c$PAtailing == "PAP"),]%>%group_by(class2,TES_count)%>%dplyr::summarise(count=n())
all7q$TES_count[which(all7q$TES_count>=3)]=">=3"
all7q=all7q%>%group_by(class2,TES_count)%>%dplyr::summarise(count=sum(count))%>%dplyr::mutate(percent=count/sum(count))
all7q$TES_count=factor(all7q$TES_count, levels=c(1,2,">=3"))
all7q$class2=factor(all7q$class2, levels=c("GENCODE_3n","Intron","Donor_site","Exon","3'UTR","Intergenic/others"))

ex1j=ggplot() + 
  labs(x="Frequency per locus", y = "% of read",title= "PAT-specific 3' end", fill=NULL) +
  scale_fill_npg()+
  geom_bar(data=all7q, mapping=aes(x=TES_count, y=percent, fill=class2), stat="identity", position_dodge(), linewidth=0.25, color="black")+
  scale_y_continuous(labels = scales::percent)+
  theme1 + theme(legend.position = c(0.75,0.7), legend.direction = "vertical")
pdf(paste0(path_fig1,"ex1j.nonPAS_TES_distribution.PATspecific.pdf"), width = 1.8, height = 1.5)
print(ex1j)
dev.off()


#=====================
#Needed
#ex1k
vennsummary=read.delim(paste0(path_fig1_data,"noPAT_PATspecific_venn.summary.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
vennsummary$group=factor(vennsummary$group, levels=c("PAT-specific","common","noPAT-specific"))
vennsummary$geneClass_T5=factor(vennsummary$geneClass_T5, levels=c("others","lncRNA","mRNA"))

vennsummary1=vennsummary%>%group_by(group)%>%dplyr::summarise(count=sum(count))

library(VennDiagram)
grid.newpage()
R9 <- draw.pairwise.venn(area1 = vennsummary1$count[which(vennsummary1$group=="PAT-specific")]+vennsummary1$count[which(vennsummary1$group=="common")], 
                         area2 = vennsummary1$count[which(vennsummary1$group=="noPAT-specific")]+vennsummary1$count[which(vennsummary1$group=="common")], 
                         cross.area = vennsummary1$count[which(vennsummary1$group=="common")],
                         euler.d = TRUE, scaled = T, inverted=F, category=c("PAT-specific","noPAT-specific"),  alpha=c(0.5,0.5), rotation.degree = 270,
                         fill = c("#3C5488FF","#7E6148FF"), lty = "blank", cex = 0.6, cat.cex = 0.7, cat.col = "black")
ex1k1 <- grid::gTree(children = R9)
pdf(paste0(path_fig1,"ex1l1.gene_content.pdf"), width = 1.2, height = 1.2)
grid.arrange(ex1k1)
dev.off()

ex1k2=ggplot(vennsummary, aes(y=count, x=geneClass_T5, fill=geneClass_T5)) + 
  labs(x=NULL, y = "Number of gene",title= "THP-1 RNA class with and without PAT", fill=NULL) +
  scale_fill_manual(values=c("#00A087FF","#4DBBD5FF","#E64B35FF"),guide=NULL)+
  facet_wrap(vars(group), scale="fixed", nrow=3)+
  geom_bar(color="black", linewidth=0.25, stat="identity")+
  coord_flip()+
  theme1+ theme(axis.text.x = element_text(color ="black", angle=0, hjust=0.5, vjust=0.5))

pdf(paste0(path_fig1,"ex1k2.gene_content.pdf"), width = 1.4, height = 1.6)
print(ex1k2)
dev.off()

#=====================
#Needed
#ex1l
roc_df1=readRDS(paste0(path_fig1_data,"roc_curve.RDS"))
auc1=0.941
auc2=0.932

ex1l=ggplot(roc_df1, aes(x = 1 - specificity, y = sensitivity, color = group)) +
  geom_line(linewidth=0.25) +
  scale_color_manual(values=c("#E64B35FF","#00A087FF")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth=0.2)+
  annotate("text", x = 0.1, y = 0.75, label = paste("FLAM-seq AUC =",auc1), size = 1.8, color = "#E64B35FF", hjust=0)+
  annotate("text", x = 0.1, y = 0.6, label = paste("THP-1 AUC =",auc2), size = 1.8, color = "#00A087FF", hjust=0)+
  labs(title = "ROC Curve of poly(A) prediction", x = "False Positive Rate (1 - Specificity)", y = "True Positive Rate (Sensitivity)", color="Ground true") +
  theme1

pdf(paste0(path_fig1,"ex1l.ROC_1_2.pdf"), width = 2, height = 1.5)
print(ex1l)
dev.off()



