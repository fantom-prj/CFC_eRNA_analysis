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

#============================================================================================
#Needed
#f3b
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)

table5$promoter_type[which(table5$promoter_type%in% c("CTCF-alone","unclassed","excluded"))]="others"
table5$promoter_type=gsub("-like","",table5$promoter_type)
table5$T4_gene_novelty=gsub("novel","Novel gene", table5$T4_gene_novelty)
table5$T4_gene_novelty=gsub("GENCODE","ENSG", table5$T4_gene_novelty)
table5$T4_gene_novelty=factor(table5$T4_gene_novelty, levels=c("Novel gene","ENSG"))

table5$transcript_class="Transcripts from\nnovel genes"
table5$transcript_class[grep("ENST",table5$model_ID)]="ENST"
table5$transcript_class[intersect(grep("ONTT",table5$model_ID),grep("ENSG",table5$T4_gene_ID))]="Novel isoforms"

piet=table5%>%group_by(T4_gene_novelty, transcript_class, promoter_type)%>%dplyr::summarise(count=n())
piet$transcript_class=factor(piet$transcript_class, levels=c("Transcripts from\nnovel genes","Novel isoforms","ENST"))
piet$promoter_type=factor(piet$promoter_type, levels=c("enhancer","promoter","others"))

piet_t=piet%>%group_by(transcript_class)%>%dplyr::summarise(count=sum(count),group="transcript_class")%>%mutate(ymin = c(0, cumsum(count)[-length(count)]), ymax = cumsum(count), xmin = 0.5, xmax = 1.4)%>%rename(category = transcript_class)
piet_g=piet%>%group_by(T4_gene_novelty)%>%dplyr::summarise(count=sum(count),group="T4_gene_novelty")%>%mutate(ymin = c(0, cumsum(count)[-length(count)]), ymax = cumsum(count), xmin = 1.41, xmax = 1.55)%>%rename(category = T4_gene_novelty)
piet_p=piet%>%group_by(transcript_class,promoter_type)%>%dplyr::summarise(count=sum(count),group="promoter_type",.groups = "drop")%>%mutate(ymin = c(0, cumsum(count)[-length(count)]), ymax = cumsum(count), xmin = 1.55, xmax = 1.69)%>%rename(category = promoter_type)
combine=rbind(piet_t, piet_g, piet_p[,c(2:8)])
combine$label1=combine$category
combine$label2=combine$count
combine$label1[which(combine$group == "promoter_type")]=NA
combine$label2[which(combine$group == "promoter_type" | combine$group == "T4_gene_novelty")]=NA

f3b1 <- ggplot(combine) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category, alpha=group), color="white", linewidth=0.2) +
  scale_fill_manual(values=c("ENSG"="lightgreen","Novel gene"="turquoise3","Novel isoforms"="lightcoral","ENST"="lightgreen","Transcripts from\nnovel genes"="violet","enhancer"="black","promoter"="grey", "others"="white"), breaks=c("enhancer","promoter"))+
  scale_alpha_manual(values=c("transcript_class"=0.5,"T4_gene_novelty"=0.8, "promoter_type"=1), guide=NULL) +
  coord_polar("y", start = 4) + 
  theme_void() +
  geom_text(aes(x = xmin+(xmax-xmin)/1.5-0.02, y = (ymin + ymax) / 2, label = label1), color = "black", size=2, lineheight = 0.8, vjust=-1) +  
  geom_text(aes(x = xmin+(xmax-xmin)/1.5-0.02, y = (ymin + ymax) / 2, label = label2), color = "darkgrey", size=1.8, lineheight = 0.8, vjust=0.1) + 
  theme(legend.position = "right", legend.key.size = unit(0.2, 'cm'), plot.title = element_text(hjust = 0.5,  color ="black"), text = element_text(size=6)) +
  labs(title = NULL, fill="promoter type")
pdf(paste0(path_fig3,"f3b1.promoter_type_pie.pdf"), width = 2.4, height = 1.7)
print(f3b1)
dev.off()

table5$T4_gene_promoter_type[which(table5$T4_gene_promoter_type%in% c("CTCF-alone","unclassed","excluded"))] <- "others"
table5$T4_gene_promoter_type <- gsub("-like","",table5$T4_gene_promoter_type)
table5g=unique(table5[,c("T4_gene_ID","T4_gene_novelty","T4_gene_promoter_type")])

pieg <- table5g%>%group_by(T4_gene_novelty, T4_gene_promoter_type)%>%dplyr::summarise(count=n())
pieg$T4_gene_promoter_type=factor(pieg$T4_gene_promoter_type, levels=c("enhancer","promoter","others"))

pieg_g <- pieg%>%group_by(T4_gene_novelty)%>%dplyr::summarise(count=sum(count),group="T4_gene_novelty")%>%mutate(ymin = c(0, cumsum(count)[-length(count)]), ymax = cumsum(count), xmin = 0.5, xmax = 1.4)%>%rename(category = T4_gene_novelty)
pieg_p <- pieg%>%group_by(T4_gene_novelty,T4_gene_promoter_type)%>%dplyr::summarise(count=sum(count),group="T4_gene_promoter_type",.groups = "drop")%>%mutate(ymin = c(0, cumsum(count)[-length(count)]), ymax = cumsum(count), xmin = 1.41, xmax = 1.55)%>%rename(category = T4_gene_promoter_type)
combine1 <- rbind( pieg_g, pieg_p[,c(2:8)])
combine1$label1 <- combine1$category
combine1$label2 <- combine1$count
combine1$label1[which(combine1$group == "T4_gene_promoter_type")] <- NA
combine1$label2[which(combine1$group == "T4_gene_promoter_type")] <- NA

f3b2 <- ggplot(combine1) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = category, alpha=group), color="white", linewidth=0.2) +
  scale_fill_manual(values=c("ENSG"="lightgreen","Novel gene"="turquoise3","Novel isoforms"="lightcoral","ENST"="lightgreen","Transcripts from\nnovel genes"="violet","enhancer"="black","promoter"="grey", "others"="white"), breaks=c("enhancer","promoter"))+
  scale_alpha_manual(values=c("T4_gene_novelty"=0.8, "T4_gene_promoter_type"=1), guide=NULL) +
  coord_polar("y", start = 4) + 
  theme_void() +
  geom_text(aes(x = xmin+(xmax-xmin)/1.5-0.02, y = (ymin + ymax) / 2, label = label1), color = "black", size=2, lineheight = 0.8, vjust=-1) +  
  geom_text(aes(x = xmin+(xmax-xmin)/1.5-0.02, y = (ymin + ymax) / 2, label = label2), color = "darkgrey", size=1.8, lineheight = 0.8, vjust=0.1) + 
  theme(legend.position = "none", legend.key.size = unit(0.2, 'cm'), plot.title = element_text(hjust = 0.5,  color ="black"), text = element_text(size=6)) +
  labs(title = NULL, fill="promoter type")
pdf(paste0(path_fig3,"f3b2.promoter_type_pie.pdf"), width = 1.2, height = 1.2)
print(f3b2)
dev.off()

#=====================
#Needed
#f3c
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)

f3c=ggplot(table5, aes(y=Coding_prob, color=factor(transcript_novelty))) + 
  labs(y="Coding probability", x ="Ranked transcripts", title="CPAT coding potential", color=NULL)+
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(labels=c("ENST","Novel isoform","Transcripts from novel genes"),values=c("chartreuse4","orangered1","purple"))+
  geom_hline(yintercept=c(0.364), linetype="dashed", color = "black", linewidth=0.25)+
  scale_x_continuous(labels = scales::percent)+
  stat_ecdf(geom = "step", linewidth=0.25)+
  theme1+theme(legend.direction = "vertical",legend.position = "bottom")
pdf(paste0(path_fig3,"f3c.CPAT_coding_potential.pdf"), width = 1.5, height = 1.7)
print(f3c)
dev.off() 

#=====================
#Needed
#f3d
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)

table5$promoter_type[which(table5$promoter_type %in% c("CTCF-alone","excluded"))]="unclassed"
lncRNA.summary=table5[which(!is.na(table5$T4_ncRNA_subclass)),]%>%group_by(T4_gene_novelty, T4_ncRNA_subclass, promoter_type)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
lncRNA.summary$T4_gene_novelty[which(lncRNA.summary$T4_gene_novelty == "GENCODE")]="ENSG"
lncRNA.summary$T4_gene_novelty[which(lncRNA.summary$T4_gene_novelty == "novel")]="Novel genes"
lncRNA.summary2=lncRNA.summary%>%group_by(T4_gene_novelty, T4_ncRNA_subclass)%>%dplyr::summarise(count=sum(count))
lncRNA.summary$promoter_type=factor(lncRNA.summary$promoter_type, levels=c("unclassed","promoter-like","enhancer-like"))

f3d=ggplot() + 
  scale_fill_manual(values=c("enhancer-like"="black","promoter-like"="grey","unclassed"="white"))+
  labs(x=NULL, y ="Number of transcripts", title="Subclass of ncRNA transcripts", fill=NULL)+
  geom_bar(data=lncRNA.summary, mapping=aes(x=T4_ncRNA_subclass, fill=promoter_type, y=count), color="black", stat = "identity", alpha=1, linewidth=0.25, width=0.8)+
  facet_grid(cols = vars(T4_gene_novelty), scale="free", space="free")+
  scale_x_discrete(limits=rev)+
  coord_flip()+
  theme1+theme(legend.position = "bottom", axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig3,"f3d.ncRNA.subclass.pdf"), width = 3.7, height = 1.5)
print(f3d)
dev.off() 

#=====================
#Needed
#f3e
gene=read.delim(paste0(path_fig3_data,"short_long_gene_TPM.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gene=gene[which(rowSums(gene[,c(3,4)])>0),]
gene$Cell=factor(gene$Cell, levels=c("iPSC","NSC","Neuron"))
gene2=gene%>%group_by(Cell)%>%summarise(r=cor.test(log10(Kallisto+0.01),log10(Bambu+0.01), method="pearson")$estimate)

f3e=ggplot()+
  geom_point_rast(data=gene, aes(x=log10(Kallisto+0.01), y=log10(Bambu+0.01), color=Cell), alpha=0.1, shape=19, size=0.1)+
  facet_grid(cols=vars(Cell))+
  scale_color_npg(guide=NULL)+
  geom_text(data=gene2, mapping=aes(x=2.5, y=4.2, label=paste0("r = ",signif(r,3))),size=2, color="black")+
  labs(title="Gene quantification between short- & long-read", y="Long-read (log10TPM)", x="Short-read (log10TPM)")+
  theme1
pdf(paste0(path_fig3,"f3e.Gene_quantification_short_long.pdf"), width = 2.4, height = 1.3)
print(f3e)
dev.off() 

#=====================
#Needed
#f3f
combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1192.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
combine2$structural_category[which(combine2$structural_category %in% c("alternative_TES", "alternative_TSS","alternative_TSS&TES"))]="full-splice-match"
combine2$structural_category[grep("ISM",combine2$structural_category)]="incomplete-splice-match"
combine2$structural_category=gsub("_","-",combine2$structural_category)

combine2a=combine2%>%group_by(structural_category,end_class)%>%dplyr::summarise(count=n())
combine2a$structural_category=factor(combine2a$structural_category, levels=c("GENCODE", "full-splice-match", "incomplete-splice-match", "novel-in-catalog", "novel-not-in-catalog"))
combine2a$end_class=factor(combine2a$end_class, levels=c("alt_TSS","alt_TSS&TES","alt_TES","n.a."))

f3f=ggplot(combine2a, aes(x=structural_category, y=count, fill=end_class)) + 
  scale_fill_manual(values=c("#E64B35FF","#3C5488FF","#4DBBD5FF","grey"))+
  labs(y ="Number of transcripts", x=NULL, title ="Switched isoforms (1,192)",  fill=NULL)+
  geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),
                 plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
pdf(paste0(path_fig3,"f3f.Switched_transcript.SQANTI3.table5.ends.pdf"), width = 1.65, height = 1.3)
print(f3f)
dev.off() 

