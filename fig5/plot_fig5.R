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
path_fig5=paste0(primary_folder,"fig5/out/")
path_fig5_data=paste0(primary_folder,"fig5/data/")

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
#f5a
junction_info_TC1=read.delim(paste0(path_fig5_data,"splicing.junction.Neuron_n_THP1.from.read.n.table5.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
junction_info_TC1$GENCODEv39[which(junction_info_TC1$GENCODEv39 == "not_supported")]="not GENCODE"
junction_info_TC1$GENCODEv39[which(junction_info_TC1$GENCODEv39 == "supported")]="GENCODE"
junction_info_TC1$Supported[which(junction_info_TC1$Supported=="Short_read")]="Short-read"
junction_info_TC1$Supported[which(junction_info_TC1$Supported=="not_supported")]="Not found"

junction_info_TC1$Supported=factor(junction_info_TC1$Supported, levels=c("GENCODE","Not found","Short-read","Both","Intropolis"))

f5a=ggplot() + 
  labs(x=NULL, y = "Number of splice junction",title= "Splice junctions from\n final transcript models", fill=NULL) +
  scale_fill_manual(values=c("black","white","#E64B35FF","darkorchid3","#4DBBD5FF"))+
  scale_alpha_manual(values=c(1,1,0.8,0.9,0.8),guide=NULL)+
  facet_grid(cols=vars(GENCODEv39), scales="free_y")+
  geom_bar(data=junction_info_TC1[which(junction_info_TC1$group == "Final transcroptome"),], mapping=aes(y=count, x=site_label, fill=Supported, alpha=Supported), linewidth=0.25, width=0.7, color = "black", stat="identity") + 
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"f5a.splicing.junction.Neuron_n_THP1.from.table5.pdf"), width = 2.3, height = 1.8)
print(f5a)
dev.off()

#=====================
#Needed
#f5c
intron_gene2=read.delim(paste0(path_fig5_data,"gene_base_SJ_eff_AI.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

intron_gene22=intron_gene2%>%group_by(occurence, gene_group)%>%dplyr::summarise(count=n(), eff_junction=median(eff_junction, na.rm = T), spliceAI_junction=median(spliceAI_junction, na.rm = T ))
intron_gene2$gene_group=factor(intron_gene2$gene_group,levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
mintron_gene2=reshape2::melt(intron_gene2, id=c(1,2,3))
mintron_gene2$group1="spliceAI"
mintron_gene2$group1[grep("eff",mintron_gene2$variable)]="efficiency"
mintron_gene2$group2=sapply(strsplit(as.character(mintron_gene2$variable),"_"),"[",2)
mintron_gene2$group1=factor(mintron_gene2$group1, levels=c("spliceAI","efficiency"))
mintron_gene2$group2=factor(mintron_gene2$group2, levels=c("donor","acceptor","junction"))
mintron_gene3=mintron_gene2%>%group_by(group1,group2,occurence, gene_group)%>%dplyr::summarise(median=median(value,na.rm=T),count=n())
k4=mintron_gene2[which(mintron_gene2$gene_group %in% c("p_ncRNA","e_ncRNA")),]%>%group_by(occurence,variable)%>%dplyr::summarise(p=wilcox.test(value ~gene_group, alternative = "two.sided")$p.value)
k4$label="***"
k4$label[which(k4$p>=0.001)]="**"
k4$label[which(k4$p>=0.01)]="*"
k4$label[which(k4$p>=0.05)]="n.s."
k5=mintron_gene2[which(mintron_gene2$gene_group %in% c("mRNA","e_ncRNA")),]%>%group_by(occurence,variable)%>%dplyr::summarise(p=wilcox.test(value ~gene_group, alternative = "two.sided")$p.value)
k5$label="***"
k5$label[which(k5$p>=0.001)]="**"
k5$label[which(k5$p>=0.01)]="*"
k5$label[which(k5$p>=0.05)]="n.s."
f5c=ggplot() + 
  labs(x=NULL, y = "weighted(splice / (splice+span))",title= "Splicing efficiency per gene region", fill="Occurrence") +
  scale_fill_manual(values=c("#E64B35FF", "#3C5488FF", "#00A087FF"), guide=NULL)+
  facet_grid(cols=vars(occurence))+
  scale_y_continuous(limits=c(0,1.15), breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=mintron_gene2[which(mintron_gene2$group1 == "efficiency" & mintron_gene2$group2=="junction"),], mapping=aes(y=value, x=gene_group, fill=occurence), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, alpha=0.7, width = 0.7, outlier.shape = NA) + 
  geom_label(data=mintron_gene3[which(mintron_gene3$group1 == "efficiency" & mintron_gene3$group2=="junction"),], mapping=aes(y=0, x=gene_group, label=count, group=occurence), position = position_dodge(width = 0.9), size=1.8, angle=25, vjust=0, hjust=0.5, linewidth = NA, alpha=0.5)+
  geom_text(data=k5[which(k5$variable == "eff_junction"),], mapping=aes(y=1.06, x=2.5, label=label),size=2.4)+
  geom_text(data=k4[which(k4$variable == "eff_junction"),], mapping=aes(y=1.14, x=2, label=label),size=2.4)+
  annotate("segment", x = 2, xend = 3, y = 1.04, yend = 1.04, linewidth = 0.25) +
  annotate("segment", x = 1, xend = 3, y = 1.12, yend = 1.12, linewidth = 0.25) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"f5c.splicing.eff.by_gene.RNAclass.pdf"), width = 2.5, height = 1.7)
print(f5c)
dev.off()

#=====================
#Needed
#f5d&e
intron_gene2=read.delim(paste0(path_fig5_data,"gene_base_SJ_eff_AI.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
table5c=read.delim(paste0(path_fig5_data,"geneID_to_ex5cluster.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)

table5c=left_join(table5c, intron_gene2, by="T4_gene_ID",copy=F)

table5c=table5c[which(table5c$CpGTATA %in% c("CGI","CGIap","CGInap","Null","TATA")),]
table5c=table5c[which(table5c$ex5cluster_class %in% c("e_ncRNA","p_ncRNA")),]
table5c=table5c[which(!is.na(table5c$occurence)),]
table5c=table5c[which(table5c$gene_group == table5c$ex5cluster_class),]
table5c$CpGTATA=factor(table5c$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))
k1=table5c%>%group_by(occurence,ex5cluster_class,CpGTATA)%>%dplyr::summarise(median=median(eff_junction, na.rm=T),median2=median(spliceAI_junction, na.rm=T),count=n())
k2=table5c[which(table5c$CpGTATA %in% c("CGIap","TATA") & table5c$ex5cluster_class == "e_ncRNA"),]%>%group_by(occurence)%>%dplyr::summarise(p=wilcox.test(eff_junction ~CpGTATA, alternative = "two.sided")$p.value)
k3=table5c[which(table5c$CpGTATA %in% c("CGInap","TATA") & table5c$ex5cluster_class == "e_ncRNA"),]%>%group_by(occurence)%>%dplyr::summarise(p=wilcox.test(eff_junction ~CpGTATA, alternative = "two.sided")$p.value)
k4=table5c[which(table5c$CpGTATA %in% c("CGI","TATA") & table5c$ex5cluster_class == "p_ncRNA"),]%>%group_by(occurence)%>%dplyr::summarise(p=wilcox.test(eff_junction ~CpGTATA, alternative = "two.sided")$p.value)
k5=table5c[which(table5c$CpGTATA %in% c("Null","TATA")),]%>%group_by(occurence,ex5cluster_class)%>%dplyr::summarise(p=wilcox.test(eff_junction ~CpGTATA, alternative = "two.sided")$p.value)

k2$label="***"
k2$label[which(k2$p>=0.001)]="**"
k2$label[which(k2$p>=0.01)]="*"
k2$label[which(k2$p>=0.05)]="n.s."
k3$label="***"
k3$label[which(k3$p>=0.001)]="**"
k3$label[which(k3$p>=0.01)]="*"
k3$label[which(k3$p>=0.05)]="n.s."
k4$label="***"
k4$label[which(k4$p>=0.001)]="**"
k4$label[which(k4$p>=0.01)]="*"
k4$label[which(k4$p>=0.05)]="n.s."
k5$label="***"
k5$label[which(k5$p>=0.001)]="**"
k5$label[which(k5$p>=0.01)]="*"
k5$label[which(k5$p>=0.05)]="n.s."

f5d=ggplot() + 
  labs(x=NULL, y = "weighted(splice / (splice+span))",title= "Splicing efficiency per e_ncRNA ex5_cluster", fill="Occurrence") +
  scale_fill_manual(values=c("#E64B35FF", "#3C5488FF", "#00A087FF"), guide=NULL)+
  facet_grid(cols=vars(occurence))+
  scale_y_continuous(limits=c(0,1.15), breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=table5c[which(table5c$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=eff_junction, x=CpGTATA, fill=occurence), position = position_dodge(width = 0.9), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, alpha=0.7, width = 0.7, outlier.shape = NA) + 
  geom_label(data=k1[which(k1$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=0, x=CpGTATA, label=signif(count,3), group=occurence), position = position_dodge(width = 0.9), size=2, angle=25, vjust=0, hjust=0.5, label.size = NA, alpha=0.5)+
  geom_text(data=k2, mapping=aes(y=1.14, x=2.5, label=label),size=2.4)+
  geom_text(data=k3, mapping=aes(y=1.06, x=3, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 4, y = 1.12, yend = 1.12, size = 0.25) +
  annotate("segment", x = 2, xend = 4, y = 1.04, yend = 1.04, size = 0.25) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"f5d.splicing.eff.by_CpGTATA.pdf"), width = 2.6, height = 1.65)
print(f5d)
dev.off()

f5e=ggplot() + 
  labs(x=NULL, y = "weighted(splice / (splice+span))",title= "Splicing efficiency per p_ncRNA ex5_cluster", fill="Occurrence") +
  scale_fill_manual(values=c("#E64B35FF", "#3C5488FF", "#00A087FF"), guide=NULL)+
  facet_grid(cols=vars(occurence))+
  scale_y_continuous(limits=c(0,1.15), breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=table5c[which(table5c$ex5cluster_class == "p_ncRNA"),], mapping=aes(y=eff_junction, x=CpGTATA, fill=occurence), position = position_dodge(width = 0.9), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, alpha=0.7, width = 0.7, outlier.shape = NA) + 
  geom_label(data=k1[which(k1$ex5cluster_class == "p_ncRNA"),], mapping=aes(y=0, x=CpGTATA, label=signif(count,3), group=occurence), position = position_dodge(width = 0.9), size=2, angle=25, vjust=0, hjust=0.5, label.size = NA, alpha=0.5)+
  geom_text(data=k4, mapping=aes(y=1.14, x=2, label=label),size=2.4)+
  geom_text(data=k5[which(k5$ex5cluster_class == "p_ncRNA"),], mapping=aes(y=1.06, x=2.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 3, y = 1.12, yend = 1.12, size = 0.25) +
  annotate("segment", x = 2, xend = 3, y = 1.04, yend = 1.04, size = 0.25) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"f5e.splicing.eff.by_CpGTATA.p_ncRNA.pdf"), width = 2.6, height = 1.65)
print(f5e)
dev.off()

#=====================
#Needed
#f5f
aRNAfold5=read.delim(paste0(path_fig5_data,"polyA_all_eRNA_ATGC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aRNAfold5=reshape2::melt(aRNAfold5, id=c(6,5))
aRNAfold5$group1=sapply(strsplit(aRNAfold5$group,"_"),"[",1)
aRNAfold5$group2=sapply(strsplit(aRNAfold5$group,"_"),"[",2)
aRNAfold5$group3=sapply(strsplit(aRNAfold5$group,"_"),"[",3)
aRNAfold5$group2[which(aRNAfold5$group2 != "mRNA")]=gsub("RNA","_ncRNA", aRNAfold5$group2[which(aRNAfold5$group2 != "mRNA")])
aRNAfold5$group1=gsub("non","non-",aRNAfold5$group1)
aRNAfold5$group1=gsub("ploy","poly",aRNAfold5$group1)
aRNAfold5$group3=gsub("less3read","< 3 reads", aRNAfold5$group3)
aRNAfold5$group3=gsub("3read",">= 3 reads", aRNAfold5$group3)

aRNAfold5$group2=factor(aRNAfold5$group2,levels=c("mRNA","p_ncRNA","e_ncRNA"))
aRNAfold5$group1=factor(aRNAfold5$group1,levels=c("polyA","non-polyA"))
aRNAfold5$group3=factor( aRNAfold5$group3, levels=c(">= 3 reads", "< 3 reads"))
aRNAfold5$variable=factor(aRNAfold5$variable, levels=c("A","U","G","C"))
aRNAfold6=aRNAfold5[which(aRNAfold5$group3 == ">= 3 reads"),]

f5f=ggplot()+
  geom_line(data=aRNAfold6, mapping=aes(x=nt, y=value, color=variable), linewidth=0.2)+
  facet_grid(rows=vars(group2), cols=vars(group1), scales="free_y")+
  labs(color=NULL, x="Stranded position from 3' end", y="Frequency", title="Sequence composition")+
  scale_color_npg()+
  coord_cartesian(xlim=c(-100,100))+
  theme1+
  theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.97,0.89))
pdf(paste0(path_fig5,"f5f.non_polyA_3_eRNA_ATGC.pdf"), width = 2.3, height = 3.4)
print(f5f)
dev.off()

#=====================
#Needed
#f5g
gather1=read.delim(paste0(path_fig5_data,"RNAfold_collapse_6groups.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)

gather1$group2=factor(gather1$group2,levels=c("mRNA","p_ncRNA","e_ncRNA"))
gather1$group1=factor(gather1$group1,levels=c("polyA","non-polyA"))
gather1$group3=factor( gather1$group3, levels=c(">= 3 reads", "< 3 reads"))

f5g=ggplot()+
  geom_line(data=gather1, aes(x=nt, y=score, color=group3), linewidth=0.2)+
  facet_grid(rows=vars(group2),cols=vars(group1))+
  scale_color_manual(values=c("black","grey"))+
  labs(title="RNA secondary structure", x="Stranded position from 3' end", y="Relative hairpin score (pairiing probability)", color=NULL)+
  coord_cartesian(xlim=c(-100,100),ylim=c(0.45,0.65))+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.85,0.7))
pdf(paste0(path_fig5,"f5g.non_polyA_eRNA_pRNA_mRNA_RNAfold.pdf"), width = 2.3, height = 3.4)
print(f5g)
dev.off()

#================
#f5h

gather1=read.delim(paste0(path_fig5_data,"RNAfold_collapse_eRNACpGTATA.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
gather1=gather1[which(gather1$CpGTATA != "Others"),]
gather1$CpGTATA=factor(gather1$CpGTATA,levels=c("CGIap","CGInap","Null","TATA"))
gather1$group3=factor( gather1$group3, levels=c(">= 3 reads", "< 3 reads"))
gather1$nt=gather1$nt-201

f5h=ggplot()+
  geom_line(data=gather1, aes(x=nt, y=score, color=group3), linewidth=0.2)+
  facet_wrap(vars(CpGTATA), ncol=4, nrow=1)+
  scale_color_manual(values=c("black","grey"))+
  labs(title="Non-poly(A) eRNA secondary structure", x="Stranded position from 3' end", y="Relative hairpin score", color=NULL)+
  coord_cartesian(xlim=c(-100,100),ylim=c(0.4,0.65))+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.85,0.55))
pdf(paste0(path_fig5,"f5h.non_polyA_eRNA_CpGTATA_RNAfold.pdf"), width = 4, height = 1.3)
print(f5h)
dev.off()


