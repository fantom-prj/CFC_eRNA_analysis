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
#ex7a
support_sets1=readRDS(paste0(path_fig5_data,"all_support_intro.summary.RDS"))
support_sets2=support_sets1%>%group_by(Supported)%>%dplyr::summarise(count=sum(count))

library(ggupset)
library(ComplexUpset)

ex7a=ggplot() + 
  labs(x=NULL, y = "Count (log10 scale)",title= "Supportive splicing junctions", fill=NULL) +
  scale_fill_npg()+
  geom_bar(data=support_sets2, mapping=aes(y=count, x=Supported), linewidth=0.25, width=0.7, color = "black", stat="identity") + 
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_x_upset()+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.box.margin = margin(t=-10, r=0, b=0, l=0))+
  guides(fill = guide_legend(ncol = 2))+
  theme_combmatrix(combmatrix.panel.point.color.fill = "black",
                   combmatrix.panel.point.size = 1,
                   combmatrix.panel.line.size = 0.5,
                   combmatrix.label.text = element_text(color ="black", size=4.5),
                   combmatrix.label.extra_spacing = 0.3,
                   combmatrix.label.make_space = FALSE)
pdf(paste0(path_fig5,"ex7a.Supported.splicing.junction.class.pdf"), width = 1.1, height = 1.7)
print(ex7a)
dev.off()

#=====================
#Needed
#ex7b
junction_info_TC_no=read.delim(paste0(path_fig5_data,"All_splicing.junction.Neuron_n_THP1.from.read.n.table5.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

junction_info_TC_no1=junction_info_TC_no[which(junction_info_TC_no$Support=="Novel"),]%>%group_by(canonical)%>%dplyr::summarise(count=length(which(max_score_per_nt >= 10 & total_count >= 3)),percent=length(which(max_score_per_nt >= 10 & total_count >= 3))/n(), total=n())
junction_info_TC_no1$label2=paste0("Passed:\nn=",junction_info_TC_no1$count,"\n",signif(junction_info_TC_no1$percent*100,2),"%")

ex7b=ggplot() + 
  labs(x="Count", y = "Maximum base-calling score", title= "Quality of novel SJs from final transcript", fill=NULL) +
  scale_color_npg()+
  facet_grid(cols=vars(canonical), scales="fixed")+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  geom_point_rast(data=junction_info_TC_no[which(junction_info_TC_no$Support=="Novel"),], mapping=aes(y=max_score_per_nt, x=total_count, color=canonical), size=0.25, alpha=0.25) + 
  geom_text(data=junction_info_TC_no1, mapping=aes(x=1000, y=10.5, label=label2), size=2, vjust=0)+
  geom_hline(yintercept=10, linewidth=0.25)+
  geom_vline(xintercept=3, linewidth=0.25)+
  theme1+theme(legend.position = "none")
pdf(paste0(path_fig5,"ex7b.score.count.splicing.junction.Neuron_n_THP1.from.read.n.table5.pdf"), width = 2.2, height = 1.6)
print(ex7b)
dev.off()

#=====================
#Needed
#ex7c
table5c=read.delim(paste0(path_fig5_data,"splicing_support_group.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

table5c$Exon=factor(table5c$Exon, levels=c("Un-spliced","Spliced"))
tt0=table5c[which(table5c$group %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),]%>%group_by(Exon,group,All_SJ_support_confidence)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
tt0$All_SJ_support_confidence[which(is.na(tt0$All_SJ_support_confidence))]="Single exon"
tt0$All_SJ_support_confidence=factor(tt0$All_SJ_support_confidence, levels=c("Single exon","Partial/not supported","All SJ supported"))
tt0$group=factor(tt0$group, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
tt0$label=paste0(signif(tt0$percent*100,3),"%")
tt0$label[which(tt0$Exon == "Un-spliced")]=NA

ex7c=ggplot() + 
  labs(x=NULL, y = "Number of transcript",title= "Splicing of RNA class", fill=NULL) +
  scale_fill_npg()+
  facet_grid(cols=vars(Exon))+
  geom_bar(data=tt0, mapping=aes(y=count, x=group, fill=All_SJ_support_confidence), linewidth=0.25, color = "black", stat="identity", width=0.85, alpha=0.7) + 
  geom_text(data=tt0, mapping=aes(y=count, x=group, group=All_SJ_support_confidence, label=label), size=1.8, color="black", position = position_stack(vjust = 0.3, reverse=F))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.25,0.75))
pdf(paste0(path_fig5,"ex7c.splicing_support_group.pdf"), width = 2.4, height = 1.7)
print(ex7c)
dev.off()

#=====================
#Needed
#ex7d
t5_intron=read.delim(paste0(path_fig5_data,"splicing.junction.RNAclass.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

t5_intron=t5_intron[which(t5_intron$group != "CTCF_ncRNA"),]
summary=t5_intron%>%group_by(group,site_label)%>%dplyr::summarise(count=n())%>%mutate(percent=count/sum(count), total=sum(count))
summary1=unique(summary[,c(1,5)])
summary1$label=paste0("n = ",signif(summary1$total/1000,3),"k")
summary1$group=factor(summary1$group, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))

ex7d=ggplot() + 
  labs(x=NULL, y = "Percentage",title= "SJ class from RNA class", fill=NULL) +
  scale_fill_npg()+
  scale_y_continuous(labels = scales::percent, limits=c(0,1.2), breaks=c(0,0.25,0.5,0.75,1))+
  geom_bar(data=summary, mapping=aes(y=percent, x=group, fill=site_label), linewidth=0.25, width=0.9, color = "black", stat="identity", alpha=0.8) + 
  geom_text(data=summary1, mapping=aes(x=group, y=1, label=label), angle=25, size=1.8, vjust=0, hjust=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position="right", legend.direction ="vertical")
pdf(paste0(path_fig5,"ex7d.splicing.junction.RNAclass.pdf"), width = 1.6, height = 1.7)
print(ex7d)
dev.off()

#=====================
#Needed
#ex7e
intron_gene2=read.delim(paste0(path_fig5_data,"gene_base_SJ_eff_AI.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

intron_gene3=intron_gene2%>%group_by(occurence)%>%dplyr::summarise(rho=cor.test(eff_junction,spliceAI_junction, method="spearman")$estimate,
                                                                   p=cor.test(eff_junction,spliceAI_junction, method="spearman")$p.value)
intron_gene3$label=paste0("rho = ",signif(intron_gene3$rho,3))
ex7e=ggplot() + 
  labs(x="Splicing efficiency", y = "spliceAI Score",title= "Correlation at gene level", fill="Occurrence") +
  scale_color_manual(values=c("#E64B35FF", "#3C5488FF", "#00A087FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,0.5))+
  geom_point_rast(data=intron_gene2, aes(x=eff_junction, y=spliceAI_junction, color=occurence), size=0.1, alpha=0.1)+
  geom_text(data=intron_gene3, aes(x=0.5, y=0.48, label=label),size=2, color="black")+
  facet_grid(cols=vars(occurence))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"ex7e.Correlation.splicing.efficiency.spliceAI.by_gene.pdf"), width = 3.3, height = 1.5)
print(ex7e)
dev.off()

#=====================
#Needed
#ex7f
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

ex7f=ggplot() + 
  labs(x=NULL, y = "spliceAI Score",title= "spliceAI score per gene region", fill="Occurrence") +
  scale_fill_manual(values=c("#E64B35FF", "#3C5488FF", "#00A087FF"), guide=NULL)+
  facet_grid(cols=vars(occurence))+
  coord_cartesian(ylim=c(0,0.55))+
  geom_boxplot(data=mintron_gene2[which(mintron_gene2$group1 == "spliceAI" & mintron_gene2$group2=="junction"),], mapping=aes(y=value, x=gene_group, fill=occurence), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, alpha=0.7, width = 0.7, outlier.shape = NA) + 
  geom_text(data=k5[which(k5$variable == "spliceAI_junction"),], mapping=aes(y=0.46, x=2.5, label=label),size=2.4)+
  geom_text(data=k4[which(k4$variable == "spliceAI_junction"),], mapping=aes(y=0.54, x=2, label=label),size=2.4)+
  annotate("segment", x = 2, xend = 3, y = 0.44, yend = 0.44, linewidth = 0.25) +
  annotate("segment", x = 1, xend = 3, y = 0.52, yend = 0.52, linewidth = 0.25) +  
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"ex7f.splicing.spliceAI.by_gene.RNAclass.pdf"), width = 2, height = 1.7)
print(ex7f)
dev.off()


#=====================
#Need 
#ex7g
both2=read.delim(paste0(path_fig5_data,"eRNA_structural_depleteion.raw.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

both3=both2%>%group_by(group1,group2, structual_depletion)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
both3$label=paste0(signif(both3$percent,3)*100,"%")
both3$label[which(both3$structual_depletion == "No")]=NA
both3$structual_depletion=factor(both3$structual_depletion, levels=c("Yes","No"))
both3$group1=factor(both3$group1, levels=c("mRNA","p_ncRNA","e_ncRNA"))
both3$group2=gsub("read","",both3$group2)
both3$group2=gsub("less","< ",both3$group2)
both3$group2[which(both3$group2=="3")]=">= 3"

ex7g=ggplot() + 
  labs(x="Read count", y = "Number of TES",title= "Non-poly(A) 3' end structural depletion", fill="depletion") +
  scale_fill_manual(values=c("black","white"))+
  facet_grid(cols=vars(group1))+
  coord_cartesian(ylim=c(0,30000))+
  geom_text(data=both3, mapping=aes(y=count+200, x=group2, label=label), size=1.8, color="black",hjust=0.5, vjust=0, position=position_stack())+
  geom_bar(data=both3, mapping=aes(y=count, x=group2, fill=structual_depletion), linewidth=0.25, width=0.7, color = "black", stat="identity") + 
  theme1+theme(legend.position = c(0.2,0.75))
pdf(paste0(path_fig5,"ex7g.eRNA_structural_depleteion.pdf"), width = 1.8, height = 1.7)
print(ex7g)
dev.off()

#=====================
#Needed
#ex7h

aRNAfold5=read.delim(paste0(path_fig5_data,"polyA_all_eRNACpGTATA_ATGC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
aRNAfold5=aRNAfold5[which(aRNAfold5$CpGTATA != "Others"),]
aRNAfold5$CpGTATA=factor(aRNAfold5$CpGTATA,levels=c("CGIap","CGInap","Null","TATA"))
aRNAfold5$group=factor(aRNAfold5$group, levels=c("A","U","G","C"))

ex7h=ggplot()+
  geom_line(data=aRNAfold5, mapping=aes(x=variable, y=value, color=group), linewidth=0.2)+
  facet_wrap(vars(CpGTATA), ncol=4, nrow=1)+
  labs(color=NULL, x="Stranded position from 3' end", y="Frequency", title="Non-poly(A) eRNA sequence composition")+
  scale_color_npg()+
  coord_cartesian(xlim=c(-100,100))+
  theme1+
  theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.97,0.89))
pdf(paste0(path_fig5,"ex7h.non_polyA_eRNA_CpGTATA_ATGC.pdf"), width = 4, height = 1.3)
print(ex7h)
dev.off()


#=====================
#Needed
#ex7i
TES_pcluster2=read.delim(paste0(path_fig5_data,"enhancer_downstream_promoter_distance_cCRE.gz"), header=T, stringsAsFactors = F, check.names = F)

TES_pcluster2$polyA = factor(TES_pcluster2$polyA, levels=c("Poly(A)","Non-poly(A)"))
TES_pcluster2$promoter_hinder=factor(TES_pcluster2$promoter_hinder, levels=c("Not Reach","Through","End"))
TES_pcluster2$label[which(TES_pcluster2$label == "CGInap_others")] = "CGInap"

ex7i=ggplot() + 
  labs(x=NULL, y = "% TES locus",title= "CGI eRNA terminataion relative\nto downstream  promoter", fill="TES location") +
  scale_y_continuous(labels = scales::percent, limits=c(0,1.05), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(polyA), scale="free_x", space="free_x")+
  geom_bar(data=TES_pcluster2[which(!is.na(TES_pcluster2$dpromoter)),], mapping=aes(y=percent, x=label, fill=promoter_hinder), linewidth=0.25, color = "black", stat="identity", width=0.85, alpha=0.7) + 
  scale_fill_manual(values=c("Not Reach"="white","End"="black","Through"="grey"), labels=c("Not Reach"="Not reaching","End"="Terminate","Through"="Read-through"))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig5,"ex7i.TESstop_upstream.pdf"), width = 2.5, height = 1.9)
print(ex7i)
dev.off()


#=====================
#Needed
#ex7j
data=read.delim(paste0(path_fig5_data,"TES_CpG_island_result.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
data=data[grep("eRNA_pN_",data$anno_region),]
data$group3=sapply(strsplit(data$anno_region,"_"),"[",3)
data$group3=factor(data$group3,levels=c("CGIap","CGInap","Null","TATA"))
data$signalID=factor(data$signalID,levels=c("n20","n40","n60","n80"))

ex7j=ggplot()+
  scale_color_npg()+
  geom_line(data = data, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_wrap(vars(group3),  nrow=1, ncol=4)+
  labs(color="CpG\nnumber", x="Stranded position from 3' end", y="% of CGI coverage", title="Distribution of CpG island near non-poly(A) eRNA TES")+
  geom_vline(xintercept=0, linewidth=0.25, linetype="dotted")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.93,0.87))
pdf(paste0(path_fig5,"ex7j.CpGisland_TES_non-polyA_eRNA.pdf"), width = 4, height = 1.3)
print(ex7j)
dev.off()

#===============================================================================
#ex7k
#MYC #subset to non-poly(A) eRNA alone
data=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_position.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)
data=data[grep("eRNA_pN_",data$anno_region),]
data$group3=sapply(strsplit(data$anno_region,"_"),"[",3)

data$group3=factor(data$group3,levels=c("CGIap","CGInap","Null","TATA"))
data$signalID=factor(data$signalID,levels=c("n50","n100","n150"))

ex7k=ggplot()+
  scale_color_npg()+
  geom_line(data = data, aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_wrap(vars(group3),  nrow=1, ncol=4)+
  labs(color="sig.\nlevel", x="Stranded position from 3' end", y="% of MYC binding peak coverage", title="Non-poly(A) eRNA MYC ChIP-seq")+
  geom_vline(xintercept=0, linewidth=0.25, linetype="dotted")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = c(0.93,0.87))
pdf(paste0(path_fig5,"ex7k.MYC_ChIP_TES_nonpolyA_eRNA.pdf"), width = 4, height = 1.3)
ex7k
dev.off()


#=====================
#Needed
#ex7l
myc_result=read.delim(paste0(path_fig5_data,"TES_ChIP_MYC_result.tsv.gz"), header=T, stringsAsFactors =F, check.names = F)

myc_result6=myc_result[which(!is.na(myc_result$exo_sensitivity_Tx)),]
myc_resulta5=myc_result6%>%group_by(TESrecur,ex5cluster_class,polyA,Exo_sensitive,MYC_TES0)%>%dplyr::summarise(count=n())

myc_result6$ex5cluster_class=factor(myc_result6$ex5cluster_class, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
myc_result6$polyA[which(myc_result6$polyA=="No")]="non-poly(A)"
myc_result6$polyA[which(myc_result6$polyA=="Yes")]="poly(A)"

myc_result9=myc_result6[which(myc_result6$ex5cluster_class != "CTCF_ncRNA"),]%>%group_by(ex5cluster_class,polyA)%>%dplyr::summarise(p=wilcox.test(exo_sensitivity_Tx ~MYC_TES0, alternative = "two.sided")$p.value, count=n())
myc_result9$label="***"
myc_result9$label[which(myc_result9$p>=0.001)]="**"
myc_result9$label[which(myc_result9$p>=0.01)]="*"
myc_result9$label[which(myc_result9$p>=0.05)]="n.s."

ex7l=ggplot() + 
  labs(x="MYC binding", y = "Exosome sensitivity",title= "MYC and exosome sensitivity") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(ex5cluster_class,polyA), scale="free_x")+
  scale_y_continuous(limits=c(0,1.15), breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=myc_result6[which(myc_result6$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA")),], mapping=aes(y=exo_sensitivity_Tx, x=MYC_TES0, fill=polyA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.5, outlier.shape = NA) + 
  geom_text(data=myc_result9[which(myc_result9$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA")),], mapping=aes(y=1.07, x=1.5, label=label),size=2.4)+
  annotate("segment", x=1,xend=2,y=1.03,yend=1.03, linewidth=0.25)+
  theme1
pdf(paste0(path_fig5,"ex7l.RNA.TES.polyA.MYC.exosome.pdf"), width = 2.5, height = 1.4)
print(ex7l)
dev.off()

#===============================================================================

