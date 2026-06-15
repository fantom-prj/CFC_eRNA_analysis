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
path_fig4=paste0(primary_folder,"fig4/out/")
path_fig4_data=paste0(primary_folder,"fig4/data/")

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
#f4a
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data1$ex5cluster_class=factor(data1$ex5cluster_class,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA", "other_ncRNA"))
length_cluster=data1%>%group_by(ex5cluster_class)%>%dplyr::summarise(length=median(read_median_length),exon=median(median_exon))
data1$ex5cluster_class=factor(data1$ex5cluster_class,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA","other_ncRNA"))

mdata2=reshape2::melt(data1[,c(1, 9,2,7)], id=c(1,2))
mdata1med=mdata2%>%group_by(variable, ex5cluster_class)%>%dplyr::summarise(median=median(value))

f4a=ggplot() + 
  labs(x = "Length (nt)",title= "Transcript length from Ex5_clusters", color=NULL) +
  scale_linetype_manual(label=c("Transcript length","Genomic range"), values=c("solid","dashed"), name=NULL)+
  scale_color_manual(values=c("grey", "#E64B35FF", "#3C5488FF", "#00A087FF"))+
  geom_density(data=mdata2[which(mdata2$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA", "other_ncRNA")),], mapping=aes(x=value, color=ex5cluster_class, linetype=variable), linewidth=0.25)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  theme1+theme(legend.position=c(0.75,0.7))
pdf(paste0(path_fig4,"f4a.end5_cluster.base.length2.pdf"), width = 1.7, height = 1.7)
print(f4a)
dev.off()

#=====================
#Needed
#f4b
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1$ex5cluster_class=factor(data1$ex5cluster_class,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA","other_ncRNA"))

exonmed=data1%>%group_by(ex5cluster_class)%>%dplyr::summarise(med=median(median_exon),atleast1=length(which(median_exon>1))/n())
exonmed$label=paste0(signif(exonmed$atleast1,3)*100,"%")
f4b=ggplot(data1[which(data1$ex5cluster_class != "CTCF_ncRNA"),], aes(x=median_exon, color=factor(ex5cluster_class))) + 
  labs(y="Ranked transcripts", x ="Number of exon",title= "Exon number of Ex5_clusters", color=NULL) +
  scale_color_manual(values=c("grey", "#E64B35FF", "#3C5488FF","#00A087FF", "#7E6148FF"))+
  coord_flip(xlim=c(0,15))+
  stat_ecdf(geom = "step", linewidth=0.25)+
  scale_y_continuous(labels = scales::percent)+
  geom_text_repel(data=exonmed[which(exonmed$ex5cluster_class != "CTCF_ncRNA"),], aes(x=1.5, y=1-atleast1, label=label), segment.size = 0.25, size=1.8, force =5, nudge_x=4, nudge_y=0.01)+
  theme1+theme(legend.position=c(0.15,0.7))
pdf(paste0(path_fig4,"f4b.n5cluster.exon.number.pdf"), width = 1.7, height = 1.7)
print(f4b)
dev.off()

#=====================
#Needed
#f4c
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data2=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","orientation","read_median_length","median_exon","median_range")]
data2=data2[which(data2$orientation != "Others"),]
data3=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","any_CpG_island","read_median_length","median_exon","median_range")]
data4=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","downstream_CpG_island","read_median_length","median_exon","median_range")]
data5=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","TATA_box","read_median_length","median_exon","median_range")]
#data6=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","conserve_up","read_median_length","median_exon","median_range")]
#data6=data6[which(data6$conserve_up != "Others"),]
data7=data1[which(data1$ex5cluster_class=="e_ncRNA"), c("ex5cluster_class","n5_string","SE_all","read_median_length","median_exon","median_range")]
data7=data7[which(data7$SE_all != "Others"),]

colnames(data2)[3]="feature"
colnames(data3)[3]="feature"
colnames(data4)[3]="feature"
colnames(data5)[3]="feature"
#colnames(data6)[3]="feature"
colnames(data7)[3]="feature"

data2$group2="1D / 2D"
data3$group2="all_CGI"
data4$group2="all_dCGI"
data5$group2="all_TATA"
#data6$group2="Conserved"
data7$group2="SE / TE"

data2=rbind(data2,data3,data4,data5,data7)
data2$group2=factor(data2$group2, levels=c("all_CGI","all_dCGI","all_TATA","1D / 2D","SE / TE"))
data2$feature=factor(data2$feature, levels=c("1D","2D","No","Yes","SE","TE"))
data4=data2%>%group_by(ex5cluster_class, group2, feature)%>%dplyr::summarise(median=median(read_median_length),count=n())
data1wilcox3=data2%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(read_median_length ~feature, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f4c=ggplot() + 
  labs(x=NULL, y = "Transcript length (bp)",title= "e_ncRNA length affected by enhancer properties") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(group2), scale="free_x")+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=data2[which(data2$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=read_median_length, x=feature, fill=group2),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=data2[which(data2$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=read_median_length, x=feature, fill=group2), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data4[which(data4$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=median, x=feature, label=signif(median,3)), size=2, angle=90, vjust=-0.5, hjust=0.5, color="black")+
  geom_text(data=data1wilcox3[which(data1wilcox3$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=2300, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 2250, yend = 2250, linewidth = 0.25)+
  theme1
pdf(paste0(path_fig4,"f4c.end5_cluster.bidirectional.SE.eRNA.length.pdf"), width = 2.6, height = 1.7)
print(f4c)
dev.off()

#=====================
#Needed
#f4d

#continue from last section
data6=data2%>%group_by(ex5cluster_class,group2, feature)%>%dplyr::summarise(median=median(median_exon), percent=length(which(median_exon >1))/n())
data6$label="> 1 exon: "

data1wilcox5=data2%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(median_exon ~feature, alternative = "two.sided")$p.value)
data1wilcox5$label="***"
data1wilcox5$label[which(data1wilcox5$p>=0.001)]="**"
data1wilcox5$label[which(data1wilcox5$p>=0.01)]="*"
data1wilcox5$label[which(data1wilcox5$p>=0.05)]="n.s."

f4d=ggplot() + 
  labs(x=NULL, y = "Number of exon",title= "e_ncRNA exon number affected by enhancer properties") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(group2), scales="free_x")+
  coord_cartesian(ylim=c(0,9.5))+
  geom_violin(data=data2[which(data2$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=median_exon, x=feature, fill=group2, alpha=feature), linewidth=0.1, color = "black", bounds=c(0,8), width = 0.9, position = position_nudge(x = +0), alpha=1) + 
  geom_text(data=data6[which(data6$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=2, x=feature, label=label), angle=90, size=2, vjust=(-0.5), hjust=0, color="grey")+
  geom_text(data=data6[which(data6$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=5.3, x=feature, label=paste0(signif(percent*100,2),"%")), angle=90, size=2, vjust=(-0.5), hjust=0, color="black")+
  geom_text(data=data1wilcox5[which(data1wilcox5$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=9.4, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 9.2, yend = 9.2, linewidth = 0.25)+
  theme1
pdf(paste0(path_fig4,"f4d.end5_cluster.bidirectional.CGI.TATA.SE.eRNA.exon.number.pdf"), width = 2.6, height = 1.7)
print(f4d)
dev.off()

#=====================
#Needed
#f4e
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data2a=data1[which(data1$CpGTATA %in% c("CGIap","Null") & data1$ex5cluster_class == "e_ncRNA"),c("CpGTATA","read_median_length","median_exon")]
data2b=data1[which(data1$CpGTATA %in% c("CGInap","Null") & data1$ex5cluster_class == "e_ncRNA"),c("CpGTATA","read_median_length","median_exon")]
data2c=data1[which(data1$CpGTATA %in% c("TATA","Null") & data1$ex5cluster_class == "e_ncRNA"),c("CpGTATA","read_median_length","median_exon")]

data2a$group2="CGIap"
data2b$group2="CGInap"
data2c$group2="TATA"

data2=rbind(data2a,data2b,data2c)

data3=data1[which(data1$CpGTATA %in% c("CGIap","CGInap","Null","TATA") & data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(CpGTATA)%>%dplyr::summarise(median=median(read_median_length, na.rm=T),count=n())
data4=data2%>%group_by(group2)%>%dplyr::summarise(p=wilcox.test(read_median_length ~CpGTATA, alternative = "two.sided")$p.value)
data4$label="***"
data4$label[which(data4$p>=0.001)]="**"
data4$label[which(data4$p>=0.01)]="*"
data4$label[which(data4$p>=0.05)]="n.s."

f4e=ggplot() + 
  labs(x=NULL, y = "Transcript length (bp)",title= "e_ncRNA length from\nregulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=data1[which(data1$CpGTATA %in% c("CGIap","CGInap","Null","TATA") & data1$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=read_median_length, x=CpGTATA, fill=CpGTATA),adjust = .75, width = .5, .width = 0, justification = -.4, point_colour = NA, alpha=0.4)+
  geom_boxplot(data=data1[which(data1$CpGTATA %in% c("CGIap","CGInap","Null","TATA") & data1$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=read_median_length, x=CpGTATA, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data3, mapping=aes(y=median, x=CpGTATA, label=signif(median,3)), angle=90, size=2, vjust=(-0.6), hjust=0.5, color="black")+
  geom_text(data=data4[which(data4$group2 == "CGIap"),], mapping=aes(y=2225, x=2, label=label),size=2.4)+
  geom_text(data=data4[which(data4$group2 == "CGInap"),], mapping=aes(y=2025, x=2.5, label=label),size=2.4)+
  geom_text(data=data4[which(data4$group2 == "TATA"),], mapping=aes(y=2425, x=3.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 3, y = 2200, yend = 2200, linewidth = 0.25)+
  annotate("segment", x = 2, xend = 3, y = 2000, yend = 2000, linewidth = 0.25)+
  annotate("segment", x = 3, xend = 4, y = 2400, yend = 2400, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4e.end5_cluster.e_ncRNA_length_CpG_TATA.pdf"), width = 1.3, height = 1.7)
print(f4e)
dev.off() 

#=====================
#Needed
#f4f

#continue from last section
data2=rbind(data2a,data2b,data2c)

edata3=data1[which(data1$CpGTATA %in% c("CGIap","CGInap","Null","TATA") & data1$ex5cluster_class == "e_ncRNA"),]%>%group_by(CpGTATA)%>%dplyr::summarise(mean=mean(median_exon), percent=length(which(median_exon !=1))/n() ,count=n())
edata3$label="> 1 exon: "
edata4=data2%>%group_by(group2)%>%dplyr::summarise(p=wilcox.test(median_exon ~CpGTATA, alternative = "two.sided")$p.value)
edata4$label="***"
edata4$label[which(edata4$p>=0.001)]="**"
edata4$label[which(edata4$p>=0.01)]="*"
edata4$label[which(edata4$p>=0.05)]="n.s."

data1$CpGTATA=factor(data1$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))

f4f=ggplot() + 
  labs(x=NULL, y = "Number of exon",title= "e_ncRNA splicing from\nregulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,10))+
  geom_violin(data=data1[which(data1$CpGTATA %in% c("CGIap","CGInap","Null","TATA") & data1$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=median_exon, x=CpGTATA, fill=CpGTATA, alpha=orientation), linewidth=0.15, color = "black", bounds=c(0,8), width = 0.9, position = position_nudge(x = +0), alpha=1) + 
  geom_text(data=edata3, mapping=aes(y=2, x=CpGTATA, label=label), angle=90, size=2, vjust=(-0.6), hjust=0, color="grey")+
  geom_text(data=edata3, mapping=aes(y=5.6, x=CpGTATA, label=paste0(signif(percent*100,2),"%")), angle=90, size=2, vjust=(-0.6), hjust=0, color="black")+
  geom_point(data=edata3, mapping=aes(y=mean, x=CpGTATA), linewidth=0.2, color="red")+
  geom_text(data=edata4[which(edata4$group2 == "CGIap"),], mapping=aes(y=8.7, x=2, label=label),size=2.4)+
  geom_text(data=edata4[which(edata4$group2 == "CGInap"),], mapping=aes(y=8, x=2.5, label=label),size=2.4)+
  geom_text(data=edata4[which(edata4$group2 == "TATA"),], mapping=aes(y=9.4, x=3.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 3, y = 8.5, yend = 8.5, linewidth = 0.25)+
  annotate("segment", x = 2, xend = 3, y = 7.8, yend = 7.8, linewidth = 0.25)+
  annotate("segment", x = 3, xend = 4, y = 9.2, yend = 9.2, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4f.end5_cluster.e_ncRNA_exon_CpG_TATA.pdf"), width = 1.3, height = 1.7)
print(f4f)
dev.off() 

#=====================
##Needed
##f4g

data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2=data1[which(data1$ex5cluster_class=="e_ncRNA" & data1$CpGTATA %in% c("CGIap","CGInap")), c("CpGTATA","n5_string","downstream_CpG_island","read_median_length","median_exon","median_range")]
data2$downstream_CpG_island=gsub("Yes","downstream",data2$downstream_CpG_island)
data2$downstream_CpG_island=gsub("No","others",data2$downstream_CpG_island)
data2$downstream_CpG_island=factor(data2$downstream_CpG_island, levels=c("downstream","others"))

data4=data2%>%group_by(CpGTATA, downstream_CpG_island)%>%dplyr::summarise(median=median(read_median_length),count=n())
data1wilcox3=data2%>%group_by(CpGTATA)%>%dplyr::summarise(p=wilcox.test(read_median_length ~downstream_CpG_island, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f4g=ggplot() + 
  labs(x=NULL, y = "Transcript length (bp)",title= "RNA length \nfrom CGI enhancer") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF"),guide=NULL)+
  facet_grid(cols=vars(CpGTATA), scale="free_x")+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=data2, mapping=aes(y=read_median_length, x=downstream_CpG_island, fill=CpGTATA),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=data2, mapping=aes(y=read_median_length, x=downstream_CpG_island, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data4, mapping=aes(y=median, x=downstream_CpG_island, label=signif(median,3)), size=2, angle=90, vjust=-0.5, hjust=0.5, color="black")+
  geom_text(data=data1wilcox3, mapping=aes(y=2300, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 2250, yend = 2250, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4g.end5_cluster.e_ncRNA_length_CpG.pdf"), width = 1.5, height = 1.7)
print(f4g)
dev.off()

#=====================
#Needed
#f4h

#continue from last section
data4=data2%>%group_by(CpGTATA, downstream_CpG_island)%>%dplyr::summarise(mean=mean(median_exon), percent=length(which(median_exon !=1))/n() ,count=n())
data4$label="> 1 exon: "
data1wilcox3=data2%>%group_by(CpGTATA)%>%dplyr::summarise(p=wilcox.test(median_exon ~downstream_CpG_island, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f4h=ggplot() + 
  labs(x=NULL, y = "Number of exon",title= "RNA splicing \nfrom CGI enhancer") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF"),guide=NULL)+
  facet_grid(cols=vars(CpGTATA), scales="free_x")+
  coord_cartesian(ylim=c(0,9.5))+
  geom_violin(data=data2, mapping=aes(y=median_exon, x=downstream_CpG_island, fill=CpGTATA), linewidth=0.1, color = "black", bounds=c(0,8), width = 0.9, position = position_nudge(x = +0), alpha=1) + 
  geom_text(data=data4, mapping=aes(y=2, x=downstream_CpG_island, label=label), angle=90, size=2, vjust=(-0.5), hjust=0, color="grey")+
  geom_text(data=data4, mapping=aes(y=6.1, x=downstream_CpG_island, label=paste0(signif(percent*100,2),"%")), angle=90, size=2, vjust=(-0.5), hjust=0, color="black")+
  geom_text(data=data1wilcox3, mapping=aes(y=9.4, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 9.2, yend = 9.2, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4h.end5_cluster.e_ncRNA_exon_CpG.pdf"), width = 1.5, height = 1.7)
print(f4h)
dev.off()

#=====================
#Needed
#f4i
d2=read.delim(paste0(path_fig4_data,"FE_CGI_dCGI_TATA.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

d2$group=factor(d2$group, levels=c("conserve_down","conserve_up","ubiquitous","SE","2D"))
d2$ex5cluster_class=factor(d2$ex5cluster_class, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
d2$enrichment=factor(d2$enrichment, levels=c("CGI","CGIap","CGInap","Null","TATA"))
d2$sig_level=factor(d2$sig_level, levels=c("ns","*","**","***"))
f4i=ggplot() + 
  labs(x=NULL, y = NULL,title= "Enrichment in ex5_cluster") +
  facet_grid(cols=vars(ex5cluster_class), scale="free", space ="free")+
  geom_point(data=d2[which(d2$ex5cluster_class %in% c("p_ncRNA","e_ncRNA")),], mapping=aes(y=group, x=enrichment, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.2,0.8,1.4,2))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4i.FE_CGI_dCGI_TATA.pdf"), width = 2.1, height = 1.7)
print(f4i)
dev.off() 

#=====================
#Needed
#f4j
data2=read.delim(paste0(path_fig4_data,"transcript.splice.and.polyA.RNA.length.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data2$group3=paste0(data2$group,"\n",data2$group2)
data2$group3=factor(data2$group3, levels=c("e_ncRNA\nSpliced","e_ncRNA\npoly(A)","p_ncRNA\nSpliced","p_ncRNA\npoly(A)"))
data2$value=factor(data2$value, levels=c("Yes","No"))
data4=data2%>%group_by(group3, value)%>%dplyr::summarise(median=median(read_median_length),count=n())
data1wilcox3=data2%>%group_by(group3)%>%dplyr::summarise(p=wilcox.test(read_median_length ~value, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f4j=ggplot() + 
  labs(x=NULL, y = "Length (bp)",title= "ncRNA transcript body features") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(group3), scale="free_x")+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=data2, mapping=aes(y=read_median_length, x=value, fill=group3),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=data2, mapping=aes(y=read_median_length, x=value, fill=group3), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data4, mapping=aes(y=median, x=value, label=signif(median,3)), size=2, angle=90, vjust=-0.6, hjust=0.5, color="black")+
  geom_text(data=data1wilcox3, mapping=aes(y=2400, x=1.5, label=label),size=2.4)+
  annotate("segment", x=1,xend=2,y=2300,yend=2300, linewidth=0.25)+
  theme1
pdf(paste0(path_fig4,"f4j.transcript.splice.and.polyA.RNA.length.pdf"), width = 2.2, height = 1.7)
print(f4j)
dev.off()

#=====================
#Needed
#f4k
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data3=data1[which(data1$CpGTATA != "Others" & data1$ex5cluster_class %in% c("p_ncRNA","e_ncRNA")),]
data3$CpGTATA=factor(data3$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))
data3$ex5cluster_class=factor(data3$ex5cluster_class, levels=c("p_ncRNA","e_ncRNA"))
q2=data3%>%group_by(ex5cluster_class, CpGTATA)%>%dplyr::summarise(median=median(polyArate), mean=mean(polyArate), percent0=length(which(polyArate==0))/n(), percent1=length(which(polyArate==1))/n(), percent_chi=length(which(polyArate>0 & polyArate<1))/n(), count=n())

#p-value
data4a=data3[which(data3$ex5cluster_class == "p_ncRNA" & data3$CpGTATA %in% c("CGI","Null")),c("ex5cluster_class","CpGTATA","polyArate")]
data4b=data3[which(data3$ex5cluster_class == "e_ncRNA" & data3$CpGTATA %in% c("CGIap","Null")),c("ex5cluster_class","CpGTATA","polyArate")]
data4c=data3[which(data3$ex5cluster_class == "e_ncRNA" & data3$CpGTATA %in% c("CGInap","Null")),c("ex5cluster_class","CpGTATA","polyArate")]
data4d=data3[which(data3$CpGTATA %in% c("TATA","Null")),c("ex5cluster_class","CpGTATA","polyArate")]

data4a$group2="CGI"
data4b$group2="CGIap"
data4c$group2="CGInap"
data4d$group2="TATA"
data4=rbind(data4a,data4b,data4c,data4d)

data5=data4%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(polyArate ~CpGTATA, alternative = "two.sided")$p.value)
data5$label="***"
data5$label[which(data5$p>=0.001)]="**"
data5$label[which(data5$p>=0.01)]="*"
data5$label[which(data5$p>=0.05)]="n.s."

f4k=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF","#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL , y ="Weighted frequency of poly(A) transcript" , title ="Rate of Poly(A) RNA",  fill=NULL)+
  scale_y_continuous(limits=c(0,1.5), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(ex5cluster_class), scales="free_y")+
  geom_boxplot(data=data3, mapping=aes(x=CpGTATA, y=polyArate, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA, alpha=0.8) + 
  geom_point(data=q2, mapping=aes(y=mean, x=CpGTATA), size=0.25, color="red")+
  geom_text(data=data5[which(data5$group2 == "CGI"),], mapping=aes(y=1.25, x=2.5, label=label),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGIap"),], mapping=aes(y=1.25, x=3, label=label),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "CGInap"),], mapping=aes(y=1.1, x=3.5, label=label),size=2.4)+
  geom_text(data=data5[which(data5$group2 == "TATA"),], mapping=aes(y=1.4, x=4.5, label=label),size=2.4)+
  annotate("segment", x=2,xend=4,y=1.2,yend=1.2, linewidth=0.25)+
  annotate("segment", x=3,xend=4,y=1.05,yend=1.05, linewidth=0.25)+
  annotate("segment", x=4,xend=5,y=1.35,yend=1.35, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4k.polyA_CGI_TATA.pdf"), width = 1.8, height = 1.7)
print(f4k)
dev.off() 

#=====================
#Needed
#f4l
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data1a=data1%>%group_by(ex5cluster_class,CpGTATA)%>%dplyr::summarise(count_chimeric=length(which(polyArate>=0.1 & polyArate<=0.9)),count_nonPolyA=length(which(polyArate<0.1)),count_total=n())
data1a$chimeric=data1a$count_chimeric/data1a$count_total
data1a$non_polyA=data1a$count_nonPolyA/data1a$count_total
data1a$polyA=1-data1a$non_polyA-data1a$chimeric

data1a=data1a[which(data1a$ex5cluster_class %in% c("p_ncRNA","e_ncRNA")),]
data1a=data1a[which(data1a$CpGTATA != "Others"),]
data1a$ex5cluster_class=factor(data1a$ex5cluster_class, levels=c("p_ncRNA","e_ncRNA"))
data1a$CpGTATA=factor(data1a$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))
data2a=reshape2::melt(data1a, id=c(1:5))
data2a$variable=factor(data2a$variable, levels=c("non_polyA","chimeric","polyA"))

f4l=ggplot() + 
  scale_fill_manual(values=c("#4DBBD5FF", "darkorchid3", "#E64B35FF"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Poly(A) from ex5_cluster",  fill=NULL)+
  facet_grid(cols=vars(ex5cluster_class), scale="free", space="free")+
  geom_col(data=data2a, mapping=aes(x=CpGTATA, y=value, fill=variable), linewidth=0.25, color = "black", width = 0.7, position = position_stack(), alpha=0.8) + 
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4l.percent.polyA_enhancer_promoter_CGI_TATA.pdf"), width = 2, height = 1.7)
print(f4l)
dev.off() 

#=====================
#Needed
#f4m
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data1a=data1[which(!is.na(data1$exo_sensitivity_Tx)),]%>%group_by(ex5cluster_class)%>%dplyr::summarise(exo_sensitivity_Tx=median(exo_sensitivity_Tx,na.rm=T), count=n(), n_read=sum(n_read))
data1$ex5cluster_class=factor(data1$ex5cluster_class, levels=c("mRNA","p_ncRNA","e_ncRNA","CTCF_ncRNA","other_ncRNA"))

f4m=ggplot() + 
  labs(x=NULL, y = "Exosome sensitivity",title= "Exosome sensitivity of\ntranscript from Ex5_cluster") +
  scale_fill_npg(guide=NULL)+
  scale_y_continuous(limits=c(0,1.15),breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=data1, mapping=aes(y=exo_sensitivity_Tx, x=ex5cluster_class, fill=ex5cluster_class), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=data1a, mapping=aes(y=exo_sensitivity_Tx, x=ex5cluster_class, label=signif(exo_sensitivity_Tx,3)), size=2, angle=90, vjust=-0.7, hjust=0.1, color="black")+
  geom_text(data=data1a, mapping=aes(y=1, x=ex5cluster_class, label=paste0("n=", count)), size=1.8, angle=25, vjust=0, hjust=0, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4m.n5_cluster_all_exo_sensitivity.pdf"), width = 1.3, height = 1.7)
print(f4m)
dev.off()

#=====================
#Needed
#f4n
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

data2=data1[which(data1$ex5cluster_class %in% c("p_ncRNA","e_ncRNA") & !is.na(data1$exo_sensitivity_Tx)),]
data2=data2[,c("ex5cluster_class", "n5_string", "exo_sensitivity_Tx", "any_CpG_island", "TATA_box", "orientation", "SE_all", "conserve_down", "ubiquitous")]
colnames(data2)[c(4:9)]=c("all_CGI","all_TATA","1D / 2D","SE / TE","Conserv.","Ubiquit.")
z1=reshape2::melt(data2,id=c(1:3))
z1=z1[which(z1$value %in% c("Yes", "No", "2D", "1D","SE","TE")),]

z1a=z1%>%group_by(ex5cluster_class, variable, value)%>%dplyr::summarise(median=median(exo_sensitivity_Tx))
z1b=z1%>%group_by(ex5cluster_class, variable)%>%dplyr::summarise(p=wilcox.test(exo_sensitivity_Tx ~value, alternative = "two.sided")$p.value)
z1b$label="***"
z1b$label[which(z1b$p>=0.001)]="**"
z1b$label[which(z1b$p>=0.01)]="*"
z1b$label[which(z1b$p>=0.05)]="n.s."

f4n=ggplot() + 
  labs(x=NULL, y = "Exosome sensitivity",title= "RNA exosome sensitivity against enhancer features") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(variable), scale="free_x")+
  coord_cartesian(ylim=c(0,1.2))+
  geom_boxplot(data=z1[which(z1$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=exo_sensitivity_Tx, x=value, fill=variable), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=z1a[which(z1a$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=median, x=value, label=signif(median,3)), size=2, angle=90, vjust=-0.7, hjust=0.5, color="black")+
  geom_text(data=z1b[which(z1b$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=1.1, x=1.5, label=label),size=2.4)+
  annotate("segment", x=1,xend=2,y=1.05,yend=1.05, linewidth=0.25)+
  theme1
pdf(paste0(path_fig4,"f4n.n5_cluster_enhancer_feature_exo_sensitivity.pdf"), width = 3, height = 1.7)
print(f4n)
dev.off()

#=====================
#Needed
#f4o
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

data2=data1[which(data1$ex5cluster_class %in% c("p_ncRNA","e_ncRNA") & !is.na(data1$exo_sensitivity_Tx)),]
data2=data2[,c("ex5cluster_class", "n5_string", "CpGTATA","exo_sensitivity_Tx")]
data2=data2[which(data2$CpGTATA != "Others"),]

data3=data2%>%group_by(ex5cluster_class,CpGTATA)%>%dplyr::summarise(median=median(exo_sensitivity_Tx, na.rm=T),count=n())

data2a=data2[which(data2$CpGTATA %in% c("CGIap","Null") & data2$ex5cluster_class == "e_ncRNA"),]
data2b=data2[which(data2$CpGTATA %in% c("CGInap","Null") & data2$ex5cluster_class == "e_ncRNA"),]
data2c=data2[which(data2$CpGTATA %in% c("TATA","Null")),]
data2d=data2[which(data2$CpGTATA %in% c("CGI","Null") & data2$ex5cluster_class == "p_ncRNA"),]
data2a$group2="CGIap"
data2b$group2="CGInap"
data2c$group2="TATA"
data2d$group2="CGI"
data4=rbind(data2a,data2b,data2c,data2d)

data4=data4%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(exo_sensitivity_Tx ~CpGTATA, alternative = "two.sided")$p.value)
data4$label="***"
data4$label[which(data4$p>=0.001)]="**"
data4$label[which(data4$p>=0.01)]="*"
data4$label[which(data4$p>=0.05)]="n.s."
data2$CpGTATA=factor(data2$CpGTATA, levels=c("CGI","CGIap","CGInap","Null","TATA"))

f4o=ggplot() + 
  labs(x=NULL, y = "Exosome sensitivity",title= "RNA exosome sensitivity\nacross regulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  scale_y_continuous(limits=c(0,1.3), breaks=c(0,0.25,0.5,0.75,1))+
  geom_boxplot(data=data2[which(data2$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=exo_sensitivity_Tx, x=CpGTATA, fill=CpGTATA), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=data3[which(data3$ex5cluster_class == "e_ncRNA"),], mapping=aes(y=median, x=CpGTATA, label=signif(median,3)), angle=90, size=2, vjust=(-0.8), hjust=0.5, color="black")+
  geom_text(data=data4[which(data4$ex5cluster_class == "e_ncRNA" & data4$group2 == "CGIap"),], mapping=aes(y=1.13, x=2, label=label),size=2.4)+
  geom_text(data=data4[which(data4$ex5cluster_class == "e_ncRNA" &data4$group2 == "CGInap"),], mapping=aes(y=1.03, x=2.5, label=label),size=2.4)+
  geom_text(data=data4[which(data4$ex5cluster_class == "e_ncRNA" &data4$group2 == "TATA"),], mapping=aes(y=1.23, x=3.5, label=label),size=2.4)+
  annotate("segment",x=1,xend=3,y=1.12,yend=1.12, linewidth=0.25)+
  annotate("segment",x=2,xend=3,y=1.02,yend=1.02, linewidth=0.25)+
  annotate("segment",x=3,xend=4,y=1.22,yend=1.22, linewidth=0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4o.end5_cluster.e_ncRNA_exoSensitivity_CpG_TATA.pdf"), width = 1.3, height = 1.7)
print(f4o)
dev.off() 

#=====================
#Needed
#f4p
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

data2=data1[which(data1$ex5cluster_class == "e_ncRNA" & !is.na(data1$exo_sensitivity_Tx)),]
data2=data2[which(data2$CpGTATA %in% c("CGIap","CGInap")), c("CpGTATA","n5_string","downstream_CpG_island","exo_sensitivity_Tx")]
data2$downstream_CpG_island=gsub("Yes","downstream",data2$downstream_CpG_island)
data2$downstream_CpG_island=gsub("No","others",data2$downstream_CpG_island)
data2$downstream_CpG_island=factor(data2$downstream_CpG_island, levels=c("downstream","others"))

data4=data2%>%group_by(CpGTATA, downstream_CpG_island)%>%dplyr::summarise(median=median(exo_sensitivity_Tx, na.rm=T),count=n())
data1wilcox3=data2%>%group_by(CpGTATA)%>%dplyr::summarise(p=wilcox.test(exo_sensitivity_Tx, na.rm=T ~downstream_CpG_island, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

f4p=ggplot() + 
  labs(x=NULL, y = "Exosome sensitivity",title= "RNA exosome sensitivity\nfrom CGI enhancer") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF"),guide=NULL)+
  facet_grid(cols=vars(CpGTATA), scale="free_x")+
  scale_y_continuous(limits=c(0,1.2), breaks=c(0,0.25,0.5,0.75,1))+
  #ggdist::stat_halfeye(data=data2, mapping=aes(y=exo_sensitivity_Tx, x=downstream_CpG_island, fill=CpGTATA),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=data2, mapping=aes(y=exo_sensitivity_Tx, x=downstream_CpG_island, fill=CpGTATA), linewidth=0.2, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=data4, mapping=aes(y=median, x=downstream_CpG_island, label=signif(median,3)), size=2, angle=90, vjust=-0.8, hjust=0.5, color="black")+
  geom_text(data=data1wilcox3, mapping=aes(y=1.03, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 1.02, yend = 1.02, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"f4p.end5_cluster.e_ncRNA_exoSensitivity_CpG.pdf"), width = 1.5, height = 1.7)
print(f4p)
dev.off()


