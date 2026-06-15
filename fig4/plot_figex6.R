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
#ex6a
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

length_cluster=data1%>%group_by(ex5cluster_class)%>%dplyr::summarise(length=median(read_median_length),exon=median(median_exon))
data1$ex5cluster_class=factor(data1$ex5cluster_class,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA","other_ncRNA"))

mdata2=reshape2::melt(data1[,c(1, 9,2,7)], id=c(1,2))
mdata1med=mdata2%>%group_by(variable, ex5cluster_class)%>%dplyr::summarise(median=median(value))

ex6a=ggplot() +
  labs(y = "Length (nt)",title= "Transcript length\nfrom Ex5_clusters", x=NULL) +
  scale_fill_manual(values=c("grey", "#E64B35FF", "#3C5488FF", "#7E6148FF","#00A087FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,4000))+
  geom_boxplot(data=mdata2[which(mdata2$variable == "read_median_length"),], mapping=aes(x=ex5cluster_class, y=value, fill=ex5cluster_class),color="black", linewidth=0.25, linetype="solid", width=0.3, alpha=1,outlier.shape = NA)+
  geom_text(data=mdata1med[which(mdata1med$variable == "read_median_length"),], mapping=aes(x=ex5cluster_class, y=median,label=signif(median,3)), size=2, angle=90, hjust=0.5, vjust=(-0.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6a.end5_cluster.base.length3.pdf"), width = 1.4, height = 1.7)
print(ex6a)
dev.off()

#=====================
#Needed
#ex6b-d
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

data1$Andersson_robust[which(data1$Andersson_robust==0)]="others"
data1$Andersson_robust[which(data1$Andersson_robust==1)]="bidirectional"

data2=data1[grep("ncRNA",data1$ex5cluster_class),]

e1=data2[which(data2$ex5cluster_class=="e_ncRNA"),c(1,9,2,6)]
e2=data2[which(data2$proximity=="distal"),c(1,9,2,6)]
e3=data2[which(data2$proximity=="distal" & data2$ATAC == "withATAC"),c(1,9,2,6)]
#e4=data2[which(data2$proximity=="distal" & data2$orientation == "2D"),c(1,9,2,6)]
e5=data2[which(data2$proximity=="distal" & data2$Andersson_robust == "bidirectional"),c(1,9,2,6)]
e6=data2[which(data2$chromatin_state_promoter_type=="enhancer-like"),c(1,9,2,6)]

e1$group2="This_study"
e2$group2="Distal_tCRE"
e3$group2="Distal_aCRE"
#e4$group2="Distal_2D"
e5$group2="Distal_2D"
e6$group2="ChromHMM"
ee=rbind(e1,e2,e3,e5,e6)
ee$group2=factor(ee$group2, levels=c("This_study","Distal_tCRE","Distal_aCRE","Distal_2D","ChromHMM"))
ee1=ee%>%group_by(group2)%>%dplyr::summarise(median=median(read_median_length),count=n())
ee1$label=paste0("n=\n",ee1$count)

ex6b=ggplot() + 
  labs(x=NULL, y = "Length (bp)",title= "ncRNA from different enhancer typing") +
  scale_fill_npg(guide=NULL)+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=ee, mapping=aes(y=read_median_length, x=group2, fill=group2),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=ee, mapping=aes(y=read_median_length, x=group2, fill=group2), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=ee1, mapping=aes(y=median, x=group2, label=signif(median,3)), size=2, angle=90, vjust=-1.1, hjust=0.5, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6b.end5_cluster.eRNA.definition.length.pdf"), width = 2.1, height = 1.7)
print(ex6b)
dev.off()

#=
ee2=ee%>%group_by(group2,ex5cluster_class)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
ee2$label=paste0(signif(ee2$percent*100,2),"%")
ee2$label[which(ee2$ex5cluster_class != "e_ncRNA")]=NA
ee2$ex5cluster_class=factor(ee2$ex5cluster_class, levels=c("e_ncRNA","p_ncRNA","CTCF_ncRNA","other_ncRNA"))

ex6c=ggplot() + 
  labs(x=NULL, y = "Number of transcript",title= "ncRNA from different enhancer typing", fill=NULL) +
  scale_fill_npg()+
  #scale_y_continuous(labels = scales::percent, limits=c(0,1.2), breaks=c(0,0.25,0.5,0.75,1))+
  geom_bar(data=ee2[which(ee2$group2 != "Distal_2D"),], mapping=aes(y=count, x=group2, fill=ex5cluster_class), width=0.7, linewidth=0.25, color = "black", stat="identity", alpha=0.8) + 
  geom_text(data=ee2[which(ee2$group2 != "Distal_2D"),], mapping=aes(y=count, x=group2, label=label, group=ex5cluster_class), size=1.8, position=position_stack(vjust=0.5), color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6c.end5_cluster.eRNA.definition.length2.pdf"), width = 2.2, height = 1.7)
print(ex6c)
dev.off()

#=
ee3=ee%>%group_by(group2)%>%dplyr::summarise(percent=length(which(median_exon>1))/n())
ee3$label=paste0(signif(ee3$percent,2) *100, "%")

ex6d=ggplot(ee, aes(x=median_exon, color=factor(group2))) + 
  labs(y="Ranked transcripts", x ="Number of exon",title= "Different enhancer typing", color=NULL) +
  scale_color_npg()+
  coord_flip(xlim=c(0,15))+
  stat_ecdf(geom = "step", linewidth=0.25)+
  scale_y_continuous(labels = scales::percent)+
  geom_text_repel(data=ee3, aes(x=1.5, y=1-percent, label=label,color=factor(group2)),segment.size = 0.25, size=1.8, force =5, nudge_x=4, nudge_y=0.01)+
  theme1+theme(legend.position=c(0.25,0.8))
pdf(paste0(path_fig4,"ex6d.end5_cluster.number_different.definition.pdf"), width = 1.9, height = 1.7)
print(ex6d)
dev.off()

#=====================
#Needed
#ex6e
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

data2=data1[which(data1$ex5cluster_class == "e_ncRNA"),]
CGIap=length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS < 2000))
CGInap=length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS >= 2000))
TATA=length(which(data2$TATA_box == "Yes"))
CGIapTATA=length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS < 2000 & data2$TATA_box == "Yes"))
CGInapTATA=length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS >= 2000 & data2$TATA_box == "Yes"))

length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS < 2000 & data2$TATA_box == "Yes" & data2$downstream_CpG_island == "Yes"))
#104
length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS < 2000 & data2$TATA_box == "No" & data2$downstream_CpG_island == "Yes"))
#3764
length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS >= 2000 & data2$TATA_box == "Yes" & data2$downstream_CpG_island == "Yes"))
#18
length(which(data2$any_CpG_island == "Yes" & data2$distance_cCRE_PLS >= 2000 & data2$TATA_box == "No" & data2$downstream_CpG_island == "Yes"))
#897
length(which(data2$CpGTATA == "Null")) #20610

library(VennDiagram)
grid.newpage()
venn.plot3 <- draw.triple.venn(area1 = CGIap, area2 = CGInap, area3= TATA, n12=0, n23=CGInapTATA, n13=CGIapTATA, n123 = 0,
                                 euler.d = TRUE, scaled = T, inverted=T, category=c("CGIap","CGInap","TATA"), fontfamily = rep("ArialMT", 7), cat.fontfamily = rep("ArialMT", 3), alpha=c(0.5,0.5,0.5),
                                 cat.just= list(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5)), fill = c("#4DBBD5FF", "#00A087FF", "#3C5488FF"), lty = "blank", cex = 0.4, cat.cex = 0.45, cat.col = "black")
pdf(paste0(path_fig4,"ex6e.e_ncRNA.venn.CpGTATA.pdf"), width = 0.75, height = 0.75)
grid.draw(venn.plot3)
dev.off()

data3=data1[which(data1$ex5cluster_class == "p_ncRNA"),]
CGI=length(which(data3$any_CpG_island == "Yes"))
CGITATA=length(which(data3$any_CpG_island == "Yes" & data3$TATA_box == "Yes"))
TATA=length(which(data3$TATA_box == "Yes"))

length(which(data3$any_CpG_island == "Yes" & data3$TATA_box == "Yes" & data3$downstream_CpG_island == "Yes"))
#166
length(which(data3$any_CpG_island == "Yes" & data3$TATA_box == "No" & data3$downstream_CpG_island == "Yes"))
#6990
length(which(data3$CpGTATA == "Null")) #4347

grid.newpage()
venn.plot3 <- draw.pairwise.venn(area1 = CGI, area2 = TATA, cross.area=CGITATA,
                               euler.d = TRUE, scaled = F, inverted=F, category=c("CGI","TATA"), fontfamily = rep("ArialMT", 3), cat.fontfamily = rep("ArialMT", 2), alpha=c(0.5,0.5),
                               cat.just= list(c(0.5, 0.5), c(0.5, 0.5)), fill = c("#4DBBD5FF", "#3C5488FF"), lty = "blank", cex = 0.4, cat.cex = 0.45, cat.col = "black")
pdf(paste0(path_fig4,"ex6e.p_ncRNA.venn.CpGTATA.pdf"), width = 0.75, height = 0.75)
grid.draw(venn.plot3)
dev.off()


#=====================
#Needed
#ex6f
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

data2=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","orientation","read_median_length","median_exon","median_range")]
data2=data2[which(data2$orientation != "Others"),]
data3=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","any_CpG_island","read_median_length","median_exon","median_range")]
data4=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","downstream_CpG_island","read_median_length","median_exon","median_range")]
data5=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","TATA_box","read_median_length","median_exon","median_range")]
data6=data1[which(data1$ex5cluster_class=="e_ncRNA" | data1$ex5cluster_class=="p_ncRNA"), c("ex5cluster_class","n5_string","conserve_up","read_median_length","median_exon","median_range")]
data6=data6[which(data6$conserve_up != "Others"),]
data7=data1[which(data1$ex5cluster_class=="e_ncRNA"), c("ex5cluster_class","n5_string","SE_all","read_median_length","median_exon","median_range")]
data7=data7[which(data7$SE_all != "Others"),]

colnames(data2)[3]="feature"
colnames(data3)[3]="feature"
colnames(data4)[3]="feature"
colnames(data5)[3]="feature"
colnames(data6)[3]="feature"
colnames(data7)[3]="feature"

data2$group2="1D / 2D"
data3$group2="all_CGI"
data4$group2="dCGI"
data5$group2="TATA box"
data6$group2="Conserved"
data7$group2="SE / TE"

data2=rbind(data2,data3,data4,data5,data6,data7)
data2$group2=factor(data2$group2, levels=c("all_CGI","dCGI","TATA box","1D / 2D","SE / TE","Conserved"))
data2$feature=factor(data2$feature, levels=c("1D","2D","No","Yes","SE","TE"))
data4=data2%>%group_by(ex5cluster_class, group2, feature)%>%dplyr::summarise(median=median(read_median_length),count=n())
data1wilcox3=data2%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(read_median_length ~feature, alternative = "two.sided")$p.value)
data1wilcox3$label="***"
data1wilcox3$label[which(data1wilcox3$p>=0.001)]="**"
data1wilcox3$label[which(data1wilcox3$p>=0.01)]="*"
data1wilcox3$label[which(data1wilcox3$p>=0.05)]="n.s."

ex6f=ggplot() + 
  labs(x=NULL, y = "Transcript length (bp)",title= "p_ncRNA length affected by promoter properties") +
  scale_fill_npg(guide=NULL)+
  facet_grid(cols=vars(group2), scale="free_x")+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=data2[which(data2$ex5cluster_class=="p_ncRNA" & data2$group2 !="Conserved"),], mapping=aes(y=read_median_length, x=feature, fill=group2),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=data2[which(data2$ex5cluster_class=="p_ncRNA" & data2$group2 !="Conserved"),], mapping=aes(y=read_median_length, x=feature, fill=group2), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data4[which(data4$ex5cluster_class=="p_ncRNA" & data4$group2 !="Conserved"),], mapping=aes(y=median, x=feature, label=signif(median,3)), size=2, angle=90, vjust=-0.7, hjust=0.5, color="black")+
  geom_text(data=data1wilcox3[which(data1wilcox3$ex5cluster_class=="p_ncRNA" & data1wilcox3$group2 !="Conserved"),], mapping=aes(y=2300, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 2250, yend = 2250, linewidth = 0.25)+
  theme1
pdf(paste0(path_fig4,"ex6f.end5_cluster.bidirectional.pRNA.length.pdf"), width = 2.8, height = 1.7)
print(ex6f)
dev.off()

#=====================
#Needed
#ex6g
#continue from last section
data6=data2%>%group_by(ex5cluster_class,group2, feature)%>%dplyr::summarise(median=median(median_exon), percent=length(which(median_exon >1))/n())
data6$label="> 1 exon: "

data1wilcox5=data2%>%group_by(ex5cluster_class,group2)%>%dplyr::summarise(p=wilcox.test(median_exon ~feature, alternative = "two.sided")$p.value)
data1wilcox5$label="***"
data1wilcox5$label[which(data1wilcox5$p>=0.001)]="**"
data1wilcox5$label[which(data1wilcox5$p>=0.01)]="*"
data1wilcox5$label[which(data1wilcox5$p>=0.05)]="n.s."

ex6g=ggplot() + 
  labs(x=NULL, y = "Number of exon",title= "p_ncRNA exon number affected by promoter properties") +
  scale_fill_npg(guide=NULL)+
  #scale_alpha_manual(values=c("bidirectional"=1, "unidirectional"=0.5), guide=NULL)+
  facet_grid(cols=vars(group2), scales="free_x")+
  coord_cartesian(ylim=c(0,9.5))+
  geom_violin(data=data2[which(data2$ex5cluster_class == "p_ncRNA" & data2$group2 !="Conserved"),], mapping=aes(y=median_exon, x=feature, fill=group2, alpha=feature), size=0.1, color = "black", bounds=c(0,8), width = 0.9, position = position_nudge(x = +0), alpha=1) + 
  geom_text(data=data6[which(data6$ex5cluster_class=="p_ncRNA" & data6$group2 !="Conserved"),], mapping=aes(y=2, x=feature, label=label), angle=90, size=2, vjust=(-0.8), hjust=0, color="grey")+
  geom_text(data=data6[which(data6$ex5cluster_class=="p_ncRNA" & data6$group2 !="Conserved"),], mapping=aes(y=5.3, x=feature, label=paste0(signif(percent*100,2),"%")), angle=90, size=2, vjust=(-0.8), hjust=0, color="black")+
  geom_text(data=data1wilcox5[which(data1wilcox5$ex5cluster_class == "p_ncRNA"& data1wilcox5$group2 !="Conserved"),], mapping=aes(y=9.4, x=1.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 2, y = 9.2, yend = 9.2, linewidth = 0.25)+
  theme1
pdf(paste0(path_fig4,"ex6g.end5_cluster.bidirectional.CGI.TATA.pRNA.exon.number.pdf"), width = 2.8, height = 1.7)
print(ex6g)
dev.off()

#=====================
#Needed
#ex6h
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

dataa=data1[which(data1$ex5cluster_class == "p_ncRNA" & data1$CpGTATA != "Others"),]
dataa$CpGTATA[which(dataa$CpGTATA=="CGI" & dataa$downstream_CpG_island == "Yes")]="dCGI"
dataa$CpGTATA[which(dataa$CpGTATA=="CGI" & dataa$downstream_CpG_island == "No")]="otherCGI"

data2a=dataa[which(dataa$CpGTATA %in% c("otherCGI","Null")),c("CpGTATA","read_median_length","median_exon")]
data2b=dataa[which(dataa$CpGTATA %in% c("dCGI","Null")),c("CpGTATA","read_median_length","median_exon")]
data2c=dataa[which(dataa$CpGTATA %in% c("TATA","Null")),c("CpGTATA","read_median_length","median_exon")]

data2a$group2="otherCGI"
data2b$group2="dCGI"
data2c$group2="TATA"
data2=rbind(data2a,data2b,data2c)

data3=dataa%>%group_by(CpGTATA)%>%dplyr::summarise(median=median(read_median_length, na.rm=T),count=n())
data4=data2%>%group_by(group2)%>%dplyr::summarise(p=wilcox.test(read_median_length ~CpGTATA, alternative = "two.sided")$p.value)
data4$label="***"
data4$label[which(data4$p>=0.001)]="**"
data4$label[which(data4$p>=0.01)]="*"
data4$label[which(data4$p>=0.05)]="n.s."

dataa$CpGTATA=factor(dataa$CpGTATA, levels=c("dCGI","otherCGI","Null","TATA"))

ex6h=ggplot() + 
  labs(x=NULL, y = "Transcript length (bp)",title= "p_ncRNA length affected\nby regulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,2500))+
  ggdist::stat_halfeye(data=dataa, mapping=aes(y=read_median_length, x=CpGTATA, fill=CpGTATA),adjust = .75, width = .5, .width = 0, justification = -.4, point_colour = NA, alpha=0.4)+
  geom_boxplot(data=dataa, mapping=aes(y=read_median_length, x=CpGTATA, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data3, mapping=aes(y=median, x=CpGTATA, label=signif(median,3)), angle=90, size=2, vjust=(-0.8), hjust=0.5, color="black")+
  geom_text(data=data4[which(data4$group2 == "dCGI"),], mapping=aes(y=2225, x=2, label=label),size=2.4)+
  geom_text(data=data4[which(data4$group2 == "otherCGI"),], mapping=aes(y=2025, x=2.5, label=label),size=2.4)+
  geom_text(data=data4[which(data4$group2 == "TATA"),], mapping=aes(y=2425, x=3.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 3, y = 2200, yend = 2200, linewidth = 0.25)+
  annotate("segment", x = 2, xend = 3, y = 2000, yend = 2000, linewidth = 0.25)+
  annotate("segment", x = 3, xend = 4, y = 2400, yend = 2400, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6h.end5_cluster.p_ncRNA_length_CpG_TATA.pdf"), width = 1.5, height = 1.7)
print(ex6h)
dev.off() 

#=====================
#Needed
#ex6i

#continue from last section
edata3=dataa%>%group_by(CpGTATA)%>%dplyr::summarise(mean=mean(median_exon), percent=length(which(median_exon !=1))/n() ,count=n())
edata3$label="> 1 exon: "
edata4=data2%>%group_by(group2)%>%dplyr::summarise(p=wilcox.test(median_exon ~CpGTATA, alternative = "two.sided")$p.value)
edata4$label="***"
edata4$label[which(edata4$p>=0.001)]="**"
edata4$label[which(edata4$p>=0.01)]="*"
edata4$label[which(edata4$p>=0.05)]="n.s."

dataa$CpGTATA=factor(dataa$CpGTATA, levels=c("dCGI","otherCGI","Null","TATA"))

ex6i=ggplot() + 
  labs(x=NULL, y = "Number of exon",title= "p_ncRNA splicing affected\nby regulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF", "#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  coord_cartesian(ylim=c(0,10))+
  geom_violin(data=dataa, mapping=aes(y=median_exon, x=CpGTATA, fill=CpGTATA, alpha=orientation), size=0.15, color = "black", bounds=c(0,8), width = 0.9, position = position_nudge(x = +0), alpha=1) + 
  geom_text(data=edata3, mapping=aes(y=2, x=CpGTATA, label=label), angle=90, size=2, vjust=(-0.8), hjust=0, color="grey")+
  geom_text(data=edata3, mapping=aes(y=5.3, x=CpGTATA, label=paste0(signif(percent*100,2),"%")), angle=90, size=2, vjust=(-0.8), hjust=0, color="black")+
  geom_point(data=edata3, mapping=aes(y=mean, x=CpGTATA), size=0.2, color="red")+
  geom_text(data=edata4[which(edata4$group2 == "dCGI"),], mapping=aes(y=8.7, x=2, label=label),size=2.4)+
  geom_text(data=edata4[which(edata4$group2 == "otherCGI"),], mapping=aes(y=8, x=2.5, label=label),size=2.4)+
  geom_text(data=edata4[which(edata4$group2 == "TATA"),], mapping=aes(y=9.4, x=3.5, label=label),size=2.4)+
  annotate("segment", x = 1, xend = 3, y = 8.5, yend = 8.5, linewidth = 0.25)+
  annotate("segment", x = 2, xend = 3, y = 7.8, yend = 7.8, linewidth = 0.25)+
  annotate("segment", x = 3, xend = 4, y = 9.2, yend = 9.2, linewidth = 0.25)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6i.end5_cluster.p_ncRNA_exon_CpG_TATA.pdf"), width = 1.5, height = 1.7)
print(ex6i)
dev.off() 

#=====================
#Needed
#ex6j
data1=read.delim(paste0(path_fig4_data,"ex5_cluster.n.transcript.base.length.relative.MAD.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/RNA_features.R

data1$group=factor(data1$group,levels=c("mRNA","p_ncRNA","e_ncRNA", "CTCF_ncRNA","other_ncRNA"))
data4=data1[which(data1$n_read>=3),]%>%group_by(group,polyArate_bin)%>%dplyr::summarise(rel_MAD=median(rel_MAD),count=n())

ex6j=ggplot() +
  labs(y = "Relative MAD (MAD/median)",title= "RNA length variation", x="Rate of poly(A) RNA", color="ex5_cluster") +
  scale_color_manual(values=c("grey", "#E64B35FF", "#3C5488FF","#00A087FF"))+
  coord_cartesian(ylim=c(0,0.5))+
  scale_x_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1))+
  geom_point(data=data4[which(data4$group != "CTCF_ncRNA"),], mapping=aes(x=polyArate_bin, y=rel_MAD, color=group), size=0.25,  alpha=1)+
  geom_line(data=data4[which(data4$group != "CTCF_ncRNA"),], mapping=aes(x=polyArate_bin, y=rel_MAD,color=group), linewidth=0.25)+
  theme1+theme(legend.position=c(0.85,0.75))
pdf(paste0(path_fig4,"ex6j.ex5_cluster.base.length.relative.MAD_polyA.bin.pdf"), width = 1.3, height = 1.7)
print(ex6j)
dev.off()

#=====================
#Needed
#ex6k
chr_table5a=read.delim(paste0(path_fig4_data,"chr_bound_iPSC_gene_count.plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

chr_table5a$group=factor(chr_table5a$group,levels=c("ENSG","ONTG","CHRG"))
chr_table5a$group2=factor(chr_table5a$group2,levels=c("p_ncRNA","e_ncRNA","other_ncRNA"))
chr_table5a$label=chr_table5a$count
chr_table5a$label[which(chr_table5a$group3=="Detected")]=NA

ex6k=ggplot(chr_table5a, aes(x=group, y=count, fill=group2, alpha=group3)) + 
  scale_fill_npg()+
  scale_alpha_manual(values=c(0.3,0.8), guide=NULL)+
  labs(x=NULL , y ="Number of gene" , title ="Genes identified from\nchromatin-bound iPSC",  fill=NULL)+
  scale_x_discrete(labels=c("Detected as\nGENCODE ncRNAs","Detected as\nONTG ncRNA","Novel ncRNA from\nthis dataset"))+
  geom_bar(linewidth=0.25, stat="identity", color="black")+
  geom_text(data=chr_table5a, mapping=aes(x=group, y=count, group=group2, label=label), size=1.8, color="black", position = position_stack(vjust = 0.4, reverse=F))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1, lineheight = 0.7))
pdf(paste0(path_fig4,"ex6k.chr_bound_iPSC_gene_count.pdf"), width = 2, height = 1.7)
print(ex6k)
dev.off() 

#=====================
#Needed
#ex6l 
data1=read.delim(paste0(path_fig4_data,"chr_bound_iPSC_length_diff.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

data1$polyA=gsub("Yes","poly(A)",data1$polyA)
data1$polyA=gsub("No","non-poly(A)",data1$polyA)
data1$polyA=factor(data1$polyA, levels=c("poly(A)","non-poly(A)"))
data1$transcript_group=factor(data1$transcript_group, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA","others"))

data2=data1%>%group_by(polyA,transcript_group)%>%dplyr::summarise(median=median(length_diff,na.rm=T), count=n())

ex6l=ggplot() +
  labs(y = "Chromatin - Total (nt)",title= "Transcript length difference\nfrom same transcript in iPSC", x=NULL) +
  scale_fill_manual(values=c("grey", "#E64B35FF", "#3C5488FF","#00A087FF"), guide=NULL)+
  facet_grid(cols=vars(polyA))+
  coord_cartesian(ylim=c(-300,300))+
  geom_boxplot(data=data1[which(data1$transcript_group != "others"),], mapping=aes(x=transcript_group, y=length_diff, fill=transcript_group),color="black", linewidth=0.25, linetype="solid", width=0.3, alpha=1,outlier.shape = NA)+
  geom_text(data=data2[which(data2$transcript_group != "others"),], mapping=aes(x=transcript_group, y=median,label=paste0(signif(median,3),", n=",count)), size=2, angle=90, hjust=0.2, vjust=-0.75)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6l.transcript_base_chr_main_iPSC_diff.pdf"), width = 2, height = 1.7)
print(ex6l)
dev.off()

#=====================
#Needed
#ex6m
exp.df2=read.delim(paste0(path_fig4_data,"Transcript_chr_total_DET_exosome_iPSC.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

exp.df2$bin[which(exp.df2$bin>5)]=5
exp.df2$bin[which(exp.df2$bin<(-3))]=(-3)
exp.df2$group2="Others"
exp.df2$group2[grep("e_",exp.df2$transcript_group)]="e_ncRNA"
exp.df2$group2[grep("p_",exp.df2$transcript_group)]="p_ncRNA"

exp.df3=exp.df2[which(exp.df2$group2!="Others" & exp.df2$logCPM>(-2)),]%>%group_by(group2, bin)%>%dplyr::summarise(count=n(), median=median(exo_sensitivity))

ex6m=ggplot()+
  geom_boxplot(data=exp.df2[which(exp.df2$group2!="Others" & exp.df2$logCPM > (-2)),], mapping=aes(x=as.factor(bin), y=exo_sensitivity), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.7, outlier.shape = NA)+
  facet_grid(cols=vars(group2))+
  scale_x_discrete(labels=c("<= -3","-2","-1","0","1","2","3","4",">= 5"))+
  labs(y="Exosome sensitivity", x="Log2FC (chromatin-bound/total)", title="Exosome sensitivty of iPSC ncRNAs")+
  geom_smooth(data = exp.df3, aes(x = as.numeric(as.factor(bin)), y = median), method = "loess", color = "blue", linewidth=0.25, se=F) +
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1, lineheight = 0.7))
pdf(paste0(path_fig4,"ex6m.chr_bound_iPSC_e_p_RNATx_exo_sensitivity.pdf"), width = 2, height = 1.7)
print(ex6m)
dev.off() 

#=====================
#Needed
#ex6n
RBP_trans_exon=read.delim(paste0(path_fig4_data,"RBP.nc.enhancer.FE.result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

need1=RBP_trans_exon$RBP[intersect(which(RBP_trans_exon$group == "coding_poential"), which(RBP_trans_exon$FE_logOR > 0.4 | RBP_trans_exon$FE_logOR < (-0.7)))]
need2=RBP_trans_exon$RBP[intersect(which(RBP_trans_exon$group == "promoter_type"), which(RBP_trans_exon$FE_logOR > 0.4 | RBP_trans_exon$FE_logOR < (-0.7)))]
need=union(need1, need2)
RBP_trans_exon1=RBP_trans_exon[which(RBP_trans_exon$RBP %in% need),]

ex6n=ggplot() + 
  labs(x=NULL, y = "Selected RBP",title= "Enrichment") +
  geom_tile(data=RBP_trans_exon1, mapping=aes(y=reorder(RBP,FE_logOR), x=group, fill=FE_logOR),linewidth=0.2, color="black")+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0, limits=c(-1.2,0.7))+ 
  scale_x_discrete(labels=c("coding vs non-coding", "promoter vs enhancer"))+
  theme1 + theme(axis.text.x = element_text(color ="black", angle=90, hjust=1, vjust=1))
pdf(paste0(path_fig4,"ex6n.RBP.nc.enhancer.FE.result.pdf"), width = 1.3, height = 3.3)
print(ex6n)
dev.off() 

#=====================
#Needed
#ex6o
data1=read.delim(paste0(path_fig4_data,"features_by_ex5cluster.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
# from ./code_n_data/Fig4_transcription_features/chromatin_bound_n_exosome.R

data2=data1[which(data1$ex5cluster_class %in% c("p_ncRNA","e_ncRNA") & !is.na(data1$exo_sensitivity_Tx)),]
data2=data2[,c("ex5cluster_class", "n5_string", "exo_sensitivity_Tx", "any_CpG_island", "downstream_CpG_island", "TATA_box", "orientation", "SE_all", "conserve_down", "ubiquitous")]
colnames(data2)[c(4:10)]=c("all_CGI","dCGI","TATA box","1D / 2D","SE / TE","Conserved","Ubiquitous")
z1=reshape2::melt(data2,id=c(1:3))
z1=z1[which(z1$value %in% c("Yes", "No", "2D", "1D","SE","TE")),]

z1a=z1%>%group_by(ex5cluster_class, variable, value)%>%dplyr::summarise(median=median(exo_sensitivity_Tx))
z1b=z1%>%group_by(ex5cluster_class, variable)%>%dplyr::summarise(p=wilcox.test(exo_sensitivity_Tx ~value, alternative = "two.sided")$p.value)
z1b$label="***"
z1b$label[which(z1b$p>=0.001)]="**"
z1b$label[which(z1b$p>=0.01)]="*"
z1b$label[which(z1b$p>=0.05)]="n.s."

ex6o=ggplot() + 
  labs(x=NULL, y = "Exosome sensitivity",title= "RNA exosome sensitivity against promoter features") +
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#8491B4FF","#91D1C2FF"),guide=NULL)+
  facet_grid(cols=vars(variable), scale="free_x")+
  coord_cartesian(ylim=c(0,1.2))+
  geom_boxplot(data=z1[which(z1$ex5cluster_class=="p_ncRNA"),], mapping=aes(y=exo_sensitivity_Tx, x=value, fill=variable), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.4, outlier.shape = NA) + 
  geom_text(data=z1a[which(z1a$ex5cluster_class=="p_ncRNA"),], mapping=aes(y=median, x=value, label=signif(median,3)), size=2, angle=90, vjust=-0.75, hjust=0.5, color="black")+
  geom_text(data=z1b[which(z1b$ex5cluster_class=="p_ncRNA"),], mapping=aes(y=1.1, x=1.5, label=label),size=2.4)+
  annotate("segment",x=1,xend=2,y=1.05,yend=1.05, linewidth=0.25)+
  theme1
pdf(paste0(path_fig4,"ex6o.n5_cluster_promoter_feature_exo_sensitivity.pdf"), width = 2.9, height = 1.7)
print(ex6o)
dev.off()

