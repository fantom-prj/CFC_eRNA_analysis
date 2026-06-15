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
path_fig2=paste0(primary_folder,"fig2/out/")
path_fig2_data=paste0(primary_folder,"fig2/data/")
path_fig3_data=paste0(primary_folder,"fig3/data/")

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
#s1a
table5=read.delim(paste0(path_fig3_data,"table5.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)

data1=separate_rows(unique(table5[,c("n5_string", "TSScluster", "CREID")]),TSScluster, sep=";")
data1a=data1%>%group_by(CREID)%>%dplyr::summarise(count=n())
data1a1=data1a%>%group_by(count)%>%dplyr::summarise(frequency=n())
data1b=unique(data1[,c(1,3)])%>%group_by(CREID)%>%dplyr::summarise(count=n())
data1b1=data1b%>%group_by(count)%>%dplyr::summarise(frequency=n())
data1c=data1[,c(1,2)]%>%group_by(n5_string)%>%dplyr::summarise(count=n())
data1c1=data1c%>%group_by(count)%>%dplyr::summarise(frequency=n())
data1a1$group="TSScluster/\ntCRE"
data1b1$group="Ex5_cluster/\ntCRE"
data1c1$group="TSScluster/\nEx5_cluster"
data1a1_p=rbind(data1a1,data1b1,data1c1)
data1a1_p=data1a1_p%>%group_by(group)%>%dplyr::mutate(percent=frequency/sum(frequency))
data1a1_p1=data1a1_p%>%group_by(group)%>%dplyr::summarise(label=sum(frequency*count)/sum(frequency))
s1a=ggplot() + 
  labs(x="Count", y = NULL ,title= "5' end annotation", fill=NULL) +
  scale_fill_manual(values=c("white","grey","black"), guide=NULL)+
  scale_y_continuous(labels = scales::percent, limits=c(0,1), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(group))+
  geom_bar(data=data1a1_p[which(data1a1_p$count<=5),], mapping=aes(y=percent, x=count, fill=group), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=data1a1_p1, mapping=aes(x=3, y=0.8, label=paste0("average=\n",signif(label,3))), size=1.8, hjust=0.5, vjust=0)+
  theme1 + theme(legend.position = "bottom")
pdf(paste0(path_fig2,"s1a.SCAFE_5n_annotation.pdf"), width = 1.9, height = 1.6)
print(s1a)
dev.off()

#=====================
#Needed
#s1b
table5=read.delim(paste0(path_fig3_data,"TS3.chimeric.194K.remove.permissive.isoform.tsv.gz"), header=T, check.names=F, stringsAsFactors=F)

data2=separate_rows(unique(table5[,c("model_ID","n5_string","TSScluster","CREID")]),TSScluster, sep=";")
data2a=data2%>%group_by(CREID)%>%dplyr::summarise(count=n())
data2a1=data2a%>%group_by(count)%>%dplyr::summarise(frequency=n())
data2b=data2%>%group_by(n5_string)%>%dplyr::summarise(count=n())
data2b1=data2b%>%group_by(count)%>%dplyr::summarise(frequency=n())
data2c=data2%>%group_by(TSScluster)%>%dplyr::summarise(count=n())
data2c1=data2c%>%group_by(count)%>%dplyr::summarise(frequency=n())
data2a1$group="Tx/\ntCRE"
data2b1$group="Tx/\nEx5_cluster"
data2c1$group="Tx/\nTSScluster"
data2a1_p=rbind(data2a1,data2b1,data2c1)
data2a1_p$group=factor(data2a1_p$group, levels=c("Tx/\nTSScluster","Tx/\nEx5_cluster","Tx/\ntCRE"))
data2a1_p=data2a1_p%>%group_by(group)%>%dplyr::mutate(percent=frequency/sum(frequency))
data2a1_p1=data2a1_p%>%group_by(group)%>%dplyr::summarise(label=sum(frequency*count)/sum(frequency))
s1b=ggplot() + 
  labs(x="Count", y = NULL ,title= "Transcript model from 5'end region", fill=NULL) +
  scale_fill_manual(values=c("white","grey","black"), guide=NULL)+
  scale_y_continuous(labels = scales::percent, limits=c(0,0.55), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(group))+
  geom_bar(data=data2a1_p[which(data2a1_p$count<=25),], mapping=aes(y=percent, x=count, fill=group), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=data2a1_p1, mapping=aes(x=12.5, y=0.4, label=paste0("average=\n",signif(label,3))), size=1.8, hjust=0.5, vjust=0)+
  theme1
pdf(paste0(path_fig2,"s1b.Tx_to_5n_annotation.pdf"), width = 1.9, height = 1.6)
print(s1b)
dev.off()

#=====================
#Needed
#s1c
#allow redundant sequence
gc2=read.delim(paste0(path_fig2_data,"end5_cluster_both_GCcontent_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
gc2=gc2[which(gc2$promoter_type != "CTCF-alone"),]
gc2$promoter_type=factor(gc2$promoter_type, levels=c("promoter-like","enhancer-like","unclassed"))
gc2$group2[which(gc2$group == "random")]="background"
gc2$group2[which(gc2$group == "anno_region")]="tested\nregion"
gc2$group2=factor(gc2$group2, levels=c("tested\nregion", "background"))

phastcon=read.delim(paste0(path_fig2_data,"end5cluster_phastconALL_plot.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
phastcon$promoter_type=factor(phastcon$promoter_type, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
phastcon=phastcon[which(phastcon$label %in% c("4way","17way","30way","100way")),]
phastcon$label=factor(phastcon$label, levels=c("4way","17way","30way","100way"))

cpg=read.delim(paste0(path_fig2_data,"end5_cluster_CpG_island_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
cpg$anno_region=factor(cpg$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
cpg=cpg[which(cpg$orientation == "Others"),]
cpg$signalID=gsub("n","min_",cpg$signalID)
cpg$signalID=factor(cpg$signalID, levels=c("min_20","min_40","min_60","min_80"))

tata=read.delim(paste0(path_fig2_data,"end5_cluster_TATA_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
tata$anno_region=factor(tata$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
tata=tata[which(tata$orientation == "Others"),]
tata$signalID=gsub("n","min_",tata$signalID)
tata$signalID=factor(tata$signalID, levels=c("min_1","min_2","min_3","min_4"))

inr=read.delim(paste0(path_fig2_data,"end5_cluster_INR_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
inr$anno_region=factor(inr$anno_region, levels=c("promoter-like","enhancer-like","unclassed", "CTCF-alone"))
inr=inr[which(inr$orientation == "Others"),]
inr$signalID=gsub("n","min_",inr$signalID)
inr$signalID=factor(inr$signalID, levels=c("min_1","min_2","min_3","min_4"))

#labs(x="Distance relative to cluster summit (nt)", y="% of ex5_clusters overlap signal peaks")
C9=ggplot(gc2, aes(x=position, y=GCprecent/100, color=group)) + 
  scale_color_npg()+
  scale_y_continuous(labels = scales::percent)+
  labs(x="tCRE summit" , y ="% of GC" , title ="GC content",  color=NULL)+
  facet_grid(rows=vars(promoter_type) ,  scales="free_y")+
  geom_line(alpha=0.75, linewidth=0.25)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C10=ggplot(phastcon[which(phastcon$promoter_type != "CTCF-alone"),], aes(x=positionV3, y=score, color=label, group=label)) + 
  scale_color_npg()+
  coord_cartesian(ylim = c(0.1, NA))+
  labs(x=NULL , y ="Mean phastCon score" , title ="Conservation",  color=NULL)+
  facet_grid(rows=vars(promoter_type) ,  scales="free_y")+
  geom_line(alpha=0.75, linewidth=0.25)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C11=ggplot()+
  scale_color_npg()+
  geom_line(data = cpg[which(cpg$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="CpG island")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C12=ggplot()+
  scale_color_npg()+
  geom_line(data = tata[which(tata$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="TATA box")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")
C13=ggplot()+
  scale_color_npg()+
  geom_line(data = inr[which(inr$anno_region != "CTCF-alone"),], aes(x=V4, y=V5, color=signalID, group=signalID),alpha=0.75, linewidth=0.25)+
  facet_grid( rows=vars(anno_region), scales="free_y")+
  labs(color=NULL, x=NULL, y=NULL, title="Initiator")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(strip.text.y = element_blank(), axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1),legend.position = "bottom", legend.direction = "vertical")

pdf(paste0(path_fig2,"s1c.end5_cluster_CpG_TATA_inia_phase.pdf"), width = 4.6, height = 3)
grid.arrange(C9,C10,C11,C12,C13, ncol = 5, widths = c(1.1,  1.1, 1, 1, 1))
dev.off()

#=====================
#Needed
#s1d
n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

n5cluster1=n5cluster[which(n5cluster$promoter_type %in% c("promoter-like","enhancer-like") & n5cluster$region1001_rep == "Yes"),]
n5cluster1=n5cluster1[which(n5cluster1$orientation %in% c("2D","1D")),]
n5cluster1$group=paste0(substr(n5cluster1$promoter_type,1,1),"_",n5cluster1$orientation)

n5cluster1$group2="Null"
n5cluster1$group2[which(n5cluster1$any_CpG_island == "Yes")]="CGI"
n5cluster1$group2[which(n5cluster1$TATA_box == "Yes")]="TATA"
n5cluster1$group2[which(n5cluster1$any_CpG_island == "Yes" & n5cluster1$TATA_box == "Yes")]="Mix"
n5cluster1$group2=factor(n5cluster1$group2, levels=c("Null","TATA","Mix","CGI"))

n5cluster2=n5cluster1%>%group_by(group, group2)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5cluster2$label=paste0(signif(n5cluster2$percent,2)*100,"%")
n5cluster2$label[which(n5cluster2$group2 %in% c("Null","Mix"))]=NA

s1d=ggplot(n5cluster2, aes(x=as.factor(group), y=percent, fill=group2)) + 
  scale_fill_manual(values=c("white","#4DBBD5FF","darkorchid3","#E64B35FF"))+
  labs(fill=NULL, x=NULL, y ="% of tCREs", title="Regulatory elements")+
  scale_y_continuous(labels = scales::percent)+
  geom_bar(stat="identity", color="black", linewidth=0.25, width=0.85,alpha=0.7)+
  geom_text(aes(label=label), size=1.8, position=position_stack(vjust=0.5), color="black", angle=0)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1), legend.position = c(0.32,0.8))
pdf(paste0(path_fig2,"s1d.CPG_TATA_1D_2D.pdf"), width = 1.25, height = 2)
print(s1d)
dev.off() 

#=====================
#Needed
#s1e
data4=read.delim(paste0(path_fig2_data,"filtered_ex5cluster_for_conservation.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
#file filtered for 1001nt overlap removed

data4a=reshape2::melt(data4[,c(7,9,3,10)], id=c(3,4))
data4a$variable=factor(data4a$variable, levels=c("TATA","CGI"))
data4b=data4a%>%group_by(group,variable, value)%>%dplyr::summarise(median=median(mean_1001, na.rm=T),count=n())
data4c=data4a%>%group_by(group,variable)%>%dplyr::summarise(p=wilcox.test(mean_1001 ~value, alternative = "two.sided")$p.value)
data4c$label="***"
data4c$label[which(data4c$p>=0.001)]="**"
data4c$label[which(data4c$p>=0.01)]="*"
data4c$label[which(data4c$p>=0.05)]="n.s."

s1e=ggplot() + 
  labs(x=NULL, y = "Mean PhastCon score",title= "Conservation from regulatory elements") +
  scale_fill_manual(values=c("#4DBBD5FF","#E64B35FF"),guide=NULL)+
  facet_grid(rows=vars(variable), cols=vars(group), scale="free_x")+
  coord_cartesian(ylim=c(0,0.95))+
  ggdist::stat_halfeye(data=data4a, mapping=aes(y=mean_1001, x=value, fill=variable),adjust = .75, width = .5, .width = 0, justification = -.4, point_colour = NA, alpha=0.4)+
  geom_boxplot(data=data4a, mapping=aes(y=mean_1001, x=value, fill=variable), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  geom_text(data=data4b, mapping=aes(y=median, x=value, label=signif(median,3)), angle=90, size=2, vjust=(-0.8), hjust=0.5, color="black")+
  geom_text(data=data4c, mapping=aes(y=0.9, x=1.5, label=label),size=2.4)+
  annotate("segment", x=1,xend=2,y=0.85,yend=0.85, linewidth=0.25)+
  theme1
pdf(paste0(path_fig2,"s1e.n5cluster_PhastCon_CpG_TATA_1D_2D.pdf"), width = 2.5, height = 2)
print(s1e)
dev.off() 

#=====================
#Needed
#s1f
n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
n5cluster1=n5cluster[which(n5cluster$representative == "Yes"),]

n5cluster1$Others="Yes"
n5cluster1$Others[which(n5cluster1$CGI == "Yes" | n5cluster1$TATA == "Yes")]="No"
n5cluster2=reshape2::melt(n5cluster1[which(n5cluster1$promoter_type == "enhancer-like"),c("distance_cCRE_PLS","CGI","CGIap","CGInap","TATA","Others")], id=1)
n5cluster2=n5cluster2[which(n5cluster2$value == "Yes"),]
n5cluster2$variable=factor(n5cluster2$variable, levels=c("CGI","CGIap","CGInap","Others","TATA"))
n5cluster3=n5cluster2%>%group_by(variable)%>%dplyr::summarise(count=n(),median=median(distance_cCRE_PLS),mean=mean(distance_cCRE_PLS))

s1f=ggplot() + 
  labs(x=NULL, y = "Minimum distance (nt)",title= "Distance of enhancer\nfrom promoter cCRE") +
  scale_fill_npg(guide=NULL)+
  ggdist::stat_halfeye(data=n5cluster2, mapping=aes(y=distance_cCRE_PLS+1, x=variable, fill=variable),adjust = .5, width = .4, .width = 0, justification = -.4, point_colour = NA, alpha=0.5)+
  geom_boxplot(data=n5cluster2, mapping=aes(y=distance_cCRE_PLS+1, x=variable, fill=variable), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.2, outlier.shape = NA) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  #geom_point(data=n5cluster3, mapping=aes(y=mean, x=variable, color=variable), linewidth=0.3, shape=19, color="red")+
  geom_text(data=n5cluster3, mapping=aes(y=median+1, x=variable, label=paste0(signif(median,3),"nt (n=",count,")")), size=2, angle=90, vjust=-0.7, hjust=0.5, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig2,"s1f.ex5_CpGTATA_PLS_distance.pdf"), width = 1.5, height = 2)
print(s1f)
dev.off() 


#=====================
#Needed
#s1g
n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
n5cluster1=n5cluster[which(n5cluster$representative == "Yes"),]

n5cluster1$group=factor(n5cluster1$group,levels=c("Un-marked","Repressed","Active","Co-marked"))
n5cluster1$express="THP-1"
n5cluster1$express[which(n5cluster1$ex5_NSC>0 | n5cluster1$ex5_Neuron>0)]="NSC or Neuron"
n5cluster1$express[which(n5cluster1$ex5_iPSC>0)]="iPSC"

#n5clusterab1=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab2=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like" & n5cluster1$CGI == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab3=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like" & n5cluster1$TATA == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab4=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like" & n5cluster1$CGIap == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab5=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like" & n5cluster1$CGInap == "Yes"),]%>%group_by(group, express)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

n5clusterab2$group2="CGI"
n5clusterab3$group2="TATA"
n5clusterab4$group2="CGIap"
n5clusterab5$group2="CGInap"

n5clusterab1=rbind(n5clusterab2,n5clusterab3,n5clusterab4,n5clusterab5)
n5clusterab1$express=factor(n5clusterab1$express, levels=c("THP-1","NSC or Neuron","iPSC"))
n5clusterab1$group2=factor(n5clusterab1$group2, levels=c("CGI","CGIap","CGInap","TATA"))
n5clusterab1$group=factor(n5clusterab1$group, levels=c("Un-marked","Repressed","Active","Co-marked"))
n5clusterab2=n5clusterab1%>%group_by(group2, group)%>%dplyr::summarise(cluster=sum(count))
n5clusterab2$label=paste0("n=",n5clusterab2$cluster)

s1g=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Expression from enhancer ex5_cluster", fill=NULL) +
  scale_fill_manual(values=c("white","grey","black"), guide = guide_legend(reverse = TRUE) )+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.35), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(rows=vars(group2))+
  geom_bar(data=n5clusterab1, mapping=aes(x=percent, y=group, fill=express), linewidth=0.2, color = "black", stat="identity", alpha=1) + 
  geom_text(data=n5clusterab2, mapping=aes(x=1.02, y=group, label=label), size=1.8, hjust=0)+
  theme1 + theme(legend.position = "right" )
pdf(paste0(path_fig2,"s1g.Expression_Repressed_enhacer.pdf"), width = 2.5, height = 2)
print(s1g)
dev.off()


#=====================
#Needed
#s1h
n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
n5cluster1=n5cluster[which(n5cluster$representative == "Yes"),]

n5cluster1$group=factor(n5cluster1$group,levels=c("Un-marked","Repressed","Active","Co-marked"))
n5clustera=reshape2::melt(n5cluster1[which(n5cluster1$promoter_type == "enhancer-like"),c("group","CGI","CGIap","CGInap","TATA")], id=1)
n5clusterab=n5clustera%>%group_by(group, variable, value)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab$value=factor(n5clusterab$value, levels=c("Yes","No"))
n5clusterab$label=paste0(n5clusterab$group,"_",n5clusterab$value)

n5cluster1$Repressed="No"
n5cluster1$Repressed[which(n5cluster1$group == "Repressed")]="Yes"
n5cluster1$Active="No"
n5cluster1$Active[which(n5cluster1$group == "Active")]="Yes"
n5cluster1$'Co-marked'="No"
n5cluster1$'Co-marked'[which(n5cluster1$group == "Co-marked")]="Yes"
n5cluster1$'Un-marked'="No"
n5cluster1$'Un-marked'[which(n5cluster1$group == "Un-marked")]="Yes"

n5clusterac=reshape2::melt(n5cluster1[which(n5cluster1$promoter_type == "enhancer-like"),c("Repressed","Active","Co-marked","Un-marked","CGI","CGIap","CGInap","TATA")], id=c(1:4))
colnames(n5clusterac)[c(5,6)]=c("variable2","value2")
n5clusterac=reshape2::melt(n5clusterac, id=c(5:6))
n5clusterac$result=paste0(n5clusterac$value,"_",n5clusterac$value2)
n5clusterac1=n5clusterac%>%group_by(variable, variable2, result)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

n5clusterad=spread(n5clusterac1[,c(1:4)],key=3, value=4)
colnames(n5clusterad)=c("group","variable","no_no","no_yes","yes_no","yes_yes")
for(i in 1:nrow(n5clusterad)){
  GSEATasting <- matrix(c(n5clusterad$no_no[i], n5clusterad$no_yes[i], n5clusterad$yes_no[i], n5clusterad$yes_yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  n5clusterad$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  n5clusterad$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
n5clusterad$label="***"
n5clusterad$label[which(n5clusterad$p.val>=0.001)]="**"
n5clusterad$label[which(n5clusterad$p.val>=0.01)]="*"
n5clusterad$label[which(n5clusterad$p.val>=0.05)]="n.s."
n5clusterad$label1=paste0("OR=",signif(n5clusterad$OR,2),n5clusterad$label)

s1h=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enrichment with iPSC chromatin state", fill="CGI/ CGIap/\nCGInap/ TATA") +
  scale_fill_manual(values=c("black","white"))+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.45), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(rows=vars(variable))+
  geom_bar(data=n5clusterab, mapping=aes(x=percent, y=group, fill=value), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=n5clusterad, mapping=aes(y=group, x=1.01, label=label1), size=1.8, color="black", hjust=0, vjust=0.5)+
  theme1+theme(legend.position=c(0.25,0.75))
pdf(paste0(path_fig2,"s1h.ex5cluster.Repressed_enhacer.pdf"), width = 2, height = 2.3)
print(s1h)
dev.off()


#=====================
#Needed
#s1i
n5cluster=read.delim(paste0(path_fig2_data,"ex5_cluster.info.p.e.se.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
n5clustera=n5cluster[which(n5cluster$representative == "Yes"),]

n5cluster1=n5clustera[,c("ex5_iPSC","K27ME3_iPS","K27Ac_iPS","promoter_type","CGI","CGIap","CGInap","TATA")]
n5cluster2=n5clustera[,c("ex5_NSC","K27ME3_NSC","K27Ac_NSC","promoter_type","CGI","CGIap","CGInap","TATA")]
n5cluster3=n5clustera[,c("ex5_Neuron","K27ME3_Neuron","K27Ac_Neuron","promoter_type","CGI","CGIap","CGInap","TATA")]
n5cluster1$cell="iPSC"
n5cluster2$cell="NSC"
n5cluster3$cell="Neuron"
colnames(n5cluster1)=c("cluster_count","K27me3","K27ac","promoter_type","CGI","CGIap","CGInap","TATA","cell")
colnames(n5cluster2)=colnames(n5cluster1)
colnames(n5cluster3)=colnames(n5cluster1)
n5cluster1=rbind(n5cluster1,n5cluster2,n5cluster3)
n5cluster1=n5cluster1[which(n5cluster1$promoter_type == "enhancer-like" & n5cluster1$cluster_count>0),]
n5cluster1$Repressed="No"
n5cluster1$Repressed[which(n5cluster1$K27me3 == "Yes" & n5cluster1$K27ac == "No")]="Yes"
n5cluster1$Active="No"
n5cluster1$Active[which(n5cluster1$K27me3 == "No" & n5cluster1$K27ac == "Yes")]="Yes"
n5cluster1$'Co-marked'="No"
n5cluster1$'Co-marked'[which(n5cluster1$K27me3 == "Yes" & n5cluster1$K27ac == "Yes")]="Yes"
n5cluster1$'Un-marked'="No"
n5cluster1$'Un-marked'[which(n5cluster1$K27me3 == "No" & n5cluster1$K27ac == "No")]="Yes"
n5cluster1$group="Un-marked"
n5cluster1$group[which(n5cluster1$K27me3=="Yes")]="Repressed"
n5cluster1$group[which(n5cluster1$K27ac=="Yes")]="Active"
n5cluster1$group[which(n5cluster1$K27me3=="Yes" & n5cluster1$K27ac=="Yes")]="Co-marked"

#for percentage
n5clustera1=reshape2::melt(n5cluster1[,c("cell","group","CGI","CGIap","CGInap","TATA")], id=c(1,2))
n5clusterab1=n5clustera1%>%group_by(cell, group, variable, value)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))
n5clusterab1$value=factor(n5clusterab1$value, levels=c("Yes","No"))
n5clusterab1$label=paste0(n5clusterab1$group,"_",n5clusterab1$value)

n5clusterak=reshape2::melt(n5cluster1[,c("cell","Repressed","Active","Co-marked","Un-marked","CGI","CGIap","CGInap","TATA")], id=c(1:5))
colnames(n5clusterak)[c(6,7)]=c("variable2","value2")
n5clusterak=reshape2::melt(n5clusterak, id=c(1,6:7))
n5clusterak$result=paste0(n5clusterak$value,"_",n5clusterak$value2)
n5clusterak1=n5clusterak%>%group_by(cell, variable, variable2, result)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count))

n5clusterak2=spread(n5clusterak1[,c(1:5)],key=4, value=5)
colnames(n5clusterak2)=c("cell","group","variable","no_no","no_yes","yes_no","yes_yes")
for(i in 1:nrow(n5clusterak2)){
  GSEATasting <- matrix(c(n5clusterak2$no_no[i], n5clusterak2$no_yes[i], n5clusterak2$yes_no[i], n5clusterak2$yes_yes[i]), nrow = 2, dimnames = list(oligo1 = c("no", "yes"), oligo2 = c("no", "yes")))
  n5clusterak2$p.val[i] = fisher.test(GSEATasting, alternative = "two.sided")$p.value
  n5clusterak2$OR[i] = fisher.test(GSEATasting, alternative = "two.sided")$estimate}
n5clusterak2$label="***"
n5clusterak2$label[which(n5clusterak2$p.val>=0.001)]="**"
n5clusterak2$label[which(n5clusterak2$p.val>=0.01)]="*"
n5clusterak2$label[which(n5clusterak2$p.val>=0.05)]=""
n5clusterak2$label1=paste0("OR=",signif(n5clusterak2$OR,2),n5clusterak2$label)
n5clusterab1$group=factor(n5clusterab1$group,levels=c("Un-marked","Repressed","Active","Co-marked"))
n5clusterab1$cell=factor(n5clusterab1$cell, levels=c("iPSC", "NSC", "Neuron"))
n5clusterak2$cell=factor(n5clusterak2$cell, levels=c("iPSC", "NSC", "Neuron"))

s1i=ggplot() + 
  labs(x=NULL, y = NULL ,title= "Enrichment with chromatin state", fill="CGI/ CGIap\nCGInap/ TATA") +
  scale_fill_manual(values=c("black","white"))+
  scale_x_continuous(labels = scales::percent, limits=c(0,1.45), breaks=c(0,0.25,0.5,0.75,1))+
  facet_grid(cols=vars(cell),rows=vars(variable))+
  geom_bar(data=n5clusterab1, mapping=aes(x=percent, y=group, fill=value), linewidth=0.25, color = "black", stat="identity", alpha=1) + 
  geom_text(data=n5clusterak2, mapping=aes(y=group, x=1.01, label=label1), size=1.8, color="black", hjust=0, vjust=0.5)+
  theme1+theme(legend.position=c(0.1,0.5))
pdf(paste0(path_fig2,"s1i.cell_type3_expressing_Repressed_enhacer_ex5.pdf"), width = 4.5, height = 2.3)
print(s1i)
dev.off()














