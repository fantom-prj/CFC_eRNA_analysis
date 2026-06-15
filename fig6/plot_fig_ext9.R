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
library(dplyr) 
library(ggrastr)

###color#####
library("ggsci")
mypal = pal_npg("nrc", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

#============
#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
path_fig6=paste0(primary_folder,"fig6/out/")
path_fig6_data=paste0(primary_folder,"fig6/data/")
setwd(path_fig6_data)

#============
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

#===============================================================================
#ext9a

data2a=read.delim(paste0(path_fig6_data,"TE_group_CpGTATA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data3a=data2a%>%group_by(ex5cluster_class,CpGTATA,group_ex5cluster)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="ex5_cluster")
data4a=data2a%>%group_by(ex5cluster_class,CpGTATA,group_exon)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="exon")
colnames(data3a)[3]="REgroup"
colnames(data4a)[3]="REgroup"
data5a=rbind(data3a, data4a)
data5a$CpGTATA=factor(data5a$CpGTATA, levels=c("CGIap","CGInap","Null","TATA","Others"))
data5a$REgroup=factor(data5a$REgroup, levels=c("DNA","LTR","LINE","SINE","Low_complexity","Satellite","Simple_repeat", "Null", "others"))
data6a=unique(data5a[,c(1,2,6,7)])

ex9a=ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15))+
  scale_fill_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "white","grey"))+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Proportion of enhancers overlapping\nwith repetitive elements",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("e_ncRNA") & data5a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=percent, fill=REgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("e_ncRNA") & data6a$CpGTATA != "Others"),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9a.Presence_of_RE_inside_ex5_cluster_exon.pdf"), width = 2.5, height = 1.7)
print(ex9a)
dev.off() 

#==================
#Needed
#ex9b
data2a=read.delim(paste0(path_fig6_data,"family_of_LTR.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data3a=data2a[which(data2a$group_ex5cluster == "LTR" & data2a$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CpGTATA,V15_ex5cluster)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="ex5_cluster")
data4a=data2a[which(data2a$group_exon == "LTR" & data2a$ex5cluster_class == "e_ncRNA"),]%>%group_by(ex5cluster_class,CpGTATA,V15_exon)%>%dplyr::summarise(count=n())%>%dplyr::mutate(percent=count/sum(count), total=sum(count), group="exon")
colnames(data3a)[c(3)]=c("REsubgroup")
colnames(data4a)[c(3)]=c("REsubgroup")
data5a=rbind(data3a, data4a)
data5a=data5a[which(data5a$CpGTATA != "Others"),]
data5a2=data5a%>%group_by(REsubgroup)%>%dplyr::summarise(count=sum(count))
notneed=data5a2$REsubgroup[which(data5a2$count<75)]
data5a$REsubgroup[which(data5a$REsubgroup %in% notneed)]="Others"
data5a$CpGTATA=factor(data5a$CpGTATA, levels=c("CGIap","CGInap","Null","TATA"))
data6a=unique(data5a[,c(1,2,6,7)])

ex9b=ggplot() + 
  scale_y_continuous(labels = scales::percent, limits=c(0, 1.15), breaks=c(0,0.25,0.5,0.75,1))+
  scale_fill_npg()+
  labs(x=NULL , y ="% of ex5_cluster" , title ="Family of LTR inside ex5_cluster \nor transcript model of enhancers",  fill=NULL)+
  facet_wrap(vars(group), ncol=2,  scales="fixed")+
  geom_bar(data=data5a[which(data5a$ex5cluster_class %in% c("e_ncRNA")),], mapping=aes(x=CpGTATA, y=percent, fill=REsubgroup), alpha=1, linewidth=0.25, stat="identity", color="black")+
  geom_text(data=data6a[which(data6a$ex5cluster_class %in% c("e_ncRNA")),], mapping=aes(x=CpGTATA, y=1, label=total), angle=25, hjust=(0), vjust=(0), nudge_x = (-0.25), size=1.8, color="black")+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9b.subfamily_of_LTR_ex5_cluster_exon.pdf"), width = 2.4, height = 1.7)
print(ex9b)
dev.off() 

#============
#ex9c
data8=read.delim(paste0(path_fig6_data,"CpGTATA_repeat_element_FE.n5_exon.tsv.gz"), header = T, stringsAsFactors = F, check.names = F)

data8$feature1=factor(data8$feature1, levels=c("TATA","Null","CGI"))
data8$sig_level=factor(data8$sig_level, levels=c("ns","*","**","***"))
data8=data8[which(data8$gene_group %in% c("p_ncRNA","other_ncRNA")),]
data8$gene_group=factor(data8$gene_group, levels=c("p_ncRNA","other_ncRNA"))

ex9c=ggplot() + 
  labs(x=NULL, y = "n5_cluster group",title= "Enrichment with LTR") +
  facet_grid(rows=vars(group))+
  geom_point(data=data8[which(data8$feature2 == "LTR"),], mapping=aes(y=feature1, x=gene_group, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high = "#DC0000FF", midpoint = 0)+ 
  scale_size_manual(values=c(0.2,1.8))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9c.all_RNA_FE_CpGTATA_LTR_ex5cluster_exon.pdf"), width = 1.6, height = 1.7)
print(ex9c)
dev.off() 

#==================
#Needed
#ex9d
exon_repeat13=read.delim(paste0(path_fig6_data,"RE_ex5_cluster.exon_coverage_nexon.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

exon_repeat13$CpGTATA[grep("CGI",exon_repeat13$CpGTATA)]="CGI"
exon_repeat13=exon_repeat13[which(exon_repeat13$CpGTATA %in% c("CGI","Null","TATA")),]
exon_repeat13=exon_repeat13[which(exon_repeat13$ex5cluster_class %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),]
exon_repeat13$ex5cluster_class=factor(exon_repeat13$ex5cluster_class, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
exon_repeat13$CpGTATA=factor(exon_repeat13$CpGTATA, levels=c("CGI","Null","TATA"))
aaaa=exon_repeat13%>%group_by(group,ex5cluster_class,CpGTATA)%>%dplyr::summarise(median=median(overlap_percent),hit_2exon=length(which(hit_exon>1)), hit_count=n())
aaaa$percent_2exon=aaaa$hit_2exon/aaaa$hit_count
aaaa$ex5cluster_class=factor(aaaa$ex5cluster_class, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))

ex9d=ggplot() + 
  scale_fill_manual(values=c("#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  labs(x=NULL, y ="Sequence coverage (%)" , title ="LTR coverage on\ntranscript model exon",  fill=NULL)+
  facet_wrap(vars(ex5cluster_class),ncol=2, nrow=2, scales="free")+
  scale_y_continuous(labels = scales::percent, breaks=c(0, 0.25,0.5,0.75,1), limits=c(0,1.2))+
  geom_boxplot(data=exon_repeat13[which(exon_repeat13$group == "LTR"),], mapping=aes(x=CpGTATA, y=overlap_percent, fill=CpGTATA), linewidth=0.25, color = "black", outlier.colour="grey25", outlier.size=3.5, notch = FALSE, width = 0.45, outlier.shape = NA, alpha=0.8) + 
  geom_text(data=aaaa[which(aaaa$group == "LTR"),], mapping=aes(y=median, x=CpGTATA, label=paste0(signif(median,3)*100,"%")),angle=90,size=1.8, hjust=0.5, vjust=(-1.8))+
  geom_text(data=aaaa[which(aaaa$group == "LTR"),], mapping=aes(y=1, x=CpGTATA, label=paste0(signif(percent_2exon,3)*100,"%")),size=1.8, hjust=0.5, vjust=0.5)+
  annotate("text", y=1.15, x=0.1, label="% >= 2 exons:", hjust=0, size=1.8)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9d.RE_ex5_cluster.exon_coverage.pdf"), width = 2.4, height = 3)
print(ex9d)
dev.off() 

#============
#ex9e
library(ComplexHeatmap)
library(circlize)
library(grid)
RLE.ex5=read.delim(paste0(path_fig6_data,"ex5_cluster_Neuron_THP1.RLE.tsv.gz"),header=T, stringsAsFactors = F, check.names = F)
data2a=read.delim(paste0(path_fig6_data,"TE_group_CpGTATA.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
data2a$CpGTATA[grep("CGI",data2a$CpGTATA)]="CGI"
data2a=data2a[which(data2a$ex5cluster_class %in% c("e_ncRNA") & data2a$ CpGTATA %in% c("CGI","Null","TATA")),]

RLE.ex5=log10(RLE.ex5+0.01)
RLE.ex5=data.frame(RLE.ex5)
RLE.ex5[RLE.ex5>3]=3

RLE.ex5t=RLE.ex5[which(rownames(RLE.ex5)%in% data2a$n5_string[which(data2a$CpGTATA=="TATA")]),]
RLE.ex5c=RLE.ex5[which(rownames(RLE.ex5)%in% data2a$n5_string[which(data2a$CpGTATA=="CGI")]),]
RLE.ex5n=RLE.ex5[which(rownames(RLE.ex5)%in% data2a$n5_string[which(data2a$CpGTATA=="Null")]),]

col_fun2 <- colorRamp2(c(min(RLE.ex5c), max(RLE.ex5c)), c("white", "#00A087FF"))
col_fun3 <- colorRamp2(c(min(RLE.ex5n), max(RLE.ex5n)), c("white", "black"))
col_fun4 <- colorRamp2(c(min(RLE.ex5t), max(RLE.ex5t)), c("white", "#3C5488FF"))

out2 <- Heatmap(as.matrix(RLE.ex5c),
                name = "value",
                col = col_fun2,
                cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE, row_dend_gp = gpar(lwd = 0.3),
                clustering_method_rows = "ward.D2",
                row_dend_width = unit(10, "mm"))
out3 <- Heatmap(as.matrix(RLE.ex5n),
                name = "value",
                col = col_fun3,
                cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE, row_dend_gp = gpar(lwd = 0.3),
                clustering_method_rows = "ward.D2",
                row_dend_width = unit(10, "mm"))
out4 <- Heatmap(as.matrix(RLE.ex5t),
                name = "value",
                col = col_fun4,
                cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
                use_raster = TRUE, row_dend_gp = gpar(lwd = 0.3),
                clustering_method_rows = "ward.D2",
                row_dend_width = unit(10, "mm"))

out2@matrix_param$height <- unit(nrow(RLE.ex5c), "null")
out3@matrix_param$height <- unit(nrow(RLE.ex5n), "null")
out4@matrix_param$height <- unit(nrow(RLE.ex5t), "null")
pdf(paste0(path_fig6,"ex9e.CpGTATA_ex5cluster_heatmap.pdf"), width = 1.8, height = 2)
draw(out2 %v% out4)
dev.off()
pdf(paste0(path_fig6,"ex9e.CpGTATA_ex5cluster_heatmap_null.pdf"), width = 1.8, height = 2)
draw(out3)
dev.off()

#============
#ex9f
meme1b=read.delim(paste0(path_fig6_data,"meme_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)
meme1b=meme1b[,-which(colnames(meme1b)%in% c("e_ncRNA_CGIap","e_ncRNA_CGInap"))]
content=read.delim(paste0(path_fig6_data,"meme_result.count.tsv"), header=T, stringsAsFactors = F, check.names = F)
content$set=gsub("RNA_","RNA.",content$set)
content$class=sapply(strsplit(content$set,"\\."),"[",1)
content$CpGTATA=sapply(strsplit(content$set,"\\."),"[",2)

dist_rows <- dist(meme1b[,-1])                  
hclust_rows <- hclust(dist_rows)         

meme1b2 <- reshape2::melt(meme1b, id=1)
meme1b2$ALT_ID=factor(meme1b2$ALT_ID, levels=rownames(meme1b)[hclust_rows$order])
meme1b2$variable=gsub("RNA_","RNA.",as.character(meme1b2$variable))
meme1b2$class=sapply(strsplit(as.character(meme1b2$variable),"\\."),"[",1)
meme1b2$CpGTATA=sapply(strsplit(as.character(meme1b2$variable),"\\."),"[",2)
meme1b2$class=factor(meme1b2$class,levels=c("p_ncRNA","e_ncRNA","other_ncRNA"))
meme1b2$CpGTATA[which(meme1b2$CpGTATA=="TATALTR")]="TATA-LTR"
meme1b2$CpGTATA=factor(meme1b2$CpGTATA,levels=c("CGI","Null","TATA","TATA-LTR"))
content=content[which(content$set %in% meme1b2$variable),]
content$CpGTATA[which(content$CpGTATA=="TATALTR")]="TATA-LTR"
content$class=factor(content$class,levels=c("p_ncRNA","e_ncRNA","other_ncRNA"))

ex9f=ggplot() +
  geom_tile_rast(data=meme1b2, aes(x=CpGTATA, y=ALT_ID, fill = value),color = "black", linewidth=0.1)+
  geom_text(data=content, aes(x=CpGTATA, y=28.5, label=paste0("n=",count)), size=1.8, angle=90, hjust=0, vjust=0)+
  facet_grid(cols=vars(class), scale="free_x", space="free_x")+
  labs(x = NULL, y = "Enriched TF motifs", fill = "ENR\nratio", title="TF motif enrichment")+
  scale_y_discrete(limits = rev, expand = expansion(mult = c(0, 0.18)))+
  scale_fill_gradientn(colors = c("white", "firebrick3", "navy"), values = scales::rescale(c(0, 15, 100)))+
  theme1 +
  theme(panel.background = element_blank(),panel.grid.major = element_blank(), 
        axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "bottom",
        text = element_text(size=6), axis.text.x = element_text(angle = 25, hjust = 1, vjust=1))
pdf(paste0(path_fig6,"ex9f.TFBS.heatmap.pdf"), width = 2.4, height = 3)
print(ex9f)
dev.off()

#===========
#ex9g
data=read.delim(paste0(path_fig6_data,"ex5_cluster_NFY_distribution_result.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data$anno_region=factor(data$anno_region, levels=c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA"))
data=data[which(data$CpGTATA %in% c("CGI","Null","TATA")),]
data$CpGTATA=factor(data$CpGTATA, levels=c("CGI","Null","TATA"))
data1=data%>%group_by(anno_region,CpGTATA)%>%dplyr::summarise(size=unique(size))
data1$label=paste0(data1$CpGTATA,": ",data1$size)
ex9g=ggplot()+
  scale_color_manual(values=c("#00A087FF", "grey", "#3C5488FF"), guide=NULL)+
  coord_cartesian(xlim=c(-1500,1500), ylim=c(0,0.35))+
  geom_line(data = data[which(data$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=V4, y=V5, color=CpGTATA, group=CpGTATA), linewidth=0.25)+
  geom_text(data=data1[which(data1$CpGTATA == "CGI" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-1500, y=0.25, label=label), vjust=0, hjust=0, size=1.8, color="#00A087FF")+
  geom_text(data=data1[which(data1$CpGTATA == "Null" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-1500, y=0.29, label=label), vjust=0, hjust=0, size=1.8, color="grey")+
  geom_text(data=data1[which(data1$CpGTATA == "TATA" & data1$anno_region %in% c("mRNA","p_ncRNA","e_ncRNA","other_ncRNA")),], aes(x=-1500, y=0.33, label=label), vjust=0, hjust=0, size=1.8, color="#3C5488FF")+
  facet_grid(cols=vars(anno_region))+
  labs(color=NULL, x="Distance from the stranded ex5_cluster summit", y="% of ex5_cluster", title="NF-Y binding with regulatory\nelements (Expressed in iPSC)")+
  scale_y_continuous(labels = scales::percent)+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9g.NFY.distribution.pdf"), width = 2.6, height = 1.7)
print(ex9g)
dev.off() 

#============
#ex9h
data3a=read.delim(paste0(path_fig6_data,"NFY_TEAD_LTR.n5_FE.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

data3a$CpGTATA=factor(data3a$CpGTATA, levels=c("CGI","Null","TATA"))
data3a$sig_level=factor(data3a$sig_level, levels=c("ns","*","**","***"))
data3a$ex5cluster_class=factor(data3a$ex5cluster_class, levels=c("other_ncRNA","e_ncRNA","p_ncRNA","mRNA"))

ex9h=ggplot() + 
  labs(x=NULL, y = NULL, title= "LTR Enrichment with\nNF-Y & TEAD binding") +
  facet_grid(cols=vars(group))+
  geom_point(data=data3a, mapping=aes(y=ex5cluster_class, x=CpGTATA, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradientn(colors = c("#3C5488FF","white","#DC0000FF", "#4b0000","#4b0000"), values = scales::rescale(c(-1, 0, 2.66,5, Inf)))+
  scale_size_manual(values=c(0.2,1.8,2.7,3.6))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9h.NFY.TEAD_LTR.n5_FE.pdf"), width = 1.9, height = 1.7)
print(ex9h)
dev.off() 


#==========
#ex9i

data4a=read.delim(paste0(path_fig6_data,"NFY_TEAD_K27_FE.tsv.gz"),  header=T, stringsAsFactors = F, check.names = F)

data4a=data4a[which(data4a$CpGTATA %in% c("TATA","Null")),]
data4a$CpGTATA=factor(data4a$CpGTATA, levels=c("TATA","Null"))
data4a$sig_level=factor(data4a$sig_level, levels=c("ns","*","**","***"))
data4a$group=factor(data4a$group, levels=c("Un-marked","Repressed","Active","Co-marked"))

ex9i=ggplot() + 
  labs(x="LTR enhancers", y =NULL, title= "Enrichment between H3K27\nmodification and TF binding") +
  facet_grid(cols=vars(group2))+
  geom_point(data=data4a[which(data4a$ex5cluster_class=="e_ncRNA"),], mapping=aes(y=group, x=CpGTATA, fill=FE_logOR, size=sig_level),shape=21)+
  scale_fill_gradient2(low = "#3C5488FF", mid = "white", high ="#DC0000FF" , midpoint = 0)+ 
  scale_size_manual(values=c(0.2,1.8,2.7,3.6))+
  theme1+theme(axis.text.x = element_text(color ="black", angle=25, hjust=1, vjust=1))
pdf(paste0(path_fig6,"ex9i.NFY.TEAD_K27_FE.pdf"), width = 1.9, height = 1.7)
print(ex9i)
dev.off() 



