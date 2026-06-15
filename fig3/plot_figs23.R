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

#======================
#Needed
#fs23
combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1279.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

colnames(combine2)[which(colnames(combine2) %in% c("iPSC_IF","NSC_IF","NRN_IF"))]=c("iPSC","NSC","Neuron")
combine2b=reshape2::melt(combine2[which(combine2$symbol %in% c("AP1S2","CA14","PORCN","PRKCZ")),c("symbol","isoform_id","iPSC","NSC","Neuron")], id=c(1,2))
combine2c=spread(combine2b, key=2, value=4)
combine2c$others=1-rowSums(combine2c[,c(3:10)],na.rm=T)
combine2c=reshape2::melt(combine2c, id=c(1,2))
combine2c=combine2c[which(!is.na(combine2c$value)),]
colnames(combine2c)[c(2,3,4)]=c("Cell_type","isoform_id","fraction")
library(forcats)

for (i in c("AP1S2","CA14","PORCN","PRKCZ")){
  fs2abcd=ggplot(combine2c[which(combine2c$symbol == i),], aes(x=Cell_type, y=fraction, fill=fct_rev(isoform_id))) + 
    scale_fill_manual(values=c("white","#E64B35FF","#3C5488FF"))+
    labs(y ="Expression fraction", x=NULL, title =i,  fill="Isoforms")+
    geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
    theme1
  pdf(paste0(path_fig3,"fs23.Switched transcript.",i,"_expression_fraction.pdf"), width = 2, height = 1.7)
  print(fs2abcd)
  dev.off() }


#======================
#Needed
#fs3
combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1279.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

colnames(combine2)[c(44:46)]=c("iPSC","NSC","Neuron")
combine2b=reshape2::melt(combine2[which(combine2$symbol %in% c("CA14","PRKCZ")),c(28,7,38,44:46)], id=c(1,2,3))
combine2b$label=paste0(combine2b$isoform_id," | ",combine2b$n5_string)
library(forcats)

for (i in c("CA14","PRKCZ")){
  fs3=ggplot(combine2b[which(combine2b$symbol == i),], aes(x=variable, y=value, fill=fct_rev(label))) + 
    scale_fill_manual(values=c("#4DBBD5FF","#00A087FF"))+
    labs(y ="ATAC (CPM)", x=NULL, title =paste0(i," ATAC"),  fill="Isoforms | ex5_cluster")+
    geom_bar(linewidth=0.25, stat="identity", color="black", alpha=0.8)+
    theme1
  pdf(paste0(path_fig3,"fs3.Switched transcript.",i,"_ATACcount.pdf"), width = 2.2, height = 1.7)
  print(fs3)
  dev.off() }

#======================
#Needed
#fs23 - ORF

library(Biostrings)
combine2=read.delim(paste0(path_fig3_data,"isoform.switch.transcript.1279.tsv.gz"), header=T, stringsAsFactors = F, check.names = F)

need=c("AP1S2","CA14","PORCN","PRKCZ")

for (i in 1:length(need)){
  combine3=combine2[which(combine2$symbol == need[i]),]
  seq1 <- AAString(combine3$ORF_seq[1])
  seq2 <- AAString(combine3$ORF_seq[2])
  align <- Biostrings::pairwiseAlignment(seq1, seq2, type = "global")
  pid(align) #82.2%
  
  a1 <- as.character(alignedPattern(align))
  a2 <- as.character(alignedSubject(align))
  v1 <- strsplit(a1, "")[[1]]
  v2 <- strsplit(a2, "")[[1]]
  mid <- ifelse(v1 == v2 & v1 != "-", "|", " ")
  df <- data.frame(pos = seq_along(v1), seq1 = v1, mid = mid, seq2 = v2, stringsAsFactors = FALSE)
  df2 <- df %>% mutate(status = case_when(seq1 == "-" | seq2 == "-" ~ "gap", seq1 == seq2 ~ "match", TRUE ~ "mismatch"))
  
  df_long <- bind_rows(
    df2 %>% transmute(pos, track = combine3$isoform_id[1], value = seq1),
    df2 %>% transmute(pos, track = "Match", value = status),
    df2 %>% transmute(pos, track = combine3$isoform_id[2], value = seq2))
  
  df_long <- df_long %>%mutate(fill_group = case_when(
    track == "Match" & value == "match" ~ "match",
    track == "Match" & value == "mismatch" ~ "mismatch",
    track == "Match" & value == "gap" ~ "gap",
    track != "Match" & value == "-" ~ "gap",
    TRUE ~ "residue"))
  
  df_long$track=factor(df_long$track, levels=c(combine3$isoform_id[1],"Match",combine3$isoform_id[2]))
  p <- ggplot(df_long, aes(x = pos, y = track, fill = fill_group)) +
    geom_tile_rast(height = 0.8, alpha=1,color = NA) +
    scale_fill_manual(values = c(match = "grey20",mismatch = "firebrick3",gap = "goldenrod2",residue = "steelblue")) +
    labs(x = "ORF Position (aa)", y = NULL, fill = NULL) +
    theme1 +theme(panel.grid = element_blank())
  pdf(paste0(path_fig3,"fs2e.ORF.alignment.",need[i],".pdf"), width = 2.1, height = 0.8)
  print(p)
  dev.off()}

#====



