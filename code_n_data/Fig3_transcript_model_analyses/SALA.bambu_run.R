# List of required packages
required_packages <- c("bambu","edgeR")


# Function to check and install missing packages
install_if_missing <- function(packages) {
  missing_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if (length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse=", "))
    devtools::install_github("GoekeLab/bambu", ref = "test_split_read_classes")
  }
}

# Ensure BiocManager is installed for Bioconductor packages
if (!requireNamespace("BiocManager", quietly=TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org", dependencies = TRUE)
}

# Install missing packages
install_if_missing(required_packages)

# load packages
suppressPackageStartupMessages(library(bambu))
library(edgeR)

#bam_path <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/data/SCAFE.step1.bam.list.txt"
#fa.file <- "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/git_folder/resources/chr17_chr18.fasta"
#please download bam files from DDBJ

#primary_folder=[primary_folder]
primary_folder="/osc-fs_home/yip/CFC_seq_paper_fig_data/"
SALA_gtf <- paste0(primary_folder,"code_n_data/SALA/Neuron_THP1_full/transcript/all_gtf_file/table5.final.partial_yes_detected.alone.gtf.gz")
results_dir <- paste0(primary_folder,"code_n_data/Fig3_transcript_model_analyses/bambu_long_t5_partialYes.ENST/")


#=====
#run bambu
print("starting running bambu...")
bam <- read.delim(bam_path, header = F, stringsAsFactors = F, check.names = F)
test.bam <- bam$V3
annotations <- prepareAnnotations(SALA_gtf)
ID_link=unique(data.frame(cbind(mcols(annotations)$GENEID, bambu:::assignGeneIds(annotations, GRangesList())$GENEID)))
ID_link=ID_link%>%group_by(X2)%>%dplyr::summarise(geneID=paste(X1,collapse=";"))
mcols(annotations)$GENEID <- bambu:::assignGeneIds(annotations, GRangesList())$GENEID
dir.create(paste0(results_dir,"/rcOut"), recursive=TRUE)
se <- bambu(reads =  test.bam, 
                    annotations = annotations, 
                    genome = fa.file, 
                    discovery = FALSE, 
                    opt.discovery = list(min.exonDistance = 0), 
                    rcOutDir = paste0(results_dir,"/rcOut"),returnDistTable=TRUE)
saveRDS(se, paste0(results_dir,"/se.rds"))
writeBambuOutput(se, results_dir)
print(paste0("finish running bambu. results located in ",results_dir))

#=====
gene_count <- read.delim(paste0(results_dir,"/counts_gene.txt.gz"), header=T)
transcript_count <- read.delim(paste0(results_dir,"/counts_transcript.txt.gz"), header=T)
transcript_CMP <- read.delim(paste0(results_dir,"/CPM_transcript.txt.gz"), header=T)
gene_count <- left_join(gene_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_count <- left_join(transcript_count, ID_link, by=c("GENEID"="X2"), copy=F)
transcript_CMP <- left_join(transcript_CMP, ID_link, by=c("GENEID"="X2"), copy=F)

write.table(gene_count, gzfile(paste0(results_dir,"/gene.count.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)
write.table(transcript_count, gzfile(paste0(results_dir,"/transcript.count.matrix.tsv.gz")), col.names=T, row.names=F, sep="\t", quote=F)



