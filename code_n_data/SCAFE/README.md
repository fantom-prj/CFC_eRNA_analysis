# SCAFE on CFC-seq and CAGE

* This section contain the codes for 4 versions of SCAFE run.
* The codes are located in each of the sub-folder
* Results are used in ext_Fig.2 & 3
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [TSS clusters and tCREs identification by SCAFE](#SCAFE)



# <a name="SCAFE"></a>TSS clusters and tCREs identification by SCAFE

The BAM files containing all alignments mapped to the hg38 genome were used as the input for SCAFE v1.01. The functions acquired from SCAFE were incorporated into SALA, where detailed workflow and explanation can be found in the GitHub (https://github.com/fantom-prj/SALA). Briefly, scafe.workflow.bk.bam_to_ctss extracted alignments with ≤ 3 nt softclip at the 5’ end and ≥ 30 matched nt, where ~70% reads passed. Approximately 30% of reads from the CFC-seq had a relatively long 5’ end unaligned region (> 3 nt), these reads were excluded from TSS cluster identification while their reads were used for quantification. The qualified alignments (~165 million) were converted into CTSS, which were grouped into TSS clusters. According to the multiple properties of each TSS cluster (percentage of CTSS with unencoded G, corrected expression, read count, summit count and flanking count), a multiple logistic regression model was trained to distinguish TSS clusters that are likely genuine. The identification of these TSS clusters was performed per replicate where cutoffs of at least 3 reads and at least 1 read with unencoded G were applied. Next, scafe.tool.cm.aggregate defines tCREs by extending and merging overlapped TSS clusters defined from the replicates. This aggregation step only merges the clusters without pooling signals to generate new TSS clusters. The TSS clusters were then filtered by scafe.tool.cm.filter according to the multiple logistic regression model. Additional commands including scafe.tool.cm.directionality and scafe.tool.cm.annotate were used to calculate bi-directional reads in each tCRE and pre-annotate the tCREs with nearby GENCODE transcript model. The codes for running SCAFE for CFC-seq in this study are listed:
```
#Convert bam into CTSS with and without unencoded G:
[1] scafe.tool.bk.bam_to_ctss
--TSS_mode=softclip \
--bamPath=[bamPath] \
--unencoded_G_upstrm_nt=3 \
--max_thread=5 \
--genome=hg38.gencode_v39 \
--max_softclip_length=3 \
--outputPrefix=[outputPrefix] \
--outDir=[outDir]

#Aggregate individual library into cell types [2-4]:
[2] scafe.tool.cm.aggregate \
--lib_list_path=[path to CTSS per library ] \
--max_thread=5 \
--genome=hg38.gencode_v39 \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/aggregate

[3] scafe.tool.cm.ctss_to_bigwig \
--genome=hg38.gencode_v39 \
--ctss_bed_path=aggregate.collapse.ctss.bed.gz \
--outputPrefix=[outputPrefix.all] \
--outDir=$baseDir/out/ctss_to_bigwig

[4] scafe.tool.cm.ctss_to_bigwig \
--genome=hg38.gencode_v39 \
--ctss_bed_path=aggregate.unencoded_G.collapse.ctss.bed.gz \
--outputPrefix=[outputPrefix.ung] \
--outDir=$baseDir/out/ctss_to_bigwig

#Identify TSS clusters, tCREs and annotate the tCREs [5-9]:
[5] scafe.tool.cm.aggregate \
--lib_list_path=[path to CTSS per library ] \
--max_thread=5 \
--genome=hg38.gencode_v39 \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/aggregate

[6] scafe.tool.cm.cluster \
--overwrite=yes \
--cluster_ctss_bed_path=aggregate.collapse.ctss.bed.gz \
--count_ctss_bed_path=aggregate.unencoded_G.collapse.ctss.bed.gz \
--min_summit_count=0 \
--min_cluster_count=1 \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/cluster

[7] scafe.tool.cm.filter \
--overwrite=yes \
--ctss_bed_path=aggregate.collapse.ctss.bed.gz \
--ung_ctss_bed_path=aggregate.unencoded_G.collapse.ctss.bed.gz \
--tssCluster_bed_path=tssCluster.bed.gz \
--genome=hg38.gencode_v39 \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/filter

[8] scafe.tool.cm.annotate \
--overwrite=yes \
--tssCluster_bed_path=tssCluster.default.filtered.bed.gz \
--tssCluster_info_path=tssCluster.log.tsv \
--min_CRE_count=3 \
--genome=hg38.gencode_v39 \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/annotate

[9] scafe.tool.cm.directionality \
--overwrite=yes \
--CRE_bed_path=CRE.coord.bed.gz \
--CRE_info_path=CRE.info.tsv.gz \
--ctss_bed_path=aggregate.collapse.ctss.bed.gz \
--outputPrefix=[outputPrefix] \
--outDir=$baseDir/out/directionality
```

