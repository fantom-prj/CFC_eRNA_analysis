# Transcriptome assembling by SALA

* SALA_Final.R -> code describing processing of SALA Final (and SALA Raw)
* Codes inside the folder ./Neuron_THP1_full showed each single step on SALA Final
* SALA_iPSchro.R -> code describing processing of SALA on chromatin-bound iPSC CFC-seq dataset 
* Codes inside the folder ./iPSC_chromatin_bound showed each single step on SALA chromatin-bound
* SALA_Default.R -> code describing processing of SALA Default
* Codes inside the folder ./Neuron_THP1_default showed each single step on SALA Default
* gtf_builder.R & gtf_builder_iPSchro.R -> associated code to generate gtf files from the transcriptome

* Results are presented in Fig.3 and ext_Fig.4
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC
* All the code described here (except lncRNA sub-classification) have been incorporated into SALA software: https://github.com/fantom-prj/SALA

# Related Methods
* [Transcript model construction using Transcript Start-site Aware Long-read Assembler (SALA)](#sala)
	* [Feature collection](#feature)
	* [Transcript model assembly and read assignment](#transcript)
	* [Initial gene annotation](#gene0)
	* [Filtering for confident transcript models](#filter)
	* [Coding potential of transcript and gene models](#CPAT)
	* [Construction of different GTF files](#GTF)
	* [SALA annotation on Chromatin-bound iPSC dataset](#chromatin)
* [Sub-classification of lncRNA transcripts and genes](#subclass)


# <a name="sala"></a>Transcript model construction using Transcript Start-site Aware Long-read Assembler (SALA)

Transcript model construction using Transcript Start-site Aware Long-read Assembler (SALA)
SALA reconstructs transcriptomes through four sequential modules: feature collection, transcript model assembly, transcript model filtering and gene annotation. The steps were documented in the wiki of SALA GitHub. The parameters used for SALA Final are the same as the “sensitive mode”, while SALA Default is from the “default mode” in GitHub. Notably, these two runs are independent and generate transcript and gene models independently. SALA Final is the primary dataset for downstream analyses and was used to build the SACAGE transcriptome of FANTOM6 Interactome by adding all GENCODE (v39) transcript models and FANTOM CAT permissive transcript models. Therefore, transcript and gene IDs from SALA Final are transferrable with SACAGE. Detailed description of SALA can be found in the GitHub wiki page (https://github.com/fantom-prj/SALA/wiki).

## <a name="feature"></a>Feature collection

SALA integrates multiple evidence to define high-confidence transcript boundaries and internal structures. 5’ end clusters and summits identified via SCAFE were utilized as primary transcription start sites (TSS). For 3’ ends, tags were collapsed into clusters using paraclu (minimum 3 tags per cluster) to define termination signal clusters and summits. Orthogonal evidence was incorporated by including SJs from sample-matched short-read RNA-seq (SJ.out.tab from STAR) and reference annotations (GENCODE v39). Splice junctions (SJs) were extracted from alignments before TranscriptClean, alongside maximum base-calling quality scores for the flanking sequences (± 3nt). For SJs without external support, we applied a threshold of a summarized base-calling score ≥ 10 and a read count ≥ 3, for technical confidence (Ext_Fig. 7b). To ensure robust boundary definition, SCAFE-defined TSS clusters were merged within a 75-nt window and resulted as extended 5’ end clusters (ex5_clusters). This effectively included GENCODE transcript model 5’ ends and linked annotated transcript models to the reads in the same ex5_clusters. Bimodal clusters (exhibiting > 10% signal across the distribution) were split at their midpoint to distinguish proximal ex5_clusters. 3’ end clusters were merged within a 150-nt window. Clusters supported by either SCAFE, paraclu cluster, or GENCODE 5’ and 3’ ends annotation were designated as confident features. In feature collection, SALA Final used the same parameter as SALA Default.
```
# Convert read bam file into bed file
[1] sh ./SALA/code/others/SALA.input.bamtobed.sh \
[bamTC_path.txt (path of bam files after transcriptclean)] \
[./SALA/input/bam_to_bed (output_directory)] \
./SALA/resources

# Prepare 3' end cluster
[2] perl ./SALA/code/SALA/3n_cluster/transcript_bed_to_end_bed_bigwig.pl \
./SALA/input/bam_to_bed/combined.bed.bgz \
./SALA/resources/chrom.sizes.tsv \
[outputPrefix] \
[./SALA/input/CTES_clusters/end3_bed_bigwig (output_directory)] \
./SALA/resources/bin/bedGraphToBigWig/bedGraphToBigWig

[3] perl ./SALA/code/SCAFEv1.0.1/scripts/scafe.tool.cm.cluster \
--overwrite=yes \
--cluster_ctss_bed_path=./SALA/input/CTES_clusters/end3_bed_bigwig/outputPrefix.end3.bed \
--count_ctss_bed_path=./SALA/input/CTES_clusters/end3_bed_bigwig/outputPrefix.end3.bed \
--min_summit_count=3 \
--min_nt_count=3 \
--min_cluster_count=5 \
--outputPrefix=[outputPrefix.CTES.s3_n3_c5] \
--outDir=./SALA/input/CTES_clusters/scafe/cluster

# Prepare splice junction
[4] perl ./SALA /code/SALA/junction_extractor/junction_extractor.pl \
--in_bam=[iPSC_rep1_run1_subset.sorted.bam(bam files before transcriptclean)] \
--chrom_size_path=./resources/chrom.sizes.tsv \
--chrom_fasta_path=./resources/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
--out_prefix=[iPSC_rep1_run1] \
--out_dir=./SALA/input/junction_extractor/output \
--max_thread=1 \
--min_nt_qual=10 \
--min_MAPQ=20 \
--samtools_bin=./resources/bin/samtools/samtools \
--bedtools_bin=./resources/bin/bedtools/bedtools \
--tabix_bin=./resources/bin/tabix/tabix \
--bgzip_bin=./resources/bin/bgzip/bgzip

[5] perl ./SALA/code/SALA/junction_extractor/junction_pool.pl \
--outDir=./SALA/input/junction_extractor/pool \
--outTag=[outputPrefix] \
--findStr=./SALA/input/junction_extractor/output/*/log/*.junct.info.tsv.gz
```

## <a name="transcript"></a>Transcript model assembly and read assignment

Each long-read was assigned a specific 5’ end cluster, a 3’ end cluster, and a set of SJ IDs. Using the default parameters, confident ex5_clusters were derived from SCAFE requiring at least 3 reads and at least one un-encoded G while confident ex3_clusters require ≥ 5 counts and ≥ 3 summit counts. If a GENCODE annotated TSS and TES is located inside the ex5_clusters and ex3_clusters, these clusters become confident independent of read count. Reads possessing both confident 5’ and 3’ boundaries, as well as all reference-matching models, were classified as complete transcript models. Reads sharing an identical triplet of features (TSS cluster, 3’ cluster, and SJ chain) were collapsed into a unique isoform. To prevent the erroneous identification of nascent or fragmented RNAs as novel isoforms, SALA utilized a hierarchical assignment strategy. Incomplete long-reads (lacking one or both confident boundaries) were first mapped to existing complete models. If an incomplete read’s structure was a subset of a complete model, it was assigned as partial support. Only remaining incomplete reads without a matching full-length scaffold were permitted to form independent models. The final coordinates for a model’s 5’ and 3’ ends were defined by the most frequent terminal positions among its constituent reads. In this step, the parameter --trnscpt_set_end_priority=commonest:summit:longest was used for SALA Final while for SALA Default, --trnscpt_set_end_priority=summit:commonest:longest was used. This setup in SALA Final enabled the observation of individual transcript 5’ ends when multiple transcript models derived from the same ex5_clusters. Notably, transcript model assignment is not affected by this parameter. The parameters used for SALA Final is shown as below:
```
[1] end5_guided_assembler_v0.1.20231102.pl
--qry_bed_bgz=[Neuron_THP1.bed.bgz]
--ref_bed_bgz=[GENCODEv39.transcript.bed.bgz]
--chrom_size_path=[hg38.gencode_v39_chrom.sizes.tsv]
--out_dir=[output]
--max_thread=3
--out_prefix=Neuron_THP1
--min_transcript_length=15
--doubtful_end_avoid_summit=yes
--min_exon_length=1
--print_trnscrptID=no
--chrom_fasta_path=[hg38.gencode_v39_genome.fa]
--min_output_qry_count=1
--trnscpt_set_end_priority=commonest:summit:longest #for SALA Final
--doubtful_end_merge_dist=150
--novel_model_prefix=ONTT
--conf_end3_merge_flank=150
--conf_end5_merge_flank=75
--conf_end5_bed_bgz=[ontCAGE.Neuron_THP1.end5.cluster.bed.bgz]
--conf_end3_bed_bgz=[ontCAGE.Neuron_THP1.end3.cluster.bed.bgz]
--min_summit_dist_split=50
--retain_no_qry_ref_bound_set=no
--doubtful_end_avoid_summit=yes
--min_size_split=100
--min_frac_split=0.2
--signal_end5_bed_bgz=[ontCAGE.Neuron_THP1.end5.signal.bed.bgz]
--signal_end3_bed_bgz=[ontCAGE.Neuron_THP1.end3.signal.bed.bgz]
--conf_end3_add_ref=yes
--conf_end5_add_ref=yes
--conf_junction_bed=[Gencode_v39.junct.bed], [long_read.hi_qual.junct.bed] ,[short_read.hi_qual.junct.bed]
--min_qry_score=0
```

## <a name="gene0"></a>Initial gene annotation

The transcript model assembly generates the SALA Raw transcript models, which include all the reference transcripts (even they are not detectable) and novel transcripts. These models were then subjected to initial gene annotation using the parameter --disable_ref_chain_bound_gene_anno=yes. This constraint is critical to prevent the erroneous fusion of distinct reference gene models based on the detection of low-level, nascent read-through transcripts. Initial gene annotation is only used for setting thresholds on transcript read count ratio per gene. Another round of gene annotation was performed on the filtered transcript models.

## <a name="filter"></a>Filtering for confident transcript models

We include different levels of filtering for benchmarking. Filtering criteria and the resulting number of transcript models are shown in Extended Figure 4a,d. In SALA Raw and SALA Final were obtained from a run using SALA sensitive. While the Raw dataset retained all the models excluding internal priming and undetected GENCODE models, the Final dataset was further subjected to different filterings according to the transcript classes. Additionally, SALA Default was also obtained using the “default” mode. <br>
Following assembly, all the transcript models generated by SALA were categorized relative to the GENCODE v39 reference into three classes: known transcripts (assigned with an ENST transcript ID), novel transcripts derived from GENCODE genes (Novel isoforms attached to an ENSG gene ID), and novel transcripts derived from novel transcriptional units (including intergenic, intronic and antisense). To minimize technical artifacts, we first excluded models associated with putative internal priming and removed GENCODE reference transcripts that showed no evidence in our dataset to generate SALA Raw dataset. For SALA Final, only the novel transcript models with their starting sites linked to the confident ex5_cluster supported by SCAFE were kept. Specifically for the novel transcripts derived from GENCODE genes, 5 supporting reads from each of the replicates and 10% transcript ratio according to complete read count from one of the cell types were required. A confident ex3_cluster supported by clustering or GENCODE was also required for this transcript group. For the transcripts from novel transcriptional units (majority are novel lncRNAs), no further filter was applied. Finally, an additional filter was applied to exclude reference transcript models without SCAFE TSS supported for consistent downstream analyses, while this filter is not part of the SALA sensitive mode. <br>
In SALA default, an additional filter was applied to the novel transcripts derived from novel transcriptional units, where at least one complete read from each replicate was required. The requirement of SCAFE TSS support from reference transcript models was voided. For each transcriptome set, transcript models were subjected to the gene annotator with the parameter “--disable_ref_chain_bound_gene_anno=no”. <br>
```
[1] perl assemble_gene_annotator_v0.1.pl \
--chrom_size_path=[hg38.gencode_v39_chrom.sizes.tsv] \
--model_bed_bgz=[table4wENST.bed12.bed.bgz] \
--model_info_gz=[table4wENST.info.gz] \
--revert_ref_model_bed_bgz=[GENCODEv39.transcript.bed.bgz] \
--ref_model_gene_link=[GENCODEv39.transcript_to_gene.tsv] \
--novel_gene_prefix=ONTG \
--disable_ref_chain_bound_gene_anno=no \
--min_ref_exon_overlap_pct=10 \
--exon_overlap_dist=-1 \
--locus_merge_dist=100000 \
--exclude_t_type=retained_intron \
--out_prefix=Neuron_THP1_T4_10percent \
--out_dir=[output] \
--max_thread=1
```

## <a name="CPAT"></a>Coding potential of transcript and gene models

The coding potential of all assembled transcripts was evaluated using the Coding Potential Assessment Tool (CPAT, v3.0.4) with default parameters. For novel transcript models, CPAT score < 0.364 were defined as non-coding RNA (ncRNA). Transcripts meeting this non-coding criterion with a length of 200 bp were further classified as long non-coding RNAs (lncRNAs). For the GENCODE annotated transcript models, we adopted the lncRNAs annotation as in GENCODE. To assign functional categories at the gene level, we utilized a majority-rule logic based on transcript abundance. A novel transcriptional unit was classified as non-coding if > 50% of its total transcriptional output (determined by complete read counts) was derived from ncRNA transcripts. Within this non-coding pool, genes containing at least one ncRNA transcript ≥ 200 bp were designated as lncRNA transcriptional units, while those exclusively producing shorter transcripts were classified as short_ncRNAs. <br>
For downstream analyses, we applied strict filtering to ensure the purity of the non-coding dataset. For novel loci, only ncRNA transcripts associated with ncRNA genes were retained. For GENCODE loci annotated as lncRNA, we included both known and novel lncRNA isoforms. Thus, novel lncRNAs discovered from GENCODE protein-coding genes were not included for downstream analyses. <br>

## <a name="GTF"></a>Construction of different GTF files

The transcriptome datasets (SALA Raw, SALA Final and SALA Default) were used to construct corresponding GTF files. These GTF files were restricted to gene, transcript and exon features. The attribute field for these features was standardized to include gene_id, transcript_id, gene_type, gene_name, transcript_type, transcript_name, gene_novelty, transcript_novelty and exon_number. Reference gene and transcript annotations were inherited from GENCODE v39. When alternative TSSs and TESs were detected from the reference genes or transcripts, their coordinate boundaries were updated in the GTF file, and the gene_novelty and transcript_novelty fields were flagged as “GENCODE_updated”. While the GTF file of SALA Raw contains all novel transcripts (without internal priming affected ones) and all reference transcripts, SALA Final and Default contains only full-length detected reference transcript models. <br>
To facilitate annotation and quantification, reference models were integrated back to the detectable transcriptomes. For instance, SALA Final with partially detectable reference transcript models were used to generate GTF files for Bambu quantification and reference of short-read RNA-seq quantification by kallisto. SALA Final with all reference transcript models were used as a reference for iPSC chromatin-bound transcriptome construction in SALA. GTF files derived from other assemblers are also provided. <br>

## <a name="chromatin"></a>SALA annotation on Chromatin-bound iPSC dataset

From the CFC-seq of iPSC chromatin-bound RNA, sequencing reads were pre-processed as described above and subjected to SALA. This run was applied to SALA sensitive mode using SALA Final and all GENCODE (v39) transcript models as reference. The parameters used and filtering conditions were the same as SALA Final.

# <a name="subclass"></a>Sub-classification of lncRNA transcripts and genes

From the CFC-seq with 5 cell-types, “SALA Final” detected 6,559 known ncRNA transcripts and identified 115,188 novel ncRNA transcripts. To characterize their genomic context relative to transcript models annotated as protein-coding and pseudogenes (GENCODE v39), we implemented a hierarchical classification system. Transcripts were assigned to one of five mutually exclusive categories using an ascending order of priority: divergent ncRNAs, sense overlap ncRNAs, sense intronic ncRNAs, other sense ncRNAs, antisense ncRNAs, anisense intronic ncRNAs, other antisense ncRNAs and intergenic ncRNAs. <br>
Divergent ncRNAs: transcripts with their TSSs located within 2 kb of a protein-coding or pseudogene transcript TSS on the opposite strand. Sense intronic ncRNAs: transcripts that intersect with protein coding or pseudogenes, with their TSS and the whole transcript lengths lying inside the intron. For those transcripts with <= 10% and > 10% exon sequence overlap with exons are classed as other sense ncRNAs and sense overlap ncRNAs respectively. Antisense ncRNAs: transcripts with >10% of their exon region overlapping with the protein-coding or pseudogenes exon on the opposite strand. Anisense intronic ncRNAs with their entire transcript region located inside the protein-coding or pseudogenes intron on the opposite strand while the remaining antisense ncRNAs are defined as other antisense ncRNAs. Intergenic ncRNAs: All remaining ncRNA transcripts that did not meet the criteria for the previous categories, residing in regions devoid of annotated gene features. <br>
For ncRNA loci containing multiple transcript isoforms, the same hierarchical priority was applied at the loci level. For example, a gene containing both sense-intronic and antisense isoforms was classified as sense-intronic due to the higher hierarchical priority of that category. This systematic approach ensured that each novel transcriptional unit was uniquely categorized based on its most structurally distinct regulatory relationship with the annotated genome. <br>

