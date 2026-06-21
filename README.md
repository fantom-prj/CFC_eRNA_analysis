# Methodology for "Genomic codes governing enhancer RNA fate"


Long-read sequencing has transformed transcriptome profiling, yet capturing full-length, non-polyadenylated transcripts like enhancer RNAs (eRNAs) remains challenging. Here, we introduce CFC-seq, combining cap-trapping and in vitro poly(A)-tailing to sequence poly(A) and non-poly(A) RNAs with precise transcription start site. Paired with our assembler, SALA, we identified 39,425 novel transcriptional units, including ~24,000 eRNAs. Our data reveal a distinct genomic code governing eRNA fate dictated by core promoter architecture. CpG-island enhancers show high chromatin connectivity but yield short, exosome-sensitive RNAs. Conversely, TATA-box enhancers systematically co-opt LTR retrotransposons to inherit structural motifs that produce long, stable, and spliced RNAs. Mechanistically, the pioneer factor NF-Y activates these viral elements to license transcription, balanced by TEAD4 activity across a dual-gear regulatory axis. Finally, non-poly(A) eRNAs terminate via exosome-associated processing at structural-depleted cleavage zones. This comprehensive annotation links enhancer sequence architecture to RNA fate, providing a new transformative framework for decoding the functional human genome.

## Method to code links

| Method                                                                                                              | Code                 |
| ------------------------------------------------------------------------------------------------------------------- | -------------------- |
| Pre-processing of the long-read ONT data                                                                            | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/preprocessing) |
| Mapping to genome and TranscriptClean                                                                               | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/preprocessing) |
| Comparative analysis of long-read 5’ end precision                                                                  | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig1_read_analyses) |
| Assessment of 3’ ends capture and filtration of potential internal priming                                          | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig1_read_analyses) |
| Computational identification and hierarchical scoring of PAS motifs                                                 | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig1_read_analyses) |
| Capture of non-poly(A) ncRNA and incomplete transcripts derived from PAT protocol                                   | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig1_read_analyses) |
| Poly(A)-tail prediction                                                                                             | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig1_read_analyses) |
| Confident splice junctions from short-read RNA-seq                                                                  | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/preprocessing) |
| TSS clusters and tCREs identification by SCAFE                                                                      | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SCAFE) |
| Promoter typing                                                                                                     | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_CRE_analysis) |
| Comparison of CFC-seq and CAGE in identification of TSSs and tCREs                                                  | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_compare_CFC_n_CAGE) |
| Mappability of repetitive elements using CFC-seq and CAGE                                                           | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_compare_CFC_n_CAGE) |
| Transcript model construction using Transcript Start-site Aware Long-read Assembler (SALA)                          | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SALAe) |
| Coding potential of transcript and gene models                                                                      | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SALA) |
| Sub-classification of lncRNA transcripts and genes                                                                  | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SALA) |
| Transcript models construction by TALON                                                                             | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/other_assemblers) |
| Transcript models construction by Isoquant                                                                          | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/other_assemblers) |
| Construction of different GTF files                                                                                 | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SALA) |
| SALA annotation on Chromatin-bound iPSC dataset                                                                     | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/SALA) |
| Analyses on short-read RNA-seq                                                                                      | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig3_transcript_model_analyses) |
| Quantification of CFC-seq samples                                                                                   | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig3_transcript_model_analyses) |
| Comparing final gene and transcript models with four reference repositories                                         | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig3_transcript_model_analyses) |
| Comparison across assemblers                                                                                        | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig3_transcript_model_analyses) |
| Isoform switching analyses                                                                                          | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig3_transcript_model_analyses) |
| Genomic properties incorporation to tCRE and ex5_cluster                                                            | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_CRE_analysis) |
| Bidirectionality of tCREs                                                                                           | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_CRE_analysis) |
| Super enhancer identification                                                                                       | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_CRE_analysis) |
| Integration of histone modification                                                                                 | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig2_CRE_analysis) |
| Grouping of ex5_clusters, transcript models and gene models into mRNA, p_ncRNA, e_ncRNA, CTCF_ncRNA and other_ncRNA | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig4_transcription_features) |
| Transcription properties analyses                                                                                   | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig4_transcription_features) |
| Exosome sensitivity                                                                                                 | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig4_transcription_features) |
| RNA binding protein interaction analysis                                                                            | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig4_transcription_features) |
| Splicing efficiency and spliceAI                                                                                    | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig5_splicing_efficiency) |
| RNA structural prediction proximal to transcript end site                                                           | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig5_TES_analyses) |
| Genomic and Regulatory Characterization of TES                                                                      | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig5_TES_analyses) |
| eQTL and GWAS SNP locations in expanded transcriptome and their effect on splicing                                  | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_GWAS_eQTL_SNP) |
| Enrichment of GWAS and eQTL SNPs in enhancers                                                                       | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_GWAS_eQTL_SNP) |
| Hi-C experimental and analytical procedures                                                                         | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_enhancer_connectivity) |
| Chromatin connectivity                                                                                              | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_enhancer_connectivity) |
| Chromatin accessibility and GC content for Hi-C bin                                                                 | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_enhancer_connectivity) |
| Prediction of enhancer-gene interactions using ABC model                                                            | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_enhancer_connectivity) |
| Repeat elements in ex5_clusters and transcript models                                                               | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_repetitive_element) |
| Transcription factor motif enrichment                                                                               | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_repetitive_element) |
| Integration of NFYB and TEAD4 ChIP-seq datasets                                                                     | [Code](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/Fig6_repetitive_element) |



| Method | Code |
|---------|------|
| Pre-processing of ONT Long-Read Data | [Link](https://github.com/fantom-prj/CFC_eRNA_analysis/tree/main/code_n_data/preprocessing) |
| Genome Mapping and TranscriptClean Processing | [Link](path/to/code) |
| Comparative Analysis of Long-Read 5′ End Precision | [Link](path/to/code) |
| 3′ End Assessment and Internal Priming Filtering | [Link](path/to/code) |
| PAS Motif Identification and Hierarchical Scoring | [Link](path/to/code) |
| Detection of Non-poly(A) ncRNAs and Incomplete PAT Transcripts | [Link](path/to/code) |
| Poly(A) Tail Prediction | [Link](path/to/code) |
| High-Confidence Splice Junction Identification | [Link](path/to/code) |
| TSS Cluster and tCRE Identification Using SCAFE | [Link](path/to/code) |
| Promoter Classification | [Link](path/to/code) |
| Comparison of CFC-seq and CAGE for TSS and tCRE Detection | [Link](path/to/code) |
| Mappability Analysis of Repetitive Elements | [Link](path/to/code) |
| Transcript Model Construction Using SALA | [Link](path/to/code) |
| Coding Potential Prediction | [Link](path/to/code) |
| lncRNA Transcript and Gene Classification | [Link](path/to/code) |
| Transcript Model Construction Using TALON | [Link](path/to/code) |
| Transcript Model Construction Using IsoQuant | [Link](path/to/code) |
| Construction of Reference and Custom GTF Files | [Link](path/to/code) |
| SALA Annotation of Chromatin-Bound iPSC Data | [Link](path/to/code) |
| Short-Read RNA-seq Analysis | [Link](path/to/code) |
| Quantification of CFC-seq Samples | [Link](path/to/code) |
| Comparison with Reference Transcriptome Repositories | [Link](path/to/code) |
| Comparison Across Transcriptome Assemblers | [Link](path/to/code) |
| Isoform Switching Analysis | [Link](path/to/code) |
| Integration of Genomic Features into tCREs and ex5 Clusters | [Link](path/to/code) |
| Analysis of tCRE Bidirectionality | [Link](path/to/code) |
| Super Enhancer Identification | [Link](path/to/code) |
| Histone Modification Integration | [Link](path/to/code) |
| Classification of ex5 Clusters, Transcripts, and Genes | [Link](path/to/code) |
| Transcriptional Property Analysis | [Link](path/to/code) |
| Exosome Sensitivity Analysis | [Link](path/to/code) |
| RNA-Binding Protein Interaction Analysis | [Link](path/to/code) |
| Splicing Efficiency and SpliceAI Analysis | [Link](path/to/code) |
| RNA Structure Prediction Near Transcript Ends | [Link](path/to/code) |
| Genomic and Regulatory Characterization of TESs | [Link](path/to/code) |
| eQTL and GWAS SNP Analysis in the Expanded Transcriptome | [Link](path/to/code) |
| Enrichment Analysis of GWAS and eQTL SNPs in Enhancers | [Link](path/to/code) |
| Hi-C Experimental and Computational Analysis | [Link](path/to/code) |
| Chromatin Connectivity Analysis | [Link](path/to/code) |
| Chromatin Accessibility and GC Content Analysis | [Link](path/to/code) |
| Enhancer–Gene Interaction Prediction Using ABC | [Link](path/to/code) |
| Repeat Element Analysis in ex5 Clusters and Transcript Models | [Link](path/to/code) |
| Transcription Factor Motif Enrichment Analysis | [Link](path/to/code) |
| Integration of NFYB and TEAD4 ChIP-seq Data | [Link](path/to/code) |


## Software and Tools

| Tool | Source |
|------|--------|
| SALA (v1.0) | https://github.com/fantom-prj/SALA |
| Tail trimmer (v1.4) | https://github.com/fantom-prj/SALA |
| Dorado (v0.2.4) | https://github.com/nanoporetech/dorado |
| Minimap2 (v2.17-r974-dirty) | https://github.com/lh3/minimap2 |
| TranscriptClean (v2.0.3) | https://github.com/mortazavilab/TranscriptClean |
| TALON (v5) | https://github.com/mortazavilab/TALON |
| SCAFE (v1.01) | https://github.com/chung-lab/SCAFE |
| STAR (v2.7.11b) | https://github.com/alexdobin/STAR |
| HOMER (v4.11) | https://github.com/javrodriguez/HOMER |
| MEME Suite (v5.5.7) | https://meme-suite.org/meme/meme-software |
| CPAT (v3.0.4) | https://github.com/liguowang/cpat |
| SQANTI3 (v5.3.0) | https://github.com/ConesaLab/SQANTI3 |
| IsoQuant (v3.4.1) | https://github.com/ablab/IsoQuant |
| bambu (v3.2.4) | https://github.com/GoekeLab/bambu |
| Bowtie2 (v2.2.6) | https://github.com/BenLangmead/bowtie2 |
| MACS2 (v2.1.0) | https://github.com/jdavisturak/MACS2-2.1.1.20160309 |
| ROSE (v1.3.1) | https://github.com/stjude/ROSE |
| RNAfold (v2.6.4) | https://github.com/ViennaRNA/ViennaRNA |
| ChromHMM (v1.24) | https://ernstlab.github.io/ChromHMM/ |
| IsoformSwitchAnalyzeR (v2.6.1) | https://github.com/kvittingseerup/IsoformSwitchAnalyzeR |
| ABC model (v1.1.2) | https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction |
| Primer-chop | https://gitlab.com/mcfrith/primer-chop |
| Paraclu | https://github.com/davetang/paraclu_prep |
| SAMtools (v1.11) | https://www.htslib.org/ |
| BEDTools (v2.30.0) | https://bedtools.readthedocs.io/en/latest/ |
| Tabix / bgzip (v1.15.1) | https://www.htslib.org/ |
| bedGraphToBigWig (v2.8) | https://github.com/ENCODE-DCC/kentUtils/tree/master/src/utils |
| bedparse (v0.2.3) | https://github.com/tleonardi/bedparse |
| R (v4.3.1) | https://cran.r-project.org/ |
| Perl (v5.26.2) | https://www.perl.org/get.html |


## Data Availability

### Core Data

| Dataset | Accession Number / Source |
|----------|----------|
| Bulk RNA-seq (iPSC, NSC, Neuron) | DRA019571 (DRR614936–DRR614941) |
| CFC-seq (iPSC, NSC, Neuron) | DRA019506 (DRR613094–DRR613110) |
| CFC-seq (THP-1, dTHP-1) | DRA019521 (DRR613160–DRR613175) |
| Chromatin-bound CFC-seq (iPSC) | DRA019506 (DRR613111, DRR613112) |

### Supportive Data

| Dataset | Accession Number / Source |
|----------|----------|
| ssCAGE (iPSC, NSC, Neuron) | DRA019567 (DRR614867–DRR614872) |
| Single-cell ATAC (iPSC, NSC, Neuron) | DRA019608 (DRR618513–DRR618515) |
| CUT&Tag (iPSC, NSC, Neuron; H3K27ac, H3K27me3, H3K4me1, H3K4me3, CTCF & control) | DRA019568 (DRR614873–DRR614905) |
| Hi-C (iPSC, NSC, Neuron) | DRA019572 (DRR614942–DRR614953) |
| iPSC TEAD4 ChIP-seq | DRA026820 (DRR959946–DRR959949) |
| iPSC ATAC-seq with TEAD4 knockdown | DRA026820 (DRR959956–DRR959961) |

### External Data

| Dataset | Accession Number / Source |
|----------|----------|
| CAP-trap PacBio long-read (WTC11 iPSC) | ENCODE: ENCSR309IKK |
| Uncapped-depleted PacBio long-read (WTC11 iPSC) | ENCODE: ENCSR507JOF |
| R2C2 ONT long-read (WTC11 iPSC) | ENCODE: ENCSR925UQZ |
| TSO ONT long-read (WTC11 iPSC) | ENCODE: ENCSR539ZXJ |
| Direct RNA ONT long-read (WTC11 iPSC) | ENCODE: ENCSR392BGY |
| MYC ChIP-seq (human ES cells) | GSM1505809 |
| NFYB ChIP-seq (WTC11 iPSC) | ENCODE: ENCSR146UIC |
| Exosome knockdown RNA-seq (iPSC) | DRA013847 |

### Public Resources

| Resource | URL |
|----------|----------|
| eQTL Catalogue | https://www.ebi.ac.uk/eqtl/ |
| CAUSALdb GWAS | http://www.mulinlab.org/causaldb/index.html |
| SCREEN cCRE database | https://screen.encodeproject.org/ |
| Intropolis splice junction resource | https://github.com/nellore/intropolis |
| Open Targets disease gene database | https://platform.opentargets.org/ |
| POSTAR3 RNA–protein interactions | http://111.198.139.65/ |
| HipSTR short tandem repeat annotations | https://github.com/HipSTR-Tool/HipSTR |
| EVLncRNAs 3.0 | https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/ |



