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

## References
1.	Yip, C. W. et al. Antisense-oligonucleotide-mediated perturbation of long non-coding RNA reveals functional features in stem cells and across cell types. Cell Rep. 41, 111893 (2022).
2.	Takahashi, H., Lassmann, T., Murata, M. & Carninci, P. 5′ end–centered expression profiling using cap-analysis gene expression and next-generation sequencing. Nat. Protoc. 7, 542–561 (2012).
3.	Takahashi, H., Nishiyori-Sueki, H., Ramilowski, J. A., Itoh, M. & Carninci, P. Low Quantity Single Strand CAGE (LQ-ssCAGE) Maps Regulatory Enhancers and Promoters. Methods Mol. Biol. Clifton NJ 2351, 67–90 (2021).
4.	Li, H. Minimap2: pairwise alignment for nucleotide sequences. Bioinforma. Oxf. Engl. 34, 3094–3100 (2018).
5.	Wyman, D. & Mortazavi, A. TranscriptClean: variant-aware correction of indels, mismatches and splice junctions in long-read transcripts. Bioinforma. Oxf. Engl. 35, 340–342 (2019).
6.	Wyman, D. et al. A Technology-Agnostic Long-Read Analysis Pipeline for Transcriptome Discovery and Quantification. http://biorxiv.org/lookup/doi/10.1101/672931 (2019) doi:10.1101/672931.
7.	Pardo-Palacios, F. J. et al. Systematic assessment of long-read RNA-seq methods for transcript identification and quantification. Nat. Methods https://doi.org/10.1038/s41592-024-02298-3 (2024) doi:10.1038/s41592-024-02298-3.
8.	Heinz, S. et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol. Cell 38, 576–589 (2010).
9.	Legnini, I., Alles, J., Karaiskos, N., Ayoub, S. & Rajewsky, N. FLAM-seq: full-length mRNA sequencing reveals principles of poly(A) tail length control. Nat. Methods 16, 879–886 (2019).
10.	Alfonso-Gonzalez, C. et al. Sites of transcription initiation drive mRNA isoform selection. Cell 186, 2438-2455.e22 (2023).
11.	Moody, J. et al. SCAFE: a software suite for analysis of transcribed cis-regulatory elements in single cells. Bioinforma. Oxf. Engl. 38, 5126–5128 (2022).
12.	Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinforma. Oxf. Engl. 29, 15–21 (2013).
13.	Willems, T. et al. Genome-wide profiling of heritable and de novo STR variations. Nat. Methods 14, 590–592 (2017).
14.	The human RNA-DNA interactome is cell type-specific and dynamic. (2026).
15.	Wang, L. et al. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model. Nucleic Acids Res. 41, e74–e74 (2013).
16.	Prjibelski, A. D. et al. Accurate isoform discovery with IsoQuant using long reads. Nat. Biotechnol. 41, 915–918 (2023).
17.	Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification. Nat. Biotechnol. 34, 525–527 (2016).
18.	Chen, Y. et al. Context-aware transcript quantification from long-read RNA-seq data with Bambu. Nat. Methods 20, 1187–1195 (2023).
19.	Zhao, W. et al. POSTAR3: an updated platform for exploring post-transcriptional regulation coordinated by RNA-binding proteins. Nucleic Acids Res. 50, D287–D294 (2022).
20.	Vitting-Seerup, K. & Sandelin, A. The Landscape of Isoform Switches in Human Cancers. Mol. Cancer Res. MCR 15, 1206–1220 (2017).
21.	Hnisz, D. et al. Super-enhancers in the control of cell identity and disease. Cell 155, 934–947 (2013).
22.	Chapuy, B. et al. Discovery and characterization of super-enhancer-associated dependencies in diffuse large B cell lymphoma. Cancer Cell 24, 777–790 (2013).
23.	Wang, Y. et al. SEdb 2.0: a comprehensive super-enhancer database of human and mouse. Nucleic Acids Res. 51, D280–D290 (2023).
24.	Yip, C. W. et al. Single cell bimodal analyses reveal the mode of activity of transcription factors on enhancers and promoters. Prep. In preparation, (2026).
25.	Lovén, J. et al. Selective inhibition of tumor oncogenes by disruption of super-enhancers. Cell 153, 320–334 (2013).
26.	Whyte, W. A. et al. Master transcription factors and mediator establish super-enhancers at key cell identity genes. Cell 153, 307–319 (2013).
27.	Ernst, J. & Kellis, M. ChromHMM: automating chromatin-state discovery and characterization. Nat. Methods 9, 215–216 (2012).
28.	Jaganathan, K. et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell 176, 535-548.e24 (2019).
29.	Lorenz, R. et al. ViennaRNA Package 2.0. Algorithms Mol. Biol. 6, 26 (2011).
30.	Kerimov, N. et al. eQTL Catalogue: A Compendium of Uniformly Processed Human Gene Expression and Splicing QTLs. http://biorxiv.org/lookup/doi/10.1101/2020.01.29.924266 (2020) doi:10.1101/2020.01.29.924266.
31.	Wang, J. et al. CAUSALdb: a database for disease/trait causal variants identified using summary statistics of genome-wide association studies. Nucleic Acids Res. gkz1026 (2019) doi:10.1093/nar/gkz1026.
32.	Mifsud, B. et al. GOTHiC, a probabilistic model to resolve complex biases and to identify real interactions in Hi-C data. PloS One 12, e0174744 (2017).
33.	Fulco, C. P. et al. Activity-by-contact model of enhancer–promoter regulation from thousands of CRISPR perturbations. Nat. Genet. 51, 1664–1669 (2019).
34.	Bailey, T. L. & Grant, C. E. SEA: Simple Enrichment Analysis of motifs. Preprint at https://doi.org/10.1101/2021.08.23.457422 (2021).
35.	The ENCODE Project Consortium et al. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699–710 (2020).
36.	Nellore, A. et al. Human splicing diversity and the extent of unannotated splice junctions across human RNA-seq samples on the Sequence Read Archive. Genome Biol. 17, 266 (2016).
37.	Koscielny, G. et al. Open Targets: a platform for therapeutic target identification and validation. Nucleic Acids Res. 45, D985–D994 (2017).
38.	Zhou, B. et al. EVLncRNAs 3.0: an updated comprehensive database for manually curated functional long non-coding RNAs validated by low-throughput experiments. Nucleic Acids Res. 52, D98–D106 (2024).
39.	Dobin, A. et al. STAR: ultrafast universal RNA-seq aligner. Bioinforma. Oxf. Engl. 29, 15–21 (2013).
40.	Pardo-Palacios, F. J. et al. SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. Nat. Methods 21, 793–797 (2024).
41.	Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie 2. Nat. Methods 9, 357–359 (2012).
42.	Zhang, Y. et al. Model-based Analysis of ChIP-Seq (MACS). Genome Biol. 9, R137 (2008).
43.	Lorenz, R. et al. ViennaRNA Package 2.0. Algorithms Mol. Biol. 6, 26 (2011).
44.	Frith, M. C. et al. A code for transcription initiation in mammalian genomes. Genome Res. 18, 1–12 (2008).
45.	Danecek, P. et al. Twelve years of SAMtools and BCFtools. GigaScience 10, giab008 (2021).
46.	Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinforma. Oxf. Engl. 26, 841–842 (2010).
47.	Li, H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics 27, 718–719 (2011).
48.	Kent, W. J. et al. The human genome browser at UCSC. Genome Res. 12, 996–1006 (2002).
49.	Leonardi, T. Bedparse: feature extraction from BED files. J. Open Source Softw. 4, 1228 (2019).
50.	R Core Team. R: A language and environment for statistical computing. (2021).



