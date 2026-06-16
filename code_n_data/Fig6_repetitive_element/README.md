# TATA box enhancer/eRNA enrichment with LTR, NFY/TEAD activating their transcription 

* repetitive.R -> main code describing all the analyses
* LTR_ex5.R -> associated code for LTR distribution across ex5_cluster
* NFY_ex5.R -> associated code for NFY binding distribution across ex5_cluster
* Results are presented in Fig.6 and ext_Fig9
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Repeat elements in ex5_clusters and transcript models](#repeat)
* [Transcription factor motif enrichment](#TF)
* [Integration of NFYB and TEAD4 ChIP-seq datasets](#NFY)


# <a name="repeat"></a>Repeat elements in ex5_clusters and transcript models

The coordinates of repeat elements in hg38 were obtained from UCSC RepeatMasker file. These repeat elements were intersected with the regions of ex5_clusters and the exons of transcript models. For ex5_clusters, a minimum overlap of 6 nt was required, and if multiple repeat elements were identified, the one with the greatest overlapping length was selected as representative. For exons, the overlap of a single feature across multiple exons of a transcript model was summed, with a minimum overlap of 200 nt required. The repeat elements intersecting at the transcript level were collapsed into the ex5_cluster level by selecting the one with the greatest overlap length. All the overlapped repeat elements that met the above criteria were retained in Supplementary Table S17.


# <a name="TF"></a>Transcription factor motif enrichment

Summits from ex5_clusters were extended 400 nt upstream and 100 nt downstream. These regions were excluded if they overlapped and involved more than one transcript class (mRNA, p_ncRNA, e_ncRNA, other_ncRNA and CTCF_ncRNA) or more than one promoter structure (CGI, and TATA box). The ex5_cluster with stronger signal was kept if their extended regions overlap. Finaly each 501 nt region contain only one ex5_cluster. These regions (n = 39,794) were grouped by transcript class and promoter structure and subjected to SEA (Simple Enrichment Analysis) from the MEME-suite (v5.5.7) to calculate the relative motif enrichment. We obtained the sequences from each region as the input sequence in fasta format and the shuffled sequences as the control sequences. The core vertebrates-non-redundant motif collection obtained JASPAR2024 was used as the scope of transcription factors. For visualization, only motifs from the TFs with CPM >1 according to CAGE quantification, from at least one of the three cell types were included. We further selected the motifs for q value < 0.01 and enrichment > 12.


# <a name="NFY"></a>Integration of NFYB and TEAD4 ChIP-seq datasets
ChIP-seq data of NFYB derived from iPSC WTC11 was obtained from ENCODE (ENCSR146UIC) where the IDR threshold peaks (n = 9759, median = 390 nt) were used. ChIP-seq of TEAD4 was performed in our iPSC (WTC11) where the merged peaks from the 2 replicates were used (n = 6986, median = 437 nt). These peaks were intersected with the ex5_cluster non-overlapped regions (summit +100 / -400 nt) described from the previous section. The ex5_clusters were grouped according to the H3K27ac/me3 signals of iPSC from the CUT&Tag results, into 4 groups: Co-mark (H3K27ac+ / H3K27me3+), Active (H3K27ac+ / H3K27me3-), Repressed (H3K27ac- / H3K27me3+) and Un-marked (H3K27ac- / H3K27me3-). Enrichment tests were performed by Fisher’s exact, using the other groups as background.<br>
ChIP-seq data of NFYB was also used to overlay with the ex5_clusters, which were restricted to the list of non-overlapped regions and active in iPSC (TSS count > 1). The cumulative coverage was determined from distinct transcript classes and promoter structures, as in Extended Data Figure 9g.










