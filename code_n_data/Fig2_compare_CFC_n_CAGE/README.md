# Comparing the TSS location, tCRE identification between CFC-seq and CAGE

* compare_CFC_CAGE.R -> main code describing all the analyses
* ssCAGE.STAR.map.sh -> associated code for ssCAGE mapping
* ont.*.sh -> associated code for extracting MAPQ from long-read bam
* SCAFE run independently from ssCAGE and CFC-seq and joint SCAFE are stored in ./SCAFE
* Results are presented in Fig.2 and ext_Fig.3
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Comparison of CFC-seq and CAGE in identification of TSSs and tCREs](TSS)
* [Mappability of repetitive elements using CFC-seq and CAGE](#repeat)



# <a name="TSS"></a>Comparison of CFC-seq and CAGE in identification of TSSs and tCREs

Sequenced reads from ssCAGE were de-multiplexed and aligned to the hg38 human genome assembly using the STAR aligner. This yielded a median depth of 60 million mapped paired-end reads across the generated libraries. Additionally, single-end alignment was also performed using the Read1 alone. Along the analysis with different MAPQ thresholds, alignments were filtered by MAPQ to obtain CTSS in the SCAFE pipeline. To assess the content of these tCREs de novo identified from CFC-seq, ssCAGE from the corresponding samples were analysed. The reads from ssCAGE were independently subjected to the SCAFE pipeline or aggregated with the CFC-seq CTSS by SCAFE (Ext_Fig. 2a). The expression level shown in Ext_Fig. 3j was derived from the aggregated tCREs derived from both CFC-seq and ssCAGE signals. Read counts per tCRE were included only if the CTSSs were positioned inside genuine TSS clusters. The counts were RLE normalized using the edgeR package.

# <a name="repeat"></a>Mappability of repetitive elements using CFC-seq and CAGE 

To determine whether the extended read length from CFC-seq improves the resolution of TSSs within or near repetitive elements, we compared the MAPQ across technologies. We utilized a set of 202,862 TSS clusters identified via joint SCAFE analysis that were detected across all three platforms: CFC-seq, single-end CAGE, and paired-end CAGE. Summits of these TSS clusters were extended 150 nt downstream and intersected with the genomic coordination of RepeatMasker obtained from UCSC. Each TSS cluster was assigned the maximum MAPQ score of its constituent reads per platform. As the MAPQ is not directly comparable between long-read and and short-read aligners, the MAPQ score from short-read datasets were revised for visualization: 0 to 0; 2 to 5; 3 to 10; and 225 to 20. <br>
To compare the evolutionary age of the repetitive elements associated to CFC-seq specific TSS clusters, we summarized the milliDiv (divergence from the consensus sequence) as provided by RepeatMasker. From the 32,365 TSS clusters overlapped with repetitive elements, we extracted the repetitive elements successfully resolved by CFC-seq but not CAGE at the TSS clusters (MAPQ ≥ 20 from CFC-seq and MAPQ < 225 from single-end CAGE). This allowed us to determine if CFC-seq's increased read length is particularly critical for resolving transcription from younger, more sequence-homologous transposable elements. <br>
Genomic coordination of human short tandem repeats (STR) were obtained from the HipSTR reference dataset (hg38.hipstr_reference.bed.gz). These STR intervals were intersected with the summits of major-strand CRE. STR groups fewer than 30 independent intersections across our tCREs dataset were excluded from downstream analysis. To identify motif overrepresentation, two tailed Fisher’s exact tests were performed for each distinct promoter class (promoter-like, enhancer-like, CTCF-alone and unclassified tCREs), using all remaining clsses as the background control set. <br>


