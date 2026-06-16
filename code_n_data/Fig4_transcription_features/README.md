# Features of RNA derived from different ex5_cluster & exosome sensitivity

* RNA_features.R -> main code describing all the analyses
* chromatin_bound_n_exosome.R -> associated code gathering information from chromatin_bound iPSC CFC-seq SALA and exosome KD
* ./exosome_sensitivity/kallisto_quantify_exosome.pl -> associated code for running kallisto
* Results are presented in Fig.4 and ext_Fig.6
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Grouping of ex5_clusters, transcript models and gene models into mRNA, p_ncRNA, e_ncRNA, CTCF_ncRNA and other_ncRNA](grouping)
* [Transcription properties analyses](#properties)
* [Exosome sensitivity](#exosome)
* [RNA binding protein interaction analysis](#RBP)


# <a name="grouping"></a>Grouping of ex5_clusters, transcript models and gene models into mRNA, p_ncRNA, e_ncRNA, CTCF_ncRNA and other_ncRNA

The ex5_clusters supported by SCAFE were divided into 4 groups according to the promoter types, which were inherited from the tCRE as promoter-like, enhancer-like, CTCF-alone with ATAC support and unclassed with ATAC support. Only the ex5_clusters producing ncRNAs (> 50% complete read count) were included while the definition of ncRNAs were described according to the CPAT results or GENCODE annotation. Additionally, ex5_clusters producing GENCODE protein-coding genes alone were grouped as a control. The ex5_clusters linked to both coding (defined by GENCODE) and non-coding (defined by GENCODE or CPAT) transcript models were excluded from the analyses. Finally, ex5_clusters linked to ncRNA transcript models belong to protein-coding genes were also excluded. This grouping was used in most of the analyses in this study. These resulted in 8,759 ex5_clusters with mRNA outputs, 16,200 promoter-like ex5_clusters with ncRNA outputs (p_ncRNA), 28,354 enhancer-like ex5_clusters with ncRNA (e_ncRNA), 190 ATAC-supported CTCF-alone ex5_clusters with ncRNA (CTCF_ncRNA) and 4,606 ATAC-supported unclassed ex5_clusters with ncRNA (other_ncRNA). For downstream analyses where other features were involved, some extra ex5_clusters were excluded. These were described in the specific analyses. Once the ex5_clusters were included, all the read derived from the clusters were included, independent of coding potential. These reads are considered as the RNA output of the ex5_clusters.
For grouping of transcript models, the identities of coding and non-coding and promoter types (promoter-like, enhancer-like, CTCF-alone with ATAC support and unclassed with ATAC support) were directly transferred from the finalized transcript information table. As described previously, only non-coding transcripts linked to non-coding genes were included into downstream analyses (Ext_Fig. 4j). From these criteria, 28,230 mRNAs, 36,723 p_ncRNAs, 64,010 e_ncRNAs, 437 CTCF_ncRNAs and 14,847 other_ncRNAs were obtained, while 49,688 transcript models were excluded. This grouping was used in summarizing if the presence of splicing and PAS correlate with transcript length, and the qualifying the libraries (Fig. 1e)
The grouping of gene models was used in comparing splicing efficiency and exosome sensitivity. The identities of coding and non-coding were derived from coding classes at gene level described in previous part. The identities of promoter types were derived from the promoter type at gene level described in the previous part. If the genes were annotated as protein-coding by GENCODE, the genes were assigned as “mRNA”. From these criteria, 14,274 mRNAs, 12,696 p_ncRNAs, 22,403 e_ncRNAs, 157 CTCF_ncRNAs and 3,787 other_ncRNAs were obtained. Gene level classification was used in splicing efficiency analyses.


# <a name="properties"></a>Transcription properties analyses

All the transcription properties including transcript length, genomic range, number of exon, splice junction and TES, are defined at the read level and grouped into ex5_clusters and transcript models. The transcript length and exon number of transcript models and those derived from ex5_clusters were determined from the median of all the reads linked to them. This information was summarized in Supplementary tables S10 and S11. The ex5_clusters were grouped into mRNAs, p_ncRNAs, e_ncRNAs, CTCF_ncRNAs and other_ncRNAs as described earlier. Ex5_clusters were also grouped for directionality (1D or 2D) and being within super enhancer cluster (SE) or not (typical enhancer, TE). This information was inherited from their linked tCREs. The e_ncRNA ex5_clusters were also grouped mutally exclusive as CGIap, CGInap, Null and TATA while p_ncRNA ex5_clusters were grouped in CGI, Null and TATA, according to as earlier described.
# <a name="exosome"></a>Exosome sensitivity
RNA-seq data of another clone of iPS cell treated with knockdown of EXOSC3 and control was annotated by the transcriptome generated in this study using pseudo alignment of Kallisto using the chromatin-bound transcriptome containing all the completely and partially detected GENCODE models, excluding ribosomal RNAs, as reference. Exosome sensitivity score was calculated at transcript and gene levels as previously described: (TPMExosome-suppressed - TPMControl) / TPMExosome-suppressed. Transcripts that were not detectable from these short-read libraries (sum of TPM from 4 libraries > 0.1) or not detectable from the both the iPSC and chromatin-bound iPSC CFC-seq libraries were excluded from the analysis. Exosome sensitivity smaller than 0 was considered as 0. To collapse exosome sensitivity at transcript and gene levels into ex5_clusters, exosome sensitivity was weighted by read counts when grouping. We performed the downstream analyses using the transcript level exosome sensitivity that one transcript only has one ex5_cluster. For single end RNA-seq with exosome-KD:
```
kallisto quant -i [kallisto_index] -o [out_dir] --single -l 200 -s 50 --rf-stranded [fq_path]
```

# <a name="RBP"></a>RNA binding protein interaction analysis

To characterize the post-transcriptional regulatory landscape of the identified transcripts, we analyzed the occupancy of 220 RNA-binding proteins (RBPs) using the human CLIPdb from POSTAR3. Binding regions were intersected with the exonic sequences of finalized transcript models (SALA Final). While 83% of exons (402544/483024) derived from promoter-like transcript models are bound by RBPs, only 26% of enhancer-derived exons (45692/162523) exhibited RBP binding. To identify specific RBP binding on eRNA exon, enrichment test was peformed for each RBP to compare proportion of peaks of a specific RBP in enhancer/promoter relative to the total RBP binding events as background. Significant enrichment of specific RBPs in enhancer exons was determined using a two-sided Fisher’s exact test, with p-values adjusted for multiple testing.








