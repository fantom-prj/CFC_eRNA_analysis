# Enhancer chromatin connectivity 

* HiC_connect.R -> main code describing all the analyses
* ABC_input_and_output.R -> associated code running ABC model
* Results are presented in Fig.6 and ext_Fig.8
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Hi-C experimental and analytical procedures](#HiC)
* [Chromatin connectivity](#connectivity)
* [Chromatin accessibility and GC content for Hi-C bin](#bin)
* [Prediction of enhancer-gene interactions using ABC model](#ABC)


# <a name="HiC"></a>Hi-C experimental and analytical procedures

The cells (iPSC, NSC and Neuron) were cultured and harvested as described above. Two replicates of each cell type were fixed by 1% formaldehyde and Hi-C library construction was performed with the Arima-HiC+ kit and the Arima Library Prep Module (Arima Genomics) according to the manufacturer’s protocols. The libraries were then subjected to Illumina NovaSeq6000 for sequencing with 150 bp paired-end mode. In order to yield a higher proportion of aligned reads, we first trim the 6 bp from the left end of each read as suggested by Arima Hi-C data analysis and then trim the 70 bp from the right end of each read using the command line “seqtk trimfq -b 6 -e 70”. The remaining 75 bp reads were used for the nf-core hic pipeline v 2.1.0 (https://nf-co.re/hic/2.1.0/). Since the two replicates of the Hi-C libraries yield consistent results in a pre Hi-C processing based on the subsampled samples. The two replicates from the same cell type were pooled together for the final nf-core hic pipeline with using the following options:
```
nextflow run nf-core/hic
--genome 'hg38'
--restriction_site '[^GATC,G^ANTC]'
--ligation_site '[GATCGATC,GANTGATC,GANTANTC,GATCANTC]'
--digestion 'arima'
--res_compartments '500000,250000,100000'
--tads_caller 'insulation,hicexplorer’
--bin_size '1000000,500000,100000,50000,25000,10000,5000,2500,1000'
```
The processed Hi-C data was transformed into valid genomic interaction pairs at 5kb, 10kb and 25kb resolutions. The statistical significance of these interactions was calculated using the Bioconductor package GOTHiC (version 1.40.0). Intra-chromosomal (cis) interactions supported by at least five read pairs and a q-value (FDR) ≤ 0.01 were deemed significant for downstream analysis. The Hi-C interaction at 5-kb resolution was used in this study, where about 10% of the valid bin pairs were found significant. 


# <a name="connectivity"></a>Chromatin connectivity

To link the Hi-C data with ncRNA ex5_clusters, 5 kb bins overlapping with their summits were extracted and grouped based on promoter types (p_ncRNA, e_ncRNA, and other_ncRNA) and regulatory elements (dCGI, uCGI, Null, and TATA box), while bins containing more than one promoter type or regulatory element were excluded from the analysis. To link mRNA to Hi-C bins, only GENCODE (v39) protein-coding transcript models detectable in the final transcriptome were included. A "target" bin was considered if it overlapped with the ex5_cluster summit of these protein-coding transcripts. The number of bin pairs grouped by the above criteria was visualized in Figure 6d. To minimize the possibility of ncRNA ex5_clusters being missed by Hi-C due to technical differences between Hi-C and CFC-seq, only ex5_clusters showing at least one significant connection to an mRNA ex5_cluster were included. However, including ncRNA ex5_clusters with zero contacts in the analysis resulted in the same overall observation.

# <a name="bin"></a>Chromatin accessibility and GC content for Hi-C bin
The chromatin accessibility of the source Hi-C bins defined in the previous section was assessed using scATAC-seq data. The scATAC-seq data were recounted into the Hi-C bins using the “CreateFragmentObject” function from the R package Signac. The resulting quantification was aggregated based on cell type. The GC content of the Hi-C bins was determined by calculating the ratio of G and C residues.

# <a name="ABC"></a>Prediction of enhancer-gene interactions using ABC model
Enhancer-gene interactions were predicted using the Activity-by-Contact (ABC) model.33 Analyses were performed independently for iPSC, NSC and Neuron using single-cell ATAC-seq (grouped by cell types), bulk H3K27ac CUT&Tag and Hi-C data. The ABC pipeline was executed using the Snakemake workflow under the hg38 genome assembly. Signals were normalized using the reference quantile-normalization distribution supplied with the ABC workflow (EnhancersQNormRef.K562.txt). Gene annotations were based on the hg38 CollapsedGeneBounds reference provided by the ABC package. <br>
Instead of generating candidate regulatory elements from MACS2 peak calling, we supplied a custom catalog of tCREs identified from CFC-seq. These tCREs were merged for overlap from opposite strand. Candidate regions include both distal enhancer-like tCREs and promoter-associated tCREs while the median merged tCRE length was 501 bp. These regions were formatted as BED intervals and substituted directly into the ABC workflow at the candidate-region generation stage. <br>
Predicted enhancer–gene interactions were extracted from the thresholded ABC output using the default score threshold of 0.02. Only interactions exceeding this threshold were considered high-confidence enhancer–gene links for downstream analyses. For comparative analyses across differentiation, ABC scores were calculated independently in iPSC, NSC and neuron samples using the same candidate regulatory element catalog. The final output was shown in Supplementary Table S16. <br>
```
snakemake -c16 --rerun-incomplete --keep-going \
  ./results/iPSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  ./results/NSC/Predictions/EnhancerPredictionsAllPutative.tsv.gz \
  ./results/Neuron/Predictions/EnhancerPredictionsAllPutative.tsv.gz
```








