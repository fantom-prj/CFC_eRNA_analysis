# Genetic variant incorporation to novel transcripts 

* GWAS_eQTL_SNP.R -> main code describing all the analyses
* spliceAI_snp.R -> associated code for running spliceAI
* Results are presented in Fig.6 and ext_Fig.8
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [eQTL and GWAS SNP locations in expanded transcriptome and their effect on splicing](#transcriptome)
* [Enrichment of GWAS and eQTL SNPs in enhancers](#enrichment)


# <a name="transcriptome"></a>eQTL and GWAS SNP locations in expanded transcriptome and their effect on splicing

The eQTL and GWAS data were obtained from eQTL Catalogue and CAUSALdb respectively. For the eQTL data, a cutoff of pip >0.3 was applied to include SNPs that are putatively casual. This dataset contains 367,989 non-redundant SNP locations. For the GWAS data, credible set version 2.0 was downloaded and converted into hg38 by crossing over the rsid with the SNP 151 dataset obtained from UCSC (n = 710,172,836). The GWAS dataset contains 917,107 non-redundant SNP locations after converting to hg38. These SNP locations were intersected with the regions of the ex5_clusters and transcript models (exons). Ex5_clusters linked to at least one GENCODE transcript model were grouped as ENST, otherwise were novel. The transcript models were intersected as they were. The SNP locations mapped to the GENCODE transcript models were removed from the novel group to reveal only the additional SNP coverage from the novel entries. <br>
The donor and acceptor sites of the finalized transcript models and all the GENCODE v39 transcript models were also subjected to the SpliceAI together with the GWAS and eQTL SNPs to identify the SNPs that affect splicing. We updated the transcript boundaries used by SpliceAI to include our new transcript model by replacing the grch38.txt file in the SpliceAI package files. After this replacement, we ran SpliceAI from the command line using the parameter -A grch38 to specify the custom gene annotation file. <br>
After collecting the SNPs that affect the predicted splicing potential, a threshold of at least delta 0.5 for either gain or loss from donor sites or acceptor sites was applied. These filtered results were maintained in Table S15. For visualization in Figure 6a-b, we further required the genomic locations of the SNPs within 2 nt from the splice junctions and collapsed the observation at the transcript level. If a SNP affecting both annotated model (ENST) and unannotated models, only the annotated one was counted as positive. <br>
```
#Running on separate GPUs the predictions for the VCFs (thus using CUDA_VISIBLE_DEVICES)
#For eQTL
CUDA_VISIBLE_DEVICES=0 nohup python3 __main__.py -I fixed_eqtl.hg38.vcf -O spliceai_out/out_eqtl.hg38.vcf -R genome.fa -A grch38 > run0.log 2>&1 &

#For GWAS
CUDA_VISIBLE_DEVICES=1 nohup python3 __main__.py -I fixed_credible_set.hg38.vcf -O spliceai_out/out_fixed_credible_set.hg38.vcf -R genome.fa -A grch38 > run1.log 2>&1 &
```

# <a name="enrichment"></a>Enrichment of GWAS and eQTL SNPs in enhancers

To divide these results according to the enhancer structures (dCGI, CGI, Null and TATA), intersection was grouped at the ex5_cluster level. Exons were merged from the transcripts that are linked to the same ex5_clusters using bedtools merge. To include the most relevant SNPs, eQTL dataset was filtered for those derived from iPSC and brain related samples while GWAS dataset was filtered for traits related to brain. In order to obtain a normalized enrichment, “common_all_20180418” dataset downloaded from NIH was used as background (n = 37,302,987, unique SNP locations = 36,363,997). The eQTL and GWAS datasets were limited to this scope of SNPs (unique SNP locations = 81,782 & 35,325 respectively) The regions per ex5_cluster were intersected with the common_all dataset, the eQTL dataset and the GWAS dataset. Exons that were shared by more than one ex5_cluster were recounted. The results were normalized as the rate of hit compared to the background SNPs or as the number of SNPs per 1kb. For the rate of hit compared to background SNPs, ex5_clusters containing no background SNPs were excluded.









