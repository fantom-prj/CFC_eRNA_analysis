# TES features for non-poly(A) RNA 

* TES_main_analysis.R -> main code describing all the analyses
* CGI_TES_.R -> associated code for distribution of CGI across TES regions
* MYC_binding_intersect.R -> associated code for intersection and distribution of MYC binding across TES regions
* ./RNAfold/output/RNAfold_run.R -> associated code for running RNAfold and collecting base-pairing probability
* Results are presented in Fig.5 and ext_Fig.7
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [RNA structural prediction proximal to transcript end site](#RNAfold)
* [Genomic and Regulatory Characterization of TES](#TESenrichment)


# <a name="RNAfold"></a>RNA structural prediction proximal to transcript end site

To assess structural features at transcription termini, the TES with the highest read count was selected from each ex3_cluster. These were categorized as poly(A) or non-poly(A) based on the presence of canonical Polyadenylation Signals (PAS) and the poly(A) classifier. We anchored our classification thresholds using an optimal cutoff derived from a THP-1 ROC curve (Ext_Fig. 1h). By linking the ex5_cluster and ex3_cluster through the transcript models, the selected TESs were then grouped by ex5_clusters (mRNA, p_ncRNA, e_ncRNA, CTCF_ncRNA and other_ncRNA) and regulatory elements. If more than one ex5_cluster groups of regulatory element groups were found, the TESs were excluded from the analyses. The TESs were further grouped into recursive (≥ 3 counts) or non-recursive according to the raw counts.
Secondary structure propensities were predicted for sequences flanking the TES (+/-200 nt from the TES). We utilized RNAfold (v2.6.4) of the ViennaRNA package using the partition function (-p) to calculate the Boltzmann-weighted equilibrium ensemble. This approach accounts for the structural diversity of eRNAs rather than relying solely on a single Minimum Free Energy conformation. Thermodynamic paired probability (P) was calculated per nucleotide per alignment from the base-pairing probability matrix (Pij) provided in the PostScript dot plots. These values were summarized into relative hairpin score for each TES group (eg, non-poly(A) e_ncRNA with 3 reads) aligned by the 3’ end by taking the average.
```
RNAfold --noPS -p [fasta_file] > [text_file]
p(i)=∑_(i≠j)▒P_ij 
```
To predict for a similar pattern of structural depletion for individual TES, mean paired probability of the RNA body (-100 to -6 nt of TES) and the structural valley (-5 to +14 nt of TES) were calculated. Predicted structural depletion was considered if the mean score for the RNA body from ≥ 0.55 and the depletion fold change (1 - valley / RNA body) ≥ 0.15. 


# <a name="TESenrichment"></a>Genomic and Regulatory Characterization of TES

To investigate the epigenetic landscape of transcription termination, we analyzed the enrichment of CGI and MYC ChIP-seq data relative to TESs. CGI annotations for the hg38 genome were retrieved from the UCSC Genome Browser. MYC ChIP-seq peak regions derived from human embryonic stem cells (GSM1505809) was lifted over to hg38 and used. Using the recursive TESs, we extended the 1-nt terminal locus 5 kb up- and down-stream and intersected these windows with CGI coordinates and MYC binding regions. For each transcript category, we calculated the percentage coverage of CGI sequences to identify specific enrichment patterns. 
For both MYC binding and GCI coordinates, we intersected the 1-nt recursive TES loci with MYC binding regions across four promoter-types (mRNA, p_ncRNA, e_ncRNA and other_ncRNA) to determine presence of co-localization. Termini were further stratified by their polyadenylation status (poly(A)+ vs. poly(A)-) based on our standardized classification criteria. 

