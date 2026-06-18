# Analyses on all the tCREs identified from CFC-seq

* fig1.analysis.R -> main code describing all the analyses
* TES_polyA.R -> associated code for poly(A) prediction and 3' ends assessment
* longread.*.sh -> processing other long-read protocols from fastq to mRNA 5' ends
* ./fantom_random_forest20240407/RunForest.Rmd -> poly(A) prediction by random forest [package available in SALA(https://github.com/fantom-prj/SALA)]
* Results are presented in ext_Fig.1
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Comparative analysis of long-read 5’ end precision](5end)
* [Assessment of 3’ ends capture and filtration of potential internal priming](#3end)
* [Computational identification and hierarchical scoring of PAS motifs](#PAS)
* [Capture of non-poly(A) ncRNA and incomplete transcripts derived from PAT protocol](#PAT)
* [Poly(A)-tail prediction](#predict)


# <a name="5end"></a>Comparative analysis of long-read 5’ end precision

To benchmark the 5’ end accuracy of CFC-seq against established methodologies, raw FASTQ datasets from five distinct long-read sequencing protocols targeting the same iPSC line (WTC11) were obtained from ENCODE LRGASP consortium. These included CAP_trap_PacBio (ENCSR309IKK), uncap_deplete_PacBio (ENCSR507JOF), TSO_ONT (ENCSR539ZXJ), R2C2_ONT (ENCSR925UQZ) and dRNA_ONT (ENCSR392BGY). All reads were aligned to hg38 by minimap2 as described above. The alignments of these libraries and our iPSC libraries were intersected with GENCODE v39 protein-coding exons using bedtools intersect. Reads overlapping with ≥ 90% of an annotated protein-coding exon were classified as mRNA and retained for downstream analysis, excluding mitochondrial transcripts. <br>
To quantify transcriptional initiation precision, the empirical 5’ end coordinates of these alignments were then intersected with three independent genomic reference datasets defining promoter and active regulatory boundaries: 1) candidate cis-regulatory elements (cCREs) from the SCREEN registry annotated as promoter-like or enhancer-like signatures, 2) open chromatin peaks derived from single-nucleus ATAC-seq data generated from the identical WTC11 iPSC line, and 3) robust CAGE-defined promoter clusters obtained from the FANTOM5 consortium (median width 300bp). Comprehensive genomic coordinates for these reference annotations are archived in Supplementary table S18. <br>


# <a name="3end"></a>Assessment of 3’ ends capture and filtration of potential internal priming

To evaluate the efficacy of in vitro poly(A)-tailing (PAT) for non-poly(A) transcript capture, we performed comparative library testing using THP-1 cells processed with or without the PAT enzyme cascade. Natural non-poly(A) control genes, including small nucleolar RNAs (snoRNAs) and replication-dependent histone gene blocks , were utilized to evaluate non-poly(A) RNA target rescue efficiency. <br>
The single nucleotide 3’ end coordinates were extracted from all alignments, and the flanking genomic sequences (± 50 nt) were retrieved from the human reference genome GRCh38 (GCA_000001405.15). To mitigate internal priming artifacts caused by oligo-dT hybridization to homopolymeric adenine stretches within transcript bodies, these flanking regions were profiled for downstream adenine composition. A transcript was flagged if adenine composition exceeded 50% within 16 nucleotides immediately downstream of the 3’ end, or exceeding 75% within the 8 nucleotides immediately downstream. Termini directly corresponding to annotated GENCODE v39 transcript models were exempt from this filtering cascade. <br>


# <a name="PAS"></a>Computational identification and hierarchical scoring of PAS motifs
For polyadenylation signal (PAS) motif identification, we performed a genome-wide search using scanMotifGenomeWide.pl from Homer (v4.11) with a position weight matrix derived from GENCODE v39 annotation. Only the loci with score > 3 were used for the analysis. The genomic locations of these PAS motifs were overlapped with the single-nucleotide locations 25 nt upstream of the 3’ ends by running bedtools closest. Presence of PAS motif with position of its last nucleotide located -5 nt to -35 nt from the 3’ ends was annotated as PAS-positive. We applied a hierarchical search strategy based on motif strength, starting from PAS motifs with score >6 ("AATAAA" "ATTAAA"), followed by PAS motifs with score > 4 ("TATAAA", "AGTAAA", "AATATA"), and by motifs with score >3 ("CATAAA", "GATAAA", "AAAAAA", "TTTAAA", "ACTAAA", "AATACA", "AATAGA", "AAGAAA", "AATAAG", "AATAAT", "AATGAA", "AATTAA", "ATTATA"). 
```
[1] scanMotifGenomeWide.pl PAS.motif hg38 -bed -keepAll -p 10 > output.bed
[2] bedtools closest -a all3n_up25.sort.bed.gz -b GENCODE_polyA_signal.motif_rng.main.FASTA.6.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS6.bed.gz
[3] bedtools closest -a all3n_up25.sort.bed.gz -b GENCODE_polyA_signal.motif_rng.main.FASTA.4.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS6.bed.gz
[4] bedtools closest -a all3n_up25.sort.bed.gz -b GENCODE_polyA_signal.motif_rng.main.FASTA.3.bed.gz -s -D a | cut -f-6,10- | gzip > all3n_up25.PAS6.bed.gz
```

# <a name="PAT"></a>Capture of non-poly(A) ncRNA and incomplete transcripts derived from PAT protocol

Following the exclusion of known non-poly(A) reference transcripts and loci affected by internal priming, the baseline transcript capture layout of non-PAT standard poly-dT protocol was compared with PAT protocol (Ext_Fig. 1g). To investigate whether the PAT protocol increases the prevalence of incomplete transcript models, protein-coding mRNAs lacking a PAS motif were classified based on genomic coordinations of their 3’ ends. Finally, PAT-specific 3’ends were isolated by subtracting the 3’ end coordinates detected in the non-PAT control. These 3’ ends were intersected with different genomic regions derived from annotations of GENCODE protein-coding transcripts.


# <a name="predict"></a>Poly(A)-tail prediction
In order to call bona fide 3’ ends for the transcripts detected by CFC-seq, we trained a random forest classifier on a curated database of 3’ ends identified by FLAM-seq and 3p-seq. Briefly, we collected all 3’ ends detected by CFC-seq and flagged those overlapping with the curated (hereon “reference”) dataset within a window of 10 nt. We then trained a random forest classifier using the ranger package (v. 0.16.0) in R 4.3.1, considering the 3’ ends found in the reference as true positives and defining a set of features to train the model on, including the presence, type and position of polyadenylation signals and the nucleotide composition of a 50 nt window upstream and downstream of each 3’ end. We then used the trained model to assign a probability score (“poly(A) score”) to each 3’ end in the CFC-seq dataset. We defined an upper threshold at which 95% of 3’ ends found in the reference and containing a polyadenylation signal were classified as true to define “poly(A)” 3’ ends (threshold > 0.34), and a lower threshold at which 95% of all 3’ ends were classified as true to flag “non-poly(A)” 3’ ends (threshold < 7.6e-5). For downstream analyses, PAS-positive 3’ ends and 3’ ends classified as true poly(A) were considered as poly(A) 3’ ends while the others were considered as non-poly(A). The read-based analysis result was transferred to the transcript-level and maintained in Table S4.

