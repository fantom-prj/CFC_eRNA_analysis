# Methodology for "Genomic codes governing enhancer RNA fate"


Long-read sequencing has transformed transcriptome profiling, yet capturing full-length, non-polyadenylated transcripts like enhancer RNAs (eRNAs) remains challenging. Here, we introduce CFC-seq, combining cap-trapping and in vitro poly(A)-tailing to sequence poly(A) and non-poly(A) RNAs with precise transcription start site. Paired with our assembler, SALA, we identified 39,425 novel transcriptional units, including ~24,000 eRNAs. Our data reveal a distinct genomic code governing eRNA fate dictated by core promoter architecture. CpG-island enhancers show high chromatin connectivity but yield short, exosome-sensitive RNAs. Conversely, TATA-box enhancers systematically co-opt LTR retrotransposons to inherit structural motifs that produce long, stable, and spliced RNAs. Mechanistically, the pioneer factor NF-Y activates these viral elements to license transcription, balanced by TEAD4 activity across a dual-gear regulatory axis. Finally, non-poly(A) eRNAs terminate via exosome-associated processing at structural-depleted cleavage zones. This comprehensive annotation links enhancer sequence architecture to RNA fate, providing a new transformative framework for decoding the functional human genome.


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



