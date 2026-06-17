# pre-processing of the CFC-seq

* isoquant.process.R ->  code describing downstream processing from Isoquant output
* talon.process.R -> code describing downstream processing from Talon output
* Results are used in ext_Fig.4
* Data is not included here, please refer to https://fantom.gsc.riken.jp/6/suppl/Yip_et_al_2026_CFC

# Related Methods
* [Pre-processing of the long-read ONT data](#preprocess)
* [Mapping to genome and TranscriptClean](#mapping)


# <a name="preprocess"></a>Pre-processing of the long-read ONT data

Basecalling was performed using Dorado (v0.2.4, Oxford Nanopore Technologies) to generate FASTQ files from raw FAST5 data files using the high-accuracy model “dna_r9.4.1_e8_sup@v3.3”. Reads were filtered for quality, retaining only those with a mean Q-score > 10 for the downstream analyses. To ensure high-quality transcript models, adapter sequences were trimmed and strands were oriented using primer chop (https://gitlab.com/mcfrith/primer-chop). Only reads containing confirmed adaptors at both ends were retained, while the relative position of the head and tail linkers was used to assign read orientation. Subsequent to adaptor removal, poly(A) tails including those basecalled as mixed adenine and guanine dinucleotides were identified, recorded and trimmed. We developed tail trimmer (v1.4), which utilizes a 20 bp sliding window to identify and trim poly(A) sequences, allowing a maximum of 5 terminal adenines to remain.
```
# Basecalling for Neuron & THP1 series: 
dorado basecaller -x cuda:all --min-qscore 10 --emit-fastq dna_r9.4.1_e8_sup@v3.3 $POD5 | gzip > $POD5/$POD5.fastq.gz
# Basecalling for iPS chromatin bound: 
dorado basecaller -x cuda:all --min-qscore 10 --emit-fastq dna_r10.4.1_e8.2_400bps_sup@v4.1.0 $POD5 | gzip > $POD5/$POD5.fastq.gz

# Pimer chop For Neuron series
primer-chop \
    -q \
    -P 12 \
    riken-yonsei-primers_minus_poly-A.fa \
    "dorado_basecaller/${LINE}.fastq.gz" \
    dorado_primer-chop_minus_poly-A
    
# Re-orentiation
sed -i '1~4 s/$/ strand=-/' dorado_primer-chop_minus_poly-A/good-rev.fq
sed -i '1~4 s/$/ strand=+/' dorado_primer-chop_minus_poly-A/good-fwd.fq

cat dorado_primer-chop_minus_poly-A/good-fwd.fq \
    dorado_primer-chop_minus_poly-A/good-rev.fq \
    > dorado_primer-chop_minus_poly-A/good.fq

# tail trimmer
tail_trimmer_v1.4.pl \
  --rescue_tail_nt=5 \
  --revise_read_prefix="$LINE" \
  --fastq_path="$LINE/dorado_primer-chop_minus_poly-A/good.fq.gz" \
  --out_dir="$LINE/dorado_primer-chop_minus_poly-A_tail-trimmer"
```

# <a name="mapping"></a>Mapping to genome and TranscriptClean
The trimmed and orientated reads were then mapped to the human reference genome (GRCh38, GCA_000001405.15), using minimap2 (v2.17-r974-dirty), guided by GENCODE (v39) annotations. Raw alignments were error-corrected using TranscriptClean (v2.0.3). Splice junction (SJ) correction was guided by a reference set comprising GENCODE v39 junctions and high confidence SJs identified short-read RNA seq of the Neuron series. At the locus-level, TranscriptClean removed 40,024 SJs and introduced 8,107 SJs while 1,654,272 (97.17%) remained unchanged. As a requirement for running TALON, an additional tag was added to the corrected alignments using the TALON (v5) talon_label_reads function which calculates the fraction of adenines in a 16-bp window within the genome immediately downstream of the last 3’-bp alignment of the read. However, internal priming was defined universally across different assemblers as described in the next section. 
```
[1] minimap2 -I 1000G -k 15 -d reference.fa.mmi reference.fa
[2] minimap2 -t [threads] -2 -ax splice -uf --MD --junc-bed gencode.v39.annotation.bed --secondary=no GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.mmi [fastq_file]> [sam_file]
[3] python TranscriptClean.py -t [threads] --sam [sam_file] --genome GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta -j gencode.v39.annotation.SJs.txt --outprefix sample --deleteTmp --tmpDir sample/tmp
[4] talon_label_reads --f pass_trim_clean_corrected.sam --ar 16 --fracA 0.5 --g GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta --t 12 --tmpDir sample/tmp --deleteTmp --o sample
```
