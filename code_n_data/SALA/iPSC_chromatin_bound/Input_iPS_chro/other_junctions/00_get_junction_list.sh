#!/bin/bash
cwd=`dirname $0`
cd $cwd

gencode_v39_junction='/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/resources/v1.0.1/resources/genome/hg38.gencode_v39/bed/intron.bed.gz'
gzip -dc $gencode_v39_junction >./gencode_v39_junct.bed

short_read_junction='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/valeria_short_read_junction/neuron_series.SJ.out.unique1.gencode_back.bed'
cp $short_read_junction ./short_read_junct.bed
