#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='ontCAGE.Neuron_THP1'

perl ../script/00_count.pl \
./output/ \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
../aggregate/pool_runs/ctss_bed/00_pool_runs.ctss_list.tsv \
/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/20231113_Neuron_THP1_para_updated/Neuron_THP1.full/bed/Neuron_THP1.S3.end5.bed.gz \
$tag
