#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='ont_ss.Neuronalone'

perl ../script/00_count_cluster.pl \
./output/ \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
../aggregate/pool_runs/ctss_bed/00_pool_runs.ctss_list.tsv \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
$tag
