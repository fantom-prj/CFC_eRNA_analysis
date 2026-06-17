#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='ontCAGE.Neuron_THP1'

perl ../script/00_count.pl \
./output/ \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
../aggregate/pool_runs/ctss_bed/00_pool_runs.ctss_list.tsv \
../aggregate/run_full/out/annotate/$tag/bed/$tag.CRE.coord.bed.gz \
$tag
