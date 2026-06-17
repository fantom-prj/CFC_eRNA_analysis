#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='ont_ss.Neuronalone'

perl ../script/00_count.pl \
./output_tCRE \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
./output_tCRE/00_pool_runs.ctss_list20250507.tsv \
../aggregate/run_full/out/annotate/$tag/bed/$tag.CRE.coord.bed.gz \
$tag
