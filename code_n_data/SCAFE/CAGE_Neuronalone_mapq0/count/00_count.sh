#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='ssCAGE_Neuron'

perl ../script/00_count.pl \
./output/ \
../aggregate/run_full/out/annotate/$tag/bed/$tag.cluster.coord.bed.gz \
../aggregate/00_20240619.CTSS_path.txt \
../aggregate/run_full/out/annotate/$tag/bed/$tag.CRE.coord.bed.gz \
$tag
