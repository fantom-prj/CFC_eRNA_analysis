#!/bin/bash
cwd=`dirname $0`
cd $cwd

cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.all/bed/ontCAGE.all.cluster.coord.bed.gz ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/run_full/out/aggregate/ontCAGE.all/bed/ontCAGE.all.aggregate.collapse.ctss.bed.gz ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/run_full/out/ctss_to_bigwig/ontCAGE.all.all/wig/ontCAGE.all.all.count.rev.bw ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/run_full/out/ctss_to_bigwig/ontCAGE.all.all/wig/ontCAGE.all.all.count.fwd.bw ./