#!/bin/bash
cwd=`dirname $0`
cd $cwd

cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/annotate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.cluster.coord.bed.gz ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/aggregate/ontCAGE.Neuron_THP1/bed/ontCAGE.Neuron_THP1.aggregate.collapse.ctss.bed.gz ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/ctss_to_bigwig/ontCAGE.Neuron_THP1.all/wig/ontCAGE.Neuron_THP1.all.count.rev.bw ./
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/Neuron_THP1_S3/ontCAGE/aggregate/run_full/out/ctss_to_bigwig/ontCAGE.Neuron_THP1.all/wig/ontCAGE.Neuron_THP1.all.count.fwd.bw ./