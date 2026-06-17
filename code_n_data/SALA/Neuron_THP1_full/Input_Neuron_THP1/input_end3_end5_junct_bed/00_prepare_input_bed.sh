#!/bin/bash
cwd=`dirname $0`
cd $cwd

gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/CTSS_clusters/ontCAGE.Neuron_THP1.aggregate.collapse.ctss.bed.gz | bgzip -c >./ontCAGE.Neuron_THP1.end5.signal.bed.bgz
gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/CTSS_clusters/ontCAGE.Neuron_THP1.cluster.coord.bed.gz | bgzip -c >./ontCAGE.Neuron_THP1.end5.cluster.bed.bgz
cat /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/CTES_clusters/end3_bed_bigwig/Neuron_THP1.end3.bed | bgzip -c >./ontCAGE.Neuron_THP1.end3.signal.bed.bgz
gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/CTES_clusters/scafe/cluster/Neuron_THP1.CTES.s3_n3_c5/bed/Neuron_THP1.CTES.s3_n3_c5.tssCluster.bed.gz | bgzip -c >./ontCAGE.Neuron_THP1.end3.cluster.bed.bgz

tabix -p bed ./ontCAGE.Neuron_THP1.end5.signal.bed.bgz
tabix -p bed ./ontCAGE.Neuron_THP1.end5.cluster.bed.bgz
tabix -p bed ./ontCAGE.Neuron_THP1.end3.signal.bed.bgz
tabix -p bed ./ontCAGE.Neuron_THP1.end3.cluster.bed.bgz

cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/junction_extractor/pool/Neuron_THP1.full.hi_qual.junct.bed ./Neuron_THP1.long_read.hi_qual.junct.bed
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/other_junctions/gencode_v39_junct.bed ./Gencode_v39.junct.bed
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/other_junctions/short_read_junct.bed ./iPSC_NSC_Neuron.short_read.hi_qual.junct.bed