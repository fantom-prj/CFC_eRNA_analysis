#!/bin/bash
cwd=`dirname $0`
cd $cwd

scafe_cluster='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/scripts/scafe.tool.cm.cluster'
scafe_ctss_to_bigwig='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/scripts/scafe.tool.cm.ctss_to_bigwig'
ONT_CAGE_neuron_CTES_path='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTES_clusters2/end3_bed_bigwig/iPSchro.end3.bed'
outputPrefix='iPSchro.CTES'

$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--count_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--min_summit_count=2 \
--min_nt_count=2 \
--min_cluster_count=3 \
--outputPrefix=$outputPrefix.s2_n3_c3 \
--outDir=../scafe/cluster

$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--count_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--min_summit_count=3 \
--min_nt_count=3 \
--min_cluster_count=5 \
--outputPrefix=$outputPrefix.s3_n3_c5 \
--outDir=../scafe/cluster

$scafe_cluster \
--overwrite=yes \
--cluster_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--count_ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--min_summit_count=5 \
--min_nt_count=5 \
--min_cluster_count=10 \
--outputPrefix=$outputPrefix.s5_n5_c10 \
--outDir=../scafe/cluster

$scafe_ctss_to_bigwig \
--ctss_bed_path=$ONT_CAGE_neuron_CTES_path \
--genome=hg38.gencode_v39 \
--outputPrefix=$outputPrefix \
--outDir=../scafe/ctss_to_bigwig/
