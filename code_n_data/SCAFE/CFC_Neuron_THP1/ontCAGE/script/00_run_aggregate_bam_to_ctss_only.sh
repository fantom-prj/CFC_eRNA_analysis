#!/bin/bash
lib_list_path=$1
#/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/00_cell_type.NSC.aggregate_list.tsv
genome=$2
#hg38.gencode_v39
baseDir=$3
#/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231110/Neuron_THP1/ontCAGE/aggregate/
outputPrefix=$4
#cell_type.NSC
scriptDir="/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/scripts"

$scriptDir/scafe.tool.cm.aggregate \
--lib_list_path=$lib_list_path \
--max_thread=5 \
--genome=$genome \
--outputPrefix=$outputPrefix \
--outDir=$baseDir/out/aggregate

$scriptDir/scafe.tool.cm.ctss_to_bigwig \
--genome=$genome \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz \
--outputPrefix=$outputPrefix.all \
--outDir=$baseDir/out/ctss_to_bigwig

$scriptDir/scafe.tool.cm.ctss_to_bigwig \
--genome=$genome \
--ctss_bed_path=$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz \
--outputPrefix=$outputPrefix.ung \
--outDir=$baseDir/out/ctss_to_bigwig

