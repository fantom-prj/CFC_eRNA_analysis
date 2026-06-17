#!/bin/bash
cwd=`dirname $0`
cd $cwd

gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTSS_clusters2/ontCAGE.all.aggregate.collapse.ctss.bed.gz | bgzip -c >./ontCAGE.iPSchro.end5.signal.bed.bgz
gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTSS_clusters2/ontCAGE.all.cluster.coord.bed.gz | bgzip -c >./ontCAGE.iPSchro.end5.cluster.bed.bgz
cat /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTES_clusters2/end3_bed_bigwig/iPSchro.end3.bed | bgzip -c >./ontCAGE.iPSchro.end3.signal.bed.bgz
gzip -dc /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/CTES_clusters2/scafe/cluster/iPSchro.CTES.s3_n3_c5/bed/iPSchro.CTES.s3_n3_c5.tssCluster.bed.gz | bgzip -c >./ontCAGE.iPSchro.end3.cluster.bed.bgz

tabix -p bed ./ontCAGE.iPSchro.end5.signal.bed.bgz
tabix -p bed ./ontCAGE.iPSchro.end5.cluster.bed.bgz
tabix -p bed ./ontCAGE.iPSchro.end3.signal.bed.bgz
tabix -p bed ./ontCAGE.iPSchro.end3.cluster.bed.bgz

cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/junction_extractor2/pool/iPSchro.full.hi_qual.junct.bed ./iPSchro.long_read.hi_qual.junct.bed
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/other_junctions/gencode_v39_junct.bed ./Gencode_v39.junct.bed
cp /analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_iPS_chro/other_junctions/short_read_junct.bed ./iPSC_NSC_Neuron.short_read.hi_qual.junct.bed