#!/bin/bash
cwd=`dirname $0`
cd $cwd
tag='Neuron_THP1'
gzip -dc ../subsample_100x_bed/*.bed.bgz | sort --parallel 10 -k1,1 -k2,2n -k6,6 | bgzip -c >./$tag.subsample_100x.bed.bgz
tabix -p bed ./$tag.subsample_100x.bed.bgz

gzip -dc ../subsample_10x_bed/*.bed.bgz | sort --parallel 10 -k1,1 -k2,2n -k6,6 | bgzip -c >./$tag.subsample_10x.bed.bgz
tabix -p bed ./$tag.subsample_10x.bed.bgz

gzip -dc ../bed/*.bed.bgz | sort --parallel 10 -k1,1 -k2,2n -k6,6 | bgzip -c >./$tag.bed.bgz
tabix -p bed ./$tag.bed.bgz
