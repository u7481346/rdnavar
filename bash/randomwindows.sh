#!/bin/bash

#Script to compare variance in rDNA copy number calculations across libraries to variation in random 45kb windows.

module load bedtools samtools
#bedtools random -l 45000 -n 1000 -g /g/data/te53/humanreference/chm13-t2t/chm13-t2t.ebv.phix.chrQ.xy.fa.fai | cut -f1,2,3 > /g/data/te53/zc4688/honours/analyses/cn/illumina_libs/random.bed

for bam in /g/data/te53/zc4688/honours/analyses/cn/illumina_libs/library/HG001/*.bam; do
    libname=$(basename "$bam"| cut -d '.' -f1)

    samtools bedcov /g/data/te53/zc4688/honours/analyses/cn/illumina_libs/random.bed "$bam" | awk '{print "random" "\t"  $4/45000}' >> "/g/data/te53/zc4688/honours/analyses/cn/illumina_libs/library/HG001/${libname}.depth.txt"
done
