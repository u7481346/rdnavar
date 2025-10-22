#!/bin/bash

base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
gc_dir="${base_directory}/gc_content"

# Make sure output directory exists
mkdir -p "$gc_dir"

> $base_directory/refmorph.stats.fa

sampleids=($(find "$base_directory" -type f -name "*rDNA.refmorph.fasta" | xargs -n 1 basename | cut -d '.' -f1))

for sampleid in "${sampleids[@]}"; do
    seqkit sliding -s 50 -W 100 ${base_directory}/${sampleid}.rDNA.refmorph.fasta | seqkit fx2tab -n -g > ${base_directory}/gc_content/${sampleid}.gc.txt
    echo "GC content calculated for $sampleid"
    seqkit stats ${base_directory}/${sampleid}.rDNA.morphs.fasta >> refmorph.stats.txt
done