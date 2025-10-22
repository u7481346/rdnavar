#!/bin/bash

base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"

sampleids=($(find "$base_directory" -type f -name "*rDNA.refmorph.fasta" | xargs -n 1 basename | cut -d '.' -f1))

for sampleid in "${sampleids[@]}"; do
    trf ${base_directory}/${sampleid}.rDNA.refmorph.fasta 2 7 7 80 10 50 2000 -h -m
done

#easier to just run on single combined fasta 
#need masked file for meme analysis
#doesnt actually need to be qsubbed at all
#trf ../new/msa/combined.fasta 2 7 7 80 10 50 2000 -h -m