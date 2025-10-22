#!/bin/bash
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=100GB
#PBS -l ncpus=24
#PBS -j oe
#PBS -N "extractseqs"
#PBS -o /g/data/te53/zc4688/honours/logs
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04

#Usage - generates commands file to be used by run_ribocop.sh (required for parallelisation)
module load pythonlib samtools parallel


base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
output_directory="$base_directory/retest"

sampleids=($(find "$base_directory" -type f -name "*ribocop.done" | xargs -n 1 basename | cut -d '.' -f1))

extractseqs(){
    
    sample=$1
    infasta=$(find /g/data/xl04/genomeprojects/referencedata/tmp/chordata/ -type f -name "${sample}*.fna.gz" | head -n 1)

    if [ ! -f "$infasta" ]; then
        echo "FASTA not found: $infasta"
        return
    fi

    if [ -e ${output_directory}/${sample}.structure.bed ]
    then
        echo "Skipping ${sample} as analysis is already done"
        return
    fi
    

    python3 /g/data/te53/zc4688/honours/scripts/extractseqs.py -o "$output_directory" -s "$sample"  --inputfasta "$infasta" -i "$base_directory"
    echo "Sequences extracted for $sample"
}

export base_directory output_directory
export -f extractseqs  # Export the function for parallel
parallel -j 6 extractseqs ::: "${sampleids[@]}"



