#!/bin/bash
#PBS -q normal
#PBS -l walltime=4:00:00
#PBS -l mem=12GB
#PBS -l ncpus=48
#PBS -j oe
#PBS -N "repeatmasker"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04


module load RepeatMasker/4.1.5 parallel

base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
output_directory="/g/data/te53/zc4688/honours/analyses/chordates/repeatmasker"

sampleids=($(find "$base_directory" -type f -name "*rDNA.refmorph.fasta" | xargs -n 1 basename | cut -d '.' -f1))

run_repeatmasker(){
    base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
    output_directory="/g/data/te53/zc4688/honours/analyses/chordates/repeatmasker"

    sampleid=$1
    input_file="${base_directory}/${sampleid}.rDNA.refmorph.fasta"

    if [ ! -f "$input_file" ]; then
        echo "File not found: $input_file"
        return
    fi

    if [ -e ${output_directory}/${sampleid}.rDNA.refmorph.fasta.out.xm ]
    then
        echo "Skipping ${sampleid} as analysis is already done"
        return
    fi

    echo "Running command for $sampleid"
    RepeatMasker -species chordates -s -xm -gff -dir $output_directory $input_file
    
    echo "Repeatmasker complete for $sampleid"

    rm -f ${output_directory}/${sampleid}.rDNA.refmorph.fasta.masked
    rm -f ${output_directory}/${sampleid}.rDNA.refmorph.fasta.cat

}

export -f run_repeatmasker  # Export the function for parallel
parallel -j 12 run_repeatmasker ::: "${sampleids[@]}"



