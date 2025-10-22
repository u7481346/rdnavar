#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=48GB
#PBS -l ncpus=48
#PBS -l jobfs=10GB
#PBS -j oe
#PBS -N "trnascan"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04


module use -a /g/data/if89/shpcroot/modules
module load quay.io/biocontainers/funannotate/1.8.15--pyhdfd78af_2 singularity parallel


indir="/g/data/te53/t2t2024/referenceresource/hprc/"
outdir="/g/data/te53/zc4688/honours/analyses/cn/trna/"
files=($(find "$indir" -type f -name "*genomic.fna"))

run_trnascan(){
    
    file=$1
    sampleid=$(basename "$file" | cut -d '.' -f1)

    if [ ! -f "$file" ]; then
        echo "File not found: $file"
        return
    fi

    if [ -e ${outdir}/${sampleid}.done ]
    then
        echo "Skipping ${sampleid} as analysis is already done"
        return
    fi

    echo "Running command for $sampleid"
    tRNAscan-SE --thread 12 -m ${outdir}/${sampleid}.summary.txt -o ${outdir}/${sampleid}.results.txt -f ${outdir}/${sampleid}.structures.txt $file
    
    echo "tRNAscan complete for $sampleid"

    touch ${outdir}/${sampleid}.done

}

export -f run_trnascan 
export indir outdir # Export the function for parallel
parallel -j 4 run_trnascan ::: "${files[@]}"



