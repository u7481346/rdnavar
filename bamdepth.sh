#!/bin/bash
#PBS -N bamdepth
#PBS -P te53
#PBS -q normal
#PBS -l storage=gdata/te53+gdata/if89
#PBS -l walltime=2:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l jobfs=400GB
#PBS -l wd
#PBS -j oe 

#copy of refgen bamdepth.sh but allows sam or bam
#usage qsub -v bam=/path/to/bam,window=windowsize,outdir=/path/to/output
#calculates depth over a window size on a bam file

set -ex

export WINDOW=${window}
bambase="$(basename "${bam}" | sed 's/\.\(bam\|sam\)$//')"

# Create BED file
samtools view -H ${bam} | perl -lne 'if ($_=~/SN:(\S+)\tLN:(\d+)/){ $c=$1;$l=$2; for ($i=0;$i<$l;$i+=$ENV{"WINDOW"}) { print "$c\t$i\t". ((($i+$ENV{"WINDOW"}) > $l) ? $l : ($i+$ENV{"WINDOW"}))  }} ' > "${PBS_JOBFS}/${bambase}.${window}.bed"

# Function to calculate depth for each BED region
calculate_depth() {
    region=$1
    samtools bedcov <(echo "$region") ${bam} | awk -v window=${WINDOW} '{print $1, $2, $3, $4/window}'
}

export -f calculate_depth

# Run depth calculation in parallel
parallel --will-cite -a "${PBS_JOBFS}/${bambase}.${window}.bed" -j ${PBS_NCPUS} calculate_depth > "${PBS_JOBFS}/${bambase}.${window}.depth.tmp.bed"
sort -k1,1 -k2,2n "${PBS_JOBFS}/${bambase}.${window}.depth.tmp.bed" > "${outdir}/${bambase}.${window}.depth.bed"
touch ${outdir}/${bambase}.${window}.depth.bed.done
