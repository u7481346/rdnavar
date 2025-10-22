#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=4GB
#PBS -j oe
#PBS -N "iqtree"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04


module load iqtree2/2.1.2
echo "input = $input"
echo "threads = $threads"
echo "prefix = $p"

#usage: qsub -v input=/g/data/te53/zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/cleaned.alignment.txt,threads=4,p=/g/data/te53/zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/eighteen -l ncpus=4 -P te53 iqtree.sh

iqtree2 -s "$input" -T "$threads" -pre "$p"