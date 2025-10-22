#!/bin/bash
#PBS -q normal
#PBS -l walltime=2:00:00
#PBS -l mem=48GB
#PBS -l jobfs=10GB
#PBS -l ncpus=12
#PBS -j oe
#PBS -N "ribocop"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04

#Script to run ribocop on single sample
module use /g/data/if89/apps/modulefiles
module load minimap2/2.28 hmmer/3.4 pythonlib/3.9.2 samtools/1.19 nci-parallel/1.0.0 htslib/1.16

#usage: qsub -v outdir=/g/data/te53/zc4688/honours/analyses/tmp/,inputfasta=/g/data/te53/t2t2024/referenceresource/hprc/GCA_018852605.2_Q100_hg002v1.0.1.pat_genomic.fna,sampleid=hg002,speciesname=homo_sapien -P te53 /g/data/te53/zc4688/honours/scripts/ribocop.singlesample.sh 


python3 /g/data/te53/zc4688/honours/scripts/ribocop/ribocop.py -o $outdir -t 12 -if $inputfasta -s  $sampleid -l /g/data/te53/zc4688/honours/data/euk.hmm  -d $PBS_JOBFS