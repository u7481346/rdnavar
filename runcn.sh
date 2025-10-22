#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=700GB
#PBS -l ncpus=192
#PBS -j oe
#PBS -N "cn"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04+gdata/xy86

#Script to run commands in parallel

module use /g/data/if89/apps/modulefiles
module load samtools seqkit jellyfish minimap2 pythonlib bbmap

tasks_per_node=4
cpupertask=12
commandsfile="/g/data/te53/zc4688/honours/scripts/command_files/cn_commands.txt"

/apps/openmpi/4.0.2/bin/mpirun -map-by ppr:$tasks_per_node:node:PE=$cpupertask /apps/nci-parallel/1.0.0/bin/nci-parallel --poll 0 --shell /bin/bash --input-file $commandsfile

#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_pb -t 4 -f /g/data/te53/zc4688/honours/tmp/hg002/pb/ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta/*.fasta.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19

#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_ill -t 24 -f /g/data/te53/zc4688/honours/tmp/hg002/ill/*.fastq.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19

#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_ont -t 4 -f /g/data/te53/zc4688/honours/tmp/hg002/ont/*.fastq.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19
