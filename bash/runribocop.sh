#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=192GB
#PBS -l jobfs=100GB
#PBS -l ncpus=48
#PBS -j oe
#PBS -N "ribocop"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04

#Script to run commands in parallel

module use /g/data/if89/apps/modulefiles
module load minimap2/2.28 hmmer/3.4 pythonlib/3.9.2 samtools/1.19 nci-parallel/1.0.0 htslib/1.16

tasks_per_node=4
cpupertask=12
commandsfile="/g/data/te53/zc4688/honours/scripts/command_files/commands_hprc.txt"

/apps/openmpi/4.0.2/bin/mpirun -map-by ppr:$tasks_per_node:node:PE=$cpupertask /apps/nci-parallel/1.0.0/bin/nci-parallel --poll 0 --shell /bin/bash --input-file $commandsfile
