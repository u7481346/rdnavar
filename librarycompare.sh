#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l jobfs=1GB
#PBS -j oe
#PBS -N "cn"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04+gdata/xy86

module load samtools/1.22 bwa-mem2/2.2.1 parallel

out_dir="/g/data/te53/zc4688/honours/analyses/cn/hg002"
#align reads to chm13

if [ ! -f "${out_dir}/${libname}.bam.bai" ]; then
    echo "No BAM index found – running alignment and indexing"
    echo "Input fastq = ${infastq}"
    echo "Output directory = ${out_dir}"
    echo "libname=${libname}"
    
    bwa-mem2 mem -t ${PBS_NCPUS} /g/data/te53/humanreference/chm13-t2t/bwa-mem2/chm13-t2t.ebv.phix.chrQ.xy.fa ${infastq} | samtools sort -@ ${PBS_NCPUS} -o ${out_dir}/${libname}.bam
    samtools index -@ ${PBS_NCPUS} ${out_dir}/${libname}.bam
else
    echo "Index already exists for ${libname} – skipping alignment"
fi


#calculate coverage and calculate median
#1kb
#window=1000 bam=${out_dir}/${libname}.bam outdir=${out_dir} bash /g/data/te53/zc4688/honours/scripts/bamdepth.sh

#sort -k4,4n "${out_dir}/${libname}.1000.depth.bed" | awk ' { a[i++]=$4; } END { print "1kb" "\t" a[int(i/2)]; }' > "${out_dir}/${libname}.depth.txt"
 
#10kb
window=10000 bam=${out_dir}/${libname}.bam outdir=${out_dir} bash /g/data/te53/zc4688/honours/scripts/bamdepth.sh
sort -k4,4n "${out_dir}/${libname}.10000.depth.bed" | awk ' { a[i++]=$4; } END { print "10kb" "\t" a[int(i/2)]; }' >> "${out_dir}/${libname}.depth.txt"

#50kb
window=50000 bam=${out_dir}/${libname}.bam outdir=${out_dir} bash /g/data/te53/zc4688/honours/scripts/bamdepth.sh
sort -k4,4n "${out_dir}/${libname}.50000.depth.bed" | awk ' { a[i++]=$4; } END { print "50kb" "\t" a[int(i/2)]; }' >> "${out_dir}/${libname}.depth.txt"


#find rdna reads
samtools bedcov /g/data/te53/zc4688/honours/analyses/cn/chm13.rdnaregions.bed "${out_dir}/${libname}.bam" | awk '{bases += $4} END {print "rdna" "\t" bases/44787}' >> "${out_dir}/${libname}.depth.txt"

touch "${out_dir}/${libname}.done"

#test is HG002_140528_D00360_0018_AH8VC6ADXX