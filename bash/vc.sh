#!/bin/bash
#PBS -q normal
#PBS -l walltime=24:00:00
#PBS -l mem=96GB
#PBS -l ncpus=24
#PBS -j oe
#PBS -N "cn"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04+gdata/xy86
#export PATH=/g/data/te53/zc4688/software/lofreq/bin:$PATH

module load samtools seqkit jellyfish minimap2 pythonlib bbmap



#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_pb -t 24 -f /g/data/te53/zc4688/honours/tmp/hg002/pb/ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/PacBio_fasta/*.fasta.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19
#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_ill -t 24 -f /g/data/te53/zc4688/honours/tmp/hg002/ill/*.fastq.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19
#python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_ont -t 24 -f /g/data/te53/zc4688/honours/tmp/hg002/ont/*.fastq.gz -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19
python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py  -o /g/data/te53/zc4688/honours/analyses/cn/hg002/ -s hg002_pb -t 24 -f /g/data/te53/zc4688/honours/tmp/hg002/pb/*.fastq -a /g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa -k 19

#lofreq call --call-indels -f /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta -o lofreq.vcf test.bam
#sed '1s/.*/>rdna/' chm13.rDNA.refmorph.fasta > chm13.rDNA.refmorph.simplified.fasta

#samtools view -H test.bam | sed 's/SN:chr21:3697875-3742662:+/SN:rdna/' > new_header.sam
#samtools reheader new_header.sam test.bam > test.renamed.bam
#bedtools genomecov -d -ibam test.bam -g /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta > coverage.txt

#lofreq call --call-indels -f /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta -o /g/data/te53/zc4688/honours/analyses/cn/result/fasta/lofreq.vcf /g/data/te53/zc4688/honours/analyses/cn/result/fasta/test.bam


#lofreq call -r rdna:0-11 -f /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.simplified.fasta -o /g/data/te53/zc4688/honours/analyses/cn/result/fasta/lofreq.vcf /g/data/te53/zc4688/honours/analyses/cn/result/fasta/small.bam



#samtools mpileup -f /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.simplified.fasta test.renamed.bam | bcftools call -mv -Oz -o rdna.vcf.gz

#bcftools mpileup --max-depth 6000 -f /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.simplified.fasta test.renamed.bam | bcftools call -mv -Ob -o calls.bcf
