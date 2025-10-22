#!/bin/bash
#PBS -q normal
#PBS -l walltime=12:00:00
#PBS -l mem=700GB
#PBS -l ncpus=192
#PBS -j oe
#PBS -N "kmers"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89+gdata/xl04+gdata/xy86


module load conda jellyfish seqtk R parallel bbmap samtools seqkit
source $CONDA_INIT_SCRIPT
conda activate /g/data/te53/condapkg/envs/emboss_v660

outdir="/g/data/te53/zc4688/honours/analyses/cn/result"
sampleids=($(find "/g/data/xy86/ont1kgp/" -type f -name "*bam" | xargs -n 1 basename | cut -d '.' -f1))

find_kmers(){
    sampleid=$1
    inbam="/g/data/xy86/ont1kgp/${sampleid}*bam"

    if [ -e ${outdir}/${sampleid}.${gene}.fa ]; then
        echo "Skipping ${sampleid} as analysis is already done"
        return
    fi

    if [ ! -e "${outdir}/${sampleid}_rdna.fa.gz" ] && [ ! -e "${outdir}/${sampleid}_rdna.fq.gz" ]; then
        echo "Generating FASTA from BAM for ${sampleid}"
        
        samtools fasta --threads 6 ${inbam} | bbduk.sh in=stdin.fa outm=${outdir}/${sampleid}_rdna.fa ref=/g/data/te53/zc4688/honours/analyses/cn/final_rdna_kmers.fa k=19 t=6 int=f zl=6 -Xmx1g
    fi

    if [ ! -e "${outdir}/${sampleid}_rdna.fa" ]; then
        echo "Converting to FA for ${sampleid}"
        if [ -e "${outdir}/${sampleid}_rdna.fq.gz" ]; then
            gunzip -c "${outdir}/${sampleid}_rdna.fq.gz"  | seqtk seq -A - > "${outdir}/${sampleid}_rdna.fa"
            rm "${outdir}/${sampleid}_rdna.fq.gz"
        elif [ -e "${outdir}/${sampleid}_rdna.fa.gz" ]; then
            gunzip -c "${outdir}/${sampleid}_rdna.fa.gz" > "${outdir}/${sampleid}_rdna.fa"
        fi
    fi

    echo "Looking for kmers in ${sampleid}"


    cat ${outdir}/${gene}kmers.patterns | parallel -j 12 "
    fuzznuc -sequence ${outdir}/${sampleid}_rdna.fa -pattern {} -pmismatch 0 -complement Y -outfile ${outdir}/${sampleid}_${gene}_{}.txt -rformat2 excel
    "

    cat ${outdir}/${sampleid}_${gene}_* | grep -v "SeqName" > ${outdir}/${sampleid}.${gene}.kmermatches.txt
    rm ${outdir}/${sampleid}_${gene}_*

    echo "Kmers found for $sampleid"
    /g/data/te53/zc4688/honours/scripts/processkmers.R ${sampleid} ${gene}

    echo "Counting kmers for $sampleid"
    seqkit subseq --bed ${outdir}/${sampleid}.${gene}.coord.bed ${outdir}/${sampleid}_rdna.fa |jellyfish count -t 10 -C -m 19 -s 3G -o ${outdir}/${sampleid}.${gene}.jf /dev/fd/0

    jellyfish dump -o ${outdir}/${sampleid}.${gene}.fa ${outdir}/${sampleid}.${gene}.jf

    if [ ! -e "${outdir}/${sampleid}_rdna.fa.gz" ]; then
        echo "Compressing FASTA for $sampleid"
        pigz -p 12 ${outdir}/${sampleid}_rdna.fa
    else
        rm ${outdir}/${sampleid}_rdna.fa
    fi

    echo "Task finished for $sampleid"
}

export -f find_kmers  # Export the function for parallel
export outdir
export gene

parallel -j 16 find_kmers ::: "${sampleids[@]}"

