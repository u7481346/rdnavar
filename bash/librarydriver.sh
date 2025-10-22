#!/bin/bash

#indir="$1"
indir="/g/data/te53/referencedata/openaccess/giab/NA12878/NIST_NA12878_HG001_HiSeq_300x/"
outdir="/g/data/te53/zc4688/honours/analyses/cn/illumina_libs/library/HG001"

#outdir="$2"
tmpdir="/g/data/te53/zc4688/honours/tmp"
libraries=($(find "$indir" -mindepth 3 -maxdepth 3 -type d -name 'Sample_U*' -exec basename {} \; | sort -u))


#per flow cell:
#for r1 in "$indir"/*_R1.fastq.gz; do
    #r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    #libname=$(basename "$r1"   | sed -E 's/^(HG[0-9]+).*_([0-9]{4}).*/\1_\2/')

    #echo $r1 $r2 $libname
    #if [ -f "${outdir}/${libname}.done" ]; then
        #echo "Task already completed for ${libname} – skipping"
        #continue
    #fi   

    #qsub -v libname="$libname",infastq="$r1 $r2" -P te53 /g/data/te53/zc4688/honours/scripts/librarycompare.sh 

#done

#per library:
for library in "${libraries[@]}"; do
    libname=$(echo "$library" | cut -f2 -d '_')
    r1files=($(find "$indir" -type f -name "$libname*R1*fastq.gz" | sort))
    r2files=($(find "$indir" -type f -name "$libname*R2*fastq.gz" | sort))


    if [ -f "${outdir}/${libname}.done" ]; then
        echo "Task already completed for ${libname} – skipping"
        continue
    fi   

    cat "${r1files[@]}" > "${tmpdir}/${libname}.R1.fastq.gz"
    cat "${r2files[@]}" > "${tmpdir}/${libname}.R2.fastq.gz"

    infastq="${tmpdir}/${libname}.R1.fastq.gz ${tmpdir}/${libname}.R2.fastq.gz"

    qsub -v libname="$libname",infastq="$infastq" -P te53 /g/data/te53/zc4688/honours/scripts/librarycompare.sh 

done


#for file in *depth.txt; do libname=$(basename $file | sed -E 's/^(HG[0-9]+).*_([0-9]{4}).*/\1_\2/'); echo $libname; genome=$(grep "50kb" $file| cut -f2); rdna=$(grep "rdna" $file| cut -f2); ratio=$(echo "scale=4; $rdna / $genome" | bc); echo $ratio; don