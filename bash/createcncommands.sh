#!/bin/bash

#Usage - generates commands file to be used by run_ribocop.sh (required for parallelisation)

#python3 /g/data/te53/zc4688/honours/ribocop/testing.py -o /g/data/te53/zc4688/honours/ribocop/all -s GCA_031021105.1_rEryReg1.hap1_genomic -l /g/data/xl04/hrp561/rdnalib/deuterostomia.18s28s.rDNA.fasta -t 12 -i /g/data/xl04/genomeprojects/referencedata/genomes/GCA_031021105.1.Erythrolamprus_reginae.txid121349/
base_directory="/g/data/xy86/ont1kgp"
output_directory="/g/data/te53/zc4688/honours/scripts/command_files/cn_commands.txt"

find "$base_directory" -type f \( -name "*.bam" -o -name "*.fna.gz" \) | while read file; do
    # Extract sampleid, txid, and species from the directory structure
    dir_name=$(dirname "$file")  # Get the directory name
    base_name=$(basename "$dir_name")  # Get the subdirectory name: sampleid.number.species.txid
    
    # Extract parts from the base name (using underscore as delimiter)
    sampleid=$(basename "$file" | cut -d '.' -f1)

    command="bash /g/data/te53/zc4688/honours/scripts/runcmd.sh 'python3 /g/data/te53/zc4688/honours/scripts/sequence_var.py -o /g/data/te53/zc4688/honours/analyses/cn/result/cn/ -a /g/data/te53/zc4688/honours/analyses/cn/kmers.fa -k 19 -s ${sampleid} -t12 -d /g/data/te53/zc4688/honours/tmp/ --inbam $file' ${sampleid}_log /g/data/te53/zc4688/honours/analyses/cn/result/${sampleid} 1"

    echo "$command " >> $output_directory 
done

echo "Commands have been written to file"

