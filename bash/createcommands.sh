#!/bin/bash

#Usage - generates commands file to be used by run_ribocop.sh (required for parallelisation)

#python3 /g/data/te53/zc4688/honours/ribocop/testing.py -o /g/data/te53/zc4688/honours/ribocop/all -s GCA_031021105.1_rEryReg1.hap1_genomic -l /g/data/xl04/hrp561/rdnalib/deuterostomia.18s28s.rDNA.fasta -t 12 -i /g/data/xl04/genomeprojects/referencedata/genomes/GCA_031021105.1.Erythrolamprus_reginae.txid121349/
base_directory="/g/data/xl04/genomeprojects/referencedata/tmp/chordata"
output_directory="/g/data/te53/zc4688/honours/scripts/command_files/commands_hmm.txt"

find "$base_directory" -type f \( -name "*.fna" -o -name "*.fna.gz" \) | while read file; do
    # Extract sampleid, txid, and species from the directory structure
    dir_name=$(dirname "$file")  # Get the directory name
    base_name=$(basename "$dir_name")  # Get the subdirectory name: sampleid.number.species.txid

    # Extract parts from the base name (using underscore as delimiter)
    sampleid=$(echo "$base_name" | cut -d '.' -f1)
    txid=$(echo "$base_name" | cut -d '.' -f4)
    species=$(echo "$base_name" | cut -d '.' -f3)

    command="bash /g/data/te53/zc4688/honours/scripts/runcmd.sh 'python3 /g/data/te53/zc4688/honours/scripts/ribocop/ribocop2.py -o /g/data/te53/zc4688/honours/analyses/chordates/new -t 12 -i $base_directory -s  $sampleid -l /g/data/te53/zc4688/honours/data/euk.hmm -x $txid -p $species' ${sampleid}_log /g/data/te53/zc4688/honours/analyses/chordates/new/${sampleid} 1"

    echo "$command " >> $output_directory 
done

echo "Commands have been written to file"

