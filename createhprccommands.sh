#!/bin/bash
base_directory="/g/data/te53/t2t2024/referenceresource/hprc/"
output_directory="/g/data/te53/zc4688/honours/scripts/command_files/commands_hprc.txt"

find "$base_directory" -type f -name "*.fna" | while read file; do
    # Extract sampleid, txid, and species from the directory structure

    # Extract parts from the base name (using underscore as delimiter)
    sampleid=$(basename $file | cut -d '.' -f1)

    command="bash /g/data/te53/zc4688/honours/scripts/runcmd.sh 'python3 /g/data/te53/zc4688/honours/scripts/ribocop/ribocop.py -o /g/data/te53/zc4688/honours/analyses/cn/hprc -t 12 -i $base_directory -s  $sampleid -if $file -l /g/data/te53/zc4688/honours/data/euk.hmm' ${sampleid}_log /g/data/te53/zc4688/honours/analyses/chordates/new/${sampleid} 1"

    echo "$command " >> $output_directory 
done

echo "Commands have been written to file"

