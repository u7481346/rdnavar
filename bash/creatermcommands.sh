#!/bin/bash

base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
output_directory="/g/data/te53/zc4688/honours/analyses/chordates/repeatmasker"
cmd_file="/g/data/te53/zc4688/honours/scripts/command_files/repeatmasker_commands.txt"

mkdir -p "$(dirname "$cmd_file")"
> "$cmd_file"

find "$base_directory" -type f -name "*.rDNA.refmorph.fasta" | while read file; do
    # Extract sampleid, txid, and species from the directory structure
    sampleid=$(basename "$file" | cut -d '.' -f1)  # Get the directory name

    echo "RepeatMasker -species chordates -s -xm -gff -dir $output_directory $file" >> "$cmd_file"
done