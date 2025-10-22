#!/bin/bash
base_directory="/g/data/te53/zc4688/honours/analyses/chordates/new"
outputfile="/g/data/te53/zc4688/honours/analyses/chordates/new/structures.txt"
> "$outputfile"
find "$base_directory" -type f \( -name "*.refmorph.fasta" \) | while read file; do
    base_name=$(basename "$file")  # Get the subdirectory name: sampleid.number.species.txid

    # Extract parts from the base name (using underscore as delimiter)
    sampleid=$(echo "$base_name" | cut -d '.' -f1)
    refsequence=$(grep ">" $file  | tr -d ">")
    if [[ $refsequence =~ (.*):([0-9]+)-([0-9]+):([-+]) ]]; then
        refseqid="${BASH_REMATCH[1]}"
        refstart="${BASH_REMATCH[2]}"
        refend="${BASH_REMATCH[3]}"
        refstrand="${BASH_REMATCH[4]}"
    else
        echo "Could not parse refsequence in $file" >&2
        continue
    fi

    primaryfile="${sampleid}.primary.txt"

    if [[ ! -f "$primaryfile" ]]; then
        echo "WARNING: primaryfile $primaryfile not found, skipping." >&2
        continue
    fi

    echo "Filtering $primaryfile for overlaps..."
    
    awk -v seqid="$refseqid" -v start="$refstart" -v end="$refend" -v sid="$sampleid" '
        $1 == seqid && $9 >= start && $10 <= end {
            print sid " " $0
        } 
    ' "$primaryfile" >> "$outputfile"
    
done
