#!/bin/bash

base_directory="/g/data/xl04/genomeprojects/referencedata/tmp/chordata"
infile="/g/data/te53/zc4688/honours/metadata/valid_ids.txt"
output_directory="/g/data/te53/zc4688/honours/metadata"
output_file="$output_directory/assembly_metadata.tsv"

echo -e "accession\torganism_name\tassembly_level\tassembly_method\tgenome_size\tassembly_type\tbioproject\ttsequencing_tech\tsubmitter\tassembly_status\tassembly_stats" > "$output_file"


while read -r id; do
    # Find the matching file (safely)
    file=$(find "$base_directory" -type f -name "${id}*_data_report.jsonl" | head -n 1)

    if [[ -f "$file" ]]; then
        sample_id=$(jq -r '.accession' "$file")
        organism_name=$(jq -r '.organism.organismName' "$file")
        assembly_level=$(jq -r '.assemblyInfo.assemblyLevel' "$file")
        assembly_method=$(jq -r '.assemblyInfo.assemblyMethod' "$file")
        genome_size=$(jq -r '.assemblyStats.totalSequenceLength' "$file")
        assembly_type=$(jq -r '.assemblyInfo.assemblyType' "$file")
        source=$(jq -r '.assemblyInfo.assembly.bioproject' "$file")
        sequencing_tech=$(jq -r '.assemblyInfo.sequencingTech' "$file")
        submitter=$(jq -r '.assemblyInfo.submitter' "$file")
        assembly_status=$(jq -r '.assemblyInfo.assemblyStatus' "$file")
        assembly_stats=$(jq -c '.assemblyStats' "$file")

        echo -e "${sample_id}\t${organism_name}\t${assembly_level}\t${assembly_method}\t${genome_size}\t${assembly_type}\t${source}\t${sequencing_tech}\t${submitter}\t${assembly_status}\t${assembly_stats}" >> "$output_file"

        echo "$sample_id successfully added to file."
    else
        echo "Warning: No file found for ID $id"
    fi
done < "$infile"
