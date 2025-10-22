import pandas as pd
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
import sys
import numpy as np
import gzip
import tarfile
import json
import glob
import re
import shutil

def morph_identification(filtered_data_df, sampleid):
    #Morph identification from filtered alignments (requires high quality 18S-28S-18S alignments)
    #Restrict to 18 and 28S for morph identification
    filtered_data_df = filtered_data_df[filtered_data_df["Type"].str.contains("18S|28S", regex=True)]


    morphs = set()
    for query_name, group in filtered_data_df.groupby("seqid"):
        group = group.sort_values(by=["Envstart"]).reset_index(drop=True)
        valid_group = set()
        i = 0
        while i < len(group):
            # start with a 18S in the plus strand, and 100bp available before this 18S on the contig
            if '18S' in group.at[i, 'Type'] and group.at[i, 'Strand'] == "+" and group.at[i, 'Envstart'] - 100 > 0:
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '28S' in group.at[i + 1, 'Type'] and group.at[i + 1, 'Strand'] == "+":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '18S' in group.at[i + 2, 'Type'] and group.at[i + 2, 'Strand'] == "+":
                        morph = group.at[i, 'seqid'], group.at[i, 'Envstart'], group.at[i + 2, 'Envend'], '+'
                        morphs.add(morph)
            # start with a 18S in the minus strand
            elif '18S' in group.at[i, 'Type'] and group.at[i, 'Strand'] == "-":
                #check if the next one is a valid 28S
                if i + 1 < len(group) and '28S' in group.at[i + 1, 'Type'] and group.at[i, 'Strand'] == "-":
                    # check if the next one is a valid 18S
                    if i + 2 < len(group) and '18S' in group.at[i + 2, 'Type'] and group.at[i + 2, 'Strand'] == "-": #had to reove bc not in filtered_gff and group.at[i + 2, 'Envend'] + 100 <= group.at[i, 'Target length']: #Note previously we also filtered to ensure the length + 100 was < total length of query but this is not part of barrnap output
                        morph = group.at[i, 'seqid'], group.at[i, 'Envstart'], group.at[i + 2, 'Envend'], '-'
                        morphs.add(morph)             

            i += 1
    # Convert morphs set to DataFrame
    if len(morphs)==0:
        print(f"No morphs pass filtering for {sampleid}.Exiting")
        exit()
    else:
        morphs_df = pd.DataFrame(list(morphs), columns=['Query sequence name', 'Start', 'End', 'Strand'])
        morphs_df["morph_number"] = [f"morph{i+1}" for i in range(len(morphs_df))]
    return morphs_df
    # return morphs

def extract_sequences(fasta_file, morphs_df, output_dir, sampleid):
    #Extracts and saves morph sequences
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "rt") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    else:
        with open(fasta_file, "rt") as handle:
            sequences = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    extracted_sequences = []
    valid_rows = []

    for _, row in morphs_df.iterrows():
        morph_id = row["morph_number"]
        seq_id = row['Query sequence name']
        start = int(row['Start'])
        end = int(row['End'])
        strand = row['Strand']

        if seq_id in sequences:
            sequence = sequences[seq_id].seq[start:end]
            if strand == '-':
                sequence = sequence.reverse_complement()
            sequence_str = str(sequence)
            if 'N' not in sequence_str:
                extracted_sequences.append((f"{morph_id}.{seq_id}:{start}-{end}:{strand}", sequence_str))
                valid_rows.append(row)

    if len(valid_rows)==0:
        print(f"No non-ambiguous morphs for {sampleid}.Exiting")
        exit()
    else:
        with open(os.path.join(output_dir, sampleid + ".rDNA.morphs.fasta"), 'w') as output_handle:
            for seq_id, sequence in extracted_sequences:
                output_handle.write(f">{seq_id}\n{sequence}\n")


def rna_builder(sampleid, inputfasta, output_dir, input_dir):
    filtered_gff = f"{input_dir}/{sampleid}.filtered.gff"

    column_names = [
        "seqid", "source", "Type", "Envstart", "Envend", "score", "Strand", "frame", "attributes"
    ]
    filtered_data_df = pd.read_csv(filtered_gff, sep="\t", header=None,names = column_names, skiprows=1)
    filtered_data_df.loc[:, "Envstart"] = filtered_data_df["Envstart"] - 1
    morphs_df = morph_identification(filtered_data_df, sampleid)
    extract_sequences(inputfasta, morphs_df, output_dir, sampleid)

    return filtered_data_df, morphs_df

def get_median_morph(sampleid, inputfasta, filtered_data_df, output_dir, morphs_df):
    #Reads morph sequences and finds median length unit which is designated as the representative morph for the genome. Uses this to find unit length and reference unit structure. 
    morph_fasta = f"{output_dir}/{sampleid}.rDNA.morphs.fasta"

    fai_file = f"{morph_fasta}.fai"

    subprocess.run(f"samtools faidx {morph_fasta}", shell=True, check=True)

    morphs = pd.read_csv(fai_file, sep='\t', header=None, names=['sequence_name', 'length', 'offset', 'linebases', 'linewidth', 'qualoffset'])

    # Sort by length
    morphs_sorted = morphs.sort_values(by='length')

    # Find the median length sequence
    median_index = len(morphs_sorted) // 2
    refrdnamorph = morphs_sorted.iloc[median_index]['sequence_name']
    refmorphlength = morphs_sorted.iloc[median_index]['length']
    min_length = morphs_sorted.iloc[0]['length']
    max_length = morphs_sorted.iloc[-1]['length']


    all_overlaps = []


    for _, morph_row in morphs_df.iterrows():
        morph_id = morph_row["morph_number"]
        seqid = morph_row["Query sequence name"]
        start = int(morph_row["Start"])
        end = int(morph_row["End"])
        strand = morph_row["Strand"]

        overlapping = filtered_data_df[
            (filtered_data_df["seqid"] == seqid) &
            (filtered_data_df["Envstart"] >= start) &
            (filtered_data_df["Envend"] <= end)
        ].copy()

        if overlapping.empty:
            continue

    # Adjust rRNA coordinates relative to morph
        if strand == "+":
            overlapping["Start"] = overlapping["Envstart"] - start
            overlapping["End"] = overlapping["Envend"] - start
        else:
            overlapping["Start"] = end - overlapping["Envend"]
            overlapping["End"] = end - overlapping["Envstart"]

        overlapping["morph_number"] = morph_id
        
        overlapping = overlapping[[
            "morph_number",
            "Start", 
            "End", 
            "Type"
        ]]

        all_overlaps.append(overlapping)
        


    if all_overlaps:
        combined = pd.concat(all_overlaps, ignore_index=True)
        output_path = os.path.join(output_dir, f"{sampleid}.structure.bed")
        combined.to_csv(output_path, sep="\t", index=False)
        print(f"Written morph structure file with {len(combined)} rRNAs from {len(morphs_df)} morphs.")
    else:
        print(f"No overlapping rRNAs found for any morph in sample {sampleid}")

    
def main():
    """Main function."""   
    parser = argparse.ArgumentParser(description="Identify and characterise rDNA units from FASTA input.")
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for results')
    parser.add_argument('-s', '--sampleid', required=True, help='Sample ID for output files')
    parser.add_argument('-i', '--input_dir', required=False, default=".", help="Input directory (default: current directory)")
    parser.add_argument('--inputfasta', required=True, help="Input FASTA file")


    args = parser.parse_args()
    

    filtered_data_df, morphs_df = rna_builder(args.sampleid, args.inputfasta, args.output_dir, args.input_dir)

    get_median_morph(args.sampleid, args.inputfasta, filtered_data_df, args.output_dir, morphs_df)



    print("Task completed. Exiting")
    
    

if __name__ == "__main__":
    main()     
