#!/bin/bash
#PBS -q normal
#PBS -l walltime=10:00:00
#PBS -l mem=6GB
#PBS -l ncpus=12
#PBS -j oe
#PBS -N "msa"
#PBS -o /g/data/te53/zc4688/honours/logs/
#PBS -l storage=gdata/te53+gdata/if89

module load seqkit/2.9.0 trimal
export PATH=$PATH:/g/data/te53/zc4688/software


#usage: qsub -v rrna=twoeight -P te53 clustermsa.sh

rrna="${rrna}"
INDIR="/g/data/te53/zc4688/honours/analyses/chordates/motif"

FASTA="$INDIR/${rrna}.fa"
CLUSTERS="$INDIR/out.${rrna}.self.blastn.mci.I20"

OUTDIR="$INDIR/${rrna}_clusters"

mkdir -p "$OUTDIR"

CLUSTER_ALIGNS=()

i=0
while IFS= read -r line; do
    cluster=($line)  # split line into array
    if [ ${#cluster[@]} -ge 10 ]; then
        ((i++))
        cluster_file="$OUTDIR/cluster_$i.ids"
        align_file="$OUTDIR/cluster_$i.aln.fasta"
        echo "${cluster[@]}" > "$cluster_file"

        # Extract sequences from the master FASTA
        seqkit grep -f <(tr ' ' '\n' < "$cluster_file") "$FASTA" > "$OUTDIR/cluster_$i.fasta"

        # Align with MUSCLE
        muscle3 -in "$OUTDIR/cluster_$i.fasta" -out "$align_file"
        CLUSTER_ALIGNS+=("$align_file")
    fi
done < "$CLUSTERS"

FINAL_ALIGN="${CLUSTER_ALIGNS[0]}"


for aln in "${CLUSTER_ALIGNS[@]:1}"; do
    NEXT_ALIGN="$OUTDIR/tmp_$$.aln.fasta"
    muscle3 -profile -in1 "$FINAL_ALIGN" -in2 "$aln" -out "$NEXT_ALIGN"
    FINAL_ALIGN="$NEXT_ALIGN"
done

mv "$FINAL_ALIGN" "$OUTDIR/final_cluster_alignment.fasta"

cat "$OUTDIR"/cluster_*.ids | tr ' ' '\n' | sort | uniq > "$OUTDIR/used.ids"

# Get unaligned sequences
seqkit grep -v -f "$OUTDIR/used.ids" "$FASTA" > "$OUTDIR/remaining_sequences.fasta"

# Add them to the final alignment

muscle3 -in1 "$OUTDIR/final_cluster_alignment.fasta" -in2 "$OUTDIR/remaining_sequences.fasta" -out "$OUTDIR/final_alignment_with_all.fasta"

# Loop over each remaining sequence
seqkit split -s 1 "$OUTDIR/remaining_sequences.fasta" -O "$OUTDIR/split_remaining"

FINAL_ALIGN="$OUTDIR/final_cluster_alignment.fasta"

for seqfile in "$OUTDIR"/split_remaining/*.fasta; do
    
    tmp_align="$OUTDIR/tmp_$$.fasta"
    muscle3 -profile -in1 "$FINAL_ALIGN" -in2 "$seqfile" -out "$tmp_align"

    # Update reference alignment
    FINAL_ALIGN="$tmp_align"
done

# Move to final output
mv "$FINAL_ALIGN" "$OUTDIR/final_alignment_with_all.fasta"

rm "$OUTDIR/remaining_sequences.fasta"
rm "$OUTDIR/used.ids"
rm -r "$OUTDIR/split_remaining"

trimal -in "$OUTDIR/final_alignment_with_all.fasta" -out "$OUTDIR/cleaned.alignment.txt" -gt 0.2