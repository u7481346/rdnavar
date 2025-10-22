import pandas as pd
import numpy as np
import subprocess
import sys
import os
import argparse

def check_output(output_dir, sampleid):
    """Check if output has already been generated for a given sample ID."""
    checkpoint_file = os.path.join(output_dir, f"{sampleid}_kmers")
    if os.path.exists(checkpoint_file + ".done"):
        print(f"Task already completed for {sampleid}. Exiting")
        sys.exit(0)
    open(checkpoint_file + ".running", 'w').close()

def bam_to_fasta(outputdir, sampleid, inbam, NCPUS):
    tmp_fasta = os.path.join(outputdir, f"{sampleid}.tmp.fasta")
    if os.path.exists(tmp_fasta):
        print(f"FASTA already present for {sampleid}")
        return
    
    print(f"Converting BAM to FASTA for {sampleid}")
    subprocess.run(f"samtools fastq --threads {NCPUS} -0 {tmp_fasta} {inbam}", check=True, shell=True)

def count_cn_kmers(sampleid, outputdir, reads, NCPUS, kmerfile, k):
    cn_kmer_file=os.path.join(outputdir, f"{sampleid}.cn.kmers.fa")
    tmp_fasta = os.path.join(outputdir, f"{sampleid}.tmp.fasta")
    if os.path.exists(cn_kmer_file):
        print(f"CN kmers alreads counted for {sampleid}")
        return

    #print(f"Counting cn kmers for {sampleid}")
    #tried with bbduk - slightly slower
    #subprocess.run(f"samtools fasta --threads {NCPUS} {inbam} | jellyfish count -t {NCPUS} -C -m 19 -s 5G --if {outputdir}/kmers.fa -o {outputdir}/{sampleid}.jf /dev/fd/0", check=True, shell=True)
    reads_str = " ".join(reads)

    try:
        print(f"Counting cn kmers for {sampleid}")

        if reads_str.endswith('.gz'):
            subprocess.run(
                f"zcat {reads_str} | jellyfish count -t {NCPUS} -C -m {k} -s 5G --if {kmerfile} -o {outputdir}/{sampleid}.jf /dev/fd/0",
                shell=True,
                check=True
            )
        else: 
            subprocess.run(
                f"jellyfish count -t {NCPUS} -C -m {k} -s 5G --if {kmerfile} -o {outputdir}/{sampleid}.jf {reads_str}",
                shell=True,
                check=True
            )
            
        print(f"Dumping cn kmers for {sampleid}")
        subprocess.run(
            f"jellyfish dump {outputdir}/{sampleid}.jf > {cn_kmer_file}",
            shell=True,
            check=True
        )
    finally:
        if os.path.exists(tmp_fasta):
            os.remove(tmp_fasta)
            print(f"Removed temporary file {tmp_fasta}")



def main():
    """Main function."""   
    parser = argparse.ArgumentParser(description="Identify and characterise rDNA units from FASTA input.")
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for results')
    parser.add_argument('-s', '--sampleid', required=True, help='Sample ID for output files')
    parser.add_argument('-t', '--ncpus', required=True, help='Threads')
    parser.add_argument('-d', '--tmp_dir', required=False, help="Temporary directory for storing modified FASTA files (default:output directory)")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-b', '--inbam', help="Input BAM file")
    group.add_argument('-f', '--infile', nargs='+', help="Input FASTQ or FASTA files (can supply multiple)")
    parser.add_argument('-a', '--kmer_file', required=True, help="Input kmer file to count. Should include both RDNA kmers and comparison kmers")
    parser.add_argument('-k', required=True, help="Kmer size")

    args = parser.parse_args()

    check_output(args.output_dir, args.sampleid)

    tmpdir = args.tmp_dir if args.tmp_dir else args.output_dir

    #Create output directory if it does not exist.
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        print(f"Directory '{args.output_dir}' created successfully.")
    except OSError as e:
        print(f"Error creating directory '{args.output_dir}': {e}")
        sys.exit(1)


    if args.inbam:
        print(f"Using BAM: {args.inbam}")
        bam_to_fasta(args.output_dir, args.sampleid, args.inbam, args.ncpus)
        tmp_fasta=os.path.join(args.output_dir, f"{args.sampleid}.tmp.fasta")
        count_cn_kmers(args.sampleid, args.output_dir, tmp_fasta, args.ncpus, args.kmer_file, args.k)
    
    if args.infile:
        print(f"Using input: {args.infile}")
        count_cn_kmers(args.sampleid, args.output_dir, args.infile, args.ncpus, args.kmer_file, args.k)




    checkpoint_file = os.path.join(args.output_dir, f"{args.sampleid}_kmers")
    
    open(checkpoint_file + ".done", 'w').close()

    os.remove(f"{checkpoint_file}.running")


    print(f"Task completed for {args.sampleid}. Exiting")
    
    

if __name__ == "__main__":
    main()     

