import subprocess
import pandas as pandas

#could add read subsampling first?
#seqtk sample - random sample 10% of reads (maybe multiple times for accuracy?)
#minimap - rdna reads
def count_kmers(reads, output_dir, k, if):
    command=f"jellyfish count -C -m {k} -s 100M -o {} --if {filtered ntsm subset} {reads}" 

def main():
    """Main function."""   
    parser = argparse.ArgumentParser(description="Calculate rDNA copy number and divergence from reads")
    parser.add_argument('-k', '--kmer_size', required=False, default=19, help='Kmer size')
    parser.add_argument('r', '--reads', required=True, help='Input reads')
    parser.add_argument('-o', '--output_dir', required=False, default='.', help="Output directory for resaults (default: current directory)")
    parser.add_argument('if', required=False)

    args = parser.parse_args()