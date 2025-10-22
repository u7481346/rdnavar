
module load RepeatMasker/4.1.5 bedtools/2.31.0 TRF/4.09.1 seqkit/2.9.0 jellyfish bwa-mem2/2.2.1

assembly="/g/data/xl04/genomeprojects/TILRUE/hifiasm/TILRUE.hap1.p_ctg.fasta"
refmorph="/g/data/xl04/genomeprojects/TILRUE/analysis/asmqc/hifiasm/ribocop/h1/TILTRUE.h1.rDNA.refmorph.fasta"
outdir="/g/data/te53/zc4688/honours/analyses/cn/TILTRUE/"

cd $outdir
RepeatMasker -species Lepidosauria -s -gff -dir $outdir $refmorph 

out=$(basename "$refmorph").out.gff
grep -v "rRNA" ${out} | tail -n +3| awk '{print $1"\t"($4-1)"\t"$5}'  > masking.bed
    
bedtools maskfasta -fi $refmorph -bed masking.bed -fo masked.fasta -mc N

trf masked.fasta 2 7 7 80 10 10 2000 -h -m


#Create kmer set
jellyfish count -m 21 -C -s 100M -o rdnakmers.jf masked.fasta.2.7.7.80.10.10.2000.mask
jellyfish dump rdnakmers.jf > rdna_dumps.fa

#filter out those occuring more than once in the rdna morph
seqkit grep -nr -p "^1$" rdna_dumps.fa -o rdna_dumps_uniq.fa



mv rdna_dumps_uniq.fa rdna_dumps.fa

seqkit fx2tab -g rdna_dumps.fa | awk '$3<65 && $3>25' | seqkit tab2fx > gcfiltered_rdna_dumps.fa
awk '/^>/ {print ">kmer_" ++i} !/^>/' gcfiltered_rdna_dumps.fa > named.fa
bwa-mem2 mem -a -t24 -T18 -h2 -k 10 $assembly named.fa > testbwa.sam


samtools view -b testbwa.sam > testbwa.bam

cut -f1,3,4 /g/data/xl04/genomeprojects/TILRUE/analysis/asmqc/hifiasm/ribocop/h1/TILTRUE.h1.rdnaregions.bed > rdnaregions.bed 


bedtools intersect -v -a testbwa.bam -b rdnaregions.bed | samtools view | cut -f1 | sort | uniq > nonunique_rdnakmers.txt
samtools view -h testbwa.bam | grep -v -F -f nonunique_rdnakmers.txt > testbwa.filtered.sam

cut -f1 testbwa.filtered.sam | sort | uniq > unique_rdnakmers.txt
seqkit grep -f unique_rdnakmers.txt named.fa  > rdna_final.fa

#classify rdna kmers and check they also only occur once in the reference morph (unmasked)
bwa-mem2 index /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/ribocop/POGVIT.rDNA.refmorph.fasta 
bwa-mem2 mem -T19 /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/ribocop/POGVIT.rDNA.refmorph.fasta rdna_final.fa -o rdna_categories.sam

#unsure if this is exactly correct.. have to check re 0based/1based. also edited to include other categories manually. i think end may be off by 1
tail -n3 /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/ribocop/POGVIT.refmorph.structure.tsv | cut -f3,9,10 | awk '{print "scaffold_40:6248-15708:-" "\t" 15708-$3 "\t" 15707-$2 "\t"  $1}' > rdna.refmorph.structure.bed
grep -v "XA" rdna_categories.sam | samtools view -b | bedtools bamtobed > rdna_categories.bed
bedtools intersect -wao -a rdna_categories.bed -b rdna.refmorph.structure.bed > test.bed

#edited the kmers overlapping the gap between an rrna and its/igs manually - if a kmer bridges both, dont count it as an rrna kmer

cut -f4 test.bed > rdna_final.txt #252 kmers
seqkit grep -f rdna_final.txt rdna_final.fa > final_rdna_kmers.fa
cut -f4,10,2 test.bed >rdna_categories.txt




jellyfish count -t 24 -C -m 21 -s 3G -o pogvit.jf /g/data/xl04/genomeprojects/Pogona_vitticeps/fasta/POGVIT.v2.1.fasta
jellyfish dump pogvit.jf > pogvit.fa
seqkit grep -nr -p "^1$" pogvit.fa -o pogvit_uniq.fa

#setdiff - only unique
seqkit sample -p 0.1 pogvit_uniq.fa| seqkit head -n 100000 > final_genome.fa


seqkit replace -p ".*" -r "genome" final_genome.fa -o renamed_final.fa

cat final_rdna_kmers.fa renamed_final.fa > kmers.fa #final kmer file


###############calculate CN based on alignments

#problem - minimap only outputs 5 secondary alignments by default. this means that not all alignments per read are output. 
#cant just subset to reference, because many reads wont align to that morph.
#if we take only the primary alignments, and calculate no. of aligned bases?
samtools view  -h -L /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/ribocop/POGVIT.rDNA.morphs.tsv -F 256,272 -o ont.subset.sam /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/depth/rearranged/POGVIT.v2.1.merged.ont.bam

cut -f4,10 ont.subset.sam | awk '{s += length($2)} END {print s}'
awk '{s += $3 - $2} END {print s}' /g/data/xl04/genomeprojects/Pogona_vitticeps/analysis/ribocop/POGVIT.rDNA.morphs.tsv

awk '!/^@/{                           
    cigar=$6;
    len=0;
    while(match(cigar, /[0-9]+[MIX=]/)) {  # exclude S
        val=substr(cigar, RSTART, RLENGTH-1);
        len += val;
        cigar=substr(cigar, RSTART+RLENGTH);
    }
    print $1, len
}' ont.subset.sam > tmp.txt