
#module load RepeatMasker/4.1.5 bedtools/2.31.0 TRF/4.09.1 seqkit/2.9.0 jellyfish bwa-mem2/2.2.1

RepeatMasker -species human -s -gff -dir /g/data/te53/zc4688/honours/analyses/cn/ /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta

grep -v "rRNA" chm13.rDNA.refmorph.fasta.out.gff | tail -n +3| awk '{print $1"\t"($4-1)"\t"$5}'  > masking.bed
    
bedtools maskfasta -fi /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta -bed masking.bed -fo masked.fasta -mc N

trf masked.fasta 2 7 7 80 10 10 2000 -h -m


#Create kmer set
jellyfish count -m 19 -C -s 100M -o rdnakmers.jf masked.fasta.2.7.7.80.10.10.2000.mask
jellyfish dump rdnakmers.jf > rdna_dumps.fa

#filter out those occuring more than once in the rdna morph
seqkit grep -nr -p "^1$" rdna_dumps.fa -o rdna_dumps_uniq.fa



mv rdna_dumps_uniq.fa rdna_dumps.fa

#seqkit fx2tab -g rdna_dumps.fa | awk '$3<65 && $3>25' | seqkit tab2fx > gcfiltered_rdna_dumps.fa
awk '/^>/ {print ">kmer_" ++i} !/^>/' rdna_dumps.fa > named.fa
bwa-mem2 mem -a -t24 -T18 -h2 -k 10 /g/data/te53/humanreference/chm13-t2t/bwa-mem2/chm13-t2t.ebv.phix.chrQ.xo.fa named.fa > testbwa.sam


samtools view -b testbwa.sam > testbwa.bam
bedtools intersect -v -a testbwa.bam -b /g/data/te53/t2t2024/referenceresource/rDNA_CHM13_NCIG_20250604.bed | samtools view | cut -f1 | sort | uniq > nonunique_rdnakmers.txt
samtools view -h testbwa.bam | grep -v -w -F -f nonunique_rdnakmers.txt > testbwa.filtered.sam

cut -f1 testbwa.filtered.sam | sort | uniq > unique_rdnakmers.txt
seqkit grep -f unique_rdnakmers.txt named.fa  > rdna_final.fa

#classify rdna kmers and check they also only occur once in the reference morph (unmasked)
bwa-mem2 index /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta 
bwa-mem2 mem -T19 /g/data/te53/zc4688/honours/analyses/chordates/new/chm13.rDNA.refmorph.fasta rdna_final.fa -o rdna_categories.sam

#unsure if this is exactly correct.. have to check re 0based/1based. also edited to include other categories manually
tail -n3 ../chordates/new/chm13.refmorph.structure.tsv | cut -f3,9,10 | awk '{print "chr21:3697875-3742662:+" "\t" $2-3697876 "\t" $3-3697875 "\t" $1}' > rdna.refmorph.structure.bed
grep -v "XA" rdna_categories.sam | samtools view -b | bedtools bamtobed > rdna_categories.bed
bedtools intersect -wao -a rdna_categories.bed -b rdna.refmorph.structure.bed > test.bed

#edited the kmers overlapping the gap between an rrna and its/igs manually - if a kmer bridges both, dont count it as an rrna kmer
awk '$11 == 19 || $10 !~ /rRNA/' test.bed > final.bed

cut -f4 final.bed > rdna_final.txt #10565
seqkit grep -f rdna_final.txt rdna_final.fa > final_rdna_kmers.fa
cut -f4,10,2 final.bed >rdna_categories.txt



#cat final_rdna_kmers.fa ntsm/data/human_sites_n10.fa > input.fa
#seqkit grep -nr -p "ref|kmer" input.fa -o tmp.fa

#subset to first 19bp - only 1 kmer per ntsm variant allele
#seqkit subseq -r 1:19 tmp.fa > tmp2.fa
#also generating chm13 single occurence kmers

zcat /g/data/te53/variantcall/referenceresource/stratification/CHM13\@all/Union/CHM13_notinalldifficultregions.bed.gz  | awk 'BEGIN {OFS="\t"} {$1="chm13#0#"$1; print}' | gzip > /g/data/te53/zc4688/honours/analyses/cn/difficult_regions.bed.gz
bedtools getfasta -fi /g/data/te53/t2t2024/referenceresource/chm13.filtered.fasta -bed /g/data/te53/zc4688/honours/analyses/cn/difficult_regions.bed.gz -fo /g/data/te53/zc4688/honours/analyses/cn/clean_regions.fa
jellyfish count -t 24 -C -m 19 -s 3G -o chm13.jf clean_regions.fa
jellyfish dump chm13.jf > chm13_kmers.fa
seqkit grep -nr -p "^1$" chm13_kmers.fa -o chm13_kmers_uniq.fa

#setdiff - only unique
seqkit sample -p 0.1 chm13_kmers_uniq.fa | seqkit head -n 100000 > final_genome.fa

#seqkit common -s -i tmp2.fa final_genome.fa | grep -v ">" > test.fa

#seqkit grep -s -f test.fa -v final_genome.fa -o final.fa

seqkit replace -p ".*" -r "genome" final_genome.fa -o renamed_final.fa

cat final_rdna_kmers.fa renamed_final.fa > kmers.fa #final kmer file


