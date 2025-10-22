#!/bin/bash

#for f in /g/data/te53/zc4688/honours/analyses/chordates/new/*.rDNA.refmorph.fasta; do awk -v fn=$(basename "$f" .rDNA.refmorph.fasta) '/^>/{$0=">"fn"_"substr($0,2)}1' "$f"; done > /g/data/te53/zc4688/honours/analyses/chordates/new/msa/combined.fasta

module load  samtools blast/2.14.1 mcl/14-137
awk -F, '{print $1 ":" $2 "-" $3}' eighteen.csv > eighteen_regions.txt
awk -F, '{print $1 ":" $2 "-" $3}' twoeight.csv > twoeight_regions.txt
awk -F, '{print $1 ":" $2 "-" $3}' fiveeight.csv > fiveeight_regions.txt

    
samtools faidx -r eighteen_regions.txt combined.fasta > eighteen.fa
samtools faidx -r twoeight_regions.txt combined.fasta > twoeight.fa
samtools faidx -r fiveeight_regions.txt combined.fasta > fiveeight.fa

makeblastdb -in eighteen.fa -dbtype nucl
blastn -query eighteen.fa -db eighteen.fa -outfmt 6 >eighteen.self.blastn.out
awk '$3>=99' eighteen.self.blastn.out >eighteen.self.blastn.abc
mcxload -abc eighteen.self.blastn.abc -o eighteen.self.blastn.mci -write-tab eighteen.self.blastn.tab
mcl eighteen.self.blastn.mci -use-tab eighteen.self.blastn.tab


makeblastdb -in twoeight.fa -dbtype nucl
blastn -query twoeight.fa -db twoeight.fa -outfmt 6 >twoeight.self.blastn.out
awk '$3>=99' twoeight.self.blastn.out >twoeight.self.blastn.abc
mcxload -abc twoeight.self.blastn.abc -o twoeight.self.blastn.mci -write-tab twoeight.self.blastn.tab
mcl twoeight.self.blastn.mci -use-tab twoeight.self.blastn.tab

makeblastdb -in fiveeight.fa -dbtype nucl
blastn -query fiveeight.fa -db fiveeight.fa -outfmt 6 >fiveeight.self.blastn.out
awk '$3>=99' fiveeight.self.blastn.out >fiveeight.self.blastn.abc
mcxload -abc fiveeight.self.blastn.abc -o fiveeight.self.blastn.mci -write-tab fiveeight.self.blastn.tab
mcl fiveeight.self.blastn.mci -use-tab fiveeight.self.blastn.tab


makeblastdb -in its1.fa -dbtype nucl
blastn -query its1.fa -db its1.fa -outfmt 6 >its1.self.blastn.out
awk '$3>=80' its1.self.blastn.out >its1.self.blastn.abc
mcxload -abc its1.self.blastn.abc -o its1.self.blastn.mci -write-tab its1.self.blastn.tab
mcl its1.self.blastn.mci -use-tab its1.self.blastn.tab

#mafft eighteen.fa > eighteen_msa.txt
#mafft twoeight.fa > twoeight_msa.txt
#mafft fiveeight.fa > fiveeight_msa.txt

#trimal -in eighteen_msa.txt -out eighteen_msa_cleaned.txt -gt 0.2
#trimal -in twoeight_msa.txt -out twoeight_msa_cleaned.txt -gt 0.2
#trimal -in fiveeight_msa.txt -out fiveeight_msa_cleaned.txt -gt 0.2







#building new hmm
#seqtk subseq twoeight_clusters/final_alignment_with_all.fasta twoeightreps.txt > twoeightreps.fa
#trimal -in twoeightreps.fa -out test.fa -noallgaps
#mv test.fa twoeightreps.fa 
#hmmbuild test.hmm twoeightreps.fa 