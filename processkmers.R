#!/usr/bin/env Rscript
.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

library(tidyverse)
library(stringr)
library(Biostrings)
library(GenomicRanges)
sampleid <- commandArgs(trailingOnly = TRUE)[1]
gene <- commandArgs(trailingOnly = TRUE)[2]


rdna_categories <- read.delim("/g/data/te53/zc4688/honours/analyses/cn/rdna_categories.txt", header = F)
colnames(rdna_categories) <- c("pos", "kmer", "category")
all_kmers <- read.delim("/g/data/te53/zc4688/honours/analyses/cn/result/kmers.fa", header = F)

all_kmers <- all_kmers %>% 
  mutate(type = str_replace(lag(V1), ">", "")) %>% 
  filter(str_detect(type, "kmer|genome|ref")) %>% 
  left_join(rdna_categories %>% dplyr::rename("type" = "kmer")) %>% 
  mutate(type = case_when(
    str_detect(type, "kmer") ~ "rdna", 
    str_detect(type, "ref") ~ "ntsm", 
    TRUE ~ type)) %>% 
  dplyr::rename("seq" = "V1") %>% 
  mutate(rc = as.character(reverseComplement(DNAStringSet(seq)))) %>% 
  mutate(category = ifelse(str_detect(category, "ITS"), "ITS", category))


file <- sprintf("/g/data/te53/zc4688/honours/analyses/cn/result/%s.%s.kmermatches.txt", sampleid, gene)

# read the file
tmp <- read.delim(file, header = FALSE)
colnames(tmp) <- c("seqID", "start", "end", "Score", "Strand", "pattern", "Mismatch")
keep <- tmp %>% group_by(seqID) %>% summarise(n = n()) %>% filter(n > 5) %>% pull(seqID)

category_to_keep <- switch(gene,
                           "eighteen" = "18S_rRNA",
                           "fiveeight" = "5_8S_rRNA",
                           "twoeight" = "28S_rRNA")


data <- all_kmers %>% filter(category == category_to_keep) %>% 
  select(seq, pos, rc)

tmp <- tmp %>% 
  filter(seqID %in% keep) %>% 
  select(seqID, pattern, start, end, Strand) %>% 
  mutate(pattern = str_replace_all(pattern, "pattern:", "")) %>% 
  left_join(data, by = c("pattern" = "seq")) %>% 
  mutate(seq_pos = pos) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc,  -pos) %>% 
  left_join(data, by = c("pattern" = "rc")) %>% 
  mutate(pattern = ifelse(match == "seq", pattern, seq), 
         pos = ifelse(match == "seq", seq_pos, pos))

if (gene == "eighteen") {
  tmp <- tmp %>%
    mutate(
      pos = pos - 99,
      start = ifelse(Strand == "+", start - pos, start - (1869 - pos - 18)),
      end = ifelse(Strand == "+", end + (1869 - pos - 18), end + pos),
      start = ifelse(start < 0, 0, start)
    )
} else if (gene == "twoeight") {
  tmp <- tmp %>%
    mutate(
      pos = pos - 4375,
      start = ifelse(Strand == "+", start - pos, start - (5065 - pos - 18)),
      end = ifelse(Strand == "+", end + (5065 - pos - 18), end + pos),
      start = ifelse(start < 0, 0, start)
    )
} else if (gene == "fiveeight") {
  tmp <- tmp %>%
    mutate(
      pos = pos - 3046,
      start = ifelse(Strand == "+", start - pos, start - (153 - pos - 18)),
      end = ifelse(Strand == "+", end + (153 - pos - 18), end + pos),
      start = ifelse(start < 0, 0, start)
    )
} else if (gene == "ITS") {
  if (pos ...){
    mutate(
      pos = pos - 3046,
      start = ifelse(Strand == "+", start - pos, start - (153 - pos - 18)),
      end = ifelse(Strand == "+", end + (153 - pos - 18), end + pos),
      start = ifelse(start < 0, 0, start)
    )
  else if (pos ...){
    tmp <- tmp %>%
    mutate(
      pos = pos - 3046,
      start = ifelse(Strand == "+", start - pos, start - (153 - pos - 18)),
      end = ifelse(Strand == "+", end + (153 - pos - 18), end + pos),
      start = ifelse(start < 0, 0, start)
    )
  }
  }
}

tmp <- GRanges(
  seqnames = tmp$seqID,
  ranges   = IRanges(start = tmp$start, end = tmp$end)
)

tmp <- GenomicRanges::reduce(tmp, min.gapwidth = 0)
tmp <- data.frame(name = as.character(seqnames(tmp)), start = start(tmp), end = end(tmp))

tmp <- tmp %>% select(name, start, end) 

write.table(
  tmp, 
  sprintf("/g/data/te53/zc4688/honours/analyses/cn/result/%s.%s.coord.bed", sampleid, gene), 
  col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE
)