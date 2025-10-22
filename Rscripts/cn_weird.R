
tmp <- read.delim("zc4688/honours/analyses/cn/result/eighteen/HG02050.eighteen.kmermatches.txt", header = F)
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

tmp2 <- tmp %>% group_by(seqID) %>% summarise(n_p = sum(Strand == "+"), n_n = sum(Strand == "-")) %>% filter(n_n > 10 & n_p > 10) %>% pull(seqID)
tmp <- tmp %>%
  mutate(
    posn = pos - 99,
    startn = ifelse(Strand == "+", start - posn, start - (1869 - posn - 18)),
    endn = ifelse(Strand == "+", end + (1869 - posn - 18), end + posn),
    startn = ifelse(startn < 0, 0, startn)
  )

ggplot(tmp %>% filter(seqID == "2ebd445e-c92e-4329-b2ba-4f229004c023")) + geom_segment(aes(x = start, xend = end, y = pos, yend = pos + 18, color = Strand)) + geom_segment(aes(x = startn, xend = endn, y = 1))



tmp %>% filter(seqID %in% tmp2) %>%ggplot() + geom_segment(aes(x = start, xend = end, y = pos, yend = pos + 18, color = Strand)) + facet_wrap(~seqID)
test <- read_table("zc4688/honours/analyses/cn/result/eighteen/HG02050.primary.txt", col_names = F, comment = "#")
ggplot(test, aes(x = X9, xend = X10, y = X12, color = X3)) + geom_segment()


############# GC
#tmp is old sample info with no gc filters on ntsm/genome kmers. new cn calculations have removed all kmers with gec < 0.25 or > 0.65

tmp %>% 
  select(sampleid, cn, Superpopulation.name) %>% 
  dplyr::rename("cn_old" = "cn") %>% 
  left_join(cn) %>% 
  mutate(diff = round(cn - cn_old)) %>% 
  left_join(chem) %>%  ggplot(aes(x = sampleid, y = diff, color = chem)) + geom_point()

#many samples have no diff, but restricting the gc only ever INCREASES cn calculations
#suggests the gc abnormal kmers are incorrectly deflating the genome coverage 




########## VCF
vcf<-read.vcfR("/g/data/te53/zc4688/honours/analyses/cn/result/fasta/test.vcf")

 
info<- vcf@fix[, "POS"]
info <- data.frame(pos = vcf@fix[, "POS"], ref = vcf@fix[, "REF"], alt = vcf@fix[, "ALT"])

info$DP<-extract.info(vcf,
                       element = "DP",
                       mask = FALSE,
                       as.numeric=TRUE)

info$AC<-extract.info(vcf,
                       element = "AC",
                       mask = FALSE,
                       as.numeric=FALSE)

info$EAS<-extract.info(vcf,
                       element = "AC_EAS_unrel",
                       mask = FALSE,
                       as.numeric=TRUE)

info$EUR<-extract.info(vcf,
                       element = "AC_EUR_unrel",
                       mask = FALSE,
                       as.numeric=TRUE)

info$SAS<-extract.info(vcf,
                       element = "AC_SAS_unrel",
                       mask = FALSE,
                       as.numeric=TRUE)
