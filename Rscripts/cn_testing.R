.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

library(tidyverse)
library(ComplexUpset)
 
library(Biostrings)
#ntsm output
test <- read.delim("zc4688/honours/analyses/cn/test.txt", comment.char = "#", header = F)

test %>% 
  mutate(ref = V4/V6, alt = V5/V7) %>% 
  select(V1, ref, alt) %>% 
  pivot_longer(cols = c("ref", "alt")) %>% 
  filter(value < 50 & value > 1) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_histogram(binwidth = 2)

tmp <- read.delim("zc4688/honours/analyses/cn/rdnanew.fa", header = F)

tmp <- tmp %>% mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% filter(!is.na(n))

#can see the heterozygous peak: *around 18)
tmp %>% filter(n < 50 & n > 5) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 1)
tmp %>% filter(n < 50 & n > 5) %>% 
  ggplot(aes(x = log(n))) + 
  geom_density()
#and rdna peaks: 4000
tmp %>% 
  filter(n > 5 & n < 10000) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 100)

#this would give cn of about 111?
tmp <- tmp %>% mutate(category = ifelse(V1 %in% rdna$V10, "rdna", "other"))
tmp %>% 
  filter(n > 5) %>% 
  ggplot(aes(x = log(n))) + 
  geom_density(aes(fill = category))

rdna <- read.delim("zc4688/honours/analyses/cn/final_sorted.sam", comment.char = "@", header = F)

#rdna plots for different parts of the unit:
rdna <- rdna %>% mutate(type = case_when(
  V4 > 100 & V4 < 1968 ~ "18S",
  V4 >= 1968 & V4 <= 4360 ~ "ITS",
  V4 > 4360 & V4 < 9435 ~ "28S",
  TRUE ~ "other"
))


tmp %>% filter(category == "rdna") %>% left_join(rdna %>% select(V10, type) %>% rename("V1" = "V10")) %>% ggplot(aes(x = n, fill = type)) + geom_density(alpha = 0.5)                   

tmp  %>% left_join(rdna %>% select(V10, type) %>% rename("V1" = "V10")) %>% ggplot(aes(x = log(n), fill = type)) + geom_density(alpha = 0.5)        

tmp %>% 
  left_join(rdna %>% select(V10, type) %>% rename("V1" = "V10")) %>% filter(n > 5) %>% 
  mutate(type = ifelse(is.na(type), "genome", type)) %>% 
  group_by(type) %>% 
  summarise(depth_median = median(n),
            depth_mode = as.numeric(names(sort(table(n), decreasing = TRUE)[1]))) %>% 
  view()

#use just 18s? 
categories <- read.delim("zc4688/honours/analyses/cn/kmers.fa", header=F) %>% 
  mutate(type = str_replace(lag(V1), ">", "")) %>% 
  filter(!str_detect(V1, ">")) %>% 
  mutate(type = case_when(
    str_detect(type, "kmer") ~ "rdna", 
    str_detect(type, "ref") ~ "ntsm", 
    TRUE ~ type))

tmp <- left_join(tmp, categories)

############################## categories

rdna_categories <- read.delim("zc4688/honours/analyses/cn/rdna_categories.txt", header = F)
colnames(rdna_categories) <- c("pos", "kmer", "category")
all_kmers <- read.delim("zc4688/honours/analyses/cn/result/kmers.fa", header = F)

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


###############################
hg01443 <- read.delim("zc4688/honours/analyses/cn/HG00729_kmers.fa", header = F)
hg01443 <- hg01443 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

hg01443 %>% filter(n > 5) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 1)

#must be an easier way but this works
hg01443 <- hg01443 %>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = category) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -category) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, category)) %>% 
  select(kmer, n, match, type, category)

ggplot(hg01443, aes(x = log(n), fill = type)) +
  geom_density() + theme_bw()

ggplot(hg01443 %>% filter(type != "rdna" & n > 5), aes(x = n, fill = type)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 100)) + theme_bw()

ggplot(hg01443 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw()

ggplot(hg01443 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density() +
  facet_wrap(~category) + theme_bw()

hg01443 %>% 
  mutate(test = ifelse(is.na(category), type, category)) %>% 
  filter(n > 5) %>% 
  group_by(test) %>% 
  summarise(depth = median(n)) %>% 
  mutate(cn = depth[which(test == "18S_rRNA")] / depth[which(test == "genome")])


###############################

hg02860 <- read.delim("zc4688/honours/analyses/cn/HG02860_kmers.fa_kmers.fa", header = F)
hg02860 <- hg02860 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

hg02860 %>% filter(n > 5) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 1) + theme_bw()

#must be an easier way but this works
hg02860 <- hg02860 %>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = V2) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -V2) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, V2)) %>% 
  select(kmer, n, match, type, category)

ggplot(hg02860, aes(x = log(n), fill = type)) +
  geom_density() + theme_bw()

ggplot(hg02860 %>% filter(type != "rdna" & n > 5), aes(x = n, fill = type)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 100)) + theme_bw()

ggplot(hg02860 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw()

ggplot(hg02860 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density() +
  facet_wrap(~category) + theme_bw()

hg02860 %>% 
  mutate(test = ifelse(is.na(category), type, category)) %>% 
  filter(n > 5) %>% 
  group_by(test) %>% 
  summarise(depth = median(n)) %>% 
  mutate(cn = depth[which(test == "18S_rRNA")] / depth[which(test == "genome")])

###############################

hg01607 <- read.delim("zc4688/honours/analyses/cn/HG00246.cn.kmers.fa", header = F)
hg01607 <- hg01607 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

hg01607 %>% filter(n > 5) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 1) + theme_bw()

#must be an easier way but this works
hg01607 <- hg01607 %>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = V2) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -V2) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, V2)) %>% 
  select(kmer, n, match, type, category)

ggplot(hg01607, aes(x = log(n), fill = type)) +
  geom_density() + theme_bw()

ggplot(hg01607 %>% filter(type != "rdna" & n > 5), aes(x = n, fill = type)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 100)) + theme_bw()

ggplot(hg01607 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw()

ggplot(hg01607 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density() +
  facet_wrap(~category) + theme_bw()

hg01607 %>% 
  mutate(test = ifelse(is.na(category), type, category)) %>% 
  filter(n > 5) %>% 
  group_by(test) %>% 
  summarise(depth = median(n)) %>% 
  mutate(cn = depth[which(test == "18S_rRNA")] / depth[which(test == "genome")])

###############################

hg03507 <- read.delim("zc4688/honours/analyses/cn/HG02820.cn.kmers.fa", header = F)
hg03507 <- hg03507 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

hg03507 %>% filter(n > 5) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(binwidth = 1) + theme_bw()

#must be an easier way but this works
hg03507 <- hg03507 %>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = V2) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -V2) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, V2)) %>% 
  select(kmer, n, match, type, category) %>% 
  filter(!is.na(type))

ggplot(hg03507, aes(x = log(n), fill = type)) +
  geom_density() + theme_bw()

ggplot(hg03507 %>% filter(type != "rdna" & n > 5), aes(x = n, fill = type)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 100)) + theme_bw()

ggplot(hg03507 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw()

ggplot(hg03507 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density() +
  facet_wrap(~category) + theme_bw()

hg03507 %>% 
  mutate(test = ifelse(is.na(category), type, category)) %>% 
  filter(n > 5) %>% 
  group_by(test) %>% 
  summarise(depth = median(n)) %>% 
  mutate(cn = depth[which(test == "18S_rRNA")] / depth[which(test == "genome")])



kmer_data <- list() 
kmer_files <- list.files("zc4688/honours/analyses/cn/result", pattern = "HG.*cn.kmers.fa$", full.names = TRUE)

for (file in kmer_files){
  sample_id <- str_replace(basename(file), "_kmers.fa", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
    filter(!is.na(n)) %>% 
    dplyr::rename("kmer" = "V1")%>%
    left_join(all_kmers, by = c("kmer" = "seq")) %>% 
    mutate(seq_type = type, seq_cat = category, seq_pos = pos) %>% 
    mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
    select(-rc, -type, -category, -pos) %>% 
    left_join(all_kmers, by = c("kmer" = "rc")) %>% 
    mutate(kmer = ifelse(match == "seq", kmer, seq), 
           type = ifelse(match == "seq", seq_type, type), 
           category = ifelse(match == "seq", seq_cat, category),
           pos = ifelse(match == "seq", seq_pos, pos)) %>% 
    select(kmer, n, match, type, category, pos) %>% 
    filter(!is.na(type)) %>% 
    mutate(sampleid = sample_id)
  kmer_data[[file]] <- data
}

kmer_data <- bind_rows(kmer_data)

chm13 <- read.delim("zc4688/honours/analyses/cn/chm13_kmers_result.fa", header = F)
chm13 <- chm13 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")%>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = category, seq_pos = pos) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -category, -pos) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, category),
         pos = ifelse(match == "seq", seq_pos, pos)) %>% 
  select(kmer, n, match, type, category, pos) %>% 
  filter(!is.na(type)) %>% 
  mutate(sampleid = "chm13")


kmer_data <- bind_rows(kmer_data, chm13)
#calculating cn based on 18S depth
cn <- kmer_data %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  filter(n > 5) %>%
  group_by(test, sampleid) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome)

tmp <- kmer_data %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  filter(n > 5) %>%
  group_by(test, sampleid) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome) %>% select(sampleid, `18S_rRNA`)

#want to look at similarity as well
test <- kmer_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  mutate(group = case_when(n < `18S_rRNA`/2 ~ "variable", TRUE ~ "common")) %>% #should be based on cn (lower for less depth samples)
  filter(group == "variable") %>% 
  mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  select(kmer, category, group, sampleid, pos) %>% 
  pivot_wider(id_cols = c("kmer", "category", "pos"), names_from = "sampleid", values_from = group, values_fill = 0)

upset(test, kmer_data %>% select(sampleid) %>% unique() %>% unlist(), 
      base_annotations = list(
        'Intersection Size'=intersection_size(
          counts=TRUE,
          mapping=aes(fill=category))
      )
)

test <- kmer_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  mutate(group = case_when(n < `18S_rRNA`/2 ~ "variable", TRUE ~ "common")) %>% #should be based on cn (lower for less depth samples)
  filter(group == "variable") %>% 
  mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  select(kmer, category, group, sampleid, pos) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name") %>% select(sampleid, Superpopulation.code)) %>% 
  group_by(Superpopulation.code, kmer, category, pos) %>% 
  summarise(number = n()) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "Superpopulation.code", values_from = number, values_fill = 0)

test$Range <- rowSums(test[, 2:6] > 0)

test$Type<-ifelse(test$Range == 5, 'Shared',
                  ifelse(test$Range %in% c(2, 3, 4,5), 'Widespread', NA))

test$Sum <- rowSums(test[, 2:6], na.rm = TRUE)
na_indices <- is.na(test$Type)

test$Type[na_indices] <- ifelse(test$Sum[na_indices] == 1, 'Private',
                                ifelse(test$Sum[na_indices] > 1, 'Population-specific', NA))

upset(test, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=brewer.pal(4, "Spectral"))+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)


kmer_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  mutate(group = case_when(n < `18S_rRNA`/2 ~ "variable", TRUE ~ "common")) %>% #should be based on cn (lower for less depth samples)
  filter(group == "variable") %>% 
  mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  select(kmer, category, group, sampleid, pos) %>% 
  group_by(kmer, category, pos) %>% 
  summarise(sum = sum(group)) %>% 
  ggplot(aes(x = pos, size = sum, color = category, y = 1)) +
  geom_jitter(height = 0.1, width = 0) +
  theme_bw()

kmer_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  group_by(kmer) %>% 
  mutate(sum = sum(n/`18S_rRNA`)) %>% 
  ggplot(aes(x = pos, y = sum, color = category)) +
  geom_point() +
  theme_bw()


#kmer distances - based on kl ratio. probaly categorise instead of present absent by the proportion of the expected number (eg 25% of the 18S mean or somethign)
tmp <- kmer_data %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  filter(n > 5) %>%
  group_by(test, sampleid) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome) %>% select(sampleid, `genome`, `18S_rRNA`) %>% 
  mutate(genome = ifelse(sampleid == "chm13", 1, genome))
kmer <- kmer_data %>% 
  filter(type == "ntsm") %>% 
  left_join(tmp) %>% 
  mutate(n = ifelse(n > genome/2, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 


# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


kdist <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))

for (i in 1:nrow(kmer)) {
  for (j in i:nrow(kmer)) {
    shared <- kmer[i,]>0 & kmer[j,]>0
    kdist[i,j] <- kdist[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
  }
}


kmer <- kmer_data %>% 
  filter(type == "rdna") %>% 
  left_join(tmp) %>% 
  mutate(n = ifelse(n > `18S_rRNA`/2, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 

sample_ids <- kmer$sampleid %>% str_extract("^HG\\d+|chm13")

# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


kdist_gene <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))

for (i in 1:nrow(kmer)) {
  for (j in i:nrow(kmer)) {
    shared <- kmer[i,]>0 & kmer[j,]>0
    kdist_gene[i,j] <- kdist_gene[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
  }
}

kmer <- kmer_data %>% 
  filter(type == "genome") %>% 
  left_join(tmp) %>% 
  mutate(n = ifelse(n > `genome`/1.5, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 


# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


kdist_genome <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))

for (i in 1:nrow(kmer)) {
  for (j in i:nrow(kmer)) {
    shared <- kmer[i,]>0 & kmer[j,]>0
    kdist_genome[i,j] <- kdist_genome[j,i] <- sum(kmer[i,][shared], kmer[j,][shared]) / (sum(kmer[i,]) + sum(kmer[j,]))
  }
}

colnames(kdist_gene) <- sample_ids
rownames(kdist_gene) <- sample_ids
colnames(kdist) <- sample_ids
rownames(kdist) <- sample_ids
colnames(kdist_genome) <- sample_ids
rownames(kdist_genome) <- sample_ids
annotations <- sampleinfo %>% filter(Sample.name %in% sample_ids) %>% select(Sample.name, Superpopulation.code) %>% column_to_rownames("Sample.name")
pheatmap::pheatmap(kdist_gene, annotation_row = annotations)
pheatmap::pheatmap(kdist, annotation_row = annotations)
pheatmap::pheatmap(kdist_genome, annotation_row = annotations)





kmer_data <- kmer_data %>% mutate(sampleid = str_extract(sampleid, "^HG\\d+|chm13"))

ggplot(kmer_data %>% left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")), 
       aes(x = log(n), fill = type)) +
  geom_density() + theme_bw() + facet_wrap(~sampleid, ) 

ggplot(kmer_data %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw() + facet_wrap(~sampleid)

ggplot(kmer_data %>% filter(type == "rdna") %>% mutate(category = ifelse(str_detect(category, "rRNA"), "gene", category)), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw() + facet_wrap(~sampleid)


sampleinfo <- read.delim("referencedata/onekg/igsr_samples.tsv")

kmer_data %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  filter(n > 5) %>%
  group_by(test, sampleid) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  write_tsv("zc4688/honours/analyses/cn/sample_info.tsv") %>% 
  ggplot(aes(x = cn, fill = Superpopulation.code)) + geom_boxplot() + theme_bw()



#this establishes similarity etc of common kmers. want to look more into the divergent kmers

eighteenkmers.sh will give all the 18S kmers per sample. look at divergence/sharing of these
#also construct 18S graphs per sample - measure difference between them
eighteen_data <- list() 
eighteen_files <- list.files("zc4688/honours/analyses/cn/result", pattern = "HG.*eighteen.fa$", full.names = TRUE)
for (file in eighteen_files){
  sample_id <- str_replace(basename(file), ".eighteen.fa", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
    filter(!is.na(n)) %>% 
    dplyr::rename("kmer" = "V1")%>%
    mutate(sampleid = sample_id) %>% 
    left_join(cn) %>% 
    filter(n > 25) 
  eighteen_data[[file]] <- data
}

eighteen_data <- bind_rows(eighteen_data)


eighteen_upset <- eighteen_data %>% 
  mutate(group = case_when(
    n < genome*5 ~ "very rare",
    n < `18S_rRNA`/2 ~ "variable", 
    TRUE ~ "common")) %>% 
  filter(group == "variable") %>% 
  mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  select(kmer, group, sampleid) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name") %>% select(sampleid, Superpopulation.code)) %>% 
  group_by(Superpopulation.code, kmer) %>% 
  summarise(number = n()) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "Superpopulation.code", values_from = number, values_fill = 0)

eighteen_upset$Range <- rowSums(eighteen_upset[, 2:6] > 0)

eighteen_upset$Type<-ifelse(eighteen_upset$Range == 5, 'Shared',
                  ifelse(eighteen_upset$Range %in% c(2, 3, 4,5), 'Widespread', NA))

eighteen_upset$Sum <- rowSums(eighteen_upset[, 2:6], na.rm = TRUE)
na_indices <- is.na(eighteen_upset$Type)

eighteen_upset$Type[na_indices] <- ifelse(eighteen_upset$Sum[na_indices] == 1, 'Private',
                                ifelse(eighteen_upset$Sum[na_indices] > 1, 'Population-specific', NA))

upset(eighteen_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=brewer.pal(4, "Spectral"))+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)




#can we compare the variable 18s positions to the conserved positions across evolution?

eighteen <- readDNAMultipleAlignment("zc4688/honours/analyses/chordates/new/msa/eighteen_clusters/final_alignment_with_all.fasta")

shannon_entropy <- function(column) {
  column <- column[column != "-"]
  number <- length(column)
  freqs <- table(column) / length(column)  # Frequency of each base
  -sum(freqs * log2(freqs), na.rm = TRUE)  # Shannon entropy formula
}

# Convert alignment to matrix
msa_matrix <- as.matrix(eighteen)

# Apply function to each column
entropy_values <- apply(msa_matrix, 2, shannon_entropy)

counts <- sapply(as.data.frame(msa_matrix), function(col) sum(col != "-"))


tmp <- data.frame(Position = 1:length(entropy_values), Score = entropy_values, count = counts)



#check chm13
chm13 <- read.delim("zc4688/honours/analyses/cn/chm13_kmers_result.fa", header = F)

chm13 <- chm13 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")%>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = category, seq_pos = pos) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -category, -pos) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, category),
         pos = ifelse(match == "seq", seq_pos, pos)) %>% 
  select(kmer, n, match, type, category, pos) %>% 
  filter(!is.na(type))

chm13 %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  filter(n > 5) %>%
  group_by(test) %>%
  summarise(depth = median(n), .groups = 'drop') 

ggplot(chm13 %>% filter(type == "rdna"), aes(x = n, fill = category)) +
  geom_density(alpha = 0.5) + theme_bw()




#de bruijn graph?
test <- read.delim("zc4688/honours/analyses/cn/sequence_var/tmp_kmers.fa", header = F)

test <- test %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

k <- nchar(test$kmer[1])
test <- test %>%
  mutate(prefix = substr(kmer, 1, k-1),
         suffix = substr(kmer, 2, k))
edges <- test %>%
  group_by(prefix, suffix) %>%
  summarise(weight = sum(n), .groups = "drop")

# Build igraph
g <- graph_from_data_frame(d = edges, directed = TRUE)
E(g)$weight <- edges$weight


#coverage - if we are capturing variaiton across the whole unit it should be reasonably constant. cov erage does drop towards end of IGS
#is this fixed by adding more reference units?
tmp <- read.delim("zc4688/honours/analyses/cn/tmp.paf", header = F)
ranges <- IRanges(start = tmp$V8, end = tmp$V9)
cov <- coverage(ranges)
cov_vec <- as.integer(cov)
cov_df <- data.frame(
  position = seq_along(cov_vec),
  depth = cov_vec
)

ggplot(cov_df, aes(x = position, y = depth)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Reference Position", y = "Coverage / Depth")
###########
tmp <- read.delim("zc4688/honours/analyses/cn/HG00257.filtered.fa", header = F)

tmp <- tmp %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")%>%
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = category, seq_pos = pos) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -category, -pos) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, category),
         pos = ifelse(match == "seq", seq_pos, pos)) %>% 
  select(kmer, n, match, type, category, pos) 

tmp %>%
  filter(n > 20 & n < 2500) %>%
  summarise(mode_n = {
    tab <- table(n)
    as.numeric(names(tab)[which.max(tab)])
  })

ggplot(tmp %>% filter(n > 15), aes(x = n)) + geom_density()

tmp <- tmp %>% mutate(
  category = case_when(
    n > 10000 ~ "super",
    n > 5000 & n < 10000 ~ "very",
    n > 3000 & n < 5000 ~ "common",
    n > 1000 & n < 3000 ~ "rare",
    n > 100 & n < 1000 ~ "very rare",
    TRUE ~ NA
  )
) %>%  filter(!is.na(category))

#the defined rdna kmers are present almost perfectly, very few missing
tmp %>% filter(type == "rdna") %>% 
  select(kmer, n, category) %>% 
  left_join(kmer_data %>% 
              filter(sampleid == "HG00729") %>% 
              dplyr::rename("n_od" = "n")) %>% 
  mutate(missing = n_od - n) %>% view


kmer_data %>% filter(sampleid == "HG00257" & type == "rdna") %>% dplyr::rename("n_od" = "n") %>% left_join(tmp %>% filter(type == "rdna") %>% select(kmer, n, category)) %>% mutate(missing = n_od - n) %>% ggplot(aes(x = pos, y = missing)) + geom_point()



####comparing all rdna kmers

hg00729 <- read.delim("zc4688/honours/analyses/cn/manual/kmers.fa", header = F)

hg00729 <- hg00729 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")
ggplot(hg00729 %>% filter(n > 15 & n < 10000), aes(x = n)) + geom_density()
tmp <- read.delim("zc4688/honours/analyses/cn/tmp2.filtered.fa", header = F)

tmp <- tmp %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")
ggplot(tmp %>% filter(n > 100 & n < 10000), aes(x = n)) + geom_histogram(binwidth = 100)

tmp <- tmp %>% filter(n > 100 & n < 10000)
hg00729 <- hg00729 %>% filter(n > 100 & n < 10000)



hg00729 <- hg00729 %>% mutate(
  category = case_when(
    n > 10000 ~ "super",
    n > 5000 & n < 10000 ~ "very",
    n > 3000 & n < 5000 ~ "common",
    n > 1000 & n < 3000 ~ "rare",
    n > 100 & n < 1000 ~ "very rare",
    TRUE ~ NA
  )
) %>%  filter(!is.na(category)) %>% 
  mutate(sampleid = "hg00729")

tmp <- read.delim("zc4688/honours/analyses/cn/HG00155_kmers.fa", header = F)

tmp <- tmp %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

tmp <- tmp %>% mutate(
  category = case_when(
    n > 10000 ~ "super",
    n > 5000 & n < 10000 ~ "very",
    n > 3000 & n < 5000 ~ "common",
    n > 1000 & n < 3000 ~ "rare",
    n > 100 & n < 1000 ~ "very rare",
    TRUE ~ NA
  )
) %>%  filter(!is.na(category)) %>% 
  mutate(sampleid = "hg02820")

tmp <- bind_rows(tmp, hg00729)

common <- tmp %>% 
  filter(category == "very rare") %>% 
  mutate(category = 1) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "sampleid", values_from = category, values_fill = 0)

upset(common, tmp %>% select(sampleid) %>% unique() %>% unlist(), 
      base_annotations = list(
        'Intersection Size'=intersection_size(
          counts=TRUE)
      )
)


#paf filtering
tmp <- read.delim("zc4688/honours/analyses/cn/HG02", header = F)
tmp <- tmp %>%
  arrange(V1, V3) %>%
  group_by(V1) %>%
  summarise(
    read_start = dplyr::first(V3),
    read_end   = dplyr::last(V4),
    rdna_start = min(V8),
    rdna_end   = max(V9), read_length = dplyr::first(V2)
  ) %>%
  ungroup() %>% view


#doesnt help with the errors at all reallt just makes more missing
tmp <- tmp %>% mutate(
  keep = case_when(
    read_length > 40000 & (rdna_end - rdna_start) < 5000 ~ "no",
    (rdna_end - rdna_start) < 500 & (read_end - read_start)/read_length < 0.8 ~ "no",
    TRUE ~ "yes"
  )
) %>% filter(keep == "yes")


tmp <- tmp %>% mutate(
  start = ifelse(read_start < rdna_start, 0, read_start),
  end = ifelse((read_length - read_end) < (45000 - rdna_end), read_length, read_end)
) %>% 
  select(V1, start, end)



#tmp is the result of just minimap filtering, tmp2 has the kmer filtering first. most of the normal frequency kmers are maintained with little to none missing. 
#there are lines of frequencies - not sure of cause?
tmp <- read.delim("zc4688/honours/analyses/cn/tmp.filtered.fa", header = F)

tmp <- tmp %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1") %>% 
  filter(n > 15)

tmp2 <- read.delim("zc4688/honours/analyses/cn/tmp3_31.filtered.fa", header = F)

tmp2 <- tmp2 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")

tmp2 <- tmp2 %>% 
  left_join(all_kmers, by = c("kmer" = "seq")) %>% 
  mutate(seq_type = type, seq_cat = category, seq_pos = pos) %>% 
  mutate(match = ifelse(!is.na(rc), "seq", "rc")) %>% 
  select(-rc, -type, -category, -pos) %>% 
  left_join(all_kmers, by = c("kmer" = "rc")) %>% 
  mutate(kmer = ifelse(match == "seq", kmer, seq), 
         type = ifelse(match == "seq", seq_type, type), 
         category = ifelse(match == "seq", seq_cat, category),
         pos = ifelse(match == "seq", seq_pos, pos)) %>% 
  select(kmer, n, match, type, category, pos) 

tmp <- tmp  %>% dplyr::rename("n_tmp" = "n") %>% left_join(tmp2) %>%  mutate(missing = n_tmp - n)
ggplot(tmp %>% filter(n_tmp > 15 & n_tmp < 10000), aes(x = n_tmp, y = n, color = missing))  + geom_point()

ggplot(tmp %>% filter(n > 15), aes(x = n)) + geom_density()


tmp %>% filter(n_tmp > 15 & n_tmp < 10000) %>% mutate(missing = ifelse(is.na(missing), n_tmp, missing))  %>% filter(missing > 0 & missing < 1000) %>% ggplot(aes(x = missing)) + geom_histogram(binwidth = 20) + coord_cartesian(xlim = c(0, 1000))




paf <- read.delim("zc4688/honours/analyses/cn/HG00653/round3.paf", header=F)
paf <- paf %>% filter(V2 > 1000)
gr <- GRanges(
  seqnames = paf$V6,
  ranges   = IRanges(start = paf$V8, end = paf$V9)
)

# coverage per target (RleList)
cov <- coverage(gr)

cov_df <- do.call(rbind, lapply(names(cov), function(seqname) {
  cov_vec <- as.integer(cov[[seqname]])
  data.frame(
    seqname = seqname,
    position = seq_along(cov_vec),
    depth = cov_vec
  )
}))

ggplot(cov_df, aes(x = position, y = depth)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Reference Position", y = "Coverage / Depth") +
  facet_wrap(~seqname)


tmp <- read.delim("zc4688/honours/analyses/cn/HG00653/test.out", header = F)
tmp <- tmp %>% 
  mutate(name = paste0(V1, "__", V2)) %>% 
  filter(V3 > 90)

gr <- GRanges(
  seqnames = tmp$name,
  ranges   = IRanges(start = tmp$V7, end = tmp$V8)
)

gr <- GenomicRanges::reduce(gr, min.gapwidth = 10)

gr <- data.frame(name = as.character(seqnames(gr)), start = start(gr), end = end(gr))
gr <- gr %>% 
  group_by(name) %>% 
  summarise(length = sum(end - start)) %>% 
  separate(name, into = c("seq1", "seq2"), sep = "__") %>% 
  pivot_wider(names_from = seq2, values_from = length) 
gr <- gr %>% column_to_rownames("seq1")


pheatmap(gr, show_colnames = F, show_rownames = F)




#18s
tmp <- read.delim("zc4688/honours/analyses/cn/HG00653/matches.txt")
keep <- tmp %>% group_by(seqID) %>% summarise(n = n()) %>% filter(n > 2) %>% pull(seqID)
eighteen <- all_kmers %>% filter(category == "18S_rRNA") %>% select(seq, pos, rc)

tmp <- tmp %>% 
  filter(seqID %in% keep) %>% 
  select(seqID, pattern, start, end) %>% 
  left_join(eighteen, by = c("pattern" = "seq")) 

tmp <- tmp %>%
  arrange(seqID, start) %>% 
  group_by(seqID) %>% 
  mutate(pos = pos - 100,
         start = start - pos,
         end = end + (1968 - pos + 19)) %>% 
  mutate(start = ifelse(start < 0, 0, start))

tmp <- GRanges(
  seqnames = tmp$seqID,
  ranges   = IRanges(start = tmp$start, end = tmp$end)
)

tmp <- GenomicRanges::reduce(tmp, min.gapwidth = 0)
tmp <- data.frame(name = as.character(seqnames(tmp)), start = start(tmp), end = end(tmp))

tmp <- tmp %>% select(name, start, end) 

write.table(tmp, "zc4688/honours/analyses/cn/HG00653/coord.bed", col.names=F, row.names=F, sep = "\t", quote=F)

hg00729 <- read.delim("zc4688/honours/analyses/cn/eighteen.fa", header = F)

hg00729 <- hg00729 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")
ggplot(hg00729 %>% filter(n > 25 & n < 10000), aes(x = n)) + geom_histogram(binwidth = 100)


hg00653 <- read.delim("zc4688/honours/analyses/cn/HG00653/eighteen.fa", header = F)

hg00653 <- hg00653 %>% 
  mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
  filter(!is.na(n)) %>% 
  dplyr::rename("kmer" = "V1")
ggplot(hg00729 %>% filter(n > 25 & n < 10000), aes(x = n)) + geom_histogram(binwidth = 100)

hg00729 <- hg00729 %>% filter(n > 25)
