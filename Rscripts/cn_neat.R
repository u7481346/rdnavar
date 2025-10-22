.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

library(tidyverse)
library(ComplexUpset)
library(pheatmap)
library(Biostrings)
library(RColorBrewer)
library(ggpubr)
library(philentropy)

############################## categories
sampleinfo <- read.delim("referencedata/onekg/igsr_samples.tsv")

rdna_categories <- read.delim("zc4688/honours/analyses/cn/rdna_categories.txt", header = F)
colnames(rdna_categories) <- c("pos", "kmer", "category")
all_kmers <- read.delim("zc4688/honours/analyses/cn/kmers.fa", header = F)

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



chem <- read_delim("zc4688/honours/analyses/cn/samplechemistry.tsv", col_names = F)
colnames(chem) <- c("sampleid", "chem")
######################
cn_data <- list() 
cn_files <- list.files("zc4688/honours/analyses/cn/result/cn", pattern = "HG.*cn.kmers.fa$", full.names = TRUE)

for (file in cn_files){
  sample_id <- str_replace(basename(file), ".cn.kmers.fa", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(n = as.numeric(str_replace(lag(V1), ">", ""))) %>% 
    filter(!is.na(n)) %>% 
    filter(n > 5) %>% 
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
  cn_data[[file]] <- data
}

cn_data <- bind_rows(cn_data)

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


cn_data <- bind_rows(cn_data, chm13)

write_tsv(cn_data, "zc4688/honours/analyses/cn/result/cn_data.tsv")
cn_data <- read_tsv("zc4688/honours/analyses/cn/result/cn_data.tsv")

#calculating cn based on 18S depth
cn <- cn_data %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  group_by(test, sampleid) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome)

cn_all <- cn_data %>% 
  group_by(type, sampleid) %>% 
  summarise(`All rDNA` = median(n), .groups="drop") %>% 
  filter(type == "rdna") %>% 
  select(-type)

cn <- left_join(cn, cn_all)
cn %>% left_join(chem) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  ggplot() + 
  geom_boxplot(aes(x = Superpopulation.code, fill = chem, y = cn))

kruskal.test(cn ~ Superpopulation.code, data = cn %>% left_join(chem) %>% 
               left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")))

cn %>% left_join(chem) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% filter(sampleid != "chm13") %>% 
  ggplot(aes(x = Superpopulation.code, fill = chem, y = cn)) + 
  geom_boxplot() + theme_bw() + scale_fill_brewer(palette = "Set1") + labs(x = "Continent", y = "Copy number", fill = "Chemistry") +
  stat_compare_means(
    aes(group = chem), 
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE
  )

ggsave("zc4688/honours/analyses/cn/figures/cn_chemistry.png")

tmp <- read_tsv("zc4688/honours/analyses/cn/41598_2020_80049_MOESM1_ESM.txt")
cn %>% 
  dplyr::rename("Sample" = "sampleid") %>% 
  left_join(tmp) %>% 
  ggplot(aes(x = cn, y = HC.18S.CN)) + 
  geom_point() + theme_bw() + labs(x = "New estimate", y = "Published value")

ggsave("zc4688/honours/analyses/cn/figures/published.comparison.png")

tmp <- tmp %>% left_join(cn %>% dplyr::rename("Sample" = "sampleid"))
cor.test(tmp %>% pull(HC.18S.CN), tmp %>% pull(cn), method = "pearson")

cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  pivot_longer(cols = c(2:9), names_to = "Variable", values_to = "value") %>% 
  filter(str_detect(Variable, "rRNA") | Variable %in% c("ITS", "IGS") | Variable == "All rDNA") %>% 
  filter(sampleid != "chm13") %>% 
  ggplot() + 
  geom_boxplot(aes(x = Variable, y = value, fill = Superpopulation.code)) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + 
  labs(fill = "Continent",
       x = NULL,
       y = "Depth")
ggsave("zc4688/honours/analyses/cn/figures/depth.png")

cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  select(sampleid, `18S_rRNA`, `28S_rRNA`, `5_8S_rRNA`, ITS, IGS, `All rDNA`, genome, Superpopulation.code) %>% 
  pivot_longer(cols = c(2:7), names_to = "Variable", values_to = "value") %>% 
  mutate(n = value/genome) %>% 
  filter(sampleid != "chm13") %>% 
  ggplot() + 
  geom_boxplot(aes(x = Variable, y = n, fill = Superpopulation.code)) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + 
  labs(fill = "Continent",
       x = NULL,
       y = "CN")
ggsave("zc4688/honours/analyses/cn/figures/cn_by_gene.png")

cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  select(sampleid, `18S_rRNA`, `28S_rRNA`, `5_8S_rRNA`, ITS, IGS, `All rDNA`, genome, Superpopulation.code) %>% 
  pivot_longer(cols = c(2:7), names_to = "Variable", values_to = "value") %>% 
  mutate(n = value/genome) %>% 
  filter(sampleid != "chm13") %>% 
  ggplot(aes(x = Superpopulation.code, y = n, fill = Variable)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set1") + 
  labs(fill = "Variable",
       x = NULL,
       y = "CN")+
  stat_compare_means(
    aes(group = Variable), 
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE
  )

tmp <- cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  select(sampleid, `18S_rRNA`, `28S_rRNA`, `5_8S_rRNA`, ITS, IGS, `All rDNA`, genome, Superpopulation.code) %>% 
  pivot_longer(cols = c(2:7), names_to = "Variable", values_to = "value") %>% 
  mutate(n = value/genome) %>% 
  filter(sampleid != "chm13") %>% mutate(group = paste0(Superpopulation.code, "_", Variable))

pairwise.wilcox.test(tmp %>% pull(n),
                     tmp %>% pull(group),
                     p.adjust.method = "BH")
ggsave("zc4688/honours/analyses/cn/figures/cn_by_gene2.png")

cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  filter(sampleid != "chm13") %>% 
  ggplot(aes(x = Sex, fill = Sex, y = cn)) + geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") + 
  labs(x = NULL, y = "CN")+
  stat_compare_means(
    aes(group = Sex), 
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE
  )

ggsave("zc4688/honours/analyses/cn/figures/cn_sex.png")


cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  filter(sampleid != "chm13") %>% 
  ggplot(aes(x = Superpopulation.code, fill = Sex, y = cn)) + geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") + 
  labs(x = "Continent", y = "CN")+
  stat_compare_means(
    aes(group = Sex), 
    method = "wilcox.test", 
    label = "p.signif", 
    hide.ns = TRUE
  )

ggsave("zc4688/honours/analyses/cn/figures/cn_sex_pop.png")



cn_upset <- cn_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  mutate(group = case_when(n < `18S_rRNA`/2 ~ "variable", TRUE ~ "common")) %>% #should be based on cn (lower for less depth samples)
  filter(group == "variable") %>% 
  mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  select(kmer, category, group, sampleid, pos) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name") %>% select(sampleid, Superpopulation.code)) %>% 
  group_by(Superpopulation.code, kmer, category, pos) %>% 
  summarise(number = n()) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "Superpopulation.code", values_from = number, values_fill = 0)

cn_upset$Range <- rowSums(cn_upset[, 2:6] > 0)

cn_upset$Type<-ifelse(cn_upset$Range == 5, 'Shared',
                  ifelse(cn_upset$Range %in% c(2, 3, 4,5), 'Widespread', NA))

cn_upset$Sum <- rowSums(cn_upset[, 2:6], na.rm = TRUE)
na_indices <- is.na(cn_upset$Type)

cn_upset$Type[na_indices] <- ifelse(cn_upset$Sum[na_indices] == 1, 'Private',
                                ifelse(cn_upset$Sum[na_indices] > 1, 'Population-specific', NA))

cn_upset <- cn_upset %>%
  mutate(Commonness = case_when(
    Sum == 1 ~ "Private",
    Sum <= 5 ~ "Rare",
    Sum <= 122 ~ "Intermediate",
    Sum > 122 ~ "Common"
  ))


#add category showing number of people in that pouplatipn
upset(cn_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Commonness)
        )+scale_y_continuous()+scale_fill_manual(values=brewer.pal(4, "RdYlBu"))+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)

ggsave('zc4688/honours/analyses/cn/figures/cn_upset_freq.png')

upset(cn_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
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

ggsave('zc4688/honours/analyses/cn/figures/cn_upset_type.png')




cn_data %>% filter(type == "rdna") %>% 
  left_join(cn) %>% 
  #mutate(group = case_when(n < `18S_rRNA`/2 ~ "variable", TRUE ~ "common")) %>% #should be based on cn (lower for less depth samples)
  #filter(group == "variable") %>% 
  #mutate(group = ifelse(group == "variable", 1, 0)) %>% 
  mutate(n = n/`18S_rRNA`) %>% 
  group_by(kmer, category, pos) %>% 
  summarise(sum = sum(n)) %>% 
  ggplot(aes(x = pos,color = category, y = sum)) +
  geom_jitter(height = 0.1, width = 0) +
  theme_bw() + 
  scale_color_brewer(palette = "Dark2")
ggsave('zc4688/honours/analyses/cn/figures/pos_freq.png')

#kmer distances - based only on presence absence
cn <- cn %>%
  mutate(genome = ifelse(sampleid == "chm13", 1, genome))
kmer <- cn_data %>% 
  filter(type == "ntsm") %>% 
  left_join(cn) %>% 
  mutate(n = ifelse(n > genome/2, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 


# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


ntsm <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))
rowsums <- rowSums(kmer)

for (i in 1:nrow(kmer)) {
  ki <- kmer[i, ]   # extract row once
  for (j in i:nrow(kmer)) {
    kj <- kmer[j, ] # extract row once
    shared <- ki > 0 & kj > 0
    ntsm[i, j] <- ntsm[j, i] <- sum(ki[shared] + kj[shared]) / (rowsums[i] + rowsums[j])
  }
}


kmer <- cn_data %>% 
  filter(type == "rdna") %>% 
  left_join(cn) %>% 
  mutate(n = ifelse(n > `18S_rRNA`/2, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 

sample_ids <- kmer$sampleid

# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


rdna <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))

rowsums <- rowSums(kmer)

for (i in 1:nrow(kmer)) {
  ki <- kmer[i, ]   # extract row once
  for (j in i:nrow(kmer)) {
    kj <- kmer[j, ] # extract row once
    shared <- ki > 0 & kj > 0
    rdna[i, j] <- rdna[j, i] <- sum(ki[shared] + kj[shared]) / (rowsums[i] + rowsums[j])
  }
}

kmer <- cn_data %>% 
  filter(type == "genome") %>% 
  left_join(cn) %>% 
  mutate(n = ifelse(n > `genome`/1.5, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0) 


# Convert to numeric matrix excluding sampleid column
kmer <- as.matrix(kmer[,-1])


genome <- matrix(nrow = nrow(kmer), ncol = nrow(kmer))

rowsums <- rowSums(kmer)

for (i in 1:nrow(kmer)) {
  ki <- kmer[i, ]   # extract row once
  for (j in i:nrow(kmer)) {
    kj <- kmer[j, ] # extract row once
    shared <- ki > 0 & kj > 0
    genome[i, j] <- genome[j, i] <- sum(ki[shared] + kj[shared]) / (rowsums[i] + rowsums[j])
  }
}

colnames(rdna) <- sample_ids
rownames(rdna) <- sample_ids
colnames(ntsm) <- sample_ids
rownames(ntsm) <- sample_ids
colnames(genome) <- sample_ids
rownames(genome) <- sample_ids
annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% sample_ids) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  column_to_rownames("Sample.name")
pheatmap::pheatmap(rdna, annotation_row = annotations)

pheatmap::pheatmap(genome, annotation_row = annotations)



###################JS
rdna <- cn_data %>%
  left_join(cn) %>% 
  filter(type == "rdna") %>% 
  mutate(n = n/genome) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- rdna / rowSums(rdna)

rdna <- distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)
rownames(rdna) <- rownames(prob_matrix)
colnames(rdna) <- rownames(prob_matrix)

genome <- cn_data %>%
  filter(type == "genome") %>% 
  left_join(cn) %>% 
  mutate(n = n/genome) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- genome / rowSums(genome)

# Compute symmetric KL divergence matrix
genome <- distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)
rownames(genome) <- rownames(prob_matrix)
colnames(genome) <- rownames(prob_matrix)



annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  left_join(chem %>% dplyr::rename("Sample.name" = "sampleid")) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code")

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   chem = c("R10" = "darkred", "R9" = "forestgreen"))
pheatmap::pheatmap(rdna %>% 
                     as.data.frame() %>% 
                     select(-chm13) %>% 
                     filter(rownames(.) != "chm13"), 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F)


ggsave("zc4688/honours/analyses/cn/figures/rdna_cn_heatmap.png")

pheatmap::pheatmap(genome %>% 
                     as.data.frame() %>% 
                     select(-chm13) %>% 
                     filter(rownames(.) != "chm13"), 
                   annotation_row = annotations,
                   annotation_colors = ann_colors,
                   annotation_col = annotations,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F)

ggsave("zc4688/honours/analyses/cn/figures/genome_heatmap.png")





cn %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")) %>% 
  write_tsv("zc4688/honours/analyses/cn/sample_info.tsv") %>% 
  ggplot(aes(x = cn, fill = Superpopulation.code)) + geom_boxplot() + theme_bw()

kruskal.test(cn ~ Superpopulation.code, data = cn %>% 
               left_join(sampleinfo %>% dplyr::rename(sampleid = "Sample.name")))

comparisons <- list( c("AFR", "AMR"), c("AFR", "EUR"), c("AFR", "SAS") , c("EAS", "EUR"), c("EAS", "SAS"))

ggboxplot(cn %>% 
            left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name")), 
          x = "Superpopulation.code", y = "cn",
               color = "Superpopulation.code", palette =c("#00AFBB", "#E7B800", "#FC4E07", "pink", "forestgreen"),
               add = "jitter") + stat_compare_means(comparisons = comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)  

pairwise.wilcox.test(cn %>% left_join(sampleinfo %>% dplyr::rename(sampleid = "Sample.name")) %>% pull(cn),
                     cn %>% left_join(sampleinfo %>% dplyr::rename(sampleid = "Sample.name")) %>% pull(Population.code),
                     p.adjust.method = "BH")

cn %>% 
  left_join(sampleinfo %>% dplyr::rename(sampleid = "Sample.name")) %>% 
  filter(!is.na(Superpopulation.code)) %>% 
  ggplot(aes(x = Superpopulation.code, y = cn, fill = Superpopulation.code)) +
  geom_boxplot() +
  geom_hline(yintercept = 223, color = "red") + 
  geom_jitter(width = 0.3, height = 0) + 
  theme_bw() + # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)  +
  scale_fill_brewer(palette = "Set2") + 
  labs(fill = "Continent", x = NULL, y = "Copy number")

ggsave("zc4688/honours/analyses/cn/figures/cn_sign.png")


cn %>% 
  left_join(sampleinfo %>% dplyr::rename(sampleid = "Sample.name")) %>% 
  filter(!is.na(Population.code)) %>% 
  mutate(
    Population.code = factor(
      Population.code,
      levels = unique(Population.code[order(Superpopulation.code, Population.code)])
    )
  ) %>%
  ggplot(aes(x = Population.code, y = cn, fill = Superpopulation.code)) +
  geom_boxplot() +
  geom_hline(yintercept = 224, color = "red") + 
  geom_jitter(width = 0.3, height = 0) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Set2") + 
  labs(fill = "Continent", x = "Population", y = "Copy number")

ggsave("zc4688/honours/analyses/cn/figures/cn_pop.png")

#this establishes similarity etc of common kmers. want to look more into the divergent kmers





#eighteenkmers.sh will give all the 18S kmers per sample. look at divergence/sharing of these
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

eighteen_data <- read_tsv("zc4688/honours/analyses/cn/result/eighteen/eighteen.kmers.tsv")
eighteen_data <- eighteen_data %>% left_join(cn) %>% filter(n > genome)
#eighteen_data <- bind_rows(eighteen_data)
#write_tsv(eighteen_data, "zc4688/honours/analyses/cn/result/eighteen_data.tsv")

eighteen_upset <- eighteen_data %>% 
  mutate(n = ifelse(n > 0, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
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

eighteen_upset <- eighteen_upset %>%
  mutate(Type = case_when(
    Type == "Shared" & Sum >=122 ~ "Common",
    TRUE ~ Type
  ))

colors <- setNames(brewer.pal(5, "Paired"), c("Common", "Shared", "Private", "Population-specific", "Widespread"))
upset(eighteen_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=colors)+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)



ggsave("zc4688/honours/analyses/cn/figures/eighteen_all_upset.png")



#KL
kmer <- eighteen_data %>% 
  mutate(n = n/genome) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- kmer / rowSums(kmer)

# Compute symmetric KL divergence matrix
kmer <- philentropy::distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)

write_tsv(as.data.frame(kmer), "zc4688/honours/analyses/cn/eighteen.js.tsv", col_names = T)
rownames(kmer) <- rownames(prob_matrix)
colnames(kmer) <- rownames(prob_matrix)


annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code")

annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  left_join(chem %>% dplyr::rename("Sample.name" = "sampleid")) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code") 

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   chem = c("R10" = "deepskyblue1", "R9" = "darkorchid"))
pheatmap::pheatmap(kmer, 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F)


ggsave("zc4688/honours/analyses/cn/figures/eighteen_heatmap.png")

######################### Eighteen R10
eighteen_upset <- eighteen_data %>% 
  left_join(chem) %>% 
  filter(chem == "R10") %>% 
  mutate(n = ifelse(n > 0, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
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

eighteen_upset <- eighteen_upset %>%
  mutate(Type = case_when(
    Type == "Shared" & Sum >=122 ~ "Common",
    TRUE ~ Type
  ))

colors <- setNames(brewer.pal(5, "Paired"), c("Common", "Shared", "Private", "Population-specific", "Widespread"))
upset(eighteen_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=colors)+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)



ggsave("zc4688/honours/analyses/cn/figures/eighteen_r10_upset.png")



#KL
kmer <- eighteen_data %>% 
  left_join(chem) %>% 
  filter(chem == "R10") %>% 
  mutate(n = n/genome) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- kmer / rowSums(kmer)

# Compute symmetric KL divergence matrix
kmer <- philentropy::distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)

rownames(kmer) <- rownames(prob_matrix)
colnames(kmer) <- rownames(prob_matrix)




annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code") 

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   chem = c("R10" = "deepskyblue1", "R9" = "darkorchid"))
pheatmap::pheatmap(kmer, 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F)


ggsave("zc4688/honours/analyses/cn/figures/eighteen_r10_heatmap.png")

#########################################
twoeight_data <- read_tsv("zc4688/honours/analyses/cn/result/twoeight/twoeight.kmers.tsv")
twoeight_data <- twoeight_data %>% left_join(cn) %>% filter(n > genome)
#eighteen_data <- bind_rows(eighteen_data)
#write_tsv(eighteen_data, "zc4688/honours/analyses/cn/result/eighteen_data.tsv")


twoeight_upset <- twoeight_data %>% 
  mutate(n = ifelse(n > 1, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name") %>% select(sampleid, Superpopulation.code)) %>% 
  group_by(Superpopulation.code, kmer) %>% 
  summarise(number = n()) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "Superpopulation.code", values_from = number, values_fill = 0)

twoeight_upset$Range <- rowSums(twoeight_upset[, 2:6] > 0)

twoeight_upset$Type<-ifelse(twoeight_upset$Range == 5, 'Shared',
                            ifelse(twoeight_upset$Range %in% c(2, 3, 4,5), 'Widespread', NA))

twoeight_upset$Sum <- rowSums(twoeight_upset[, 2:6], na.rm = TRUE)
na_indices <- is.na(twoeight_upset$Type)

twoeight_upset$Type[na_indices] <- ifelse(twoeight_upset$Sum[na_indices] == 1, 'Private',
                                          ifelse(twoeight_upset$Sum[na_indices] > 1, 'Population-specific', NA))

twoeight_upset <- twoeight_upset %>%
  mutate(Type = case_when(
    Type == "Shared" & Sum >=122 ~ "Common",
    TRUE ~ Type
  ))

colors <- setNames(brewer.pal(5, "Paired"), c("Common", "Shared", "Private", "Population-specific", "Widespread"))
upset(twoeight_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=colors)+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)



ggsave("zc4688/honours/analyses/cn/figures/twoeight_all_upset.png")


#KL
twoeight_js <- twoeight_data %>% 
  mutate(n = n/genome) %>% 
  left_join(chem) %>% 
  filter(chem == "R10") %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- twoeight_js / rowSums(twoeight_js)

# Compute symmetric KL divergence matrix
twoeight_js <- philentropy::distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)

rownames(twoeight_js) <- rownames(prob_matrix)
colnames(twoeight_js) <- rownames(prob_matrix)
#write_tsv(as.data.frame(twoeight_js), "zc4688/honours/analyses/cn/twoeight.js.tsv", col_names = T)


annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  left_join(chem %>% dplyr::rename("Sample.name" = "sampleid")) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code") 

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   chem = c("R10" = "deepskyblue1", "R9" = "darkorchid"))

pheatmap::pheatmap(kmer, 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F)



ggsave("zc4688/honours/analyses/cn/figures/twoeight_heatmap.png")
ggplot(eighteen_data %>% filter(n > genome*5) %>% mutate(n = n/genome), aes(x = n)) + geom_histogram() + facet_wrap(~sampleid)


################fiveeight
fiveeight_data <- read_tsv("zc4688/honours/analyses/cn/result/fiveeight/fiveeight.kmers.tsv")
fiveeight_data <- fiveeight_data %>% left_join(cn) %>% filter(n > genome)
#eighteen_data <- bind_rows(eighteen_data)
#write_tsv(eighteen_data, "zc4688/honours/analyses/cn/result/eighteen_data.tsv")



fiveeight_upset <- fiveeight_data %>% 
  mutate(n = ifelse(n >0, 1, 0)) %>% 
  select(kmer, n, sampleid) %>% 
  left_join(sampleinfo %>% dplyr::rename("sampleid" = "Sample.name") %>% select(sampleid, Superpopulation.code)) %>% 
  group_by(Superpopulation.code, kmer) %>% 
  summarise(number = n()) %>% 
  pivot_wider(id_cols = c("kmer"), names_from = "Superpopulation.code", values_from = number, values_fill = 0)

fiveeight_upset$Range <- rowSums(fiveeight_upset[, 2:6] > 0)

fiveeight_upset$Type<-ifelse(fiveeight_upset$Range == 5, 'Shared',
                            ifelse(fiveeight_upset$Range %in% c(2, 3, 4,5), 'Widespread', NA))

fiveeight_upset$Sum <- rowSums(fiveeight_upset[, 2:6], na.rm = TRUE)
na_indices <- is.na(fiveeight_upset$Type)

fiveeight_upset$Type[na_indices] <- ifelse(fiveeight_upset$Sum[na_indices] == 1, 'Private',
                                          ifelse(fiveeight_upset$Sum[na_indices] > 1, 'Population-specific', NA))

fiveeight_upset <- fiveeight_upset %>%
  mutate(Type = case_when(
    Type == "Shared" & Sum >=122 ~ "Common",
    TRUE ~ Type
  ))

colors <- setNames(brewer.pal(5, "Paired"), c("Common", "Shared", "Private", "Population-specific", "Widespread"))
upset(fiveeight_upset, sampleinfo %>% select(Superpopulation.code) %>% unique() %>% unlist(), 
      base_annotations=list(
        'Intersection Size'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=Type)
        )+scale_y_continuous()+scale_fill_manual(values=colors)+labs(fill=NULL)+theme(legend.position="bottom")
      ),
      width_ratio=0.1,
      themes=upset_default_themes(text=element_text(size=15)),
      set_sizes=FALSE
)



ggsave("zc4688/honours/analyses/cn/figures/fiveeight_all_upset.png")


#KL
fiveeight_js <- fiveeight_data %>% 
  mutate(n = n/genome) %>% 
  left_join(chem) %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- fiveeight_js / rowSums(fiveeight_js)

# Compute symmetric KL divergence matrix
#fiveeight_js <- philentropy::distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)

#rownames(fiveeight_js) <- rownames(prob_matrix)
#colnames(fiveeight_js) <- rownames(prob_matrix)
#write_tsv(as.data.frame(fiveeight_js), "zc4688/honours/analyses/cn/fiveeight.js.tsv", col_names = T)

kmer <- read_tsv("zc4688/honours/analyses/cn/fiveeight.js.tsv", col_names = T)


rownames(kmer) <- colnames(kmer)
annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code) %>% 
  left_join(chem %>% dplyr::rename("Sample.name" = "sampleid")) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code") 

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   chem = c("R10" = "deepskyblue1", "R9" = "darkorchid"))
pheatmap::pheatmap(fiveeight_js, 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = F,
                   show_colnames = F,
                   annotation_names_row = F,
                   annotation_names_col = F)






########## rRNAs

rrna <- bind_rows(eighteen_data %>% 
                    mutate(kmer = paste0(kmer, "_eighteen")), 
                  twoeight_data %>% 
                    mutate(kmer = paste0(kmer, "_twoeight")), 
                  fiveeight_data %>% 
                    mutate(kmer = paste0(kmer, "_fiveeight")))

rrna <- rrna %>% 
  select(sampleid, n, kmer, genome) %>% 
  mutate(n = n/genome) %>% 
  left_join(chem) %>% 
  filter(chem == "R10") %>% 
  select(kmer, n, sampleid) %>% 
  pivot_wider(id_cols = sampleid, names_from = kmer, values_from = n, values_fill = 0)  %>% 
  column_to_rownames("sampleid")
prob_matrix <- rrna / rowSums(rrna)

# Compute symmetric KL divergence matrix
rrna <- philentropy::distance(prob_matrix, method = "jensen-shannon", test.na = FALSE)

rownames(rrna) <- rownames(prob_matrix)
colnames(rrna) <- rownames(prob_matrix)
#write_tsv(as.data.frame(fiveeight_js), "zc4688/honours/analyses/cn/fiveeight.js.tsv", col_names = T)






annotations <- 
  sampleinfo %>% 
  filter(Sample.name %in% (cn_data %>% pull(sampleid) %>% unique())) %>% 
  select(Sample.name, Superpopulation.code, Population.code) %>% 
  column_to_rownames("Sample.name") %>% 
  dplyr::rename("Continent" = "Superpopulation.code")

ann_colors <- list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent)),
                   Population.code = setNames(c(
                     "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                     "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                     "#a55194", "#393b79", "#637939", "#8c6d31", "#843c39"), unique(annotations$Population.code)))


pheatmap::pheatmap(rrna, 
                   annotation_row = annotations,
                   annotation_col = annotations,
                   annotation_colors = ann_colors,
                   show_rownames = T,
                   show_colnames = F,
                   annotation_names_row = F,
                   fontsize = 5)


.######### compare graphs

#de bruijn graph?
tmp1 <- eighteen_data %>% filter(sampleid == "HG00116")
tmp2 <- eighteen_data %>% filter(sampleid == "HG01607")

k <- 19
tmp1 <- tmp1 %>%
  mutate(prefix = substr(kmer, 1, k-1),
         suffix = substr(kmer, 2, k))
edges <- tmp1 %>%
  group_by(prefix, suffix) %>%
  summarise(weight = sum(n), .groups = "drop")

tmp2 <- tmp2 %>%
  mutate(prefix = substr(kmer, 1, k-1),
         suffix = substr(kmer, 2, k))
edges2 <- tmp2 %>%
  group_by(prefix, suffix) %>%
  summarise(weight = sum(n), .groups = "drop")

# Build igraph
g1 <- graph_from_data_frame(d = edges, directed = TRUE)
E(g1)$weight <- edges$weight

g2 <- graph_from_data_frame(d = edges2, directed = TRUE)
E(g2)$weight <- edges2$weight

K <- CalculateWLKernel(list(g1, g2), h = 2)  # h = # iterations
distance <- sqrt(K[1,1] + K[2,2] - 2*K[1,2])





HG02134.primary.txt
