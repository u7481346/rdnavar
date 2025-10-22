.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(Biostrings)
rdna_categories <- read.delim("zc4688/honours/analyses/cn/pogvit/rdna_categories.txt", header = F)
colnames(rdna_categories) <- c("pos", "kmer", "category")
all_kmers <- read.delim("zc4688/honours/analyses/cn/pogvit/kmers.fa", header = F)

all_kmers <- all_kmers %>% 
  mutate(type = str_replace(lag(V1), ">", "")) %>% 
  filter(str_detect(type, "kmer|genome|ref")) %>% 
  left_join(rdna_categories %>% dplyr::rename("type" = "kmer")) %>% 
  mutate(type = case_when(
    str_detect(type, "kmer") ~ "rdna", 
    TRUE ~ type)) %>% 
  dplyr::rename("seq" = "V1") %>% 
  mutate(rc = as.character(reverseComplement(DNAStringSet(seq)))) %>% 
  mutate(category = ifelse(str_detect(category, "ITS"), "ITS", category))



pacbio <- read.delim("zc4688/honours/analyses/cn/pogvit/POGVIT_pb.cn.kmers.fa", header = F) %>% 
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
    mutate(tech = "PB")

pacbio %>%
  group_by(type, tech) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = depth) %>%
  mutate(cn = rdna / genome)


illumina <- read.delim("zc4688/honours/analyses/cn/pogvit/POGVIT_ill.cn.kmers.fa", header = F) %>% 
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
  mutate(tech = "Illumina")

illumina %>%
  group_by(type, tech) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = depth) %>%
  mutate(cn = rdna / genome)

ont <- read.delim("zc4688/honours/analyses/cn/pogvit/POGVIT_ont.cn.kmers.fa", header = F) %>% 
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
  mutate(tech = "ONT")

ont %>%
  group_by(type, tech) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = depth) %>%
  mutate(cn = rdna / genome)

all <- bind_rows(pacbio, illumina, ont)

all %>%
  group_by(type, tech) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = type, values_from = depth) %>%
  mutate(cn = case_when(
    tech == "Illumina" ~ rdna/42,
    tech == "ONT" ~ rdna/33,
    TRUE ~ rdna/39
  ))


all %>% filter(n < 500) %>% ggplot(aes(x = n, fill = type)) + geom_density() + facet_wrap(~tech)

all %>% filter(n < 500) %>% mutate(d = case_when(tech == "Illumina" ~ 42, tech == "ONT" ~ 33, TRUE ~ 39)) %>% ggplot() + geom_density(aes(x = n, fill = type)) + geom_vline(aes(xintercept = d)) + facet_wrap(~tech)


all %>% filter(n < 1000000) %>% mutate(d = case_when(tech == "Illumina" ~ 42, tech == "ONT" ~ 33, TRUE ~ 39)) %>% ggplot() + geom_density(aes(x = n, fill = ifelse(type != "rdna", type, category)), alpha = 0.5) + geom_vline(aes(xintercept = d)) + facet_wrap(~tech, scales="free_x") + theme_bw()
all %>% 
  filter(n < 1000000 & type != "genome") %>%
  mutate(n = case_when(tech == "Illumina" ~ n/42, tech == "ONT" ~ n/33, TRUE ~ n/39)) %>% 
  ggplot() + geom_density(aes(x = n, y =..count.., fill = category), alpha = 0.5, stat = "density") + 
  #geom_vline(xintercept = )
  theme_bw() +facet_wrap(~tech) + 
  labs(x = "Copy number", y = "Frequency", fill = "Category") + scale_fill_brewer(palette = "Set1")

ggsave("zc4688/honours/analyses/cn/figures/pogvit.density.png")
cn <- all %>%
  mutate(test = ifelse(is.na(category), type, category)) %>%
  group_by(test, tech) %>%
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA` / genome)
cn %>% 
  select(tech, `18S_rRNA`, `28S_rRNA`, `5_8S_rRNA`, ITS, IGS,  genome) %>% 
  pivot_longer(cols = c(2:6), names_to = "Variable", values_to = "value") %>% 
  mutate(n = case_when(
    tech == "Illumina" ~ value/42,
    tech == "ONT" ~ value/33,
    TRUE ~ value/39
    )) %>% view %>% 
  ggplot(aes(x = tech, color = Variable, y = n)) + geom_point()


pacbio <- pacbio %>% filter(type == "rdna") %>% dplyr::rename("PacBio" = "n") %>% mutate(PacBio = PacBio/39)%>% select(kmer, category, pos, PacBio)
illumina <- illumina %>% filter(type == "rdna") %>% dplyr::rename("Illumina" = "n") %>% mutate(Illumina = Illumina/42)%>% select(kmer, category, pos, Illumina)
ont <- ont %>% filter(type == "rdna") %>% dplyr::rename("ONT" = "n") %>% mutate(ONT = ONT/33)%>% select(kmer, category, pos, ONT)

tmp <- left_join(pacbio, illumina)
tmp <- left_join(tmp, ont)

plot(log1p(tmp$n_pb), log1p(tmp$n_ill)); cor(tmp$n_pb, tmp$n_ill, method="spearman")
[1] 0.08839049

plot(log1p(tmp$n_pb), log1p(tmp$n_ont)); cor(tmp$n_pb, tmp$n_ont, method="spearman")
[1] 0.4872049

plot(log1p(tmp$n_ill), log1p(tmp$n_ont)); cor(tmp$n_ill, tmp$n_ont, method="spearman")
[1] 0.1569039


ggplot(tmp, aes(x = ONT, y = Illumina)) + geom_point() + theme_bw() 
ggsave("zc4688/honours/analyses/cn/figures/pogvit.ont.ill.png")
ggplot(tmp, aes(x = PacBio, y = Illumina)) + geom_point() + theme_bw()
ggsave("zc4688/honours/analyses/cn/figures/pogvit.pb.ill.png")
ggplot(tmp, aes(x = n_pb, y = n_ont)) + geom_point() + theme_bw()
ggsave("zc4688/honours/analyses/cn/figures/pogvit.pb.ont.png")

tmp <- tmp %>% select(PacBio, Illumina, ONT)
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cornflowerblue", ...)
}
pairs(tmp, diag.panel = panel.hist, upper.panel = NULL)

#gc
test <- read.delim("zc4688/honours/analyses/cn/pogvit/gc.txt", header = F)
test %>% separate(V1, into = c("a", "b", "c", "d","e", "start", "end")) %>%
  ggplot(aes(x = as.numeric(start), y = as.numeric(V2))) + geom_line()



#hg002
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


pacbio <- read.delim("zc4688/honours/analyses/cn/hg002/hg002_pb.cn.kmers.fa", header = F) %>% 
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
  mutate(tech = "PB")

ont <- read.delim("zc4688/honours/analyses/cn/hg002/hg002_ont.cn.kmers.fa", header = F) %>% 
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
  mutate(tech = "ONT")

illumina <- read.delim("zc4688/honours/analyses/cn/hg002/hg002_ill.cn.kmers.fa", header = F) %>% 
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
  mutate(tech = "Illumina")

all <- bind_rows(pacbio, illumina, ont)

all %>%mutate(test = ifelse(is.na(category), type, category)) %>%
  group_by(test, tech) %>%  
  summarise(depth = median(n), .groups = 'drop') %>%
  pivot_wider(names_from = test, values_from = depth) %>%
  mutate(cn = `18S_rRNA`/genome)


all %>% 
  filter(n < 1000000 & type == "rdna") %>%
  mutate(n = case_when(tech == "Illumina" ~ n/28, tech == "ONT" ~ n/18, TRUE ~ n/19)) %>% 
  ggplot() + geom_density(aes(x = n, y =after_stat(count), fill = category), alpha = 0.5, stat = "density") + 
  theme_bw() +facet_wrap(~tech) + 
  labs(x = "Copy number", y = "Frequency", fill = "Category") + scale_fill_brewer(palette = "Set1")

all %>% 
  filter(n < 1000000 & type == "rdna") %>%
  mutate(n = case_when(tech == "Illumina" ~ n/28, tech == "ONT" ~ n/18, TRUE ~ n/19)) %>% 
  ggplot() + geom_density(aes(x = n, y =after_stat(count), fill = category), alpha = 0.5, stat = "density") + 
  theme_bw() +facet_wrap(~tech) + 
  labs(x = "Copy number", y = "Density", fill = "Category") + scale_fill_brewer(palette = "Set1")
ggsave("zc4688/honours/analyses/cn/figures/hg002.rdna.png")

all %>% 
  filter(n < 10000 & type != "ntsm") %>%
  ggplot() + geom_density(aes(x = log(n), fill = type), alpha = 0.5, stat = "density") + 
  theme_bw() +facet_wrap(~tech) + 
  labs(x = "Log(count)", y = "Frequency", fill = "Category")

ggsave("zc4688/honours/analyses/cn/figures/hg002.cn.png")

all %>% 
  filter(n < 100000 & type == "rdna") %>%
  mutate(n = case_when(tech == "Illumina" ~ n/28, tech == "ONT" ~ n/18, TRUE ~ n/19)) %>% 
  ggplot() + geom_boxplot(aes(y = n, x =tech, fill = category)) + 
  theme_bw()  + 
  labs(y = "Copy number", fill = "Category", x = NULL) + scale_fill_brewer(palette = "Set1")
ggsave("zc4688/honours/analyses/cn/figures/hg002.boxplot.png")



pacbio <- pacbio %>% filter(type == "rdna") %>% dplyr::rename("PacBio" = "n") %>% mutate(PacBio = PacBio/19)%>% select(kmer, category, pos, PacBio)
illumina <- illumina %>% filter(type == "rdna") %>% dplyr::rename("Illumina" = "n") %>% mutate(Illumina = Illumina/28)%>% select(kmer, category, pos, Illumina)
ont <- ont %>% filter(type == "rdna") %>% dplyr::rename("ONT" = "n") %>% mutate(ONT = ONT/18)%>% select(kmer, category, pos, ONT)

tmp <- left_join(pacbio, illumina)
tmp <- left_join(tmp, ont)



ggplot(tmp, aes(x = n_ont, y = n_ill, color = category)) + geom_point() + theme_bw()  + scale_color_brewer(palette = "Set1")
ggsave("zc4688/honours/analyses/cn/figures/hg002.ont.ill.png")
ggplot(tmp, aes(x = n_pb, y = n_ill, color = category)) + geom_point() + theme_bw()+ scale_color_brewer(palette = "Set1")
ggsave("zc4688/honours/analyses/cn/figures/hg002.pb.ill.png")
ggplot(tmp, aes(x = n_pb, y = n_ont, color = category)) + geom_point() + theme_bw()+ scale_color_brewer(palette = "Set1")

ggsave("zc4688/honours/analyses/cn/figures/hg002.pb.ont.png")

plot(log1p(tmp$n_pb), log1p(tmp$n_ill)); cor(tmp$n_pb, tmp$n_ill, method="spearman")
[1] 0.4651369

plot(log1p(tmp$n_pb), log1p(tmp$n_ont)); cor(tmp$n_pb, tmp$n_ont, method="spearman")
[1]0.4097751
plot(log1p(tmp$n_ill), log1p(tmp$n_ont)); cor(tmp$n_ill, tmp$n_ont, method="spearman")
[1] 0.2046107

tmp <- tmp %>% select(PacBio, Illumina, ONT)
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cornflowerblue", ...)
}
pairs(tmp, diag.panel = panel.hist, upper.panel = NULL)




################## flow cell comparison
libdata <- list() 
libfiles <- list.files("zc4688/honours/analyses/cn/illumina_libs/flowcell/", pattern = "HG002.*depth.txt$", full.names = TRUE)

for (file in libfiles){
  libname <- str_replace(basename(file), ".depth.txt", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(V1 = ifelse(str_detect(V1, "random"), paste0("random", row_number()),
                       V1)) %>% 
    column_to_rownames("V1") %>% 
    dplyr::rename(!!libname := "V2")
  libdata[[file]] <- data
}

libdata <- bind_cols(libdata)

libdata <- libdata %>% rownames_to_column("test") %>% 
  pivot_longer(-test, names_to = "lib") %>% 
  pivot_wider(names_from = test, values_from = value) %>% 
  mutate(across(-c(lib, `50kb`), ~ .x / `50kb`))  


variances <- libdata %>%
  summarise(across(where(is.numeric), var))

variances %>% select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% pivot_longer(everything()) %>% ggplot(aes(x = value)) + geom_boxplot()

libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% view %>% 
  ggplot(aes(x = cv)) + 
  geom_histogram(binwidth=0.005, fill= "cornflowerblue") + 
  geom_density() + 
  geom_vline(xintercept = 0.04240057, linetype = "dashed") + 
  theme_bw() +
  labs(x = "CV", y = "Count")

ggsave("zc4688/honours/analyses/cn/figures/hg002.flowcell.cv.png")

tmp <- libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% 
  pull(cv)
t.test(tmp, mu = 0.04240057, alternative = "less")
#p-value < 2.2e-16

tmp <- libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`, lib))
correlations <- cor(tmp$rdna, tmp[ , setdiff(names(tmp), "rdna")], use = "pairwise.complete.obs")

################## hg001
hg001 <- list() 
libfiles <- list.files("zc4688/honours/analyses/cn/illumina_libs/flowcell", pattern = "HG001.*depth.txt$", full.names = TRUE)
libfiles <- libfiles[libfiles != "zc4688/honours/analyses/cn/illumina_libs/flowcell/HG001_1312.depth.txt"]
for (file in libfiles){
  libname <- str_replace(basename(file), ".depth.txt", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(V1 = ifelse(str_detect(V1, "random"), paste0("random", row_number()),
                       V1)) %>% 
    column_to_rownames("V1") %>% 
    dplyr::rename(!!libname := "V2")
  hg001[[file]] <- data
}


hg001 <- bind_cols(hg001)

hg001 <- hg001 %>% rownames_to_column("test") %>% 
  pivot_longer(-test, names_to = "lib") %>% 
  pivot_wider(names_from = test, values_from = value) %>% 
  mutate(across(-c(lib, `50kb`), ~ .x / `50kb`))  


variances <- hg001 %>%
  summarise(across(where(is.numeric), var))

variances %>% select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% pivot_longer(everything()) %>% ggplot(aes(x = value)) + geom_boxplot()

hg001 %>% 
  select(-c(`1kb`, `10kb`, `50kb`)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% 
  view %>% 
  ggplot(aes(x = cv)) + 
  geom_histogram(binwidth=0.01, fill= "cornflowerblue") + 
  geom_density() +
  geom_vline(xintercept = 0.04283992,linetype = "dashed") +
  theme_bw() + coord_cartesian(xlim = c(0,0.5))+
  labs(x = "CV", y = "Count")

ggsave("zc4688/honours/analyses/cn/figures/hg001.flowcell.cv.png")
tmp <- hg001 %>% 
  select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% 
  pull(cv)
#p-value = 0.7125


t.test(tmp, mu = 0.04283992, alternative = "less")
ggplot() + geom_density(data = libdata, aes(x = rdna, y = after_stat(count), fill = "HG002"), alpha = 0.8) +
  geom_density(data = hg001, aes(x = rdna, y = after_stat(count), fill = "HG001"), alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg002.hg001.flowcell.cn.png")

ggplot() + geom_density(data = libdata, aes(x = rdna, y = after_stat(count)), fill = "#E41A1C", alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg002.flowcell.cn.png")

ggplot()  +
  geom_density(data = hg001, aes(x = rdna, y = after_stat(count)), fill = "#377EB8", alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg001.flowcell.cn.png")


hg001 <- hg001 %>% 
  select(-c(`1kb`, `10kb`, `50kb`, lib))
hg001 <- cor(hg001$rdna, hg001[ , setdiff(names(hg001), "rdna")], use = "pairwise.complete.obs") %>% as.data.frame() %>%  pivot_longer(everything())

hg001 %>% dplyr::rename("hg001" = "value") %>% left_join(correlations %>% as.data.frame() %>% pivot_longer(everything(),values_to = "hg002"))  %>% 
  view %>% 
  ggplot(aes(x = hg001, y = hg002)) + geom_point()




################## library comparison
libdata <- list() 
libfiles <- list.files("zc4688/honours/analyses/cn/illumina_libs/library/", pattern = "*depth.txt$", full.names = TRUE)

for (file in libfiles){
  libname <- str_replace(basename(file), ".depth.txt", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(V1 = ifelse(str_detect(V1, "random"), paste0("random", row_number()),
                       V1)) %>% 
    column_to_rownames("V1") %>% 
    dplyr::rename(!!libname := "V2")
  libdata[[file]] <- data
}

libdata <- bind_cols(libdata)

libdata <- libdata %>% rownames_to_column("test") %>% 
  pivot_longer(-test, names_to = "lib") %>% 
  pivot_wider(names_from = test, values_from = value) %>% 
  mutate(across(-c(lib, `50kb`), ~ .x / `50kb`))  


variances <- libdata %>%
  summarise(across(where(is.numeric), var))

variances %>% select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% pivot_longer(everything()) %>% ggplot(aes(x = value)) + geom_boxplot()

libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% view %>% 
  ggplot(aes(x = cv)) + 
  geom_histogram(binwidth=0.005, fill= "cornflowerblue") + 
  geom_density() + 
  geom_vline(xintercept = 0.012615513, linetype = "dashed") + 
  theme_bw() +
  labs(x = "CV", y = "Count") + coord_cartesian(xlim = c(0,0.5))

ggsave("zc4688/honours/analyses/cn/figures/hg002.lib.cv.png")
tmp <- libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% 
  pull(cv)
t.test(tmp, mu = 0.01261551, alternative = "less")
tmp <- libdata %>% 
  select(-c(`1kb`, `10kb`, `50kb`, lib))
correlations <- cor(tmp$rdna, tmp[ , setdiff(names(tmp), "rdna")], use = "pairwise.complete.obs")


#compare same dna vs different dna

################## hg001
hg001 <- list() 
libfiles <- list.files("zc4688/honours/analyses/cn/illumina_libs/library/HG001/", pattern = "*depth.txt$", full.names = TRUE)

for (file in libfiles){
  libname <- str_replace(basename(file), ".depth.txt", "")
  data <- read.delim(file, header = F)
  data <- data %>% 
    mutate(V1 = ifelse(str_detect(V1, "random"), paste0("random", row_number()),
                       V1)) %>% 
    column_to_rownames("V1") %>% 
    dplyr::rename(!!libname := "V2")
  hg001[[file]] <- data
}


hg001 <- bind_cols(hg001)

hg001 <- hg001 %>% rownames_to_column("test") %>% 
  pivot_longer(-test, names_to = "lib") %>% 
  pivot_wider(names_from = test, values_from = value) %>% 
  mutate(across(-c(lib, `50kb`), ~ .x / `50kb`))  


variances <- hg001 %>%
  summarise(across(where(is.numeric), var))

variances %>% select(-c(`1kb`, `10kb`, `50kb`, rdna)) %>% pivot_longer(everything()) %>% ggplot(aes(x = value)) + geom_boxplot()

hg001 %>% 
  select(-c(`1kb`, `10kb`, `50kb`)) %>% 
  pivot_longer(-lib) %>% 
  group_by(name) %>% 
  summarise(mean = mean(value), SD = sd(value), cv = SD/mean) %>% 
  view %>% 
  ggplot(aes(x = cv)) + 
  geom_histogram(binwidth=0.01, fill= "cornflowerblue") + 
  geom_density() +
  geom_vline(xintercept = 0.01968755,linetype = "dashed") +
  theme_bw() + coord_cartesian(xlim = c(0,0.5))+
  labs(x = "CV", y = "Count")

ggsave("zc4688/honours/analyses/cn/figures/hg001.cv.png")

ggplot() + geom_density(data = libdata, aes(x = rdna, y = after_stat(count), fill = "HG002"), alpha = 0.8) +
  geom_density(data = hg001, aes(x = rdna, y = after_stat(count), fill = "HG001"), alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg002.hg001.lib.cn.png")

ggplot() + geom_density(data = libdata, aes(x = rdna, y = after_stat(count)), fill = "#E41A1C", alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg002.lib.cn.png")

ggplot()  +
  geom_density(data = hg001, aes(x = rdna, y = after_stat(count)), fill = "#377EB8", alpha = 0.8) +
  theme_bw() +
  labs(x = "Copy number", y = "Frequency", fill = "Sample") + scale_fill_brewer(palette = "Set1")


ggsave("zc4688/honours/analyses/cn/figures/hg001.lib.cn.png")


hg001 <- hg001 %>% 
  select(-c(`1kb`, `10kb`, `50kb`, lib))
hg001 <- cor(hg001$rdna, hg001[ , setdiff(names(hg001), "rdna")], use = "pairwise.complete.obs") %>% as.data.frame() %>%  pivot_longer(everything())

hg001 %>% dplyr::rename("hg001" = "value") %>% left_join(correlations %>% as.data.frame() %>% pivot_longer(everything(),values_to = "hg002"))  %>% 
  view %>% 
  ggplot(aes(x = hg001, y = hg002)) + geom_point()

#######gc
tmp <- read.delim("zc4688/honours/analyses/cn/chm13.refmorph.gc.txt", header = F)
test <- read.delim("zc4688/honours/analyses/cn/pogvit/gc.txt", header = F)
tmp <- tmp %>% 
  separate(V1, into = c("a", "b", "c", "d", "start", "end")) %>% 
  select(start, V2) %>% mutate(s = "CHM13")
test <- test %>% 
  separate(V1, into = c("a", "b", "c", "d", "e", "start", "end")) %>% 
  select(start, V2) %>% mutate(s = "P. Vitticeps")


bind_rows(tmp, test) %>% 
  ggplot(aes(x = as.numeric(start), y = V2)) + geom_line() + facet_wrap(~s) + theme_bw() +
  labs(x = "Position in rDNA unit", y = "GC content (%)")
ggsave("zc4688/honours/analyses/cn/figures/chm13.pogvit.gc.png")

