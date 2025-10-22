.libPaths(c("/g/data/if89/apps/Rlib/4.3.1/", .libPaths()))

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(DBI)
library(RSQLite)


###################### trna
trnadata <- list() 
trnafiles <- list.files("zc4688/honours/analyses/cn/trna/", pattern = "*results.txt$", full.names = TRUE)

for (file in trnafiles){
  if (file.info(file)$size > 0) {
    sampleid <- str_replace(basename(file), ".results.txt", "")
    
    data <- read.delim(file, header = F, skip=3)
    data <- data %>% 
      mutate(sampleid = sampleid)
    trnadata[[file]] <- data
  }
}

trnadata <- bind_rows(trnadata)
colnames(trnadata) <- c("contig", "trna", "start", "end", "AA", "codon", "istart", "iend", "score", "pseudo", "sampleid")
trnadata <- trnadata %>% mutate(name = paste0(AA, "_", codon))
###################### metadata
population <- read.delim("t2t2024/referenceresource/ethnic_info/ethnic_info_cont_151124.txt", sep = "\t", header = FALSE)

population <- population %>% 
  dplyr::rename(Population = V4, Continent = V5, donorID = V2) %>% 
  mutate(V1 = str_replace_all(V1, "\\.1", ""), V1 = str_replace_all(V1, "\\.2", "")) %>% 
  select(V1, Continent) %>% 
  filter(V1 %in% trnadata$sampleid) %>% 
  column_to_rownames("V1")

# Path to your SQLite database
asmdb <- "/g/data/te53/t2t2024/db/asm.db"

# Connect to the database
con <- dbConnect(RSQLite::SQLite(), asmdb)

# Import the mapping table
mapping_df <- dbGetQuery(con, "SELECT * FROM seq")

# Close the connection
dbDisconnect(con)

#####################

trnadata <- trnadata %>% 
  dplyr::rename("seqid" = "contig") %>% 
  mutate(seqid = str_replace_all(seqid, " ", "")) %>% 
  left_join(mapping_df)


trnadata %>% 
  group_by(sampleid) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = n)) + geom_density(aes(y = after_stat(count))) +
  theme_bw() +
  labs(x = "Number of tRNAs", y = "Count")
ggsave("zc4688/honours/analyses/cn/figures/n_trnas.png")

trnadata %>% 
  group_by(sampleid) %>% 
  filter(score > 50) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = n)) + geom_density(aes(y = after_stat(count))) +
  theme_bw() +
  labs(x = "Number of high confidence tRNAs", y = "Count")
ggsave("zc4688/honours/analyses/cn/figures/n_hc_trnas.png")

trnadata %>% 
  group_by(sampleid, name) %>% 
  filter(score > 50) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = n)) + geom_histogram(binwidth = 1) +
  theme_bw() +
  labs(x = "tRNA CN", y = "Count")
ggsave("zc4688/honours/analyses/cn/figures/trna_cn.png")


trnadata %>% 
  ggplot(aes(x = score)) + geom_density(aes(y = after_stat(count))) +
  theme_bw() +
  labs(x = "tRNA score (bits)", y = "Count")
ggsave("zc4688/honours/analyses/cn/figures/trna_score.png")


annotations <- trnadata %>% 
  select(sampleid, superpopulation) %>% 
  distinct() %>% 
  filter(!is.na(superpopulation)) %>% 
  column_to_rownames("sampleid")


tmp <- trnadata %>% 
  filter(score > 50) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, AA) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = AA, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid")


pheatmap(tmp, cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$superpopulation))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  filter(score < 50) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, AA) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = AA, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid")


pheatmap(tmp, cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  filter(score > 50) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, name) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = name, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid")


pheatmap(tmp, cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  filter(score > 50) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, codon) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = codon, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid")


pheatmap(tmp, cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(Continent = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$Continent))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  filter(score > 50) %>% 
  mutate(name = paste0(name, "_", communitylabel)) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, name) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = name, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid") %>% 
  select(where(~ n_distinct(.) > 1))


pheatmap(tmp,  cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(superpopulation = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$superpopulation))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  filter(score > 50) %>% 
  #filter(!str_detect(sampleid, "GCA_0188526")) %>% 
  group_by(sampleid, communitylabel) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = communitylabel, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid") 


pheatmap(tmp,  cluster_cols = F, annotation_row = annotations,
         annotation_colors = list(superpopulation = setNames(brewer.pal(n = 5, name = "Set2"), unique(annotations$superpopulation))),
         show_rownames = F,
         annotation_names_row = F)

tmp <- trnadata %>% 
  group_by(sampleid, trna) %>% 
  summarise(n = n()) %>% 
  pivot_wider(names_from = trna, values_from = n, values_fill = 0) %>% 
  column_to_rownames("sampleid")

trnadata %>% group_by(as.character(trna), sampleid) %>% summarise(n = n()) %>% 
  group_by(`as.character(trna)`) %>%
  summarise(mean = mean(n), SD = sd(n), cv = SD/mean) %>% 
  ggplot() + geom_bar(aes(x = `as.character(trna)`, y = mean))


trnadata %>% group_by(as.character(trna), sampleid) %>% summarise(n = n()) %>% 
  group_by(`as.character(trna)`) %>%
  summarise(mean = mean(n), SD = sd(n), cv = SD/mean) %>% 
  ggplot() + geom_bar(aes(x = `as.character(trna)`, y = cv), stat = "")
