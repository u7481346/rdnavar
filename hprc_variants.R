

var_data <- list() 
var_files <- list.files("zc4688/honours/analyses/cn/hprc", pattern = "*var.txt$", full.names = TRUE)

for (file in var_files){
  sample_id <- str_replace(basename(file), ".var.txt", "")
  data <- read.delim(file, header = F) %>% 
    mutate(sampleid = sample_id)
  var_data[[file]] <- data
}


var_data <- bind_rows(var_data)

colnames(var_data) <- c("b", "a", "refstart", "refend", "c", "d", "ref", "alt", "morph", "morphstart", "morphend", "strand", "sampleid")


tmp <- var_data %>% mutate(type = case_when(
  refstart > 99 & refstart < 1969 ~ "18S",
  refstart >= 1969 & refstart <= 3045 ~ "ITS1",
  refstart > 3045 & refstart < 3200 ~ "5.8S",
  refstart >= 3200 & refstart <= 4374 ~ "ITS2",
  refstart > 4374 & refstart < 9441 ~ "28S",
  TRUE ~ "IGS"))%>% 
  mutate(category = ifelse(type %in% c("18S", "5.8S", "28S"), "rRNA", "other")) %>% 
  group_by(sampleid, morph, category) %>% summarise(n = n()) %>% 
  pivot_wider(names_from = category, values_from = n) %>% 
  mutate(ratio = rRNA/other)


tmp <- var_data %>% filter(refstart < 10000) %>% group_by(sampleid, refstart) %>% summarise(n = n())

total_units <- var_data %>%
  distinct(sampleid, morph) %>%
  dplyr::count(sampleid, name = "n_units")

tmp <- tmp %>% left_join(total_units) %>% 
  mutate(freq = n/n_units) %>%
  filter(n > 1) %>% 
  #mutate(variant_id = paste(refstart, sep = "_")) %>%
  mutate(variant_id = refstart) %>% 
  ungroup() %>% 
  select(sampleid, variant_id, freq) %>%
  pivot_wider(names_from = variant_id, values_from = freq, values_fill = 0)

mat <- as.data.frame(tmp)
rownames(mat) <- mat$sampleid
mat <- mat[,-1]

res.pca <- PCA(mat, scale.unit = TRUE, graph = FALSE)

library(DBI)
library(RSQLite)
# Path to your SQLite database
asmdb <- "/g/data/te53/t2t2024/db/asm.db"
# Connect to the database
con <- dbConnect(RSQLite::SQLite(), asmdb)
# Import the mapping table
mapping_df <- dbGetQuery(con, "SELECT * FROM seq")
# Close the connection
dbDisconnect(con)

mapping_df <- mapping_df %>% filter(str_detect(asmid, "GCA")) %>%  separate(asmid, into = c("sampleid", "no"), sep = "\\.") %>% select(sampleid, superpopulation) %>% distinct(sampleid, superpopulation)

mapping_df <- mapping_df %>%
  filter(sampleid %in% rownames(mat)) %>%
  arrange(match(sampleid, rownames(mat)))


# Plot PCA
fviz_pca_ind(res.pca, geom = "point", repel = TRUE,
             habillage = factor(mapping_df$superpopulation))



library(FactoMineR)




tmp <- var_data %>% filter(refstart > 99 & refstart < 1969 |
                             refstart > 3045 & refstart < 3200 |
                             refstart > 4374 & refstart < 9441) %>% 
  group_by(sampleid, refstart) %>% summarise(n = n())



tmp <- tmp %>% left_join(total_units) %>% 
  mutate(freq = n/n_units) %>%
  filter(n > 1) %>% 
  #mutate(variant_id = paste(refstart, sep = "_")) %>%
  mutate(variant_id = refstart) %>% 
  ungroup() %>% 
  select(sampleid, variant_id, freq) %>%
  pivot_wider(names_from = variant_id, values_from = freq, values_fill = 0)

mat <- as.data.frame(tmp)
rownames(mat) <- mat$sampleid
mat <- mat[,-1]

res.pca <- PCA(mat, scale.unit = TRUE, graph = FALSE)


# Plot PCA
fviz_pca_ind(res.pca, geom = "point", repel = TRUE,
             habillage = factor(mapping_df$superpopulation))


tmp <- var_data %>% 
  mutate(type = case_when(ref == "-" & nchar(alt)  > 5 ~ "Larger insertion (>=5 bp)",
                          ref == "-" & nchar(alt)  <= 5 ~ "Small insertion (<5 bp)",
                          alt == "-" & nchar(ref)  > 5 ~ "Larger deletion (>=5 bp)",
                          alt == "-" & nchar(ref)  <= 5 ~ "Small deletion (<5 bp)",
                          TRUE ~ "SNP"))


ggplot(tmp, aes(x = refstart, fill = type)) + geom_histogram(binwidth = 200) + theme_bw() +
  geom_rect(aes(xmin= 99, xmax = 1969, ymin = -600, ymax = -100), fill = "black") +
  geom_rect(aes(xmin =  3045, xmax = 3200, ymin = -600, ymax = -100), fill = "black") +
  geom_rect(aes(xmin = 4374, xmax = 9441, ymin = -600, ymax = -100), fill = "black") +
  scale_fill_brewer(palette = "Paired") +
  labs(y = "Count", x = "Position in rDNA unit", fill = "Variant type")
ggsave("zc4688/honours/analyses/cn/figures/hprc.pos.variants.png")


var_data %>%filter(refstart > 99 & refstart < 1969 |
                     refstart > 3045 & refstart < 3200 |
                     refstart > 4374 & refstart < 9441) %>% 
  group_by(sampleid, morph) %>% 
  summarise(n = n()) %>% 
  left_join(mapping_df) %>% 
  ggplot(aes(x = n, fill = superpopulation)) + geom_histogram(binwidth = 20)

var_data %>% group_by(sampleid, morph) %>% 
  summarise(n = n()) %>% 
  left_join(mapping_df) %>% 
  ggplot(aes(x = n, fill = superpopulation)) + geom_histogram(binwidth = 20) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  labs(x = "Number of variants", y = "Number of morphs", fill = "Continent")
ggsave("zc4688/honours/analyses/cn/figures/hprc.n.variants.png")

library(mclust)
tmp <- var_data %>% filter(str_detect(morph, "morph")) %>% 
  group_by(sampleid, morph) %>% 
  summarise(n = n())
model <- densityMclust(tmp$n)
summary(model)
plot(model, what = "density")
tmp$class <- model$classification


ggplot(tmp %>% left_join(mapping_df), aes(x = sampleid, y = n, color = factor(class))) +
  geom_jitter(width = 0.2, alpha = 0.7, size = 1) +
  facet_wrap(~ superpopulation, scales = "free_x") +
  theme_bw() +
  labs(color = "Variant class", y = "Variants per morph", x = "Haplotype") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


ggplot(tmp, aes(x = n, fill = class)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(x = "Variants per morph", y = "Density", fill = "Mclust class",
       title = "Variant count distribution (4-component mixture model)")



vars <- var_data %>% mutate(morph_id = paste0(sampleid, "_", morph)) %>% 
  mutate(var_id = paste0(refstart, "_", alt)) %>% 
  select(morph_id, var_id) %>%
  distinct() %>% 
  mutate(freq = 1) %>% 
  pivot_wider(names_from = var_id, values_from = freq, values_fill = 0)

mat <- as.data.frame(vars)
rownames(mat) <- mat$morph_id
mat <- mat[,-1]
mat <- mat[, colSums(mat) > nrow(mat) * 0.05]

res.pca <- PCA(mat, scale.unit = TRUE, graph = FALSE)
hab <- var_data %>% mutate(morph_id = paste0(sampleid, "_", morph)) %>% 
  select(morph_id, sampleid, morph) %>% 
  distinct() %>% 
  left_join(mapping_df) %>% 
  left_join(tmp)

fviz_pca_ind(res.pca, geom = "point", repel = TRUE, habillage = factor(hab$class))

mat_class <- mat %>% 
  mutate(morph_id = rownames(mat)) %>%
  left_join(hab %>% select(morph_id, class))

# Summarize by class: proportion of morphs with each variant
mat_class <- mat_class %>%
  group_by(class) %>%
  summarise(across(-morph_id, ~mean(.x))) %>%  # proportion with variant
  ungroup() %>% 
  filter(!is.na(class))
mat_class <- mat_class %>%
  pivot_longer(-class, names_to = "variant", values_to = "freq") %>%
  group_by(variant) %>%
  filter(max(freq) > 0.8 & min(freq) < 0.2) %>%
  ungroup()

mat_class %>% separate(variant, into = c("pos", "alt")) %>% filter(!is.na(class)) %>% ggplot(aes(x = as.numeric(pos), y = freq, color = class)) + geom_point() + coord_cartesian(xlim = c(0, 45000))

trf_files <- list.files("/g/data/te53/zc4688/honours/analyses/cn/hprc", 
                        pattern = ".trf.out", 
                        full.names = TRUE, 
                        recursive = TRUE)
trf <- list()
for (i in 1:length(trf_files)) {
  x <- read_table(trf_files[i], col_names = F)
  x <- x %>%
    mutate(`Sample ID` = str_replace(basename(trf_files[i]), ".trf.out", ""))
  if (nrow(x) > 0) {
    trf[[length(trf) + 1]] <- x  # Only append if not empty
  }
}

trf <- bind_rows(trf)
trf <- trf %>% mutate(morph_id = paste0(`Sample ID`, "_", X1))

trf <- trf %>% filter(!is.na(X10))

trf_ranges <- GRanges(seqnames = trf$morph_id, ranges = IRanges(start = trf$X2, end = trf$X3))
trf_ranges <- GenomicRanges::reduce(trf_ranges, min.gapwidth = 10)
trf_ranges <- data.frame(`Sample ID` = as.character(seqnames(trf_ranges)), start = start(trf_ranges), end = end(trf_ranges))
colnames(trf_ranges) <- c("morph_id", "start", "end")

#not just overall burden
trf_ranges %>% 
  group_by(morph_id) %>% 
  summarise(len = sum(end - start), .groups = "drop") %>% 
  filter(len < 50000)  %>% left_join(tmp) %>% view %>% 
  ggplot(aes(x = len, fill = factor(class))) + geom_density(alpha = 0.5) 



#some but not all class 4 have long stretch straight after 28S
#nothing obvious from the segment plots - could it be tandem repeats instead?


rm <- read.table("zc4688/honours/analyses/cn/hprc/GCA_018505835.rDNA.morphs.fasta.out.gff",header = F, sep = "\t")

rm <- bind_rows(rm) %>% 
  separate(V9, into = c("id", "motif", "class", "length"), sep = ";")  %>% 
  mutate(type = str_replace(motif, "Target Motif:", ""), class = str_replace(class, "RepeatType ClassFamily:", "")) %>% 
  dplyr::rename("start" = "V4", "end" = "V5", "strand" = "V7") %>% 
  filter(!(class == "rRNA" | str_detect(class, "Simple") | str_detect(class, "Satellite") | str_detect(class, "Low_complexity")))


#also not just te burden 
rm %>% 
  group_by(V1) %>% dplyr::rename("morph" = "V1") %>% 
  summarise(len = sum(end - start), .groups = "drop") %>% 
  filter(len < 50000)  %>% left_join(tmp %>% filter(sampleid == "GCA_018505835") %>% dplyr::rename("morphclass" = "class")) %>% view %>% 
  ggplot(aes(x = len, fill = factor(morphclass))) + geom_histogram()
pos <- var_data %>% group_by(refstart) %>% summarise(n = n()) %>% view
pos <- pos %>% filter(n > 3000) %>% pull(refstart)

var_data %>% filter(refstart %in% pos) %>% 
  mutate(type = case_when(ref == "-" & nchar(alt)  > 10 ~ "big insertion",
                          ref == "-" & nchar(alt)  <= 10 ~ "little insertion",
                          alt == "-" & nchar(ref)  > 10 ~ "big deletion",
                          alt == "-" & nchar(ref)  <= 10 ~ "little deletion",
                          TRUE ~ "SNP")) %>% 
  ggplot(aes(x = as.factor(refstart), fill = type)) + geom_bar() 

pos <- var_data %>% group_by(alt, refstart) %>% summarise(n = n()) %>% 
  mutate(freq = n/3981) %>% 
  ungroup() %>% 
  group_by(refstart) %>% 
  summarise(t = sum(freq)) %>% filter(t > 0.2) %>% 
  filter(refstart > 99 & refstart < 1969 |
           refstart > 3045 & refstart < 3200 |
           refstart > 4374 & refstart < 9441) %>% pull(refstart)


var_data %>% filter(refstart %in% pos) %>% 
  mutate(type = case_when(ref == "-" & nchar(alt)  > 10 ~ "big insertion",
                          ref == "-" & nchar(alt)  <= 10 ~ "little insertion",
                          alt == "-" & nchar(ref)  > 10 ~ "big deletion",
                          alt == "-" & nchar(ref)  <= 10 ~ "little deletion",
                          TRUE ~ "SNP")) %>% 
  group_by(alt, sampleid, refstart, type) %>% summarise(n = n()) %>% 
  left_join(total_units) %>% 
  mutate(freq = n/n_units) %>%
  ggplot(aes(x = as.factor(refstart), fill = type, y = freq)) + geom_bar(stat = "identity") 


var_data %>% filter(refstart %in% pos) %>% 
  mutate(type = case_when(ref == "-" & nchar(alt)  > 5 ~ "Larger insertion (>=5 bp)",
                          ref == "-" & nchar(alt)  <= 5 ~ "Small insertion (<5 bp)",
                          alt == "-" & nchar(ref)  > 5 ~ "Larger deletion (>=5 bp)",
                          alt == "-" & nchar(ref)  <= 5 ~ "Small deletion (<5 bp)",
                          TRUE ~ "SNP")) %>% 
  group_by(refstart, type) %>% summarise(n = n()) %>% 
  mutate(freq = n/3981) %>%
  ggplot(aes(x = as.factor(refstart), fill = type, y = freq)) + geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  labs(x = "Position in rDNA unit",
       y = "Frequency",
       fill = "Variant type")

ggsave("zc4688/honours/analyses/cn/figures/hprc.variants.png")
