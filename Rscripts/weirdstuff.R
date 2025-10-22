

# Mixes of unit sizes
#Often there is a lot of variation in unit sizes in a single genome. In this example, 16599, there are 35 morphs. min size and median size are both around 16000 but Q3 is 241836 (ie there are quite a few massive morphs)
#closer inspection shows this seem to be the case (not a filtering error or somethig):
tmp <- read_table("zc4688/honours/analyses/chordates/new/GCA_039720435.primary.txt", comment = "#", col_names=F)
tmp %>% select(X1, X3, 5, X6, X9, X10) %>% filter(X1 == "JBCFXG010000032.1") %>% filter(X3 == "18S_rRNA") %>% arrange(X9) %>% mutate(size = X9 - lag(X10)) %>% view

#23 of the morphs are very close in size at 16500. the remaining 12 range between 190000 and 490000 (can have both sizes on single contig)
#inspection of these longer 'units' shows these seem to be old fragments of units with veryt very long IGS. some are right next to a 'true' array. almost all are missing the middle 1500 bp of the 28S.

#in others it seems to be due to presence of fragments of old units in the IGS - eg GCA_016271365
tmp <- read_table("zc4688/honours/analyses/chordates/new/GCA_016271365.primary.txt", comment = "#", col_names=F)
tmp %>% select(X1, X3, 5, X6, X9, X10) %>% filter(X1 == "CM057456.1") %>%  view

#other times the units seem to just be actually massive:
tmp <- read_table("zc4688/honours/analyses/chordates/new/GCA_964289735.primary.txt", comment = "#", col_names=F)
tmp %>% select(X1, X3, 5, X6, X9, X10) %>% filter(X1 == "CAXYMO010000414.1") %>% arrange(X9) %>%  mutate(gap = X9 - lag(X10)) %>% view
#same with 
tmp <- read_table("zc4688/honours/analyses/chordates/new/GCA_965151615.primary.txt", comment = "#", col_names=F)
tmp %>% select(X1, X3, 5, X6, X9, X10) %>% filter(X1 == "CBCUAN010000271.1") %>% arrange(X9) %>%  mutate(gap = X9 - lag(X10)) %>% view

#question is - some of these large ones are likely to be errors, but some of them seem to be real - how to distinguish? how to validate?



############################## GC content
gcdirectory <- "/g/data/te53/zc4688/honours/analyses/chordates/new/gc_content"

gc_files <- list.files(gcdirectory, pattern = "\\gc.txt$", full.names = TRUE)


# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in gc_files) {
  data <- read.delim(file, header = F, sep="\t")
  samplename <- str_replace(basename(file), ".gc.txt", "")
  data <- data %>% 
    mutate(`Sample ID` = samplename) %>% 
    mutate(startpos = as.numeric(str_extract(V1, "(?<=_sliding:)\\d+"))) %>% 
    filter(startpos < 15000) %>% 
    select(-V1)
  data_list[[file]] <- data
}

gc_content <- bind_rows(data_list)
gc_content %>% left_join(structure %>% 
                           group_by(`Sample ID`, Type) %>% 
                           filter((End - Start) == max(End - Start)) %>% 
                           arrange(`Sample ID`, Type, Start) %>% 
                           slice_head(n = 1) %>% filter(Type == "28S_rRNA") %>% 
                           select(`Sample ID`, Start, End)) %>% 
  filter(startpos > Start & startpos < End) %>% 
  group_by(`Sample ID`) %>% 
  summarise(gc = median(V2)) %>% 
  left_join(classifications) %>% 
  filter(group %in% classtree$tip.label) %>% 
  ggplot() + geom_boxplot(aes(x = group, y = gc, fill = group), notch = T) + scale_fill_manual(values = colours) + theme_bw() + labs(x = element_blank(), y = "GC content (%)", fill = "Taxonomic group") + theme(legend.position = "none")

#mammals have highest 28S gc content - roughly corresponds to legnth
tmp <- gc_content %>% left_join(structure %>% 
                           group_by(`Sample ID`, Type) %>% 
                           filter((End - Start) == max(End - Start)) %>% 
                           arrange(`Sample ID`, Type, Start) %>% 
                           slice_head(n = 1) %>% filter(Type == "28S_rRNA") %>% 
                           select(`Sample ID`, Start, End)) %>% 
  filter(startpos > Start & startpos < End) %>% 
  group_by(`Sample ID`) %>% 
  summarise(gc = median(V2), length = median(End - Start)) %>% 
  left_join(classifications) %>% 
   mutate(prop = gc/length)

ggplot(tmp) + geom_point(aes(x = length, y = gc, color = group))

gc_content %>% left_join(structure %>% 
                           group_by(`Sample ID`, Type) %>% 
                           filter((End - Start) == max(End - Start)) %>% 
                           arrange(`Sample ID`, Type, Start) %>% 
                           slice_head(n = 1) %>% filter(Type == "28S_rRNA") %>% 
                           select(`Sample ID`, Start, End)) %>% 
  filter(startpos > Start & startpos < End) %>% left_join(classifications) %>% 
  filter(group == "Mammalia") %>% filter(`Sample ID` == "hg002") %>% 
  ggplot(aes(x = startpos, y = V2)) + geom_line()

tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/twoeight_gc.txt", header = F)

tmp <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  mutate(start = str_extract(V1, "(?<=sliding:)\\d+")) %>% 
  left_join(classifications)

tmp <- tmp %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())
  

p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p + new_scale_color() + geom_facet(
  data = tmp,
  mapping = aes(x = as.numeric(start), color = V2),
  geom = geom_point,       # Horizontal bar chart
  panel = "Position in 288S"
) +
  scale_y_discrete()+scale_color_viridis_b(option = "viridis", n.breaks=10) + labs(color = "GC content (%)")+ coord_cartesian(clip = 'off')+ xlim_expand(c(0, 25), "Tree")


ggsave("zc4688/honours/analyses/chordates/figures/twoeight_gc.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/twoeight_gc.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = n, fill = group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colours) + 
  theme_bw() + 
  labs(x = NULL, y = "Number of windows with >80% GC content", fill = "Taxonomic Group") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 

ggsave("zc4688/honours/analyses/chordates/figures/gc.length.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.length.png", height = 6, width = 10)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(n ~ group, p.adj.method = "BH") %>% view

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = max, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "Maximum GC % in 28S", fill = "Taxonomic Group")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 
ggsave("zc4688/honours/analyses/chordates/figures/gc.max.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.max.png", height = 6, width = 10)

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(max ~ group, p.adj.method = "BH") %>% view


classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  left_join(
    structuremetadata
  ) %>%
  mutate(remain = `28S_length` - n*50) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = remain, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "28S length without GC rich regions", fill = "Taxonomic Group") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6)) 
ggsave("zc4688/honours/analyses/chordates/figures/gc.removed.pdf", height = 6, width = 10)
ggsave("zc4688/honours/analyses/chordates/figures/gc.removed.png", height = 6, width = 10)


tmp <- read.delim("zc4688/honours/analyses/chordates/motif/its1_gc.txt", header = F)

tmp <- tmp %>% 
  mutate(`Sample ID` = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  mutate(start = str_extract(V1, "(?<=sliding:)\\d+")) %>% 
  left_join(classifications)

tmp <- tmp %>% 
  dplyr::rename("label" = "Species") %>% 
  relocate(label, .before = everything())


p <- ggtree(tree) %<+% classifications +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p + new_scale_color() + geom_facet(
  data = tmp,
  mapping = aes(x = as.numeric(start), color = V2),
  geom = geom_point,       # Horizontal bar chart
  panel = "Position in 288S"
) +
  scale_y_discrete()+scale_color_viridis_b(option = "viridis", n.breaks=10) + labs(color = "GC content (%)")+ coord_cartesian(clip = 'off')+ xlim_expand(c(0, 25), "Tree")


ggsave("zc4688/honours/analyses/chordates/figures/its1_gc.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1_gc.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = n, fill = group)) + 
  geom_boxplot() + 
  scale_fill_manual(values = colours) + 
  theme_bw() + 
  labs(x = NULL, y = "Number of windows with >80% GC content", fill = "Taxonomic Group")
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.length.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.length.png", height = 20, width = 15)

classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(n ~ group, p.adj.method = "BH") %>% view

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = max, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "Maximum GC % in ITS1", fill = "Taxonomic Group")
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.max.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.max.png", height = 20, width = 15)

classifications %>% 
  dplyr::rename("label" = "Species") %>% 
  left_join(tmp %>% group_by(label) %>% summarise(max = max(V2))) %>%  
  filter(group %in% common_class) %>% 
  pairwise_wilcox_test(max ~ group, p.adj.method = "BH") %>% view


classifications %>% dplyr::rename("label" = "Species") %>% 
  left_join(
    tmp %>%
      filter(V2 > 80) %>%
      count(label, name = "n"),  # count how many rows per label
    by = "label"
  ) %>%
  mutate(n = replace_na(n, 0)) %>% 
  left_join(
    structuremetadata
  ) %>%
  mutate(remain = ITS1 - n*50) %>% 
  filter(group %in% common_class) %>% 
  ggplot(aes(x = group, y = remain, fill = group)) + 
  geom_boxplot() + scale_fill_manual(values = colours) + 
  theme_bw() + labs(x = NULL, y = "ITS1 length without GC rich regions", fill = "Taxonomic Group") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.removed.pdf", height = 20, width = 15)
ggsave("zc4688/honours/analyses/chordates/figures/its1.gc.removed.png", height = 20, width = 15)



#################### SLOTHFISH
tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/giraffablast.out", header = F)
tmp %>% 
  mutate(`Sample ID` = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% 
  filter(`Sample ID` != "GCF_015220235") %>% 
  filter(V3 > 90) %>% view
  ggplot() + geom_boxplot(aes(x = group, y = V4))
                                                                                                                                                                                                     + )
tmp %>% 
  mutate(`Sample ID` = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% 
  filter(`Sample ID` != "GCF_015220235") %>% 
  filter(V3 > 90) %>% ggplot(aes(x = V7, xend = V8, y = V9, yend = V10)) + geom_segment(aes(color = group)) + facet_wrap(~group)

#seems like misassembly of contamination. the rRNAs cluster with mamalian IGS in the iqtree output but only actinoperygii has any significant matches in the igs

#################### Weird bat weird frog
#together on trees. more similar to each other, not to others in same clade. 

tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/antrozousblast.out", header = F)
tmp %>% 
  mutate(`Sample ID` = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% filter(Species == "Leptobrachium leishanense") %>% ggplot(aes(x = V9, xend = V10, y = V7, yend = V8))+ geom_segment()

############### creating new Hmm
tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/out.twoeight.self.blastn.mci.I20", header = F)
tmp %>% 
  rownames_to_column() %>% 
  pivot_longer(names_to = "test", values_to = "sample", cols = -rowname) %>% 
  select(-test) %>% 
  mutate(`Sample ID` = str_extract(sample, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% group_by(group) %>% 
  mutate(rowname = as.numeric(rowname)) %>%  
  arrange(group, rowname) %>% 
  filter(!is.na(group)) %>% view
  slice(c(1, round(n()/2), n())) %>% 
  ungroup() %>% 
  select(sample) %>% 
  write.table("zc4688/honours/analyses/chordates/new/msa/twoeightreps.txt", row.names = F, col.names = F, quote = F)
  
  
  
  
  #############################Mammalgaps
  #inspect breaks
  primarystructure %>% filter(type == "28S_rRNA") %>% select(label, hmmstart, hmmend, Start, End, `Sample ID`, group) %>% group_by(label) %>% arrange(label, Start)%>% filter((End - Start) > 100) %>% ggplot() + geom_jitter(aes(x = hmmstart, y = 0, color = group), height = 1) + theme_bw() + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
  
  
  #######################################what are the breaks?
  mammalgaps <- primarystructure %>% 
    filter(group == "Mammalia" & type == "28S_rRNA") %>% 
    group_by(label) %>% mutate(n = n()) %>% 
    filter(n > 1) %>% 
    arrange(label, Start) %>% 
    mutate(gap = Start - lag(End), gapstart = lag(End), gapend = Start) %>% 
    filter(gap > 500 & gap < 800) %>% 
    mutate(name = paste(`Sample ID`, refmorph, sep = "_")) %>% 
    ungroup()  %>% select(name, gapstart, gapend, `Sample ID`) 
  
  write.table(mammalgaps, "zc4688/honours/analyses/chordates/mammalgaps.csv", sep = ",", quote = F, row.names = F, col.names = F)
  
  #from blast - they are highly similar
  tmp <- read.delim("zc4688/honours/analyses/chordates/mammalgaps.self.blastn.out", header = F)
  tmp <- tmp %>% 
    mutate(sample1 = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002"), sample2 = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
    group_by(sample1, sample2) %>% slice_max(order_by = V8 - V7, with_ties = FALSE) %>% select(sample1, sample2, V3) %>% 
    pivot_wider(names_from = sample2, values_from = V3, values_fill = 0) %>% column_to_rownames("sample1")
  tmp <- tmp[intersect(rownames(tmp), colnames(tmp)), intersect(rownames(tmp), colnames(tmp))]
  pheatmap(tmp, cluster_rows = F, cluster_cols = F)
  
  tmp <- read.delim("zc4688/honours/analyses/chordates/mammalgaps.self.blastn.out", header = F)
  tmp <- tmp %>% 
    mutate(sample1 = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002"), sample2 = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% left_join(mammalgaps %>% dplyr::rename("sample1" = "Sample ID")) %>% filter(V3 > 75) %>% 
    group_by(sample1, sample2) %>% mutate(length = V8 - V7, prop = length/(gapend - gapstart)) %>% slice_max(order_by = prop, with_ties = FALSE)  %>% select(sample1, sample2, prop) %>% 
    pivot_wider(names_from = sample2, values_from = prop, values_fill = 0) %>% column_to_rownames("sample1")
  tmp <- tmp[intersect(rownames(tmp), colnames(tmp)), intersect(rownames(tmp), colnames(tmp))]
  pheatmap(tmp, cluster_rows = F, cluster_cols = F)
  
  
  tmp <- read.delim("zc4688/honours/analyses/chordates/mammalgaps.self.blastn.out", header = F)
  tmp <- tmp %>% 
    mutate(sample1 = str_extract(V1, "^GCA_\\d+|^GCF_\\d+|hg002"), sample2 = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
    left_join(classifications %>% dplyr::rename("sample2" = "Sample ID") %>% select(sample2, Species)) %>% 
    dplyr::rename("species2" = "Species") %>% 
    select(-sample2) %>% 
    left_join(classifications %>% dplyr::rename("sample1" = "Sample ID")) %>% dplyr::rename("label" = "Species") %>% 
    left_join(mammalgaps %>% dplyr::rename("sample1" = "Sample ID")) %>% 
    select(-sample1) %>% 
    filter(V3 > 75) %>% 
    group_by(label, species2) %>% 
    mutate(length = V8 - V7, prop = length/(gapend - gapstart)) %>% 
    slice_max(order_by = prop, with_ties = FALSE)  %>% 
    select(label, species2, prop) %>% 
    pivot_wider(names_from = species2, values_from = prop, values_fill = 0) %>% 
    relocate(label, .before = everything())
  
  mammaltree <- ape::drop.tip(tree, setdiff(tree$tip.label, tmp$label))
  
  tmp <- tmp[, c(1, match(mammaltree$tip.label, colnames(tmp)[-1]) + 1)]
  
  p <- ggtree(mammaltree) + 
    geom_tiplab(size=2) + 
    theme_tree2()
  
  gheatmap(p, tmp,
           colnames=TRUE, legend_title="similarity",
           low = "blue", high = "red")
  
  
  #from trf - they are repeat heavy but not completely repeats (around 50%)
  trfgaps <- read.delim("zc4688/honours/analyses/chordates/mammalgaps.fa.2.7.7.80.10.10.12.dat", header = F)
  trfgaps <- trfgaps %>% filter(!str_detect(V1, "Parameters|Tandem|Gary|Program|Boston|Version")) %>% mutate(sequence = ifelse(str_detect(V1, "Sequence"), V1, NA)) %>% fill(sequence) %>% mutate(sequence = str_replace_all(sequence, "Sequence: ", "")) %>% filter(!str_detect(V1, "Sequence")) %>% separate(V1, into = c("start", "end"), sep = " ") %>% mutate(start = as.numeric(start), end = as.numeric(end))
  trfgaps <- GRanges(seqnames = trfgaps$sequence, ranges = IRanges(start = trfgaps$start, end = trfgaps$end))
  trfgaps <- GenomicRanges::reduce(trfgaps, min.gapwidth = 10)
  trfgaps <- data.frame(sequence = as.character(seqnames(trfgaps)), start = start(trfgaps), end = end(trfgaps))
  trfgaps %>% mutate(`Sample ID` =str_extract(sequence, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% group_by(`Sample ID`) %>% summarise(repeatlength = sum(end - start)) %>% left_join(mammalgaps) %>% mutate(prop = repeatlength/gap) %>% view
  #TOMORROW: CALCULATE LENGTH AS PROPORTION OF GAP LENGTH. NEED TO EXTRACT GAPLENGTH FROM MAMMALGAPS THEN PLOT
  
  
  
  
  
  
  
  
  
  
  ########################### 28S clusters
  tmp <- read.delim("zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/all_clusters.tsv", header = F)
  
  tmp <- tmp %>% mutate('Sample ID' = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% left_join(lengths)
  
  tmp <- tmp %>% relocate(Species, .before = everything()) %>% dplyr::rename("label" = "Species")
  
  
  p + geom_facet(
    data = plotdata,
    mapping = aes(x = ifelse(`Unit length` < 60000, `Unit length`/1000, 60), fill = group),
    geom = geom_col,       # Horizontal bar chart
    panel = "Unit length (Kb)",
    width = 0.6 
  ) + geom_facet(data = tmp, mapping = aes(x = 0, color = V1), geom = geom_point,       # Horizontal bar chart
                 panel = "Unit length (Kb)") +
    scale_y_discrete() +  
    theme_tree2() + xlim_expand(c(0, 23), "Tree") + xlim_expand(c(0, 60), "Unit length (Kb)") +
    scale_fill_manual(values = colours) + labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 
  
  #mammals all cluster together very consistently - even those with the break. seems that it is just more extended in some mammals above probably an arbitrary hmm threshold
tmp <- classifications %>% filter(group == "Mammalia") %>% select(`Sample ID`) %>% unlist()
msa <- readDNAStringSet("zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/final_alignment_with_all.fasta")
sample_id_from_name <- sub("^(GCA_\\d+|GCF_\\d+|hg002).*", "\\1", names(msa))
subset_fasta <- msa[sample_id_from_name %in% tmp]
writeXStringSet(subset_fasta, "zc4688/honours/analyses/chordates/new/msa/twoeight_clusters/mammaltwoeight.fa")





#looking further at the gaps
# aligned all 28S sequences to the yeast 28S - aim to determine whether the insertions are big/small, where they are etc.

#as expected, mammals have an abundance of large gaps
tmp %>% 
  select(V1, V2, V7, V8, V9, V10) %>% 
  group_by(V2) %>% 
  arrange(V2, V9) %>% 
  mutate(yeastgap = V7 - lag(V8), chorgap = V9 - lag(V10))  %>% 
  mutate(percent = abs(chorgap - yeastgap)/yeastgap) %>% 
  mutate(`Sample ID` = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% 
  mutate(diff = (chorgap - yeastgap)) %>% 
  filter(diff > 100) %>% 
  view %>% 
  ggplot(aes(x = diff)) + 
  geom_density(aes(fill = group), alpha = 0.6) + 
  scale_fill_manual(values = colours) + theme_bw()


#These gaps correspond 2 to main insertion segments. the first is NOT the big mammalian specific one, its earlier in the sequence. the second is the massive one in mammals
tmp %>% 
  select(V1, V2, V7, V8, V9, V10) %>% 
  group_by(V2) %>% 
  arrange(V2, V9) %>% 
  mutate(yeastgap = V7 - lag(V8), chorgap = V9 - lag(V10))  %>%
  mutate(percent = abs(chorgap - yeastgap)/yeastgap) %>% 
  mutate(`Sample ID` = str_extract(V2, "^GCA_\\d+|^GCF_\\d+|hg002")) %>% 
  left_join(classifications) %>% 
  mutate(diff = (chorgap - yeastgap)) %>% 
  view %>% filter(diff > 200) %>% 
  ggplot(aes(x = V7, fill = group)) + 
  geom_histogram() + 
  scale_fill_manual(values = colours) + 
  theme_bw()


####################Inspect fails

allspecies <-read.tree("/g/data/te53/zc4688/honours/trees/allspecies.phy")
allspecies$edge.length <- NULL

allspeciestree <- combined_df %>% 
  mutate(success = ifelse(is.na(`Errors.Alignment filtering`) & is.na(`Errors.Morph identification`) & is.na(`Errors.Non-ambiguous morph identification`) & is.na(`Errors.Primary alignments`), "yes", "no")) %>% 
  select(`Sample_details.Species`, success) %>% 
  dplyr::rename("label" = "Sample_details.Species")

allspecies$tip.label <- str_replace_all(allspecies$tip.label, " ", "_")
allspecies$tip.label <- str_replace_all(allspecies$tip.label, "'", "")

common_organisms <- intersect(allspeciestree$label, allspecies$tip.label)

allspecies <- ape::drop.tip(allspecies, setdiff(allspecies$tip.label, common_organisms))

allspeciestree <- allspeciestree %>% filter(label %in% common_organisms) %>% left_join(classifications %>% dplyr::rename("label" = "Species") %>% mutate(label = str_replace_all(label, " ", "_")))

ggtree(allspecies) %<+% allspeciestree +
  geom_tiplab(size=1, aes(color = success)) +
  scale_color_manual(values = c("yes" = "black", "no" = "red")) + 
  new_scale_color() + 
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group", na.translate = FALSE)


rm(allspecies)
rm(allspeciestree)



bam <- BamFile("zc4688/honours/analyses/cn/HG00729.test.bam")

gal <- readGAlignments(bam)

# Extract CIGAR strings
cigars <- cigar(gal)
names(cigars) <- mcols(gal)$qname
ops_list <- explodeCigarOps(cigars)   # list of operations per read
lens_list <- explodeCigarOpLengths(cigars)  # lengths of each operation
ref_len <- max(end(gal))

# Initialize matrix
mat <- matrix(NA, nrow = length(gal), ncol = ref_len)

for(i in seq_along(gal)) {
  ref_pos <- start(gal[i]) - 1
  ops <- ops_list[[i]]
  lens <- lens_list[[i]]
  
  for(j in seq_along(ops)) {
    op <- ops[j]
    len <- lens[j]
    
    if(op == "M") {
      mat[i, (ref_pos + 1):(ref_pos + len)] <- 0  # match/mismatch
      ref_pos <- ref_pos + len
    } else if(op == "I") {
      mat[i, ref_pos] <- 2  # insertion at current reference
      # ref_pos not advanced
    } else if(op == "D") {
      mat[i, (ref_pos + 1):(ref_pos + len)] <- 1  # deletion
      ref_pos <- ref_pos + len
    }
  }
}


library(reshape2)

df <- melt(mat)
colnames(df) <- c("Read", "RefPos", "VarType")

ggplot(df, aes(x = RefPos, y = Read, fill = factor(VarType))) +
  geom_tile() +
  scale_fill_manual(values = c("0"="grey", "1"="red", "2" = "blue"), 
                    labels = c("Match/Mismatch", "Deletion", "Insertion")) +
  theme_bw() +
  labs(x = "Reference position", y = "Read", fill = "Variant type") +
  theme(axis.text.y = element_blank())



################################correlations
moremetadata <- read.delim("zc4688/honours/metadata/assembly_metadata.tsv")
metadata <- left_join(metadata, moremetadata %>% separate(accession, into = c("Sample ID", "number"), sep = "\\."))
summary(lm(`Unit length` ~ `28S_length`, data = metadata))

summary(lm(`Unit length` ~ `genome_size`, data = metadata))

tmp <- metadata %>% separate(assembly_stats, sep = ",", into = c("contigL50", "contigN50", "gcCount", "gcPercent", "other"))  %>% mutate(gcPercent = as.numeric(str_replace(gcPercent, "gcPercent:", "")))

ggplot(metadata, aes(x = `ITS2`, y = `Unit length`, color = group)) + 
  geom_point() + 
  scale_color_manual(values = colours) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 120000))

tes <- read.delim("zc4688/honours/metadata/tes.txt")
tmp <- metadata %>% select(`Unit length`, `group`, genome_size, `28S_length`, `18S_length`, IGS_length, ITS1, ITS2, assembly_stats)%>% separate(assembly_stats, sep = ",", into = c("contigL50", "contigN50", "gcCount", "gcPercent", "other"))  %>% mutate(gcPercent = as.numeric(str_replace(gcPercent, "gcPercent:", "")))
predictor_cols <- setdiff(names(tmp), "Unit length")

map(predictor_cols, ~ {
  formula_str <- paste0("`Unit length` ~ `", .x, "`")
  formula <- as.formula(formula_str)
  model <- lm(formula, data = tmp, na.action = "na.omit")
  summary(model)$coefficients[2, ]  # extract slope info only (not intercept)
}) %>% bind_rows(.id = "predictor") %>%
  dplyr::rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    t_value = `t value`,
    p_value = `Pr(>|t|)`
  ) %>% view


age <- read.delim("zc4688/honours/metadata/anage_data.txt")
age <- age %>% mutate(Species = paste(Genus, Species, sep = " ")) %>% 
  left_join(metadata) %>% filter(!is.na(`Unit length`))

ggplot(age) + geom_point(aes(x = log(Metabolic.rate..W.), y = `Unit length`, color = group)) + coord_cartesian(ylim = c(0, 60000)) + scale_color_manual(values = colours) + theme_bw()
#doesnt seem to much/anything of significance here



age <- read.csv("zc4688/honours/metadata/observations.csv")
age %>% group_by(species) %>% summarise(brain = max(`brain size`), mass = max(`body mass`)) %>% left_join(metadata %>% rename("species" = "Species")) %>% filter(!is.na(Unit.length)) %>% ggplot(aes(x = log(mass), y = Unit.length)) + geom_point(aes(color = group))
