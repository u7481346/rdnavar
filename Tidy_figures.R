library(ggtree)
library(tidyverse)
library(ape)
library(jsonlite)
library(GenomicRanges)
############################### Figure 1: length tree
directory <- "/g/data/te53/zc4688/honours/ribocop/results/chordata_filter2"

# List all JSON files in the directory
json_files <- list.files(directory, pattern = "\\.json$", full.names = TRUE)

# Initialize an empty list to store the data frames
data_list <- list()

# Loop through each file and read the JSON into a data frame
for (file in json_files) {
  tryCatch({
    # Read the file as plain text
    json_text <- paste(readLines(file, warn = FALSE), collapse = "")
    
    # Replace invalid `NaN` with `null`
    json_text <- gsub("\\bNaN\\b", "null", json_text)
    
    # Parse the JSON
    data <- fromJSON(json_text, flatten = TRUE)
    
    # Ensure data is a data frame with one row
    if (!is.data.frame(data)) {
      data <- as.data.frame(t(unlist(data)), stringsAsFactors = FALSE)
    }
    
    data_list[[file]] <- data
  })
}


# Combine all data frames into one large data frame
combined_df <- bind_rows(data_list)

combined_df <- combined_df %>% 
  filter(is.na(`Errors.rDNA identification`) & is.na(`Errors.Morph identification`)) %>% 
  filter(!is.na(`rDNA_details.Number of primary alignments`)) %>% 
  dplyr::select(`Sample_details.Sample ID`, `Sample_details.Species`, `Sample_details.TaxID`, 
                `rDNA_details.Unit length`, `rDNA_details.Total length`, `rDNA_details.Number of contigs`,
                `rDNA_details.Minimum unit length`, `rDNA_details.Maximum unit length`) %>% 
  dplyr::rename(
    `Sample ID` = `Sample_details.Sample ID`,
    Species = `Sample_details.Species`,
    TaxID = `Sample_details.TaxID`,
    `Unit length` = `rDNA_details.Unit length`,
    `Total length` = `rDNA_details.Total length`,
    `Contigs` = `rDNA_details.Number of contigs`,
    `Minimum unit length` = `rDNA_details.Minimum unit length`,
    `Maximum unit length` = `rDNA_details.Maximum unit length`
  ) %>% 
  mutate(`Total length` = as.numeric(`Total length`)) %>% 
  mutate(`Unit length` = as.numeric(`Unit length`)) %>% 
  mutate(`Minimum unut length` = as.numeric(`Minimum unit length`)) %>% 
  mutate(`Maximum unit length` = as.numeric(`Maximum unit length`)) %>% 
  mutate(CN = round(`Total length`/`Unit length`, digits=2))


lengths <- combined_df

lengths <- lengths %>% 
  filter(!is.na(`Unit length`))

tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict
plotdata <- lengths %>% select(Species, `Unit length`, `CN`, `Total length`, `Minimum unit length`, `Maximum unit length`, `Sample ID`, Contigs) %>% 
  group_by(Species) %>% dplyr::summarize(across(everything(), ~ dplyr::first(.))) %>% 
  filter(`Unit length` < 120000)

test <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

test <- test %>% mutate(Species = str_replace_all(Species, " ", "_"), Replacement = str_replace_all(Replacement, " ", "_") )

plotdata <- left_join(plotdata, test) 

plotdata <- plotdata %>% 
  mutate(label = ifelse(is.na(Replacement), Species, Replacement))

metadata <- plotdata
tree$tip.label <- str_replace_all(tree$tip.label, "'", "")

#If there are multiple, take the first (GCA) - only 2 species actually have different results across


common_organisms <- intersect(plotdata$label, tree$tip.label)


plotdata <- plotdata %>%
  filter(label %in% common_organisms)



#classifications <- classification(plotdata$label, db="ncbi")
#classifications <- cbind(classifications)
#classifications <- classifications %>% 
  #select(kingdom, phylum, subphylum, superclass, class, subclass, superorder, order, suborder, family, subfamily, genus, species, query)
#write.table(classifications, "zc4688/honours/ribocop/classifications.tsv", sep = "\t", quote = F, row.names = F)
classifications <- read.table("zc4688/honours/ribocop/classifications.tsv", sep = "\t", header = TRUE)




classifications <- classifications %>% 
  dplyr::select(class, order, species) %>% 
  dplyr::rename("label" = "species") 

classifications <- classifications %>% 
  group_by(class, order, label) %>% summarize(across(everything(), ~ first(.))) 
plotdata <- plotdata %>% mutate(label = str_replace_all(label, "_", " "))
plotdata <- left_join(plotdata, classifications, by="label")






tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))


tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")

plotdata <- plotdata %>% mutate(class = ifelse(is.na(class), "Other", class)) %>% relocate(label, .before = everything())

############New class

plotdata <- plotdata %>%
  mutate(
    class = case_when(
      label %in% c("Arripis georgianus", "Cyprinodon nevadensis mionectes", "Gymnocypris eckloni scoliostomus") ~ "Actinopteri",
      label %in% c("Bubalus carabanensis", "Canis lupus dingo", "Diceros bicornis minor", "Elephas maximus indicus", 
                   "Hippopotamus amphibius kiboko", "Lagenorhynchus acutus", "Ovis ammon polii", 
                   "Perognathus longimembris pacificus", "Gorilla gorilla gorilla") ~ "Mammalia",
      label %in% c("Elgaria multicarinata webbii") ~ "Lepidosauria",
      label %in% c("Melospiza melodia melodia", "Motacilla alba alba") ~ "Aves",
      TRUE ~ class  # Keep existing values if no match
    ),
    order = case_when(
      label %in% c("Chrysemys picta bellii", "Malaclemys terrapin pileata") ~ "Testudines",
      label %in% c("Latimeria chalumnae") ~ "Coelacanthiformes",
      TRUE ~ order  # Keep existing values if no match
    ),
    class = case_when(
      order %in% c("Coelacanthiformes", "Ceratodontiformes") ~ "Sarcopterygii",
      order %in% c("Crocodylia", "Testudines") ~ "Testudines and Crocodylia",
      TRUE ~ class
    ),
    group = case_when(
      class %in% c("Actinopteri", "Cladistia") ~ "Actinopterygii",
      class %in% c("Ascidiacea", "Appendicularia") ~"Tunicata",
      TRUE ~ class
    )
  )

classifications <- plotdata %>% select(label, group)
groups <- unique(plotdata$group)
colours <- setNames(brewer.pal(n = 11, name = "Paired"), groups)

p <- ggtree(tree) + 
  geom_tiplab(size=1) 



p <- p + geom_facet(
  data = plotdata,
  mapping = aes(x = `Unit length`/1000, fill = group),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
)  +
  geom_facet(
    data = plotdata,
    mapping = aes(x = log(CN, base = 10), fill = group),
    geom = geom_col,      
    panel = "log10(CN)",
    width = 0.6
  ) +
  geom_facet(
    data = plotdata,
    mapping = aes(x = `Total length`/1000000, fill = group),
    geom = geom_col,      
    panel = "Total length (Mb)",
    width = 0.6
  )  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree")+scale_fill_manual(values = colours) + labs(fill = "Taxonomic Group") +coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(3, 1, 1, 1, 1))

ggsave("zc4688/honours/ribocop/results/figures/tree_current.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/tree_current.pdf", height = 20, width = 15)

##############################################################  Figure 2: alignment tree (blast vs library)
blast <- read_delim("/g/data/te53/zc4688/honours/ribocop/results/chordata_filter2/alignments/blast/libalignment.tsv", col_names = F)

colnames(blast) <- c("queryseq", "libsequence", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart",  "send", "evalue", "bitscore")

species <- blast %>% pull(queryseq) %>% unique


ranges_list <- list()


for (sample in species){
  tmp <- blast %>% filter(queryseq == sample)
  for (i in c("_18S", "_28S")){
    alignments <- tmp %>% filter(str_detect(libsequence, i)) %>% filter(length > 500)
    ranges <- IRanges(start = alignments$qstart, end = alignments$qend)
    reduced_ranges <- reduce(ranges, min.gapwidth = 50)
    
    start_pos <- start(reduced_ranges)
    end_pos <- end(reduced_ranges)
    
    if (!is.list(ranges_list[[sample]])) {
      ranges_list[[sample]] <- list()
    }
    
    ranges_list[[sample]][[i]] <- data.frame(unit = i, start = start_pos, end = end_pos, samplename = sample)
  }
  
}

all_ranges <- bind_rows(lapply(ranges_list, bind_rows))

all_ranges <- all_ranges %>%
  group_by(samplename, unit) %>%
  mutate(first_alignment = row_number() == 1,
         length = end - start,
         longest_alignment = length == max(length)) %>% 
  ungroup()


ggplot(all_ranges) +
  geom_segment(aes(
    y = as.numeric(factor(samplename)), 
    yend = as.numeric(factor(samplename)), 
    x = start, 
    xend = end, 
    color = ifelse(longest_alignment, unit, paste0(unit, "_secondary"))  # Color longest alignment, grey for others
  )) +
  scale_color_manual(values = c("_18S" = "lightpink", "_28S" = "skyblue", "_18S_secondary" = "red2", "_28S_secondary" = "darkblue" )) + 
  labs(x = "Morph position", y = "Species", color = "Unit") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



itstreetmp <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (4).nwk") 


itsplot <- all_ranges %>%
  select(-first_alignment) 
itsplot <- itsplot %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", samplename)) %>% 
  left_join(lengths %>%  select(`Sample ID`, Species)) %>% 
  mutate(label = str_replace_all(Species, "_", " ")) %>% 
  left_join(classifications, by = "label") %>% 
  select(-Species, -samplename, -`Sample ID`) %>% 
  select(label, start, end, unit, longest_alignment)


itstreetmp$tip.label <- str_replace_all(itstreetmp$tip.label, "_", " ")

common_organisms <- intersect(itsplot$label, itstreetmp$tip.label)





itsplot <- itsplot %>%
  filter(label %in% common_organisms)



itstreetmp$edge.length <- NULL
itstreetmp <- ape::drop.tip(itstreetmp, setdiff(itstreetmp$tip.label, common_organisms))




p <- ggtree(itstreetmp) + 
  geom_tiplab(size=1) 


#replacement data frame error is due to order of columns! just need 'label' to be first

p <- p + geom_facet(
  data = itsplot,
  mapping = aes(x = start, 
                xend = end, color = ifelse(longest_alignment, unit, paste0(unit, "_secondary"))),
  geom = geom_segment,       
  panel = "Morph Position",
  width = 0.6
) +
  scale_y_discrete()+
  scale_color_manual(values = c("_18S" = "lightpink", "_28S" = "skyblue", "_18S_secondary" = "red2", "_28S_secondary" = "darkblue" )) +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position = "right") + labs(color = "Unit") + coord_cartesian(clip = 'off') 

facet_widths(p, widths = c(1, 2))

ggsave("zc4688/honours/ribocop/results/figures/itstree.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/itstree.pdf", height = 20, width = 15)



############################### Figure 3: BLAST boxplots

classtree <- read.tree("/g/data/te53/zc4688/honours/ribocop/Chordates_class.nwk")
all_wide <- all_ranges %>% filter(longest_alignment == TRUE) %>% select(-first_alignment, -length) %>% pivot_wider(
  names_from = unit,
  values_from = c(start, end),
  names_glue = "{unit}_{.value}"
) %>% 
  mutate(itslength = `_28S_start` - `_18S_end`,
         `18S_length` = `_18S_end` - `_18S_start`,
         `28S_length` = `_28S_end` - `_28S_start`)
boxplotdata <- all_wide %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", samplename)) %>% 
  left_join(lengths) %>% 
  mutate(label = str_replace_all(Species, "_", " ")) %>% 
  left_join(classifications, by = "label") %>%
  dplyr::rename(`ITS Length` = itslength, `18S Length` = `18S_length`, `28S Length` = `28S_length`, `Unit Length` = `Unit length`) %>% 
  mutate(`IGS Length` = `Unit Length` - `18S Length` - `28S Length` - `ITS Length`) %>% 
  pivot_longer(cols = c(`ITS Length`, `18S Length`, `28S Length`, `Unit Length`, `CN`, `IGS Length`), names_to = "Variable", values_to = "Value") %>% 
  filter(Variable != "CN") %>% #stop here for old plot
  select(-label) %>% 
  dplyr::rename(label = class) %>% 
  select(label, Variable, Value)



#note ascidiacea is excluded here
common_class <- intersect(boxplotdata$label, classtree$tip.label)





boxplotdata <- boxplotdata %>%
  filter(label %in% common_class)



classtree$edge.length <- NULL
classtree <- ape::drop.tip(classtree, setdiff(classtree$tip.label, common_class))




p <- ggtree(classtree) + p + new_scale_color()+ geom_facet(
  data = boxplotdata %>% filter(Variable == "Unit Length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable, color = Variable),
  geom = geom_boxplot,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = boxplotdata %>% filter(Variable == "Unit Length"), 
  mapping = aes(xintercept = median(Value)/1000), 
  geom = geom_vline,
  panel = "Unit Length (Kb)", linetype = "dashed"
  ) + geom_facet(
    data = boxplotdata %>% filter(Variable == "IGS Length"),
    mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable, color = Variable),
    geom = geom_boxplot,       
    panel = "IGS Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  ) + geom_facet(
    data = boxplotdata %>% filter(Variable == "IGS Length"), 
    mapping = aes(xintercept = median(Value)/1000), 
    geom = geom_vline,
    panel = "IGS Length (Kb)", linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "18S Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "18S Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "18S Length"),
  mapping = aes(xintercept = median(Value)/1000),
  geom = geom_vline, 
  panel = "18S Length (Kb)", 
  linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "28S Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "28S Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "28S Length"),
  geom = geom_vline, aes(xintercept = median(Value)/1000), panel = "28S Length (Kb)", linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "ITS Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "ITS Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
    data = boxplotdata %>% filter(Variable == "ITS Length"),
    geom = geom_vline, 
    mapping = aes(xintercept = median(Value)/1000), panel = "ITS Length (Kb)", linetype = "dashed"
    )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position =  "none")  + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1")

  geom_tiplab(size=3) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours)


p <- p + new_scale_color() +geom_facet(
  data = boxplotdata %>% filter(Variable == "Unit Length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable, color = Variable),
  geom = geom_boxplot,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = boxplotdata %>% filter(Variable != "Unit Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "Subunit Lengths (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position = "right") + labs(color = NULL, fill = NULL) + scale_fill_brewer(palette = "Set1")+ scale_color_brewer(palette = "Set1") +
  guides(color = "none")

facet_widths(p, widths = c(1, 1, 1))

ggsave("zc4688/honours/ribocop/results/figures/boxplottree.png", height = 8, width = 12)
ggsave("zc4688/honours/ribocop/results/figures/boxplottree.pdf", height = 8, width = 12)

p <- ggtree(classtree) + 
  geom_tiplab(size=3) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours)


#replacement data frame error is due to order of columns! just need 'label' to be first
#need to annotate classes somehow

p <- p + new_scale_color()+ geom_facet(
  data = boxplotdata %>% filter(Variable == "Unit Length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable, color = Variable),
  geom = geom_boxplot,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = boxplotdata %>% filter(Variable == "Unit Length"), 
  mapping = aes(xintercept = median(Value)/1000), 
  geom = geom_vline,
  panel = "Unit Length (Kb)", linetype = "dashed"
  ) + geom_facet(
    data = boxplotdata %>% filter(Variable == "IGS Length"),
    mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable, color = Variable),
    geom = geom_boxplot,       
    panel = "IGS Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  ) + geom_facet(
    data = boxplotdata %>% filter(Variable == "IGS Length"), 
    mapping = aes(xintercept = median(Value)/1000), 
    geom = geom_vline,
    panel = "IGS Length (Kb)", linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "18S Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "18S Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "18S Length"),
  mapping = aes(xintercept = median(Value)/1000),
  geom = geom_vline, 
  panel = "18S Length (Kb)", 
  linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "28S Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "28S Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "28S Length"),
  geom = geom_vline, aes(xintercept = median(Value)/1000), panel = "28S Length (Kb)", linetype = "dashed"
  )+
  geom_facet(
  data = boxplotdata %>% filter(Variable == "ITS Length"),
  mapping = aes(x = Value/1000, fill = Variable, color = Variable, group = interaction(label, Variable)),
  geom = geom_boxplot,       
  panel = "ITS Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
    data = boxplotdata %>% filter(Variable == "ITS Length"),
    geom = geom_vline, 
    mapping = aes(xintercept = median(Value)/1000), panel = "ITS Length (Kb)", linetype = "dashed"
    )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position =  "none")  + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1")

facet_widths(p, widths = c(1, 1, 1, 1, 1, 1))
ggsave("zc4688/honours/ribocop/results/figures/boxplottree_separate.png", height = 8, width = 12)
ggsave("zc4688/honours/ribocop/results/figures/boxplottree_separate.pdf", height = 8, width = 12)


ggplot(boxplotdata %>% filter(Variable != "CN"), aes(y = Value/1000, fill = label)) +
  geom_boxplot(notch = TRUE, outliers = FALSE, coef = 2.5) +
  facet_wrap(~ Variable, scales = "free_y")  + 
  theme_bw() + 
  labs(fill = "Class", y = "Length (Kb)")+
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
                     strip.background = element_blank(),strip.placement = "outside") + 
  scale_fill_manual(values = colours)

ggsave("zc4688/honours/ribocop/results/figures/boxplot.png", height = 8, width = 12)



########################### Figure 4: Original barrnap
barrnap <- read.delim("zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/barrnap", header = F)


barrnap <- barrnap %>% 
  separate(V9, into = c("Name", "Other"), sep = ";") %>% 
  separate(Name, into = c("Name", "Unit"), sep = "=") %>% 
  separate(Other, into = c("Other", "Completeness"), sep = "\\(") %>% 
  mutate(Completeness = str_replace_all(Completeness, "\\)", "")) %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", V1)) %>% 
  left_join(lengths) %>% 
  mutate(label = str_replace_all(Species, "_", " ")) %>% 
  left_join(classifications, by = "label") 

ggplot(barrnap) +
  geom_segment(aes(x = V4, xend = V5, y = V1, yend = V1, color = ifelse(is.na(Completeness), Unit, Completeness))) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
  labs(y = "Species", color = "Unit", x = "Length (bp)") +
  scale_color_manual(values = c("partial" = "grey10", "18S_rRNA" = "pink", "28S_rRNA" = "lightblue", "5_8S_rRNA" = "lightgreen", "5S_rRNA" = "purple")) 



barrnapplot <- barrnap %>% select(Species, Unit, V4, V5, Completeness) %>% 
  dplyr::rename("label" = "Species") 
tree$tip.label <- str_replace_all(tree$tip.label, " ", "_")



common_organisms <- intersect(barrnapplot$label, tree$tip.label)


barrnapplot <- barrnapplot %>%
  filter(label %in% common_organisms)


classifications <- read.table("zc4688/honours/ribocop/classifications.tsv", sep = "\t", header = TRUE)




classifications <- classifications %>% 
  dplyr::select(class, order, species, family) %>% 
  dplyr::rename("label" = "species") 

classifications <- classifications %>% 
  group_by(class, order, label, family) %>% summarize(across(everything(), ~ first(.))) 
barrnapplot <- barrnapplot %>% mutate(label = str_replace_all(label, "_", " "))
barrnapplot <- left_join(barrnapplot, classifications, by="label")

#pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")


barrnapplot <- barrnapplot %>% mutate(class = ifelse(is.na(class), "Other", class))
p <- ggtree(tree) + 
  geom_tiplab(size=1) 



p <- p + geom_facet(
  data = barrnapplot,
  mapping = aes(x = V4, xend = V5, color = ifelse(is.na(Completeness), Unit, Completeness)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  width = 0.6 
)  +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree")+labs(color = "Unit") +
  scale_color_manual(values = c("partial" = "grey10", "18S_rRNA" = "pink", "28S_rRNA" = "lightblue", "5_8S_rRNA" = "lightgreen", "5S_rRNA" = "purple"))





########################### Figure 5: Barrnap with edits (marked Palign for every unit)
#Most recently - edited filtering slightly and added new columns. This made little differences except for removing some very low HMM palign alignments (3 total)
barrnap <- read.delim("zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/barrnap_newpalign.gff", header = F)
tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict

#colnames(barrnap) <- c("Sample", "Version", "rRNA", "Start", "End", "Envstart", "Envend", "Hmmstart", "Hmmend", "Score", "Strand", "Other", "Product")
colnames(barrnap) <- c("Sample", "Version", "Start", "End", "Envstart", "Envend", "Hmmstart", "Hmmend", "Score", "Bitscore", "Hmmpalign", "Targetpalign", "Strand", "Other", "Product")

barrnap <- barrnap %>% 
  filter(!is.na(Start)) %>% 
  separate(Product, into = c("Name", "Other", "Note"), sep = ";") %>% 
  separate(Name, into = c("Name", "Unit"), sep = "=") %>% 
  #separate(Other, into = c("Other", "Completeness"), sep = "\\(") %>% 
 # separate(Note, into = c("Note", "Palign"), sep = "=") %>% 
  #mutate(Completeness = str_replace_all(Completeness, "\\)", "")) %>% 
  select(-Note, -Name, -Other, -Version, -Strand) %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", Sample)) %>% 
  left_join(lengths) 



barrnap <- barrnap %>% group_by(Species) %>%
  filter(`Sample ID` == first(`Sample ID`)) %>%
  ungroup() %>% view


barrnapplot <- barrnap %>% relocate(Species, .before = everything()) %>% rename(label = Species)

barrnapmetadata <- barrnapplot %>% 
  group_by(Unit, label, `Sample ID`) %>% 
  mutate(Hmmpalign = ifelse(Hmmpalign > 1, 1, Hmmpalign)) %>% 
  mutate(Hmmpalign = round(Hmmpalign, digits = 2)) %>% 
  filter(Hmmpalign == max((Hmmpalign))) %>% 
  arrange(Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(label, `Sample ID`, Unit, Envstart, Envend, `Unit length`) %>% pivot_wider(
    names_from = Unit,
    values_from = c(Envstart, Envend),
    names_glue = "{Unit}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`) %>% 
  select(ITS1, ITS2, `18S_length`, `28S_length`, `IGS_length`, label, `Unit length`, `Sample ID`)

common_organisms <- intersect(barrnapplot$label, tree$tip.label)


barrnapplot <- barrnapplot %>%
  filter(label %in% common_organisms)


#classifications <- read.table("zc4688/honours/ribocop/classifications.tsv", sep = "\t", header = TRUE)




#classifications <- classifications %>% 
 # dplyr::select(class, order, species, family) %>% 
  #dplyr::rename("label" = "species") 

#classifications <- classifications %>% 
 # group_by(class, order, label, family) %>% summarize(across(everything(), ~ first(.))) 
barrnapplot <- barrnapplot %>% mutate(label = str_replace_all(label, "_", " "))
barrnapplot <- left_join(barrnapplot, classifications, by="label")

#pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")




barrnapplot <- barrnapplot %>% mutate(Hmmpalign = Hmmpalign * 100)




tmp <- barrnapplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  + 
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.pdf", height = 20, width = 15)

#########################################separate vs broken alignbemnts - required for boxplots and for msa
test <- barrnapplot %>% 
  arrange(Sample, Envstart) %>% 
  group_by(label, Sample, Unit) %>% 
  mutate(query_gap = Envstart - lag(Envend, default = first(Envstart)), 
         ref_gap = Hmmstart - lag(Hmmend, default = first(Hmmstart)), 
         diff = query_gap - ref_gap, 
         type = case_when((ref_gap < -100 & query_gap > 100) | query_gap > 5000 ~ "separate",
                          ref_gap == 0 & query_gap == 0 ~ "single", 
                          TRUE ~ "gap")) %>% 
  select(label, Sample, Unit, Envstart, Envend, type, `Unit length`, group) %>% 
  mutate(
    merged_Envstart = if_else(type == "gap", lag(Envstart, default = first(Envstart)), Envstart),
    merged_Envend = Envend)  %>% 
  filter(!lead(type, default = "none") == "gap") %>%
  ungroup() %>% 
  filter(Unit != "5S_rRNA") %>% 
  mutate(Envstart = merged_Envstart,
         Envend = merged_Envend)

#inspect the merging
p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")
p <- p + new_scale_color() + geom_facet(
  data = test,
  mapping = aes(x = Envstart, xend = Envend, color = type),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  + 
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')
p
ggsave("zc4688/honours/ribocop/results/figures/gaptree.png", height = 20, width = 15)

msadf <- test %>% 
  group_by(Sample, Unit) %>% 
  filter((Envstart - Envend) == max(Envstart - Envend)) %>% 
  arrange(Sample, Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(Sample, Unit, Envstart, Envend)

eighteen <- msadf %>% filter(Unit == "18S_rRNA")
write.table(eighteen, "zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/msa/eighteen.tsv", sep = ",", quote = F, row.names = F, col.names = F)
twoeight <- msadf %>% filter(Unit == "28S_rRNA")
write.table(twoeight, "zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/msa/twoeight.tsv", sep = ",", quote = F, row.names = F, col.names = F)

twoeight_filtered <- msadf %>% filter(Unit == "28S_rRNA" & (Envend - Envstart) < 5500)
write.table(twoeight_filtered, "zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/msa/twoeight_filtered.tsv", sep = ",", quote = F, row.names = F, col.names = F)


fiveeight <- msadf %>% filter(Unit == "5_8S_rRNA")
write.table(fiveeight, "zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/msa/fiveeight.tsv", sep = ",", quote = F, row.names = F, col.names = F)

sequences <- readDNAMultipleAlignment("zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/msa/twoeight_msa.txt")
conservation_scores <- msaConservationScore(sequences, BLOSUM62)
human <- DNAStringSet(sequences)["GCF_000001405_NC_000021.9:8392567-8436777:+:99-1967"]
human <- DNAStringSet(sequences)["GCF_000001405_NC_000021.9:8392567-8436777:+:4361-9412"]
human <- strsplit(as.character(human), "")[[1]]

df <- data.frame(Position = 1:length(conservation_scores), Score = conservation_scores, 
                 Sequence = human)

ggplot(df) +
  geom_smooth(aes(x = Position, y = Score / 1000000)) +
  geom_point(aes(x = Position, y = 0, color = Sequence), size = 2) +  # Add sequence points
  theme_bw() +
  labs(title = NULL, x = "Position", y = "Conservation Score") +
  scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue"))

ggplot(df %>% filter(Sequence != "-")) +
  geom_point(aes(x = Position, y = Score / 1000000), size = 0.5) +
  geom_point(aes(x = Position, y = 0, color = Sequence), size = 2) +  # Add sequence points
  theme_bw() +
  labs(title = NULL, x = "Position", y = "Conservation Score (x1000000)") +
  scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue"))

ggplot() +
  geom_point(data = df %>% filter(Sequence != "-" & Score > 3000000), aes(x = Position, y = Score / 1000000), size = 0.5, color = "red") +geom_point(data = df %>% filter(Sequence != "-" & Score < 3000000 & Score > 0), aes(x = Position, y = Score / 1000000), size = 0.5, color = "orange") +geom_point(data = df %>% filter(Sequence != "-" & Score < 0), aes(x = Position, y = Score / 1000000), size = 0.5, color = "green") +
  geom_point(data = df %>% filter(Sequence != "-"), aes(x = Position, y = 0, color = Sequence), size = 2) +  # Add sequence points
  theme_bw() +
  labs(title = NULL, x = "Position", y = "Conservation Score (x1000000)") +
  scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue"))

df %>% mutate(Score = as.numeric(Score), category = case_when(Score > 3000000 ~ "conserved", Score <= 3000000 & Score > 0 ~"medium", Score <=0 ~ "not", TRUE ~ "NA")) %>%  ggplot(aes(x = Position, fill = category)) +
  geom_density(alpha = 0.5) +   geom_point(aes(x = Position, y = 0, color = Sequence), size = 2) +      scale_fill_manual(values = c("conserved" = "red", "medium" = "orange", "not" = "green")) +scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue")) +
  theme_bw() +
  labs(title = NULL, 
       x = "Position", 
       y = "Density") +
  theme(legend.title = element_blank())
df %>% mutate(Score = as.numeric(Score), category = case_when(Score > 3000000 ~ "conserved", Score <= 3000000 & Score > 0 ~"medium", Score <=0 ~ "not", TRUE ~ "NA")) %>%  ggplot(aes(x = Position, fill = category)) +
  geom_density(aes(y = after_stat(density * n/40434)), alpha = 0.5, bw = 250) +   geom_rug(aes(x = Position, y = 0, color = Sequence), size = 2, sides = "b") +      scale_fill_manual(values = c("conserved" = "red", "medium" = "orange", "not" = "green")) +scale_color_manual(values = c("A" = "blue", "C" = "purple", "G" = "purple", "T" = "blue")) +
  theme_bw() +
  labs(title = NULL, 
       x = "Position", 
       y = "Density") +
  theme(legend.title = element_blank())

barnapboxplot <- test %>% 
  group_by(Sample, Unit) %>% 
  filter((Envstart - Envend) == max(Envstart - Envend)) %>% 
  arrange(Sample, Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(Sample, Unit, Envstart, Envend, group, `Unit length`) %>% pivot_wider(
    names_from = Unit,
    values_from = c(Envstart, Envend),
    names_glue = "{Unit}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`)
################################ Figure 6: Barrnap boxplots
classtree <- read.tree("/g/data/te53/zc4688/honours/ribocop/Chordates_class.nwk")

#tnis section no longer used as per new gap/separate
barnapboxplot <- barrnapplot %>% 
  group_by(Sample, Unit) %>% 
  filter(Hmmpalign == max(Hmmpalign)) %>% 
  arrange(Sample, Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(Sample, Unit, Envstart, Envend, group, `Unit length`) %>% pivot_wider(
  names_from = Unit,
  values_from = c(Envstart, Envend),
  names_glue = "{Unit}_{.value}"
) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`)
         
barnapboxplot <- barnapboxplot %>% 
  pivot_longer(cols = c(`ITS1`, ITS2, `18S_length`, `28S_length`, `IGS_length`, `Unit length`), names_to = "Variable", values_to = "Value") %>% 
  dplyr::rename(label = group) %>% 
  select(label, Variable, Value)


ggplot(barnapboxplot, aes(y = Value, fill = label)) +
  geom_boxplot(coef = 2.5, outliers = FALSE, notch = TRUE) +
  facet_wrap(~ Variable, scales = "free_y")  + theme_bw() + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  scale_fill_manual(values = colours)
common_class <- intersect(barnapboxplot$label, classtree$tip.label)



barnapboxplot <- barnapboxplot %>%
  filter(label %in% common_class)



classtree$edge.length <- NULL
classtree <- ape::drop.tip(classtree, setdiff(classtree$tip.label, common_class))




p <- ggtree(classtree) + 
  geom_tiplab(size=3) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")


#replacement data frame error is due to order of columns! just need 'label' to be first
#need to annotate classes somehow

p <- p + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position = "none") + labs(color = NULL, fill = NULL) + scale_fill_brewer(palette = "Set1") 

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/ribocop/results/figures/itsboxplot.pdf", height = 8, width = 12)



p <- ggtree(classtree) + 
  geom_tiplab(size=3, offset = 0.5) +
  geom_tippoint(aes(color = label), show.legend = FALSE, size = 1.5) +
  scale_color_manual(values = colours) + theme(legend.position = "none")


#replacement data frame error is due to order of columns! just need 'label' to be first
#need to annotate classes somehow

p <- p + new_scale_color() + geom_facet(
  data = barnapboxplot %>% filter(Variable == "Unit length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "Unit Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS1"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS1 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) + geom_facet(
  data = barnapboxplot %>% filter(Variable == "ITS2"),
  mapping = aes(x = Value, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "ITS2 Length (bp)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +  geom_facet(
  data = barnapboxplot %>% filter(Variable == "IGS_length"),
  mapping = aes(x = Value/1000, group = interaction(label, Variable), fill = Variable),
  geom = geom_boxplot,       
  panel = "IGS Length (Kb)",
  notch = TRUE, outliers = FALSE, coef = 2.5
) +
  geom_facet(
    data = barnapboxplot %>% filter(Variable == "18S_length"),
    mapping = aes(x = Value/1000, fill = Variable, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "18S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  geom_facet(
    data = barnapboxplot %>% filter(Variable == "28S_length"),
    mapping = aes(x = Value/1000, fill = Variable, group = interaction(label, Variable)),
    geom = geom_boxplot,       
    panel = "28S Length (Kb)",
    notch = TRUE, outliers = FALSE, coef = 2.5
  )+
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 20), "Tree") +
  theme(legend.position =  "none")  + scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
  coord_cartesian(clip = "off")

facet_widths(p, widths = c(1, 1, 1))
ggsave("zc4688/honours/ribocop/results/figures/combinedboxplots.pdf", height = 8, width = 12)

barnapboxplot %>% 
  group_by(label, Variable) %>% 
  summarise(Median = median(Value, na.rm = TRUE), Min = min(Value, na.rm = TRUE),  Max = max(Value, na.rm = TRUE), Species = n()) %>% 
  view




unfinished <- bind_rows(data_list)

unfinished <- unfinished %>% 
  filter(is.na(`rDNA_details.Unit length`)) 
unfinished <- unfinished %>% rename("query" = "Sample_details.Species")
test <- classification(unfinished$Sample_details.Species, db="ncbi")

test <- cbind(test)
test <- test %>% 
  +    select(kingdom, phylum, subphylum, superclass, class, subclass, superorder, order, suborder, family, subfamily, genus, species, query)
write.table(test, "zc4688/honours/ribocop/unfinished.tsv", sep = "\t", quote = F, row.names = F)

test <- read.table("zc4688/honours/ribocop/unfinished.tsv", sep = "\t", header = TRUE)

test <- left_join(test, unfinished, by = "query")


test %>% select(class, query, `Errors.rDNA identification`, `Errors.Morph identification`) %>% 
  mutate(other = ifelse(is.na(`Errors.rDNA identification`) & is.na(`Errors.Morph identification`), "ambig", NA)) %>% 
  pivot_longer(cols = c("Errors.rDNA identification", "Errors.Morph identification", "other"), values_to = "blah", names_to = "error", values_drop_na = TRUE) %>% 
  select(-blah) %>% group_by(class, error) %>% summarise(n = n()) %>% view




##############Metadata table



metadata <- metadata %>% 
         left_join(barrnapmetadata)

metadata <- metadata %>% 
  rename("Species" = "label") %>% rename("rDNA contigs" = "Contigs")

write.table(metadata, "zc4688/honours/ribocop/metadata.tsv", sep = "\t", quote = F, row.names = F)

##############barrnap comparison

test <- read.delim("zc4688/honours/ribocop/results/chordata_filter2/morph_fasta/test.txt", header = F)

colnames(test) <- c("Sample", "Version", "rRNA", "Start", "End", "Envstart", "Envend", "Hmmstart", "Hmmend", "Score", "Strand", "Other", "Product")
test <- test %>% 
  filter(!is.na(Start)) %>% 
  separate(Product, into = c("Name", "Other", "Note"), sep = ";") %>% 
  separate(Name, into = c("Name", "Unit"), sep = "=") %>% 
  separate(Other, into = c("Other", "Completeness"), sep = "\\(") %>% 
  separate(Note, into = c("Note", "Palign"), sep = "=") %>% 
  mutate(Completeness = str_replace_all(Completeness, "\\)", "")) %>% 
  select(-Note, -Name, -Other, -Version, -rRNA, -Strand) 



test %>% filter(Sample == "CM066328.1" & Start > 137550003 & Start < 137588669) %>% ggplot() + geom_segment(aes(x = Envstart, xend = Envend, y = Sample, yend = Sample, color = Unit))



######################################## previous barrnap method
########################### Figure 2: Barrnap with edits (marked Palign for every unit)
barrnap <- read.delim("zc4688/honours/ribocop/results/main/morph_fasta/barrnap_newpalign.gff", header = F)
tree <- read.tree("/g/data/te53/zc4688/honours/ribocop/species (5).nwk") #ORIGINAL tree before filtering, just for checking. most current results are in refiltered directory (with current filtering), all has more because filtering was lest strict

treenames <- as.data.frame(tree$tip.label)
colnames(treenames) <- c("Timetreename")
replacementnames <- read.delim("zc4688/honours/ribocop/speciestimetree.tsv", header = T, sep = ",")

replacementnames <- replacementnames %>% mutate(Species = str_replace_all(Species, " ", "_"), Timetreename = str_replace_all(Replacement, " ", "_") ) %>% select(-Replacement)

treenames <- treenames %>% left_join(replacementnames) %>% mutate(Species = ifelse(is.na(Species), Timetreename, Species)) %>% select(Species)
tree$tip.label <- as.character(treenames$Species)

colnames(barrnap) <- c("Sample", "Version", "Start", "End", "Envstart", "Envend", "Hmmstart", "Hmmend", "Score", "Bitscore", "Hmmpalign", "Targetpalign", "Strand", "Other", "Product")

barrnap <- barrnap %>% 
  filter(!is.na(Start)) %>% 
  separate(Product, into = c("Name", "Other", "Note"), sep = ";") %>% 
  separate(Name, into = c("Name", "Unit"), sep = "=") %>% 
  #separate(Other, into = c("Other", "Completeness"), sep = "\\(") %>% 
  # separate(Note, into = c("Note", "Palign"), sep = "=") %>% 
  #mutate(Completeness = str_replace_all(Completeness, "\\)", "")) %>% 
  select(-Note, -Name, -Other, -Version, -Strand) %>% 
  mutate(`Sample ID` = sub("^([^_]*_[^_]*)_.*$", "\\1", Sample)) %>% 
  left_join(lengths) 


barrnap <- barrnap %>% group_by(Species) %>%
  filter(`Sample ID` == last(`Sample ID`)) %>%
  ungroup()

barrnap <- barrnap %>% 
  dplyr::rename("label" = "Species")

barrnapplot <- barrnap %>% relocate(label, .before = everything())

barrnapmetadata <- barrnapplot %>% 
  group_by(Unit, label, `Sample ID`) %>% 
  mutate(Hmmpalign = ifelse(Hmmpalign > 1, 1, Hmmpalign)) %>% 
  mutate(Hmmpalign = round(Hmmpalign, digits = 2)) %>% 
  filter(Hmmpalign == max((Hmmpalign))) %>% 
  arrange(Unit, Envstart) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(label, `Sample ID`, Unit, Envstart, Envend, `Unit length`) %>% pivot_wider(
    names_from = Unit,
    values_from = c(Envstart, Envend),
    names_glue = "{Unit}_{.value}"
  ) %>% 
  mutate(ITS1 = `5_8S_rRNA_Envstart` - `18S_rRNA_Envend`,
         `ITS2` = `28S_rRNA_Envstart` - `5_8S_rRNA_Envend`,
         `18S_length` = `18S_rRNA_Envend` - `18S_rRNA_Envstart`,
         `28S_length` = `28S_rRNA_Envend` - `28S_rRNA_Envstart`,
         `IGS_length` = `Unit length` - `28S_rRNA_Envend`) %>% 
  select(ITS1, ITS2, `18S_length`, `28S_length`, `IGS_length`, label, `Unit length`, `Sample ID`)

common_organisms <- intersect(barrnapplot$label, tree$tip.label)


barrnapplot <- barrnapplot %>%
  filter(label %in% common_organisms)





barrnapplot <- barrnapplot %>% mutate(label = str_replace_all(label, "_", " "))
barrnapplot <- left_join(barrnapplot, classifications, by="label")

#pruned_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$edge.length <- NULL
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, common_organisms))

tree$tip.label <- str_replace_all(tree$tip.label, "_", " ")




barrnapplot <- barrnapplot %>% mutate(Hmmpalign = Hmmpalign * 100)




tmp <- barrnapplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
)  +  geom_facet(
  data = barrnapplot %>% filter(`Unit length` < 35000),
  mapping = aes(x = `Unit length`, xend = 35000),
  geom = geom_segment,       # Horizontal bar chart
  panel = " ",
  color = "grey80",
  alpha = 0.5
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree.pdf", height = 20, width = 15)

tmp <- barrnapplot %>% select(label, group)

p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

p <- p + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = " "
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1)) + coord_cartesian(clip = 'off')

p

ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_nolength.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_nolength.pdf", height = 20, width = 15)


p <- ggtree(tree) %<+% tmp +
  geom_tiplab(size=1) +
  geom_tippoint(aes(color = group), size = 1.5) +
  scale_color_manual(values = colours, name = "Taxonomic Group")

test <- barrnapplot %>% select(label, `Unit length`, group) %>% rename("testgroup" = "group")
p <- p  + new_scale_color() + geom_facet(
  data = barrnapplot,
  mapping = aes(x = Envstart, xend = Envend, color = Unit, alpha = as.numeric(Hmmpalign)),
  geom = geom_segment,       # Horizontal bar chart
  panel = "Unit Structure"
)  + geom_facet(
  data = test,
  mapping = aes(x = `Unit length`/1000, fill = testgroup),
  geom = geom_col,       # Horizontal bar chart
  panel = "Unit length (Kb)",
  width = 0.6 
) +
  scale_y_discrete() +  
  theme_tree2() + xlim_expand(c(0, 40), "Tree")+labs(alpha = "% HMM Match") + scale_fill_manual(values = colours) + scale_alpha(range = c(0.1, 1))+
  scale_color_manual(values = c("18S_rRNA" = "#F75F86", "28S_rRNA" = "cornflowerblue", "5_8S_rRNA" = "#50C878", "5S_rRNA" = "purple"), name = "Unit")+
  guides(
    color = guide_legend(order = 1),
    fill = guide_none()) + coord_cartesian(clip = 'off')

facet_widths(p, widths = c(3, 3, 1))


ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_lengths.png", height = 20, width = 15)
ggsave("zc4688/honours/ribocop/results/figures/barrnaptree_lengths.pdf", height = 20, width = 15)

############################### Separate vs broken alignments
test <- barrnapplot %>% 
  arrange(Sample, Start) %>% 
  group_by(label, Sample, Unit) %>% 
  mutate(query_gap = Envstart - lag(Envend, default = first(Envstart)), 
         ref_gap = Hmmstart - lag(Hmmend, default = first(Hmmstart)), 
         diff = query_gap - ref_gap, 
         type = case_when((ref_gap < -100 & query_gap > 100) | query_gap > 5000 ~ "separate",
                          ref_gap == 0 & query_gap == 0 ~ "single", 
                          TRUE ~ "gap")) %>% 
  select(label, Sample, Unit, Envstart, Envend, type, `Unit length`, group) %>% 
  mutate(
    merged_Envstart = if_else(type == "gap", lag(Envstart, default = first(Envstart)), Envstart),
    merged_Envend = Envend)  %>% 
  filter(!lead(type, default = "none") == "gap") %>%
  ungroup() %>% 
  filter(Unit != "5S_rRNA") %>% 
  mutate(Envstart = merged_Envstart,
         Envend = merged_Envend)
