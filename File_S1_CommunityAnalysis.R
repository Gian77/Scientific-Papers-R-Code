# ************ DATA ANALYSIS UPARSE ********************************** ----------------------------
# Project name: Tuber magantum project Giorgio
# Manuscript:   
# Authors:      
# Affiliation:  University of Perugia - Michigan State University
# Journal:     
# Date:         May 24, 2019
# ******************************************************************** ----------------------------

# WORKING ENVIRONMENT SETUP -----------------------------------------------------------------------
library(styler) # hilight the code and press CTRL+SHIFT+A to style the code.
options(scipen = 9999) #to use decimals
options(max.print=100000000) # to print more lines on the display
# rm(list= ls()) # remove all objects

# loading required packages -----------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(Biostrings)
library(ape)
library(ggplot2)

library(vegan); packageVersion("vegan")
library(dplyr)
library(indicspecies)

library(multcompView)
library(ggpubr)
library(grid)
library(gridExtra)

# Color Palettes ----------------------------------------------------------------------------------
palette_site = c("#332288", "#88CCEE", "#117733")
pie(rep(1, length(palette_site)), labels = sprintf("%d (%s)",
     seq_along(palette_site),palette_site), col = palette_site)

palette_year = c("#FF899d", "#663F05") #,  "#E69F00", 
pie(rep(1, length(palette_year)), labels = sprintf("%d (%s)",
    seq_along(palette_year),palette_year), col = palette_year)


# ******************************************************************** ----------------------------
# IMPORTING DATASETS ------------------------------------------------------------------------------
# A) uparse ITS R1 --------------------------------------------------------------------------------
# Importing taxonomies at 0.6 confidence ----------------------------------------------------------
taxonomy_ITS06 <-
  read.delim(
    "taxonomies/consensus_taxonomy_ITS06.txt",
    header = TRUE,
    row.names = 1,
    sep = "\t")

head(taxonomy_ITS06)

# Importing taxonomies at 0.8 confidence ----------------------------------------------------------
taxonomy_ITS08 <-
  read.delim(
    "taxonomies/consensus_taxonomy_ITS08.txt",
    header = TRUE,
    row.names = 1, 
    sep = "\t")

head(taxonomy_ITS08)

# check for identical ordering
identical(rownames(taxonomy_ITS06), rownames(taxonomy_ITS08))
# sample_order <- match(rownames(taxonomy_ITS08), rownames(taxonomy_ITS06))
# taxonomy_ITS06 <- taxonomy_ITS06[,sample_order]

taxonomy_ITS08$Kingdom_06 <- taxonomy_ITS06$Kingdom
head(taxonomy_ITS08)
dim(taxonomy_ITS08)

levels(taxonomy_ITS08$Kingdom)

nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom!="Fungi",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Choanoflagellozoa",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Metazoa",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Stramenopila",])
nrow(taxonomy_ITS08[taxonomy_ITS08$Kingdom=="Rhizaria",])

# removing non-fungal OTUs
subset(taxonomy_ITS08, taxonomy_ITS08$Kingdom_06 == "Fungi") -> taxonomy_ITS08_filt
dim(taxonomy_ITS08_filt)
str(taxonomy_ITS08_filt)

# checking filtering results
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom != "Fungi", ])
nrow(taxonomy_ITS08_filt[taxonomy_ITS08_filt$Kingdom_06 != "Fungi", ])

# classify all OTUs at 60% identity to Fungi
taxonomy_ITS08_filt$Kingdom <- taxonomy_ITS08_filt$Kingdom_06

# formatting ranks
taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, 1:9]
taxonomy_ITS08_filt$OTU_ID <- rownames(taxonomy_ITS08_filt)

taxonomy_ITS08_filt <- taxonomy_ITS08_filt[, c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species",
  "OTU_ID",
  "Isolate",
  "Isolate_percent_id"
)]

head(taxonomy_ITS08_filt)

# importing other tables 
otus_ITS_uparse_R1 <-
  read.delim(
    "otu_tables/otu_table_ITS_UPARSE_R1.txt",
    header = TRUE,
    row.names = 1,)

metadata_ITS_uparse_R1 <-
  read.delim(
    "metadata/map_new_ITS.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t")

otus_seq_ITS_uparse_R1 <-
  readDNAStringSet(
    "fasta/otus_ITS_R1.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE)

physeq_ITS_uparse <- phyloseq(otu_table(otus_ITS_uparse_R1, taxa_are_rows = TRUE),
                                     sample_data(metadata_ITS_uparse_R1),
                                     tax_table(as.matrix(taxonomy_ITS08_filt)),
                                     otus_seq_ITS_uparse_R1) 

physeq_ITS_uparse
str(physeq_ITS_uparse)
head(sample_data(physeq_ITS_uparse))
tax_table(physeq_ITS_uparse)[tax_table(physeq_ITS_uparse)==""]<- NA
head(tax_table(physeq_ITS_uparse))

# checking the phyloseq object
sort(unique(as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom)) # everything is Fungi
nrow(as.data.frame(tax_table(physeq_ITS_uparse))[as.data.frame(tax_table(physeq_ITS_uparse))$Kingdom!="Fungi",])

# ******************************************************************** ----------------------------
# B) uparse 16S PAIRED ----------------------------------------------------------------------------

# # Importing RDP taxonomy ------------------------------------------------------------------------
# rdp_taxonomy_16s <-
#   read.delim(
#     "taxonomies/rdp_taxonomy_16s.txt",
#     header = TRUE,
#     row.names = 1)
# 
# head(rdp_taxonomy_16s)
# 
# ifelse(taxonomy_16s_uparse_RDP$D_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Domain), NA) -> Kingdom
# ifelse(taxonomy_16s_uparse_RDP$P_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Phylum), NA) -> Phylum
# ifelse(taxonomy_16s_uparse_RDP$C_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Class), NA) -> Class
# ifelse(taxonomy_16s_uparse_RDP$O_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Order), NA) -> Order
# ifelse(taxonomy_16s_uparse_RDP$F_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Family), NA) -> Family
# ifelse(taxonomy_16s_uparse_RDP$G_Score>=0.8, paste(taxonomy_16s_uparse_RDP$Genus), NA) -> Genus
# taxonomy <- cbind(Kingdom, Phylum, Class, Order, Family, Genus)
# rownames(taxonomy) <- rownames(taxonomy_16s_uparse_RDP)
# taxonomy_phy_16s_uparse_RDP <- tax_table(as.matrix(taxonomy))
# head(taxonomy_phy_16s_uparse_RDP)

# Importing CONSTAX SILVA taxonomy ----------------------------------------------------------------
silva_taxonomy_16s <-
  read.delim("taxonomies/consensus_taxonomy_16s.txt",
             header = TRUE,
             row.names = 1)

head(silva_taxonomy_16s)

# cleaning taxonomy labels ------------------------------------------------------------------------
colnames(silva_taxonomy_16s) <-
  c("Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "OTU_ID",
    "Isolate",
    "Isolate_percent_id")

silva_taxonomy_16s$OTU_ID <- rownames(silva_taxonomy_16s)

silva_taxonomy_16s[, "Kingdom"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Kingdom"]))
silva_taxonomy_16s[, "Phylum"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Phylum"]))
silva_taxonomy_16s[, "Class"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Class"]))
silva_taxonomy_16s[, "Order"] <- as.factor(gsub("_1", "",silva_taxonomy_16s[, "Order"]))
silva_taxonomy_16s[, "Family"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Family"]))
silva_taxonomy_16s[, "Genus"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Genus"]))
silva_taxonomy_16s[, "Species"] <- as.factor(gsub("_1", "", silva_taxonomy_16s[, "Species"]))

head(silva_taxonomy_16s)
# silva_taxonomy_16s[1:50,]
# endsWith("aceae", as.character(silva_taxonomy_16s$Family))

any(silva_taxonomy_16s$Kingdom == "Chloroplast" | silva_taxonomy_16s$Kingdom == "Mitochondria")
any(silva_taxonomy_16s$Phylum == "Chloroplast" | silva_taxonomy_16s$Phylum == "Mitochondria")
any(silva_taxonomy_16s$Class == "Chloroplast" | silva_taxonomy_16s$Class == "Mitochondria")
any(silva_taxonomy_16s$Family == "Chloroplast" | silva_taxonomy_16s$Family == "Mitochondria")
any(silva_taxonomy_16s$Genus == "Chloroplast" | silva_taxonomy_16s$Genus == "Mitochondria")

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Family == "Chloroplast")
silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Family == "Mitochondria")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Family != "Mitochondria")

# Check for unclassified OTUs
any(silva_taxonomy_16s$Kingdom == "")

silva_taxonomy_16s %>% dplyr::filter(silva_taxonomy_16s$Kingdom == "")
silva_taxonomy_16s <- subset(silva_taxonomy_16s, Kingdom != "")

# importing other tables
otus_16s_uparse <-
  read.delim("otu_tables/otu_table_16s_UPARSE.txt",
             row.names = 1,
             header = TRUE)

metadata_16s_uparse <-
  read.delim(
    "metadata/map_new_16s.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t"
  )

otus_seq_16s_uparse <-
  readDNAStringSet(
    "fasta/otus_16s.fasta",
    format = "fasta",
    seek.first.rec = TRUE,
    use.names = TRUE
  )


physeq_16s_uparse <- phyloseq(otu_table(otus_16s_uparse, taxa_are_rows = TRUE), 
                              sample_data(metadata_16s_uparse),
                              tax_table(as.matrix(silva_taxonomy_16s)),
                              otus_seq_16s_uparse) 

physeq_16s_uparse
head(sample_data(physeq_16s_uparse))
tax_table(physeq_16s_uparse)[tax_table(physeq_16s_uparse)==""]<- NA
head(tax_table(physeq_16s_uparse))

# ******************************************************************** ----------------------------
# CREATIN NEW PHYLOSEQ OBJETCS --------------------------------------------------------------------
physeq_ITS_uparse -> physeq_fungi_uparse
physeq_16s_uparse -> physeq_prok_uparse

# REMOVING SAMPLES non included in this study -----------------------------------------------------
sample_data(physeq_fungi_uparse)
physeq_fungi_uparse
sample_data(physeq_prok_uparse)
physeq_prok_uparse

saveRDS(physeq_fungi_uparse, "physeq_fungi_uparse_prefiltering.rds")
saveRDS(physeq_prok_uparse, "physeq_prok_uparse_prefiltering.rds")


physeq_fungi_uparse <-
  subset_samples(
    physeq_fungi_uparse,
    Description != "2_100_2015" &
      Description != "2_T_2015" &
      Description != "3_100_2015" &
      Description != "3_T_2015" &
      Description != "4_100_2015" &
      Description != "4_T_2015" &
      Description != "5_100_2015" &
      Description != "5_T_2015" &
      Description != "6_100_2015" &
      Description != "6_T_2015"  &
      Description != "7_100_2015" &
      Description != "7_T_2015" &
      Description != "8_100_2015" &
      Description != "8_T_2015"
  )
otu_table(physeq_fungi_uparse) <-
  otu_table(physeq_fungi_uparse)[which(rowSums(otu_table(physeq_fungi_uparse)) > 0), ]

physeq_fungi_uparse

physeq_prok_uparse <-
  subset_samples(
    physeq_prok_uparse,
    Description != "2_100_2015" &
      Description != "2_T_2015" &
      Description != "3_100_2015" &
      Description != "3_T_2015" &
      Description != "4_100_2015" &
      Description != "4_T_2015" &
      Description != "5_100_2015" &
      Description != "5_T_2015" &
      Description != "6_100_2015" &
      Description != "6_T_2015"  &
      Description != "7_100_2015" &
      Description != "7_T_2015" &
      Description != "8_100_2015" &
      Description != "8_T_2015"
  )
otu_table(physeq_prok_uparse) <-
  otu_table(physeq_prok_uparse)[which(rowSums(otu_table(physeq_prok_uparse)) > 0), ]

physeq_prok_uparse

# ******************************************************************** ----------------------------
# FILTERING OUT CONTAMINANTS ----------------------------------------------------------------------
library(decontam)

# detecting contaminants Fungi
sample_data(physeq_fungi_uparse)$is.neg <-
  sample_data(physeq_fungi_uparse)$Site == "control"

contam_prev_fungi <-
  isContaminant(
    physeq_fungi_uparse,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_fungi$contaminant)

# detecting contaminants Prokaryotes
sample_data(physeq_prok_uparse)$is.neg <-
  sample_data(physeq_prok_uparse)$Site == "control"
contam_prev_prokaryote <-
  isContaminant(
    physeq_prok_uparse,
    method = "prevalence",
    neg = "is.neg",
    threshold = 0.5)

table(contam_prev_prokaryote$contaminant)

# INSPECTING LIBRARY SIZES ------------------------------------------------------------------------
# Before removing controls
df_fungi_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_fungi_uparse)))
df_fungi_uparse$LibrarySize <- sample_sums(physeq_fungi_uparse)
df_fungi_uparse <-
  df_fungi_uparse[order(df_fungi_uparse$LibrarySize), ]
df_fungi_uparse$Index <-
  seq(nrow(df_fungi_uparse)) # sample numbering
df_fungi_uparse

# reorder the factor
df_fungi_uparse$Site <- factor(
  df_fungi_uparse$Site,
  levels = c(
    "Citta_di_Castello",
    "Ripabianca",
    "San_Giovanni_dAsso",
    "control"))

df_prok_uparse <-
  as.data.frame(as.matrix(sample_data(physeq_prok_uparse)))
df_prok_uparse$LibrarySize <- sample_sums(physeq_prok_uparse)
df_prok_uparse <- df_prok_uparse[order(df_prok_uparse$LibrarySize), ]
df_prok_uparse$Index <- seq(nrow(df_prok_uparse)) # sample numbering
df_prok_uparse

df_prok_uparse$Site <- factor(
  df_prok_uparse$Site,
  levels = c(
    "Citta_di_Castello",
    "Ripabianca",
    "San_Giovanni_dAsso",
    "control"))

# *** FIGURES S1 - Distribution of Sample Libraries -----------------------------------------------
ggarrange(
ggplot(data=df_fungi_uparse, aes(x=Index, y=LibrarySize, color=is.neg)) +
  geom_point(alpha =0.7) +
  labs(title="Fungi", x="Sample number", y="Read number") + 
  theme_classic() +
      scale_colour_manual("Negative control", values = c("grey", "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)),
ggplot(data=df_prok_uparse, aes(x=Index, y=LibrarySize, color=is.neg)) +
  geom_point(alpha =0.7) +
  labs(title="Prokaryotes", x="Sample number", y="Read number") + 
  theme_classic() +
      scale_colour_manual("Negative control", values = c("grey", "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size = 8)),
ggplot(df_fungi_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
  geom_histogram(color = "indianred", fill = "indianred", binwidth = 1000) +
  theme_classic() +
  #facet_grid(~Treatment, scales = "free_x", space="free_x") +
  labs(title=NULL, x="Read number", y="Sample number") +
  theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)),
ggplot(df_prok_uparse, aes(x = LibrarySize)) + # Histogram of sample read counts
  geom_histogram(color = "indianred", fill = "indianred", binwidth = 1000) +
  theme_classic() +
  #facet_grid(~Treatment, scales = "free_x", space="free_x") +
  labs(title=NULL, x="Read number", y="Sample number") + 
  theme(axis.text.x = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5)),
labels = c("A", "B", "C", "D"),
widths = c(1,1,1,1),
align = "v", ncol = 2, nrow = 2,
common.legend = TRUE,
legend = "bottom") -> lib_size_all

lib_size_all


# scale_colour_manual(values = c(paletteCB3_fungi, "red"), 
#                     labels = c(Citta_di_Castello="Citta' di Castello", 
#                                San_Giovanni_dAsso = "San Giovanni", 
#                                Ripabiance = "Ripabianca",
#                                control = "control")) +

# removing contaminants form the phyloseq object - Fungi ------------------------------------------
remove_taxa = function(physeq, badTaxa) {
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

contaminants_fungi = rownames(subset(contam_prev_fungi, contaminant %in% c("TRUE")))
physeq_fungi_uparse_clean <-
  remove_taxa(physeq_fungi_uparse, contaminants_fungi)
physeq_fungi_uparse_clean
sample_data(physeq_fungi_uparse_clean)

# Now removing controls samples
physeq_fungi_uparse_clean <-
  subset_samples(physeq_fungi_uparse_clean, Type %in% c("truffle", "soil"))
otu_table(physeq_fungi_uparse_clean) <-
  otu_table(physeq_fungi_uparse_clean)[which(rowSums(otu_table(physeq_fungi_uparse_clean)) > 0), ]
physeq_fungi_uparse_clean
levels(sample_data(physeq_fungi_uparse_clean)$Site)

any(taxa_sums(physeq_fungi_uparse_clean) == 0)
any(sample_sums(physeq_fungi_uparse_clean) == 0)
sort(sample_sums(physeq_fungi_uparse_clean))

# removing contaminants form the phyloseq object - Prokaryotes ------------------------------------
contaminants_prokaryote = rownames(subset(contam_prev_prokaryote, contaminant %in%
                                            c("TRUE")))
length(contaminants_prokaryote)

physeq_prok_uparse_clean <-
  remove_taxa(physeq_prok_uparse, contaminants_prokaryote)
physeq_prok_uparse_clean

# Now removing controls samples
physeq_prok_uparse_clean <-
  subset_samples(physeq_prok_uparse_clean, Type %in% c("truffle", "soil"))
otu_table(physeq_prok_uparse_clean) <-
  otu_table(physeq_prok_uparse_clean)[which(rowSums(otu_table(physeq_prok_uparse_clean)) > 0), ]
physeq_prok_uparse_clean
levels(sample_data(physeq_prok_uparse_clean)$Site)

any(taxa_sums(physeq_prok_uparse_clean) == 0)
any(sample_sums(physeq_prok_uparse_clean) == 0)
sort(sample_sums(physeq_prok_uparse_clean))

# FILTERING DATASTES following quality  -----------------------------------------------------------

# Oliver et al. 2015, PCR errors and tag switching
# Barberan et al. 2012, removing OTUs that appear in less than x samples
# Lindhal et al. 2013, tag switching - that's a good  one!

physeq_fungi_uparse_qc <- physeq_fungi_uparse_clean
otu_table(physeq_fungi_uparse_qc) <-
  otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) >= 10), ] ### PCR errors
otu_table(physeq_fungi_uparse_qc)[otu_table(physeq_fungi_uparse_qc) < 4] <- 0 ### tag switching
otu_table(physeq_fungi_uparse_qc) <-
  otu_table(physeq_fungi_uparse_qc)[which(rowSums(otu_table(physeq_fungi_uparse_qc)) > 0), ]

physeq_fungi_uparse_qc

physeq_prok_uparse_qc <- physeq_prok_uparse_clean
otu_table(physeq_prok_uparse_qc) <-
  otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) >= 10), ]
otu_table(physeq_prok_uparse_qc)[otu_table(physeq_prok_uparse_qc) < 4] <- 0 ### tag switching
otu_table(physeq_prok_uparse_qc) <-
  otu_table(physeq_prok_uparse_qc)[which(rowSums(otu_table(physeq_prok_uparse_qc)) > 0), ]

physeq_prok_uparse_qc

# Change OTU names --------------------------------------------------------------------------------
ChangeName <- function(physeq, org){
  otu_names<- rownames(otu_table(physeq))
  if (startsWith(otu_names, org)[1] == TRUE) {
    stop("OTUs are already renamed")
  } else {
    otu_names<- paste(org, otu_names, sep="")
    taxa_names(physeq) <- paste(otu_names)
    return(physeq)
  }
}

ChangeName(physeq_fungi_uparse_qc, "F") -> physeq_fungi_uparse_qc
ChangeName(physeq_prok_uparse_qc, "B") -> physeq_prok_uparse_qc

# EXTRACT LAS TAXONOMIC LEVEL ---------------------------------------------------------------------
# thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R

source("../R_functions/ReformatTaxonomy.R")

head(tax_table(physeq_fungi_uparse_qc))
ReformatTaxonomy(physeq_fungi_uparse_qc) -> physeq_fungi_uparse_filt
head(tax_table(physeq_fungi_uparse_filt))

head(tax_table(physeq_prok_uparse_qc))
ReformatTaxonomy(physeq_prok_uparse_qc) -> physeq_prok_uparse_filt
head(tax_table(physeq_prok_uparse_filt))


# RAREFACTION CURVES ------------------------------------------------------------------------------
# Sample40 has zero reads that's why the warning!
# converting to table objetcs ---------------------------------------------------------------------
sample_sums(physeq_fungi_uparse_filt) >0 
subset_samples(physeq_fungi_uparse_filt, sample_sums(physeq_fungi_uparse_filt) >0) -> physeq_fungi_uparse_filt

otu_fungi_up <- as.data.frame(otu_table(physeq_fungi_uparse_filt))
taxa_fungi_up <- as.data.frame(as.matrix(tax_table(physeq_fungi_uparse_filt)))
metadata_fungi_up <- as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))
dim(otu_fungi_up)
metadata_fungi_up

sample_sums(physeq_prok_uparse_filt) >0 
subset_samples(physeq_prok_uparse_filt, sample_sums(physeq_prok_uparse_filt) >0) -> physeq_prok_uparse_filt

otu_prok_up<- as.data.frame(otu_table(physeq_prok_uparse_filt))
taxa_prok_up <- as.data.frame(as.matrix(tax_table(physeq_prok_uparse_filt)))
metadata_prok_up <- as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))
dim(metadata_prok_up)
metadata_prok_up

# *** FIGURE S2 - Rarefaction curves ******--------------------------------------------------------

# creating and empty character 
curve_colors_fungi <- vector(mode="character", length=ncol(as.data.frame(otu_table(physeq_fungi_uparse_filt))))
curve_colors_fungi[as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))$Type=="truffle"] <- "gold"
curve_colors_fungi[as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt)))$Type=="soil"] <- "brown"

curve_colors_prok <- vector(mode="character", length=ncol(as.data.frame(otu_table(physeq_prok_uparse_filt))))
curve_colors_prok[as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))$Type=="truffle"] <- "gold"
curve_colors_prok[as.data.frame(as.matrix(sample_data(physeq_prok_uparse_filt)))$Type=="soil"] <- "brown"


par(mar = c(3.5, 3.5, 2, 0) + 0.1,mgp = c(2.4, 0.8, 0),las = 0,mfrow = c(1, 2))
# Change the panel layout to 2 x 1
rarecurve(t(otu_fungi_up), col = curve_colors_fungi, label = FALSE, step = 1000,
          main="Fungi", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("topright", legend=c("Truffle", "Soil"),
       col=c("gold", "brown"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
# if you rarefy you can add at what depth you rarefy here 
#abline(v = min(sample_sums(physeq_fungi_uparse_qc)), col="red", lwd=3, lty=2)
rarecurve(t(otu_prok_up), col = curve_colors_prok, label = FALSE,  step = 1000,
          main="Prokaryotes", ylab = "Number of OTUs", xlab = "Number of DNA reads")
legend("topright", legend=c("Truffle", "Soil"),
       col=c("gold", "brown"), lty=2, cex=0.8, box.lty=0,lwd=3,bty = "n") 
#abline(v = min(sample_sums(physeq_prok_uparse_qc)), col="red", lwd=3, lty=2)
dev.off()

# ******************************************************************** ----------------------------
# NEW DATASETS REDY TO GO -------------------------------------------------------------------------
saveRDS(physeq_fungi_uparse_filt, "physeq_fungi_uparse.rds")
saveRDS(physeq_prok_uparse_filt, "physeq_prok_uparse.rds")


# ******************************************************************** ----------------------------
# ******************************************************************** ----------------------------
# ******************************************************************** ----------------------------
# TESITNG CSS TRANSFORMATION ----------------------------------------------------------------------
library(metagenomeSeq)

CSSNorm <-function(dataframe){
  require(metagenomeSeq)
  dataframe %>% 
    phyloseq_to_metagenomeSeq() -> physeq_CSS
  p_biom <-cumNormStatFast(physeq_CSS)
  biom_quant <-cumNorm(physeq_CSS, p=p_biom)
  physeq_CSS <- MRcounts(biom_quant, norm=T)
  physeq_mSeq <- dataframe
  otu_table(physeq_mSeq) <- otu_table(physeq_CSS, taxa_are_rows=TRUE)
  return(physeq_mSeq)
}


# LET"S TRY USING METAGENOMESED INSTESD OF RAREFACTION
# MetagenomeSeq normalization - Gaussian model ----------------------------------------------------
CSSNorm(physeq_fungi_uparse_filt) -> physeq_fungi_uparse_mSeq
head(otu_table(physeq_fungi_uparse_mSeq))

CSSNorm(physeq_prok_uparse_filt) -> physeq_prok_uparse_mSeq
head(otu_table(physeq_prok_uparse))


# >>> ALPHA DIVERSITY --------------------------------------------------------
library(multcompView)
library(ggpubr)

# creating an additional variable for filtering
ExtrDataFrame <- function(physeq, var){
sample_data(physeq)$LocYear <- paste(
                               sample_data(physeq)$Site,
                               sample_data(physeq)$Year,
                               sep = "_")
if (var == "Site"){
physeq %>% subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                                       "San_Giovanni_dAsso_2015",
                                       "Ripabianca_2014")) -> physeq_filt
}else{
physeq %>% subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                                        "Citta_di_Castello_2015")) -> physeq_filt
}
otu_table(physeq_filt) <- otu_table(physeq_filt)[which(rowSums(otu_table(physeq_filt)) > 0),] 
# adding aplha metrics to metadata file
df_alpha <- as.data.frame(as.matrix(sample_data(physeq_filt))) 
df_alpha$Observed <- specnumber(as.data.frame(otu_table(physeq_filt)), MARGIN = 2)
df_alpha$Shannon <- diversity(as.data.frame(otu_table(physeq_filt)), index="shannon", MARGIN = 2)
# re-order factor levels
df_alpha$Distance <- factor(df_alpha$Distance,levels=c("Inside", "0", "40", "100", "200"))
return(df_alpha)
}


ExtrDataFrame(physeq_fungi_uparse_filt, "Site") -> df_alpha_fungi
df_alpha_fungi
ExtrDataFrame(physeq_prok_uparse_filt, "Site") -> df_alpha_prok
df_alpha_prok

ExtrDataFrame(physeq_fungi_uparse_filt, "Year") -> df_alpha_fungi_yr
df_alpha_fungi_yr
ExtrDataFrame(physeq_prok_uparse_filt, "Year") -> df_alpha_prok_yr
df_alpha_prok_yr


# count variables levels
df_alpha_fungi %>% group_by(Site) %>% count(vars = Distance)
df_alpha_fungi %>% group_by(Distance) %>% count(vars = Site)
df_alpha_prok %>% group_by(Site) %>% count(vars = Distance)
df_alpha_prok %>% group_by(Distance) %>% count(vars = Site)


# comparing Distance across sites using Wilcox.test
# Since we have done many independent tests (one for each taxon), 
# we should correct for multiple comparisions. However if we use p.adj,
# even even Inside form 0 is not significant 
CompSampl <- function(dataframe, formula, level, var){
	require(ggpubr)
  if (var == "Site"){
    compare_means(formula, data = dataframe[dataframe$Site==level,],
                  p.adjust.method = "bonferroni") -> test_CC 
        test_CC <- as.data.frame(test_CC)[,c(2,3,4)] # for p.adj 4 to 5
        test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
          colnames(test_CC2) <- c("group1", "group2", "p") # change p to p.adj
    rbind(test_CC, test_CC2) -> test_all
        as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
        data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
          res_CC$sample <- rownames(res_CC)
          res_CC %>% slice(match(c("Inside", "0","40","100","200"), sample)) -> res_CC
return(res_CC)
  }else{
    compare_means(formula, data = dataframe[dataframe$Distance==level,],
                  p.adjust.method = "bonferroni") -> test_CC 
    test_CC <- as.data.frame(test_CC)[,c(2,3,4)] 
    test_CC2 <- data.frame(test_CC[,2], test_CC[,1], test_CC[,3])
    colnames(test_CC2) <- c("group1", "group2", "p") 
    rbind(test_CC, test_CC2) -> test_all
    as.dist(xtabs(test_all[, 3] ~ (test_all[, 2] + test_all[, 1])), diag = TRUE) -> dist_CC
    data.frame(multcompLetters(dist_CC)['Letters']) -> res_CC
    res_CC$sample <- rownames(res_CC)
    res_CC %>% slice(match(c("2014", "2015"), sample)) -> res_CC
    return(res_CC)
  }
}


# fungi -------------------------------------------------------------------------------------------
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "Citta_di_Castello", "Site") -> Wtest_fungi_Cast_rich
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "Ripabianca", "Site") -> Wtest_fungi_Rip_rich
CompSampl(df_alpha_fungi, formula(Observed ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_fungi_Gio_rich

my_lables_fungi_rich <-c(as.character(Wtest_fungi_Cast_rich$Letters), 
                        as.character(Wtest_fungi_Rip_rich$Letters),
                        as.character(Wtest_fungi_Gio_rich$Letters))
my_lables_fungi_rich

CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "Citta_di_Castello", "Site") -> Wtest_fungi_Cast_shan
CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "Ripabianca", "Site") -> Wtest_fungi_Rip_shan
CompSampl(df_alpha_fungi, formula(Shannon ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_fungi_Gio_shan

my_lables_fungi_shan <-c(as.character(Wtest_fungi_Cast_shan$Letters), 
                         as.character(Wtest_fungi_Rip_shan$Letters),
                         as.character(Wtest_fungi_Gio_shan$Letters))
my_lables_fungi_shan




# prok --------------------------------------------------------------------------------------------
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "Citta_di_Castello", "Site") -> Wtest_prok_Cast_rich
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "Ripabianca", "Site") -> Wtest_prok_Rip_rich
CompSampl(df_alpha_prok, formula(Observed ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_prok_Gio_rich

my_lables_prok_rich <-c(as.character(Wtest_prok_Cast_rich$Letters), 
                        as.character(Wtest_prok_Rip_rich$Letters),
                        as.character(Wtest_prok_Gio_rich$Letters))
my_lables_prok_rich


CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "Citta_di_Castello", "Site") -> Wtest_prok_Cast_shan
CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "Ripabianca", "Site") -> Wtest_prok_Rip_shan
CompSampl(df_alpha_prok, formula(Shannon ~ Distance), "San_Giovanni_dAsso", "Site") -> Wtest_prok_Gio_shan

my_lables_prok_shan <-c(as.character(Wtest_prok_Cast_shan$Letters), 
                        as.character(Wtest_prok_Rip_shan$Letters),
                        as.character(Wtest_prok_Gio_shan$Letters))
my_lables_prok_shan

# Alpha diversity metric plus anova letters ------------------------------------------------------
# library(agricolae)
# library(car)

# metadata_fungi_up_mSeq %>%subset(Site=="Citta_di_Castello") -> richness_fungi_CC
# aov(Observed ~ Distance, data=richness_fungi_CC) -> aov_richness_fungi_CC
# HSD.test(aov_richness_fungi_CC, "Distance", group=TRUE) -> HSD_richness_fungi_CC
# HSD_richness_fungi_CC
# 
# # test assumtions and visualize q-q- plot 
# # correlation between a given sample and the normal distribution
# shapiro.test(richness_fungi_CC$Observed)
# ggqqplot(richness_fungi_CC$Observed)
# leveneTest(richness_fungi_CC$Observed ~ richness_fungi_CC$Distance, center=mean)
# # Can't use anova, data are not normal and eteroschedastic!
# #rep(rev(HSD_richness_fungi_CC$groups$groups), times=3)


# Plot Richness -----------------------------------------------------------------------
PlotRich <- function(dataframe, var, palette, my_labels){
#calculating where to put the letters
max(dataframe[,var] + 0.1* max(dataframe[,var])) -> labels_y
ggplot(dataframe, aes(x = Distance, y = get(var), color = Site)) +
  geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
      stat_summary(fun=mean, geom="point", shape=18, size=1.5, color="red", fill="red") +
      stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=2.2, color="black") +
  facet_grid(~Site) +
   #         facet_grid(~Site, scales = "free_x", 
   #          labeller = as_labeller(c(Citta_di_Castello = "Site 1",
   #                                   Ripabianca = "Site 2",
   #                                   San_Giovanni_dAsso = "Site 3"))) +
  theme_classic() +
  #scale_colour_manual(values = palette, 
  #                    labels = c("Site 1", "Site 2", "Site 3")) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
  theme(strip.text.x = element_text(size = 8)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
  theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 0.9, vjust = 1)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
  grids(linetype = "dashed") +
  theme(legend.position="bottom") -> rich_plot 
return(rich_plot)
}

PlotRich(df_alpha_fungi, "Observed", palette_site, my_lables_fungi_rich) + 
  labs(x=NULL, y="Richness", title = "Fungi") -> plot_rich_fungi_rich
plot_rich_fungi_rich
PlotRich(df_alpha_prok, "Observed", palette_site, my_lables_prok_rich) + 
  labs(x=NULL, y=NULL, title = "Prokaryotes")  -> plot_rich_prok_rich
plot_rich_prok_rich


PlotRich(df_alpha_fungi, "Shannon", palette_site, my_lables_fungi_rich) + 
  labs(x="Distance", y="Shannon", title = NULL) -> plot_rich_fungi_shan
plot_rich_fungi_shan
PlotRich(df_alpha_prok,"Shannon", palette_site, my_lables_prok_rich) + 
  labs(x="Distance", y=NULL, title = NULL)  -> plot_rich_prok_shan
plot_rich_prok_shan



ggarrange(
df_alpha_fungi %>%
  subset(Site%in%c("Citta_di_Castello")) %>%
  ggplot(aes(Distance, Observed)) +
  geom_boxplot(),

df_alpha_fungi %>%
  subset(Site%in%c("Ripabianca")) %>%
  ggplot(aes(Distance, Observed)) +
  geom_boxplot(),

df_alpha_fungi %>%
  subset(Site%in%c("San_Giovanni_dAsso")) %>%
  ggplot(aes(Distance, Observed)) +
  geom_boxplot(),
ncol = 3)




# *** FIGURE 1A - Richness -------------------------------------------------------
Fig_1A <- ggarrange(plot_rich_fungi_rich, 
                    plot_rich_prok_rich,
                    plot_rich_fungi_shan, 
                    plot_rich_prok_shan,
                    labels = c("A", "B", "C", "D"),
                    widths = c(1,1,1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 2,
                    common.legend = TRUE,
                    legend = "none")

Fig_1A

title1=text_grob("Site",size=12, face=2)
grid.arrange(Fig_1A, top = title1) -> Fig_1A_alpha_site
Fig_1A_alpha_site

# We can combine these plots with ggarrange. Before that, Girogio, 
# try to make the shannon graphs too, or you can use evennes in vegan
# ALSO, VERY ALARMING, it looks like that the bacteria were more diverse
# in the truffle rather than in the soil. How can be that possible?
# let's check read number, maybe we should try to rarefy and re rurn the code
# on rarefied data.

# now the plotting between years
Calclabel <- function(df, formula){
  CompSampl(df, formula(Observed ~ Year), "Inside", "Year") -> pval_Inside
  CompSampl(df, formula(Observed ~ Year), "0", "Year") -> pval_0
  CompSampl(df, formula(Observed ~ Year), "40", "Year") -> pval_40
  CompSampl(df, formula(Observed ~ Year), "100", "Year") -> pval_100
  CompSampl(df, formula(Observed ~ Year), "200", "Year") -> pval_200
  my_lables <-c(as.character(pval_Inside$Letters), 
                as.character(pval_0$Letters),
                as.character(pval_40$Letters),
                as.character(pval_100$Letters),
                as.character(pval_200$Letters))
  return(my_lables)
}


Calclabel(df_alpha_fungi_yr, formula(Observed ~ Year)) -> Wtest_fungi_year_rich
Calclabel(df_alpha_fungi_yr, formula(Shannon ~ Year)) -> Wtest_fungi_year_shan

Calclabel(df_alpha_prok_yr, formula(Observed ~ Year)) -> Wtest_prok_year_rich
Calclabel(df_alpha_prok_yr, formula(Shannon ~ Year)) -> Wtest_prok_year_shan



#Plotting
PlotRichY <- function(dataframe, var, palette, my_labels){
  #calculating where to put the letters
  max(dataframe[,var] + 0.1* max(dataframe[,var])) -> labels_y
  ggplot(dataframe, aes(x = Year, y = get(var), color = Year)) +
    geom_boxplot(outlier.colour="black", outlier.size=1.2, outlier.shape = 8, lwd=0.5, width = 0.7) +
    stat_summary(fun=mean, geom="point", shape=18, size=1.5, color="red", fill="red") +
    stat_summary(geom = 'text', label = my_labels, fun= max, aes(y = labels_y), size=2.2, color="black") +
    facet_grid(~factor(dataframe$Distance, levels = c("Inside", "0", "40", "100", "200"))) +
    theme_classic() +
    scale_colour_manual(values = palette, 
                        labels = c("Inside", "0", "40", "100", "200")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) +
    theme(strip.text.x = element_text(size = 8)) +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_blank(), legend.text = element_text(size = 8)) +
    theme(axis.text.x = element_text(angle = 45, size = 7, hjust = 0.9, vjust = 1)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) + 
    grids(linetype = "dashed") +
    theme(legend.position="none") -> rich_plot 
  return(rich_plot)
}

PlotRichY(df_alpha_fungi_yr, "Observed", palette_year, Wtest_fungi_year_rich) + 
  labs(x=NULL, y="Richness", title = "Fungi") 

# **** FIGURE 2A richness year ----------------------------------------------------------------- 
Fig_1B <- ggarrange(PlotRichY(df_alpha_fungi_yr, "Observed", palette_year, Wtest_fungi_year_rich) + 
                      labs(x=NULL, y="Richness", title = "Fungi"),
                    PlotRichY(df_alpha_prok_yr, "Observed", palette_year, Wtest_prok_year_rich) + 
                      labs(x=NULL, y=NULL, title = "Prokaryotes"),
                    PlotRichY(df_alpha_fungi_yr, "Shannon", palette_year, Wtest_fungi_year_shan) + 
                      labs(x="Distance", y="Shannon", title = NULL),
                    PlotRichY(df_alpha_prok_yr,"Shannon", palette_year, Wtest_prok_year_shan) + 
                      labs(x="Distance", y=NULL, title = NULL),
                    labels = c("A", "B", "C", "D"),
                    widths = c(1,1,1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 2,
                    common.legend = TRUE,
                    legend = "none")

Fig_1B

title2=text_grob("Year",size=12, face=2)
grid.arrange(Fig_1B, top = title2) -> Fig_1B_year
Fig_1B_year


# *****************************************************************************************--------
# >>> BETA DIVERSITY ------------------------------------------------------------------------------
library(ggpubr)

# >>> PERMANOVA -----------------------------------------------------------------------------------
count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Site)

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Distance)

ExtractOTU <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "San_Giovanni_dAsso_2015",
                                               "Ripabianca_2014")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "San_Giovanni_dAsso_2015",
                                               "Ripabianca_2014")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}

ExtractOTU(physeq_fungi_uparse_mSeq) -> otu_fungi_uparse_mSeq
ExtractMETA(physeq_fungi_uparse_mSeq) -> metadata_fungi_uparse_mSeq

ExtractOTU(physeq_prok_uparse_mSeq) -> otu_prok_uparse_mSeq
ExtractMETA(physeq_prok_uparse_mSeq) -> metadata_prok_uparse_mSeq

# test the effect Site 
adonis_fungi <- adonis(t(otu_fungi_uparse_mSeq) ~ Site,
                       metadata_fungi_uparse_mSeq, 
                       strata=metadata_fungi_uparse_mSeq$Distance, 
                       permutations=9999)
adonis_fungi

# test for the effect of Distance 
adonis_fungi1 <- adonis(t(otu_fungi_uparse_mSeq) ~ Distance,  
                        metadata_fungi_uparse_mSeq, 
                        strata=metadata_fungi_uparse_mSeq$Site, 
                        permutations=9999)
adonis_fungi1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_fungi2 <- adonis(t(otu_fungi_uparse_mSeq) ~ Distance * Site,  
                        metadata_fungi_uparse_mSeq, 
                        permutations=9999)
adonis_fungi2

adonis_fungi3 <- adonis(t(otu_fungi_uparse_mSeq) ~ Site * Distance,
                        metadata_fungi_uparse_mSeq, 
                        permutations=9999)
adonis_fungi3
round(p.adjust(adonis_fungi3$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0027 1.0000 1.0000     NA     NA


# *** TABLE 1 - PERMANOVA ------------------------------------------------------------------------

# Fungi
#                   Df  SumsOfSqs  MeanSqs  F.Model      R2       Pr(>F)    
#Site           2    2.2735    1.13677  3.02143   0.12183     0.0001 ***
#Distance           3    1.1452    0.38173  1.01462   0.06137     0.3989    
#Site:Distance  6    2.0745    0.34574  0.91896   0.11116     0.9213    
#Residuals         35    13.1682   0.37623            0.70564           
#Total             46    18.6614                      1.00000 


# repeat this for prokaryotes and please develop the following: 
# what if we test with the other year for Citta' di Castello?
# what about the effect of the year in the dataset that contain only 
# Citta di Castello?
# We should also try one including the truffles.

# test the effect Site 
adonis_prok <- adonis(t(otu_prok_uparse_mSeq) ~ Site,
                      metadata_prok_uparse_mSeq, 
                      strata=metadata_prok_uparse_mSeq$Distance, 
                      permutations=9999)
adonis_prok

# test for the effect of Distance 
adonis_prok1 <- adonis(t(otu_prok_uparse_mSeq) ~ Distance,  
                       metadata_prok_uparse_mSeq, 
                       strata=metadata_prok_uparse_mSeq$Site, 
                       permutations=9999)
adonis_prok1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_prok2 <- adonis(t(otu_prok_uparse_mSeq) ~ Distance * Site,  
                       metadata_prok_uparse_mSeq, 
                       permutations=9999)
adonis_prok2

adonis_prok3 <- adonis(t(otu_prok_uparse_mSeq) ~ Site * Distance,
                       metadata_prok_uparse_mSeq, 
                       permutations=9999)
adonis_prok3
round(p.adjust(adonis_prok3$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0003 0.0891 1.0000     NA     NA

# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Site           2    1.5820 0.79099  7.7525 0.24844 0.0001 ***
#   Distance       3    0.5178 0.17260  1.6916 0.08132 0.0077 ** 
#   Site:Distance  6    0.5949 0.09914  0.9717 0.09342 0.5378    
# Residuals     36    3.6731 0.10203         0.57683           
# Total         47    6.3677                 1.00000      
# 


# Now testing the years --------------------------------------------------------------------------
count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "soil" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "Citta_di_Castello_2015"))))), Year)

ExtractOTU_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "Citta_di_Castello_2015")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014",
                                               "Citta_di_Castello_2015")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_y(physeq_fungi_uparse_mSeq) -> otu_fungi_uparse_mSeq_y
ExtractMETA_y(physeq_fungi_uparse_mSeq) -> metadata_fungi_uparse_mSeq_y

ExtractOTU_y(physeq_prok_uparse_mSeq) -> otu_prok_uparse_mSeq_y
ExtractMETA_y(physeq_prok_uparse_mSeq) -> metadata_prok_uparse_mSeq_y

# Fungi ------------------------------------------------------------
# test the effect Site 
adonis_fungi <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Year,
                       metadata_fungi_uparse_mSeq_y, 
                       strata=metadata_fungi_uparse_mSeq_y$Distance, 
                       permutations=9999)
adonis_fungi

# test for the effect of Distance 
adonis_fungi1 <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Distance,  
                        metadata_fungi_uparse_mSeq_y, 
                        strata=metadata_fungi_uparse_mSeq_y$Site, 
                        permutations=9999)
adonis_fungi1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_fungi2 <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Distance * Year,  
                        metadata_fungi_uparse_mSeq_y, 
                        permutations=9999)
adonis_fungi2

adonis_fungi3_y <- adonis(t(otu_fungi_uparse_mSeq_y) ~ Year * Distance,
                        metadata_fungi_uparse_mSeq_y, 
                        permutations=9999)
adonis_fungi3_y
round(p.adjust(adonis_fungi3_y$aov.tab$`Pr(>F)`, "bonferroni"), 4)
# [1] 0.0033 1.0000 1.0000     NA     NA

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Year           1    0.7042 0.70416 1.76779 0.05688 0.0009 ***
#   Distance       3    1.0710 0.35699 0.89622 0.08650 0.8611    
# Year:Distance  3    1.0455 0.34850 0.87491 0.08445 0.9075    
# Residuals     24    9.5598 0.39833         0.77217           
# Total         31   12.3804                 1.00000  


# PRokaryotes -----------------------------------------------------------
# test the effect Site 
adonis_prok <- adonis(t(otu_prok_uparse_mSeq_y) ~ Year,
                      metadata_prok_uparse_mSeq_y, 
                      strata=metadata_prok_uparse_mSeq_y$Distance, 
                      permutations=9999)
adonis_prok

# test for the effect of Distance 
adonis_prok1 <- adonis(t(otu_prok_uparse_mSeq_y) ~ Distance,  
                       metadata_prok_uparse_mSeq_y, 
                       strata=metadata_prok_uparse_mSeq_y$Site, 
                       permutations=9999)
adonis_prok1

# test for interactions - in these case I do no include the strata
# since Distance which is at the beginning of the model is treated as a random
# effect 
adonis_prok2 <- adonis(t(otu_prok_uparse_mSeq_y) ~ Distance * Year,  
                       metadata_prok_uparse_mSeq_y, 
                       permutations=9999)
adonis_prok2

adonis_prok3_y <- adonis(t(otu_prok_uparse_mSeq_y) ~ Year * Distance,
                       metadata_prok_uparse_mSeq_y, 
                       permutations=9999)
adonis_prok3_y
round(p.adjust(adonis_prok3_y$aov.tab$`Pr(>F)`, "bonferroni"), 4)
#[1] 0.0012 0.0912 1.0000     NA     NA

#                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Year           1    0.4289 0.42890  3.2038 0.09383 0.0004 ***
#   Distance       3    0.5744 0.19147  1.4303 0.12567 0.0304 *  
#   Year:Distance  3    0.3547 0.11824  0.8832 0.07760 0.7023    
# Residuals     24    3.2130 0.13387         0.70290           
# Total         31    4.5710                 1.00000  


# BETADISPER - analysis of betadispersion ---------------------------------------------------------

# Gio', I think you have the code for this, right?
# It needs to be run on both fungi and bacteria on Site
# and Distance. You could also plot the distance to centroids
# as boxplots, it is nice. 

# Fungi -------------------------------------------------------------------------------------------
sample_data(physeq_fungi_uparse_mSeq)$LocYear <- paste(sample_data(physeq_fungi_uparse_mSeq)$Site,
                                                      sample_data(physeq_fungi_uparse_mSeq)$Year,
                                                      sep = "_")

sample_data(physeq_fungi_uparse_mSeq)$Distance <- factor(sample_data(physeq_fungi_uparse_mSeq)$Distance, 
                                                         levels=c("Inside", "0", "40", "100", "200"))
physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_mSeq_site

physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_mSeq_year


# Prokaryotes -------------------------------------------------------------------------------------
sample_data(physeq_prok_uparse_mSeq)$LocYear <- paste(sample_data(physeq_prok_uparse_mSeq)$Site,
                                                      sample_data(physeq_prok_uparse_mSeq)$Year,
                                                      sep = "_")

sample_data(physeq_prok_uparse_mSeq)$Distance <- factor(sample_data(physeq_prok_uparse_mSeq)$Distance, 
                                                        levels=c("Inside", "0", "40", "100", "200"))
physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_site

physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="soil" & LocYear%in%c("Citta_di_Castello_2014","Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_year

physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_mSeq_truf


PlotOrdin <-function(dataframe, ord, var){
  ord <- plot_ordination(dataframe, ord, color=var, shape="Distance") + 
    geom_point(size=1.2, alpha=0.9, aes(shape=Distance)) +
    scale_shape_manual(values = c(15, 16, 17, 18)) + #  labels=c("","", ...)
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 8, face = "bold")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") +
    guides(color=guide_legend(nrow=2, order = 2), shape=guide_legend(nrow=2, order = 1) )
  return(ord)
}

PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_site, "Site") + 
  scale_colour_manual("Site", values=palette_site, 
          labels = c("Site 1", "Site 2", "Site 3")) + 
  labs(title = "Fungi") +
  annotate("text", Inf, -Inf, 
           label = expression(paste("Site ", italic(R) ^ 2,"= 7.2%***"), parse=TRUE), 
           size = 3, hjust = 1.2, vjust = -0.5) +
  annotate("text", -Inf, -Inf, 
           label = expression(paste("Depth ", italic(R) ^ 2,"= 7.2%***"), parse=TRUE),
           size = 3, hjust = -0.2, vjust = -0.5)


PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_year, "Year") + 
  scale_colour_manual("Year", values=palette_year, 
                      labels = c("2014", "2015")) + 
  labs(title = "Fungi")

PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_truf, "Site") + 
  scale_colour_manual("Site", values=paletteCB3_prok, 
                      labels = c("Citta di Castello'", "San Giovanni", "Ripabianca")) + 
  labs(title = "Prokaryotes")


# *** FIGURE 2 - Ordination (JUST AN EXAMPLE) -----------------------------------------------------
library(grid)
library(gridExtra)

plots1 <- ggarrange(PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_site, "Site") + 
                      scale_colour_manual("Site", values=palette_site, 
                                          labels = c("Site 1", "Site 2", "Site 3")) + 
                      labs(title = "Fungi") +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Site ", italic(R) ^ 2,"= 12.2% **"), parse=TRUE),
                               size = 2.7, hjust = -0.05, vjust =  1),
                    PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_site, "Site") + 
                      scale_colour_manual("Site", values=palette_site, 
                                          labels = c("Site 1", "Site 2", "Site 3")) + 
                    labs(title = "Prokaryotes") +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Site ", italic(R) ^ 2,"= 24.8% ***"), parse=TRUE),
                               size = 2.7, hjust = -0.05, vjust = 1) +
                      annotate("text", -Inf, Inf, 
                               label = expression(paste("Distance ", italic(R) ^ 2,"= 8.1% **"), parse=TRUE),
                               size = 2.7, hjust = -0.04, vjust = 2.5),
                    labels = c("A", "B"),
                    widths = c(1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

title1=text_grob("Site",size=12, face=2)
grid.arrange(plots1, top = title1)

plots2 <- ggarrange(PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_year, "Year") + 
                      scale_colour_manual("Site", values=palette_year, 
                                          labels = c("2014", "2015")) + 
                      labs(title = "Fungi") +
                      annotate("text", Inf, -Inf, 
                               label = expression(paste("Year ", italic(R) ^ 2,"= 5.7% **"), parse=TRUE), 
                               size = 2.7, hjust = 1.2, vjust = -0.5) , 
                    PlotOrdin(physeq_prok_uparse_mSeq, pcoa_prok_mSeq_year, "Year") + 
                      scale_colour_manual("Site", values=palette_year, 
                                          labels = c("2014", "2015")) + 
                      labs(title = "Prokaryotes") +
                      annotate("text", Inf, -Inf, 
                               label = expression(paste("Year ", italic(R) ^ 2,"= 9.4% **"), parse=TRUE), 
                               size = 2.7, hjust = 1.2, vjust = -0.5) ,
                    labels = c("C", "D"),
                    widths = c(1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

title2=text_grob("Year",size=12, face=2)
grid.arrange(plots2, top = title2)


Fig_2 = ggarrange(grid.arrange(plots1, top = title1),
                  grid.arrange(plots2, top = title2),
                   widths = c(1,1),
                   align = "hv", 
                     ncol = 1, 
                     nrow = 2)

Fig_2

# I removed the Truffles as you see but I am not sure if if is a good idea. Compasing the
# soils form 0 to 200 is interesting but maybe it is interesiting also looking into how different
# the soils are form the truffles. 

# DA FARE
# Lo script lo hai sull'altro file con UNOISE. Potresti provare a vedere se veien meglio.
# inoltre bosognerebbe farla anche con l'atro sito non incluso nei 3 qui.
# Anche la comparazione 2014 e 2014 nello stesso sito va fatta.

# TEST on ordinations ----------------------------------------------------------------------------
library(microbiome)

physeq_fungi_uparse_mSeq %>%
  transform("log10", target = "OTU") %>%
  subset_samples(LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_mSeq_test

PlotOrdin(physeq_fungi_uparse_mSeq, pcoa_fungi_mSeq_test, "Site") + 
  scale_colour_manual("Site", values=palette_site, 
                      labels = c("Citta' di Castello", "San Giovanni D'Asso", "Ripabianca")) +
  scale_shape_manual(values = c(0, 1, 2, 5, 6)) 




# test the effect Site 
physeq_fungi_uparse_mSeq %>%
  transform("hellinger", target = "OTU") %>%
  subset_samples(LocYear%in%c("Citta_di_Castello_2014","San_Giovanni_dAsso_2015","Ripabianca_2014")) -> fungi_test

otu <- as.data.frame(otu_table(fungi_test))
meta <- as(sample_data(fungi_test), "data.frame")

adonis(t(otu) ~ Type * Site * Distance, meta,permutations=9999)

# *******************************************************************************************-------
# ANALYZING TRUFFLES MICROBIOME --------------------------------------------------------------------

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "truffle" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "San_Giovanni_dAsso_2015",
          "Ripabianca_2014"))))), Site)

count(as.data.frame(as.matrix(
  sample_data(
    physeq_fungi_uparse_mSeq %>%
      subset_samples(
        Type == "truffle" & LocYear %in% c(
          "Citta_di_Castello_2014",
          "Citta_di_Castello_2015"))))), Year)


ExtractOTU_t <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                               "San_Giovanni_dAsso_2015",
                                               "Ripabianca_2014")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_t <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "San_Giovanni_dAsso_2015",
                                                  "Ripabianca_2014")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_t(physeq_fungi_uparse_mSeq) -> otu_fungi_truf
ExtractMETA_t(physeq_fungi_uparse_mSeq) -> metadata_fungi_truf

ExtractOTU_t(physeq_prok_uparse_mSeq) -> otu_prok_truf
ExtractMETA_t(physeq_prok_uparse_mSeq) -> metadata_prok_truf



# test the effect Site 
adonis_fungi_t <- adonis(t(otu_fungi_truf) ~ Site,
                         metadata_fungi_truf, 
                         permutations=9999)
adonis_fungi_t

adonis_prok_t <- adonis(t(otu_prok_truf) ~ Site,
                         metadata_prok_truf, 
                         permutations=9999)
adonis_prok_t


ExtractOTU_t_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "Citta_di_Castello_2015")) -> Physeq_filt
  otu <- as.data.frame(otu_table(Physeq_filt))
  return(otu)
}

ExtractMETA_t_y <- function(physeq){
  physeq %>%
    subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                  "Citta_di_Castello_2015")) -> Physeq_filt
  meta <- as(sample_data(Physeq_filt), "data.frame")
  return(meta)
}


ExtractOTU_t_y(physeq_fungi_uparse_mSeq) -> otu_fungi_truf_year
ExtractMETA_t_y(physeq_fungi_uparse_mSeq) -> metadata_fungi_truf_year

ExtractOTU_t_y(physeq_prok_uparse_mSeq) -> otu_prok_truf_year
ExtractMETA_t_y(physeq_prok_uparse_mSeq) -> metadata_prok_truf_year


# test the effect Year 
adonis_fungi_t_year <- adonis(t(otu_fungi_truf_year) ~ Year,
                         metadata_fungi_truf_year, 
                         permutations=9999)
adonis_fungi_t_year

adonis_prok_t_year <- adonis(t(otu_prok_truf_year) ~ Year,
                              metadata_prok_truf_year, 
                              permutations=9999)
adonis_prok_t_year


# Test the homogenity of group variances
permdisp_fungi_truf_year <- 
  betadisper(vegdist(t(otu_fungi_truf_year), method="bray"), metadata_fungi_truf_year$Year) 
anova(permdisp_fungi_truf_year, permutations = 9999)


permdisp_prok_truf_year <- 
  betadisper(vegdist(t(otu_prok_truf_year), method="bray"), metadata_prok_truf_year$Year) 
anova(permdisp_prok_truf_year, permutations = 9999)



# Plotting
PlotOrdinTruf <-function(dataframe, ord, var){
  ord <- plot_ordination(dataframe, ord, color=var) + 
    geom_point(size=1.2, alpha=0.9) +
    scale_shape_manual(values = c(15, 16, 17, 18)) + #  labels=c("","", ...)
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
    theme(strip.text.x = element_text(size = 8, face = "bold")) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 8,hjust = 0.5, vjust = 0.5)) +
    theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
    theme(axis.title = element_text(angle = 0, size = 7, face = "bold")) +
    grids(linetype = "dashed") +
    theme(legend.position="bottom") +
    guides(color=guide_legend(nrow=2, order = 2), shape=guide_legend(nrow=2, order = 1) )
  return(ord)
}



physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "San_Giovanni_dAsso_2015",
                                                "Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_truf_site


physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "San_Giovanni_dAsso_2015",
                                                "Ripabianca_2014")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_truf_site



physeq_fungi_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_fungi_truf_year


physeq_prok_uparse_mSeq %>%
  subset_samples(Type=="truffle" & LocYear%in%c("Citta_di_Castello_2014",
                                                "Citta_di_Castello_2015")) %>%
  ordinate(method ="PCoA", distance="bray") -> pcoa_prok_truf_year



PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_year, "Year") +
  scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
  labs(title = "Fungi")

PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_year, "Year") +
  scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
  labs(title = "Fungi")



# FIGURE suppl. XXX - truffle geography and year -------------------------------------------------

plots_truf_site <- ggarrange(
  PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_site, "Site") +
    scale_colour_manual("Site", values=palette_site, 
                        labels = c("Site 1", "Site 2", "Site 3")) +
    labs(title = "Fungi") +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Site ", italic(R) ^ 2,"= 6.8% ns"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust =  1),
  PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_site, "Site") +
    scale_colour_manual("Site", values=palette_site, 
                        labels = c("Site 1", "Site 2", "Site 3")) +
    labs(title = "Prokaryotes") +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Site ", italic(R) ^ 2,"= 19.9% ns"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust = 1),
                    labels = c("A", "B"),
                    widths = c(1,1),
                    align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

plots_truf_site


plots_truf_year <- 
  ggarrange(
  PlotOrdinTruf(physeq_fungi_uparse_mSeq, pcoa_fungi_truf_year, "Year") +
    scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) +
    labs(title = "Fungi") +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Year ", italic(R) ^ 2,"= 30.0% *"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust =  1) +
    annotate("text", -Inf, Inf, 
             label = expression(paste("Dispersion ", italic(p),"= 0.038"), parse=TRUE),
             size = 2.3, hjust = -0.05, vjust =  2.5),
  PlotOrdinTruf(physeq_prok_uparse_mSeq, pcoa_prok_truf_year, "Year") +
    scale_colour_manual("Year", values=palette_year, labels = c("2014", "2015")) + 
    labs(title = "Prokaryotes") +
    annotate("text", Inf, -Inf, 
             label = expression(paste("Year ", italic(R) ^ 2,"= 27.9% *"), parse=TRUE),
             size = 2.3, hjust = 1, vjust = -1.5) +
    annotate("text", Inf, -Inf, 
             label = expression(paste("Dispersion ", italic(p),"= 0.036"), parse=TRUE),
             size = 2.3, hjust = 1, vjust = -0.5),
  labels = c("C", "D"),
  widths = c(1,1),
  align = "hv", 
  ncol = 2, 
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom")


plots_truf_year


Fig_3_suppl = ggarrange(
                  grid.arrange(plots_truf_site, top = title1),
                  grid.arrange(plots_truf_year, top = title2),
                  widths = c(1,1),
                  align = "hv", 
                  ncol = 1, 
                  nrow = 2)

Fig_3_suppl



# VENN --------------------------------------------------------------------------------------------
# try making some venns. I think you have my code, right?
# Just simple ones, comparing truffle vs. 0 for example.
# and I would also do truffle OTUs vs. truffle Isolates.
# samll example below.

library(limma)
source("../R_functions/my_venn_diag.R")

taxa_fungi_up[taxa_fungi_up$Genus%in%"Tuber", ]
taxa_fungi_up[!(rownames(taxa_fungi_up) %in% c("FOTU_1","FOTU_3","FOTU_809","FOTU_3308",
                                               "FOTU_1844","FOTU_1987","FOTU_2741","FOTU_2998")), ] 
physeq_fungi_uparse_filt -> physeq_fungi_uparse_ven

tax_table(physeq_fungi_uparse_ven) <- tax_table(as.matrix(
  taxa_fungi_up[!(rownames(taxa_fungi_up) %in% c("FOTU_1","FOTU_3","FOTU_809","FOTU_3308",
                                                 "FOTU_1844","FOTU_1987","FOTU_2741","FOTU_2998")), ]))

tax_table(physeq_fungi_uparse_ven) <-
  tax_table(physeq_fungi_uparse_ven)[,!(colnames(tax_table(physeq_fungi_uparse_ven)) %in% c("OTU_ID"))]

# Fungi
physeq_fungi_uparse_ven %>% 
  tax_glom(taxrank="Genus") %>% #how many shared genera
  subset_samples(Distance%in%c("Inside","0")) %>%
  merge_samples("Type") -> physeq_fungi_0inside

otu_table(physeq_fungi_0inside) <- otu_table(physeq_fungi_0inside)[,which(
  colSums(otu_table(physeq_fungi_0inside)) > 0)] 

otu_fungi_0inside <- as.data.frame(t(otu_table(physeq_fungi_0inside)))
counts_fungi_0inside <- vennCounts(otu_fungi_0inside, include="both")
counts_fungi_0inside

venn_Gian(counts_fungi_0inside,
          cex=c(1),
          circle.col =c("brown", "gold"),
          mar = c(1,1,1,1),
          lwd = 2, main="Fungi")


# Prokaryotes
physeq_prok_uparse_filt -> physeq_prok_uparse_ven
tax_table(physeq_prok_uparse_ven) <-
  tax_table(physeq_prok_uparse_ven)[,!(colnames(tax_table(physeq_prok_uparse_ven)) %in% c("OTU_ID"))]

physeq_prok_uparse_ven %>% 
  tax_glom(taxrank="Genus") %>%
  subset_samples(Distance%in%c("Inside","0")) %>%
  merge_samples("Type") -> physeq_prok_0inside

otu_table(physeq_prok_0inside) <- otu_table(physeq_prok_0inside)[, which(
  colSums(otu_table(physeq_prok_0inside)) > 0)] 

otu_prok_0inside <- as.data.frame(t(otu_table(physeq_prok_0inside)))
counts_prok_0inside <- vennCounts(otu_prok_0inside, include="both")
counts_prok_0inside

venn_Gian(counts_prok_0inside,
          cex=c(1),
          circle.col =c("brown", "gold"),
          mar = c(1,1,1,1),
          lwd = 2, main="Prokaryotes")

# Extract OTU names for shared and truffle only compartment
library(reshape2)

head(otu_fungi_0inside)
subset(otu_fungi_0inside, truffle>0) -> fungi_core
fungi_core$niche <- ifelse(fungi_core$soil>0, "shared", "truffle")
left_join(rownames_to_column(fungi_core), 
          rownames_to_column(taxa_fungi_up), by = ("rowname" = "rowname")) -> fungi_core
fungi_core$rel_soil <- fungi_core$soil/sum(fungi_core$soil)*100
fungi_core$rel_truffle <- fungi_core$truffle/sum(fungi_core$truffle)*100
head(fungi_core)
melt(fungi_core[, c(11,17:18)]) -> melt_fungi_core


head(otu_prok_0inside)
subset(otu_prok_0inside, truffle>0) -> prok_core
prok_core$niche <- ifelse(prok_core$soil>0, "shared", "truffle")
left_join(rownames_to_column(prok_core), 
          rownames_to_column(taxa_prok_up), by = ("rowname" = "rowname")) -> prok_core
prok_core$rel_soil <- prok_core$soil/sum(prok_core$soil)*100
prok_core$rel_truffle <- prok_core$truffle/sum(prok_core$truffle)*100

# selecting top 30 bacteria
prok_core$sum <- prok_core$rel_soil + prok_core$rel_truffle
prok_core[order(-prok_core$sum),][1:30,] -> prok_core

# correctin taxonomy 
refseq(physeq_prok_uparse_filt)[prok_core$rowname,] 

write.dna(refseq(refseq(physeq_prok_uparse_filt)[prok_core$rowname,]), 
          format="fasta", 
          colsep="",
          file="bacteria_core30.fasta")

#write.csv(prok_core, "prok_core.csv")
read.csv("prok_core_new.csv", sep = ",", header = TRUE, row.names = 1) -> prok_core
head(prok_core)
melt(prok_core[, c(11,17:18)]) -> melt_prok_core

melt_fungi_core
range(melt_fungi_core$value)
melt_fungi_core$value <- ifelse(melt_fungi_core$value==0, NA, paste(melt_fungi_core$value))
melt_fungi_core$value <- as.numeric(as.character(melt_fungi_core$value))

melt_prok_core
range(melt_prok_core$value)
melt_prok_core$value <- ifelse(melt_prok_core$value==0, NA, paste(melt_prok_core$value))
melt_prok_core$value <- as.numeric(as.character(melt_prok_core$value))


#https://stackoverflow.com/questions/23919064/figure-output-with-same-scale-size-between-plots-but-different-axis-length


# Plot Baloons 
PlotBaloon <- function(df){
plotb <- ggplot(df, aes(x = variable, y = reorder(Genus, desc(Genus)))) +
  geom_point(aes(size=value, fill=value, color=value), shape=16)+
  #scale_size_area("Abundance %", max_size = 20) +
  #scale_size("Abundance %", range = c(1,20)) +
  scale_x_discrete(labels= c(rel_soil = "Soil", rel_truffle="Truffle")) +
  scale_colour_gradientn(colours = c("#d73027","darkgrey","#4575b4")) +
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 9, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, face ="italic", size = 7, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(angle = 90, size = 9, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="bottom") +
  guides(color = FALSE, fill = FALSE) +
  grids(linetype = "dashed") 
return(plotb)
}

PlotBaloon(melt_fungi_core) +
  labs(title = "Fungi", x=NULL, y="Genus")+
  scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))

PlotBaloon(melt_prok_core) + 
  labs(title = "Prokaryotes",  x=NULL, y="Genus") +
  scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))

get_legend(PlotBaloon(melt_fungi_core) +
             labs(title = "Fungi", x=NULL, y="Genus")+
             scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40))) -> legend_baloon
as_ggplot(legend_baloon)

write.csv(melt_fungi_core, "melt_fungi_core.csv")
write.csv(melt_prok_core, "melt_prok_core.csv")

# FIGURE 5 - core taxa truffle/soil -----------------------------------------------
plots_baloon <- ggarrange(PlotBaloon(melt_fungi_core) +
                            labs(title = "Fungi", x=NULL, y="Genus")+
                            scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)),
                          PlotBaloon(melt_prok_core) + 
                            labs(title = "Prokaryotes",  x=NULL, y="Genus") +
                            scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)),
                    labels = c("A", "B"),
                    widths = c(1,1),
                    #align = "hv", 
                    ncol = 2, 
                    nrow = 1,
                    common.legend = TRUE,
                    legend = "bottom")

plots_baloon

# plotting with cowplot to customize the plot sizes 
library(cowplot)

ggdraw() +
  draw_plot(PlotBaloon(melt_fungi_core) +
              labs(title = "Fungi", x=NULL, y="Genus")+
              scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)) +
              theme(legend.position="none"), 
            x = 0, y = 0.24, width = 0.45, height = 0.6) +
  draw_plot(PlotBaloon(melt_prok_core) + 
              labs(title = "Prokaryotes",  x=NULL, y="Genus") +
              scale_size_continuous("Abundance %",range = c(0,15), breaks = c(1, 10, 20, 30, 40)) +
              theme(legend.position="none"), 
            x = 0.45, y = 0.1, width = 0.5, height = 0.74) +
  draw_plot(as_ggplot(legend_baloon), 
            x = 0, y = 0, width = 1, height = 0.1) +
  draw_plot_label(c("A", "B"), x = c(0, 0.45), y = c(0.84, 0.84), size = 12)



# *******************************************-------------------------
# *******************************************-------------------------
# *******************************************-------------------------







ggplot(melt_fungi_core, aes(x =variable, y = Genus)) +
  geom_point(aes(size=value, fill=value, color=value), shape=16)+
  scale_size_area("Abundance %", max_size = 20) +
  #scale_size("Abundance %", range = c(1,15)) +
  scale_x_discrete(labels= c(rel_soil = "Soil", rel_truffle="Truffle")) +
  scale_colour_gradientn(colours = c("#d73027","darkgrey","#4575b4")) +
  theme_classic() +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, face = "bold", size = 8, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, face ="italic", size = 7, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_text(angle = 90, size = 7, face = "bold", hjust = 0.5, vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position="bottom") +
  grids(linetype = "dashed") +
  guides(color = FALSE, fill = FALSE,
         size = guide_legend(override.aes = list(size = c(0.1, 2, 4, 6, 8)), nrow  = 2))

# >>> BARPLOTS OF CORE TAXA ------------------------------------------------------------------------------------
# I am using original count data 
sample_data(physeq_fungi_uparse_filt)$LocYear <- paste(sample_data(physeq_fungi_uparse_filt)$Site,
                                                       sample_data(physeq_fungi_uparse_filt)$Year,
                                                       sep = "_")

physeq_fungi_uparse_filt %>%
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                              "San_Giovanni_dAsso_2015",
                              "Ripabianca_2014")) %>%
  tax_glom(taxrank = "Genus") %>%
  merge_samples("Site") %>% 
  psmelt() %>%                                         
  #filter(Abundance > 0.01) %>%                         
  arrange(Phylum) -> otu_fungi_phylum

head(otu_fungi_phylum)
dim(otu_fungi_phylum)
levels(otu_fungi_phylum$Phylum)
levels(otu_fungi_phylum$Distance)

rownames(fungi_core)
subset(otu_fungi_phylum, OTU %in% rownames(fungi_core)) -> subset_fungi

# palette
unlined <-c("#C1B4C0","#D191B4","#E8A0C1","#F9C0B0","#F7C789",
            "#E8E990","#C7BC7B","#8FB26A","#717B71","#77A9C2",
            "#33A0DAD4","#1D5A91","#33197BAD","#F39B3A")


# plotting 
PlotBarGraph <-function(df){ 
barplot <- ggplot(df, aes(x = Sample, y = Abundance, fill = Genus)) + 
  theme_classic() +
  geom_bar(stat = "identity") +
  labs(title="", x="Samples", y = "Relative Abundance") +
  scale_fill_manual(values = unlined) +
  #scale_fill_manual(values = palette_fiji) +
  facet_grid(~Site, scales = "free_x",
             space="free_x", 
             labeller = as_labeller(c(Citta_di_Castello = "Citta' di\nCastello",
                                      Ripabianca = "Ripabianca",
                                      San_Giovanni_dAsso = "San Giovanni\nD'Asso"))) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")
return(barplot)
}


PlotBarGraph(fungi_core)


# plotting heatmaps function -----------------------------------------------------------------------
PlotHeat <- function(df, rank){
  plot_heat <- ggplot(df, aes(x=Sample, y=get(rank), fill=sqrt(Abundance))) +
    geom_tile(aes(height = 0.92, width = 0.92))+ #colour="white",size=0.1
    scale_y_discrete(expand=c(0,0),) +
    facet_grid(~Site, scales = "free_x", space="free_x") + 
    theme_classic() +
    theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
    theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7, color = "black")) +
    theme(strip.text.x = element_text(size = 9, face = "bold")) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) + 
    theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5, color = "black")) +
    theme(axis.text.y = element_text(angle = 0, face = "italic", size = 6, hjust = 1, vjust = 0.5, colour = "black")) +
    theme(axis.title = element_text(angle = 0, size = 9, face = "bold")) +
    theme(plot.background = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0.2, "lines"))
  return(plot_heat)
}


PlotHeat(fungi_core, "Genus")

















# BARPLOTS ----------------------------------------------------------------------------------------
library(evoPalette)
launch_evo_palette()
palette_box()

yummy_orn <- c("#C69DC3","#F57E4C","#FDD00E","#6F9358","#24548EBA","#22FEBDA1")
q_ms <- c("#DDDDDD","#CCA937","#DE8E1A","#DF5327","#91AE42","#526F80")


# I would do two barplots. One for fungi and one for bacteria
# probably at class level. Then I would do another one showing the
# first 20 genera. and an additional one that shows only Tuber spp.
# those can eba smaall and stackenon over the other.
# See example below.

# >>> BARPLOTS ------------------------------------------------------------------------------------
# I am using original count data 
sample_data(physeq_fungi_uparse_filt)$LocYear <- paste(sample_data(physeq_fungi_uparse_filt)$Site,
                                                       sample_data(physeq_fungi_uparse_filt)$Year,
                                                       sep = "_")

physeq_fungi_uparse_filt %>%
  subset_samples(LocYear%in%c("Citta_di_Castello_2014",
                              "San_Giovanni_dAsso_2015",
                              "Ripabianca_2014")) %>%
  tax_glom(taxrank = "Phylum") %>%                     
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>%                                         
  #filter(Abundance > 0.01) %>%                         
  arrange(Phylum) -> otu_fungi_phylum

head(otu_fungi_phylum)
levels(otu_fungi_phylum$Phylum)
levels(otu_fungi_phylum$Distance)

# otu_fungi_phylum$Distance <- factor(otu_fungi_phylum$Distance,
#                                     levels=c("Inside", "0", "40", "100", "200"))

barplot_fungi <- ggplot(otu_fungi_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  theme_classic() +
  geom_bar(stat = "identity") +
  labs(title="", x="Samples", y = "Relative Abundance") +
  scale_fill_manual(values = unlined) +
  #scale_fill_manual(values = palette_fiji) +
  facet_grid(~Site, scales = "free_x",
             space="free_x", 
             labeller = as_labeller(c(Citta_di_Castello = "Citta' di\nCastello",
                                      Ripabianca = "Ripabianca",
                                      San_Giovanni_dAsso = "San Giovanni\nD'Asso"))) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5)) +
  #theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")

barplot_fungi



# Couple of things. You should reorder the samples according to Distance,
# like I did for the richness plots. You have a column called Sample, you have
# to recode it according to the right order Inside, 0, 40, 100, 200.
# Also, I want to figure out how to draw lines at the x-axis to group samples
# as "Inside", "0", "40", "100", "200" fro each facet.





# >>> HEATMAP TREES -------------------------------------------------------------------------------
# https://grunwaldlab.github.io/metacoder_documentation/example.html
# https://grunwaldlab.github.io/metacoder_documentation/publication--06--tara_oceans.html
# https://grunwaldlab.github.io/metacoder_documentation/workshop--07--diversity_stats.html

library(metacoder)
library(tibble)

physeq_fungi_uparse_filt %>% 
  subset_samples(Distance%in%c("Inside","0")) %>%
  tax_glom(taxrank="Genus") -> physeq_fungi_uparse_filt_gen

otu_table(physeq_fungi_uparse_filt_gen) <- otu_table(physeq_fungi_uparse_filt_gen)[which(
                  rowSums(otu_table(physeq_fungi_uparse_filt_gen)) > 0),] 

physeq_fungi_uparse_filt_gen
head(tax_table(physeq_fungi_uparse_filt_gen))
head(sample_data(physeq_fungi_uparse_filt_gen))


# otu table
otu_fungi_gen <- as.data.frame(otu_table(physeq_fungi_uparse_filt_gen))
otu_fungi_gen$sample_id <- rownames(otu_fungi_gen)
head(otu_fungi_gen)

metadata_fungi_gen <- as.data.frame(as.matrix(sample_data(physeq_fungi_uparse_filt_gen)))
metadata_fungi_gen$sample_id <- rownames(metadata_fungi_gen)
metadata_fungi_gen
str(metadata_fungi_gen)

# remove NA
tax_table(physeq_fungi_uparse_filt_gen)[is.na(tax_table(physeq_fungi_uparse_filt_gen))] <- ""

taxa_fungi_gen <- as.matrix(tax_table(physeq_fungi_uparse_filt_gen))
taxa_fungi_gen[,c(2:7)] -> taxa_fungi_gen
as.data.frame(taxa_fungi_gen) -> taxa_fungi_gen
head(taxa_fungi_gen)

# Reformatting taxonomy
taxa_fungi_gen[, "Genus"] <- gsub(" sp.", "", taxa_fungi_gen[, "Genus"])
cols <- c("Kingdom","Phylum","Class","Order","Family","Genus")
taxa_fungi_gen$Taxonomy <- do.call(paste, c(taxa_fungi_gen[cols], sep=";"))
head(taxa_fungi_gen)
head(taxa_fungi_gen)

identical(colnames(otu_fungi_gen), metadata_fungi_gen$sample_id)
sample_order <- match(metadata_fungi_gen$sample_id, colnames(otu_fungi_gen))
otu_fungi_gen <- otu_fungi_gen[,sample_order]

identical(rownames(taxa_fungi_gen), rownames(otu_fungi_gen))
otu_fungi_gen$Taxonomy <- taxa_fungi_gen$Taxonomy
head(otu_fungi_gen)

# transform a data.frame into a tibble ------------------------------------------------------------
tibbl_fungi <- tibble::as_tibble(otu_fungi_gen)
tibbl_fungi

tibbl_metadata_fungi <- tibble::as_tibble(metadata_fungi_gen)
tibbl_metadata_fungi

identical(colnames(tibbl_fungi[1:31]), tibbl_metadata_fungi$sample_id)
# sample_order <- match(tibbl_metadata_fungi$sample_id, colnames(tibbl_fungi))
# tibbl_fungi <- tibbl_fungi[,sample_order]

# create metacoder object 
heat_tree_ITS <- parse_tax_data(tibbl_fungi,
                                class_cols = "Taxonomy",
                                class_sep = ";")
print(heat_tree_ITS)
heat_tree_ITS$data$tax_data

#data_test <- tibbl_fungi[1:500,1:14]

# Removing low-abundance counts ------------------------------------------------------------------
# We need to thing about this very well. Maybe better usinf the
# non-normalized data!
heat_tree_ITS$data$tax_data <- zero_low_counts(heat_tree_ITS, "tax_data", min_count = 50)
no_reads <- rowSums(heat_tree_ITS$data$tax_data[, tibbl_metadata_fungi$sample_id]) == 0
sum(no_reads)

# removing OTUs with no reads left 
heat_tree_ITS <- filter_obs(heat_tree_ITS, data  = "tax_data", ! no_reads, drop_taxa = TRUE)
print(heat_tree_ITS)

# calculate taxa proprostion ----------------------------------------------------------------------
heat_tree_ITS$data$tax_data <- calc_obs_props(heat_tree_ITS, "tax_data")

# abundance per-taxon -----------------------------------------------------------------------------
heat_tree_ITS$data$tax_abund <- calc_taxon_abund(heat_tree_ITS, "tax_data", cols = tibbl_metadata_fungi$sample_id)

# number of samples have reads for each taxon -----------------------------------------------------
heat_tree_ITS$data$tax_occ <- calc_n_samples(heat_tree_ITS, "tax_abund", 
                                             cols = tibbl_metadata_fungi$sample_id,
                                             groups = tibbl_metadata_fungi$Type)

print(heat_tree_ITS$data$tax_occ)
# I think the column order match but we should verify! 
# Maybe adding the cols parameter would be better.

# adding taxon abundance
heat_tree_ITS$data$type_abund <- calc_group_mean(heat_tree_ITS, "tax_abund",
                                       cols = tibbl_metadata_fungi$sample_id,
                                       groups = tibbl_metadata_fungi$Type)

print(heat_tree_ITS$data$type_abund)


# plotting the heatmap tree
set.seed(1)
heat_tree(heat_tree_ITS, 
          node_label = heat_tree_ITS$taxon_names(),
          node_size = heat_tree_ITS$n_obs(),
          node_color = heat_tree_ITS$data$tax_occ$truffle, 
          overlap_avoidance = 0.01,
          node_label_size_range = c(0.008, 0.05), 
          tree_label_size = as.numeric(0.6),
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", # The primary layout algorithm
          initial_layout = "reingold-tilford") # The layout algorithm that initializes node Sites


set.seed(2)
heat_tree_ITS %>%
  #taxa::filter_taxa(truffle > 0.001) %>% # taxa:: needed because of phyloseq::filter_taxa
     heat_tree(node_label = heat_tree_ITS$taxon_names(),
               node_size = heat_tree_ITS$data$type_abund$truffle, 
               node_color = heat_tree_ITS$data$type_abund$truffle, 
               layout = "da", initial_layout = "re", 
               title = "Taxa inside truffles")



# Comparing groups --------------------------------------------------------------------------------
heat_tree_ITS$data$diff_table <- compare_groups(heat_tree_ITS, dataset = "tax_abund",
                                                cols = tibbl_metadata_fungi$sample_id, # What columns of sample data to use
                                                groups = tibbl_metadata_fungi$Type) # What category each sample is assigned to
print(heat_tree_ITS$data$diff_table)

heat_tree_ITS$data$diff_table$wilcox_p_value <- p.adjust(heat_tree_ITS$data$diff_table$wilcox_p_value, method = "fdr")

# The most useful statistic for plotting is the log of ratio of median abundances
# in the two groups, since it is centered on 0 and is symmetric (e.g., a value of -4
# is the same magnitude as 4). Lets set any differences that are not significant to 0
# so all differences shown on the plot are significant.

heat_tree_ITS$data$diff_table$log2_median_ratio[heat_tree_ITS$data$diff_table$wilcox_p_value > 0.05] <- 0
heat_tree_ITS$data$diff_table$wilcox_p_value[is.nan(heat_tree_ITS$data$diff_table$wilcox_p_value)] <- 0

heat_tree_matrix(heat_tree_ITS,
                 data = "diff_table",
                 node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                 node_label = taxon_names,
                 node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                 node_color_range = diverging_palette(), # The built-in palette for diverging data
                 node_color_trans = "linear", # The default is scaled by circle area
                 #node_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 #edge_color_interval = c(-3, 3), # The range of `log2_median_ratio` to display
                 node_size_axis_label = "Number of OTUs",
                 node_color_axis_label = "Log2 ratio median proportions",
                 layout = "davidson-harel", # The primary layout algorithm
                 initial_layout = "reingold-tilford") # The layout algorithm that initializes node Sites

# Saves the plot as a pdf file  output_file = "differential_heat_tree.pdf"

# This doesn't really work well. I will try to fix it later. You focus on the other graphs.

# FORM THE CORAL DATASET --------------------------------------------------------------------------
#>>> RANDOM FOREST --------------------------------------------------------------------------------
library(randomForest)
library(plyr) # for the "arrange" function
library(rfUtilities) # to test model significance
library(caret) # to get leave-one-out cross-validation and also contains the nearZeroVar function 
library(compositions)

# *********************************************************************----------------------------
# **** FUNGI *****---------------------------------------------------------------------------------
sort(colSums(otu_fungi))
sort(rowSums(otu_fungi))

# prune OTUs with non-zero values below a certain threshold ---------------------------------------
source("../R_functions/remove_rare.R")

# I removed OTUs that have non-zero values in <= 10% of samples
# that is OTUs present in at least 90% of the samples
otu_fungi_nonrare <- remove_rare(table = otu_fungi, cutoff_pro = 0.08)
dim(otu_fungi_nonrare)
head(otu_fungi_nonrare)

# Prune samples with super low read number
sort(colSums(otu_fungi_nonrare))
otu_fungi_nonrare[colSums(otu_fungi_nonrare) > 1000] -> otu_fungi_filt
otu_fungi_filt[rowSums(otu_fungi_filt) > 0,] -> otu_fungi_filt
dim(otu_fungi_filt)
sort(colSums(otu_fungi_filt))

# DATA TRANSFORMATION -----------------------------------------------------------------------------
# centred log-ratio (CLR) transformation for microbiome sequencing
# data - used for compositional data like otu_tables
library("compositions")

otus_fungi_clr <- clr(otu_fungi_filt)
dim(otus_fungi_clr)
head(otus_fungi_clr)
range(otus_fungi_clr)

# > CASSIFYING BY STATUS --------------------------------------------------------------------------
#metadata_fungi[colnames(otus_symb_clr),] -> meta_fungi_filt

otus_fungi_clr_Status <- data.frame(t(otus_fungi_clr))
otus_fungi_clr_Status$Sample <- rownames(otus_fungi_clr_Status)
otus_fungi_clr_Status$Status <-
  metadata_fungi[rownames(otus_fungi_clr_Status), "Status"]
head(otus_fungi_clr_Status)

# try tuning the model first ----------------------------------------------------------------------
set.seed(30)
sqrt(nrow(otus_fungi_clr))

bestmtry_fungi_St <-
  tuneRF(
    x = otus_fungi_clr_Status[, 1:(ncol(otus_fungi_clr_Status) - 2)],
    y = otus_fungi_clr_Status$Status,
    stepFactor = 1.5,
    mtryStart = 26,
    improve = 0.0001,
    ntree = 501,
    nodesize = 1,
    doBest = TRUE,
    importance = TRUE
  )

bestmtry_fungi_St

set.seed(31)
RF_fungi_St <-
  randomForest(
    x = otus_fungi_clr_Status[, 1:(ncol(otus_fungi_clr_Status) - 2)],
    y = otus_fungi_clr_Status$Status,
    ntree = 501,
    mtry = 12,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_St
plot(RF_fungi_St)

# removing unim-ortant variables using MeanDecreaseGini -------------------------------------------
head(importance(RF_fungi_St))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_fungi_St) 

# ASSESSING MODEL FIT -----------------------------------------------------------------------------
# Now check to see if the random forest is actually big enough... Up to a point, the more 
# trees in the forest, the better. You can tell when you've made enough when the OOB
# no longer improves.

oob_error_fungi_st <- data.frame(
  Trees = rep(1:nrow(RF_fungi_St$err.rate), times = 3),
  Status = rep(
    c("OOB", "bleached", "unbleached"),
    each = nrow(RF_fungi_St$err.rate)
  ),
  Error = c(
    RF_fungi_St$err.rate[, "OOB"],
    RF_fungi_St$err.rate[, "bleached"],
    RF_fungi_St$err.rate[, "unbleached"]
  )
)

oob_error_fungi_st$Status <- factor(oob_error_fungi_st$Status,
                                    levels = c("bleached", "unbleached", "OOB"))

p_oob_error_fungi_st <- ggplot(data=oob_error_fungi_st, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 48.44%", size=3) +
  geom_line(aes(color=Status)) + 
  scale_color_manual(values = c(paletteCB2, "red")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), legend.text = element_text(size =8)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_fungi_st

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS ------------------------------------------------------
set.seed(33)
perm_RF_fungi_St <-
  rf.significance(
    x = RF_fungi_St,
    xdata = otus_fungi_clr_Status[, 1:(ncol(otus_fungi_clr_Status) - 2)],
    nperm = 999,
    nmtry = 26,
    ntree = 501
  )
perm_RF_fungi_St

# PLOTTING THE RESULTS ----------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each 
# other. Start by converting the proximity matrix into a distance matrix.
dist_fungi_st <- as.dist(1 - RF_fungi_St$proximity)
mds_fungi_st <- cmdscale(dist_fungi_st, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_fungi_st <- round(mds_fungi_st$eig / sum(mds_fungi_st$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_fungi_st <- mds_fungi_st$points
mds_fungi_st <- data.frame(
  Sample = rownames(values_fungi_st),
  X = values_fungi_st[, 1],
  Y = values_fungi_st[, 2],
  Status = otus_fungi_clr_Status$Status
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

p_nmds_fungi_st = ggplot(data=mds_fungi_st, aes(x=X, y=Y, label=Sample, color=Status)) + 
  geom_point(size=1.2) +
  #geom_text_repel(label=rownames(mds_data), size = 2, min.segment.length = 0.5) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_fungi_st[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_fungi_st[2], "%", sep="")) +
  ggtitle("Fungi - Status") +
  scale_color_manual(values = paletteCB2) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 48.44%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8)) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  #guides(color=guide_legend(ncol=1)) +
  theme(legend.position="bottom")

p_nmds_fungi_st

# Identifying Important Features ------------------------------------------------------------------
imp_fungi_st <- as.data.frame(RF_fungi_St$importance)
imp_fungi_st$features <- rownames(imp_fungi_st)
imp_fungi_st <- arrange(imp_fungi_st, desc(MeanDecreaseAccuracy))
head(imp_fungi_st)
dim(imp_fungi_st)

taxa_fungi[rownames(imp_fungi_st), ] -> taxa_fungi_filt
dim(taxa_fungi_filt)
head(taxa_fungi_filt)

identical(rownames(taxa_fungi_filt), rownames(imp_fungi_st))
order_fungi_st <-
  match(rownames(taxa_fungi_filt), rownames(imp_fungi_st))
imp_fungi_st <- imp_fungi_st[order_fungi_st,]

imp_fungi_st$Taxonomy <- taxa_fungi_filt$Taxonomy
imp_fungi_st <-
  imp_fungi_st[order(imp_fungi_st$MeanDecreaseAccuracy, decreasing = TRUE),]
head(imp_fungi_st)

p_bar_fungi_st = ggplot(data=imp_fungi_st[1:20,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  labs(title = "Top Fungi OTUs to\nclassify coral Status", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.03) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 90, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_fungi_st

# *********************************************************************----------------------------
# > CASSIFYING BY SPECIES -------------------------------------------------------------------------
otus_fungi_clr_species <- data.frame(t(otus_fungi_clr))
otus_fungi_clr_species$Sample <- rownames(otus_fungi_clr_species)

metadata_fungi$Species <- metadata_fungi$Species %>%
  recode_factor(
    "Acropora_millepora" = "Acropora millepora",
    "Montipora_digitata" = "Montipora digitata",
    "Pocillopora_damicornis" = "Pocillopora damicornis",
    "Porites_cylindrica" = "Porites cylindrica"
  )

otus_fungi_clr_species$Species <-
  metadata_fungi[rownames(otus_fungi_clr_species), "Species"]
head(otus_fungi_clr_species)
dim(otus_fungi_clr_species)

# try tuning the model first ----------------------------------------------------------------------
set.seed(40)
sqrt(nrow(otus_fungi_clr))

bestmtry_fungi_sp <-
  tuneRF(
    x = otus_fungi_clr_species[, 1:(ncol(otus_fungi_clr_species) - 2)],
    y = otus_fungi_clr_species$Species,
    stepFactor = 1.5,
    mtryStart = 26,
    improve = 0.0001,
    ntree = 501,
    nodesize = 1,
    doBest = TRUE,
    importance = TRUE
  )

bestmtry_fungi_sp

set.seed(41)
RF_fungi_sp <-
  randomForest(
    x = otus_fungi_clr_species[, 1:(ncol(otus_fungi_clr_species) - 2)],
    y = otus_fungi_clr_species$Species,
    ntree = 501,
    mtry = 12,
    importance = TRUE,
    proximity = TRUE
  )

RF_fungi_sp
plot(RF_fungi_sp)

# removing unim-ortant variables using MeanDecreaseGini -------------------------------------------
head(importance(RF_fungi_sp))

##PlotMean Decrease Gini and Mean Decrease Accuracy 
varImpPlot(RF_fungi_sp) 

# ASSESSING MODEL FIT -----------------------------------------------------------------------------
oob_error_fungi_sp <- data.frame(
  Trees = rep(1:nrow(RF_fungi_sp$err.rate), times = 5),
  Species = rep(
    c(
      "OOB",
      "Acropora millepora",
      "Montipora digitata",
      "Pocillopora damicornis",
      "Porites cylindrica"
    ),
    each = nrow(RF_fungi_sp$err.rate)
  ),
  Error = c(
    RF_fungi_sp$err.rate[, "OOB"],
    RF_fungi_sp$err.rate[, "Acropora millepora"],
    RF_fungi_sp$err.rate[, "Montipora digitata"],
    RF_fungi_sp$err.rate[, "Pocillopora damicornis"],
    RF_fungi_sp$err.rate[, "Porites cylindrica"]
  )
)

oob_error_fungi_sp$Species <- factor(
  oob_error_fungi_sp$Species,
  levels = c(
    "Acropora millepora",
    "Montipora digitata",
    "Pocillopora damicornis",
    "Porites cylindrica",
    "OOB"
  )
)

p_oob_error_fungi_sp <- ggplot(data=oob_error_fungi_sp, aes(x=Trees, y=Error)) +
  ggtitle("Model Errors") +
  theme_classic() +
  annotate("text", x=250, y=0.9, label="OOB estimate: 20.31%", size=3) +
  geom_line(aes(color=Species)) + 
  scale_color_manual(values = c(paletteCB4, "black"), # to make legend text italic by level
                     labels = c(expression(italic("Acropora millepora")),
                                expression(italic("Montipora digitata")),
                                expression(italic("Pocillopora damicornis")),
                                expression(italic("Porites cylindrica")),"OOB")) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 9, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0)) +
  theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_oob_error_fungi_sp

# ASSESS MODEL PERFORMANCE WITH PERMUTATIONS ------------------------------------------------------
set.seed(43)
perm_RF_fungi_sp <-
  rf.significance(
    x = RF_fungi_sp,
    xdata = otus_fungi_clr_species[, 1:(ncol(otus_fungi_clr_species) - 2)],
    nperm = 999,
    nmtry = 12,
    ntree = 501
  )
perm_RF_fungi_sp

# PLOTTING THE RESULTS-----------------------------------------------------------------------------
# Now let's create an MDS-plot to show how the samples are related to each
# other. Start by converting the proximity matrix into a distance matrix.
dist_fungi_sp <- as.dist(1 - RF_fungi_sp$proximity)
mds_fungi_sp <- cmdscale(dist_fungi_sp, eig = TRUE, x.ret = TRUE)

## calculate the percentage of variation that each MDS axis accounts for...
var_fungi_sp <-
  round(mds_fungi_sp$eig / sum(mds_fungi_sp$eig) * 100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
values_fungi_sp <- mds_fungi_sp$points
mds_fungi_sp <- data.frame(
  Sample = rownames(values_fungi_sp),
  X = values_fungi_sp[, 1],
  Y = values_fungi_sp[, 2],
  Species = otus_fungi_clr_species$Species
)

# Annotation position and hjust, vjust use to modify it
# xpos ypos        annotateText        hjust    vjust
# 1 -Inf -Inf Bottom Left (h0,v0)        0        0
# 2 -Inf  Inf    Top Left (h0,v1)        0        1
# 3  Inf -Inf  Bottom Right h1,v0        1        0
# 4  Inf  Inf     Top Right h1,v1        1        1

p_nmds_fungi_sp = ggplot(data=mds_fungi_sp, aes(x=X, y=Y, label=Sample, color=Species)) + 
  geom_point(size=1.2) +
  theme_classic() +
  xlab(paste("MDS1 - ", var_fungi_sp[1], "%", sep="")) +
  ylab(paste("MDS2 - ", var_fungi_sp[2], "%", sep="")) +
  ggtitle("Fungi - Species") +
  scale_color_manual(values = paletteCB4) +
  annotate("text", -Inf, -Inf,  hjust = - 0.2, vjust = - 1.2, 
           label="OOB = 20.31%", size=2.5) +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) + # to adjust decimals
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), 
        legend.text = element_text(size = 7.5, hjust = 0, face = "italic")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 0, size = 7,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 7, hjust = 0.5, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="bottom")

p_nmds_fungi_sp

# Identifying Important Features ------------------------------------------------------------------
imp_fungi_sp <- as.data.frame(RF_fungi_sp$importance)
imp_fungi_sp$features <- rownames(imp_fungi_sp)
imp_fungi_sp <- arrange(imp_fungi_sp, desc(MeanDecreaseAccuracy))
head(imp_fungi_sp)

identical(rownames(taxa_fungi_filt), rownames(imp_fungi_sp))
order_fungi_sp <-
  match(rownames(taxa_fungi_filt), rownames(imp_fungi_sp))
imp_fungi_sp <- imp_fungi_sp[order_fungi_sp, ]

imp_fungi_sp$Taxonomy <- taxa_fungi_filt$Taxonomy
imp_fungi_sp <-
  imp_fungi_sp[order(imp_fungi_sp$MeanDecreaseAccuracy, decreasing = TRUE), ]
head(imp_fungi_sp)

p_bar_fungi_sp = ggplot(data=imp_fungi_sp[1:20,]) + 
  geom_bar(aes(x= reorder(Taxonomy, -MeanDecreaseAccuracy),
               y= MeanDecreaseAccuracy), color="indianred", fill="indianred",stat="identity") +
  labs(title = "Top Fungi OTUs to\nclassify coral Species", 
       x= "OTU", y= "Mean Decrease\nin Accuracy") +
  theme_classic() +
  coord_flip() +
  #ylim(0, 0.06) +
  scale_y_continuous(breaks=c(0.00, 0.01, 0.02, 0.03), limits = c(0, 0.03)) +
  theme(legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.3, "cm")) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) + 
  theme(axis.text.x = element_text(angle = 45, size = 7 ,hjust = 0.5, vjust = 0.5)) +
  theme(axis.text.y = element_text(angle = 0, size = 5, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) + 
  guides(color=guide_legend(ncol=2)) +
  theme(legend.position="none") #theme(legend.position=c(0.8, 0.8))

p_bar_fungi_sp
# I can explain you the logic behind or you can look at the expanation 
# I included in Reid's paper below, it's pretty straightforward 
# https://www.frontiersin.org/articles/10.3389/fmicb.2020.01116/full


# INDICATOR TAXA ----------------------------------------------------------------------------------
# Indicator species multipatt() analysis or ANCOM or LEFSE, you decide.
# We can look at this later on, though...




# ***************************************************************************----------------------



# NEW ***************************************_------------------------------------------------------
heat_tree_ITS <- parse_tax_data(tib_ITS,
                                class_cols = "Taxonomy",
                                class_sep = c(";"))


names(heat_tree_ITS$data) <- "otu_count"
print(heat_tree_ITS)

sample_names_its <- intersect(colnames(heat_tree_ITS$data$otu_count), tib_meta_ITS$sample_id)
sample_names_its <- tib_meta_ITS[match(sample_names_its, tib_meta_ITS$sample_id), ]
sample_names_its

heat_tree_ITS$data$otu_count <- heat_tree_ITS$data$otu_count[, c("taxon_id", "OTU_ID", sample_names_its)]
heat_tree_ITS$data$otu_count

head(taxon_names(heat_tree_ITS), n = 20)

heat_tree_ITS$data$otu_prop <- calc_obs_props(heat_tree_ITS,
                                              data = "otu_count",
                                              cols = tib_meta_ITS$sample_id)
heat_tree_ITS$data$otu_prop 


heat_tree_ITS$data$tax_prop <- calc_taxon_abund(heat_tree_ITS, 
                                                data = "otu_prop",
                                                cols = tib_meta_ITS$sample_id)
heat_tree_ITS$data$tax_prop


heat_tree_ITS


# GARBAGE CODE TO REMOVE --------------------------------------
heat_tree_ITS$taxon_names() -> taxa_names_its
taxa_names_its <- as.data.frame(taxa_names_its)
taxa_names_its

heat_tree_ITS$n_obs() -> n_obs_its
n_obs_its <- as.data.frame(n_obs_its)
n_obs_its

heat_tree_ITS$data$tax_occ_type -> tax_occ_its_type
tax_occ_its_type <- as.data.frame(tax_occ_its_type)
rownames(tax_occ_its_type) <- tax_occ_its_type$taxon_id
tax_occ_its_type

heat_tree_ITS$data$tax_occ_mark -> tax_occ_mark
tax_occ_mark <- as.data.frame(tax_occ_mark)
rownames(tax_occ_mark) <- tax_occ_mark$taxon_id
tax_occ_mark

heat_tree_ITS$data$tax_prop-> tax_prop_its
tax_prop_its <- as.data.frame(tax_prop_its)
rownames(tax_prop_its) <- tax_prop_its$taxon_id
tax_prop_its

tax_prop_its$taxa_abund <-
  rowSums(tax_prop_its[, 2:ncol(tax_prop_its)]) / rowSums(tax_prop_its[, 2:ncol(tax_prop_its)])[1] *
  100

tax_prop_its

identical(rownames(taxa_names_its), rownames(n_obs_its))
identical(rownames(n_obs_its), rownames(tax_occ_its_type))
identical(rownames(tax_occ_its_type), rownames(tax_occ_mark))
identical(rownames(tax_occ_mark), rownames(tax_prop_its))

dplyr::bind_cols(list(taxa_names_its, n_obs_its, tax_occ_its_type, tax_occ_mark, tax_prop_its)) -> df_its_heat
df_its_heat <- as.data.frame(df_its_heat[, c(1,2,4,5, 7:10,84)])
head(df_its_heat)
(df_its_heat$'in'* 100)/36






taxa_names_its -> df_its_heat
df_its_heat$n_obs_its <- n_obs_its$n_obs_its
df_its_heat









identical(rownames(df_its_heat), rownames(tax_prop_its))
sample_order_its <- match(rownames(tax_prop_its), rownames(df_its_heat))
df_its_heat <- df_its_heat[sample_order_its, ]

df_its_heat$taxa_abund <- tax_prop_its$taxa_abund
arrange(df_its_heat, desc(n_obs_its)) -> df_its_heat
df_its_heat


colnames(df_its_heat) <- c("OTU_number", "Taxon", "Taxon_id", "In", "Out")
df_its_heat <- df_its_heat[, c(3,2,1,4,5)]
head(df_its_heat)









heat_tree_ITS$n_obs() -> n_obs_its
n_obs_its <- as.data.frame(n_obs_its)
n_obs_its

heat_tree_ITS$taxon_names() -> taxa_names_its
taxa_names_its <- as.data.frame(taxa_names_its)
taxa_names_its

identical(rownames(n_obs_its), rownames(taxa_names_its))
cbind(n_obs_its, taxa_names_its) -> df_its_heat
df_its_heat

heat_tree_ITS$data$tax_occ -> tax_occ_its
tax_occ_its <- as.data.frame(tax_occ_its)
rownames(tax_occ_its) <- tax_occ_its$taxon_id
tax_occ_its

identical(rownames(df_its_heat), rownames(tax_occ_its))
cbind(df_its_heat, tax_occ_its) -> df_its_heat

colnames(df_its_heat) <- c("OTU_number", "Taxon", "Taxon_id", "In", "Out")
df_its_heat <- df_its_heat[, c(3,2,1,4,5)]
head(df_its_heat)





# ----------------------------

taxonomy_ITS_uparse_R1
rownames(taxonomy_ITS_uparse_R1) <- paste("F", rownames(taxonomy_ITS_uparse_R1), sep="")
taxonomy_ITS_uparse_R1[rownames(old_taxa_its), ] -> new_taxa_ITS
identical(rownames(new_taxa_ITS), rownames(old_taxa_its))

new_taxa_ITS -> taxa_ITS

nrow(taxa_ITS[taxa_ITS$Kingdom=="Plantae",])
taxa_ITS[taxa_ITS$Kingdom=="Plantae",]
otu_ITS[rownames(otu_ITS)=="FOTU_173",]
metadata_ITS[rownames(metadata_ITS)=="sample74",]



#taxa_16s -> old_taxa_16s

taxonomy_16s_uparse_R1
rownames(taxonomy_16s_uparse_R1) <- paste("B", rownames(taxonomy_16s_uparse_R1), sep="")
taxonomy_16s_uparse_R1[rownames(taxa_16s), ] -> new_taxa_16s
identical(rownames(new_taxa_16s), rownames(taxa_16s))
new_taxa_16s -> taxa_16s





