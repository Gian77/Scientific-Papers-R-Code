
## New script for barplots and relative abundances


biom_ITS_symb = import_biom("mycorrhizal_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_symb) <- map_ITS
colnames(tax_table(biom_ITS_symb)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_symb <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_symb <- merge_phyloseq(biom_ITS_symb, otus_rep_ITS_symb)
biom_ITS_symb
sample_data(biom_ITS_symb)
tax_table(biom_ITS_symb)
refseq(biom_ITS_symb)


biom_ITS_symb_gen = tax_glom(biom_ITS_symb, "Genus")
biom_ITS_symb_gen
biom_ITS_symb_gen_merg = merge_samples(biom_ITS_symb_gen, "Site1")
biom_ITS_symb_gen_merg
sample_data(biom_ITS_symb_gen_merg)
otu_table(biom_ITS_symb_gen_merg)

sample_data(biom_ITS_symb_gen_merg)$Site1 <- levels(sample_data(biom_ITS_symb_gen)$Site1)
biom_ITS_symb_gen_merg_ab = transform_sample_counts(biom_ITS_symb_gen_merg, function(x) 100 * x/sum(x))
biom_ITS_symb_gen_merg_ab
otu_table(biom_ITS_symb_gen_merg_ab)
tax_table(biom_ITS_symb_gen_merg_ab)
sum(otu_table(biom_ITS_symb_gen_merg_ab))

source("../plot_ordered_bar_byGian.R")
p1 <- plot_ordered_bar(biom_ITS_symb_gen_merg_ab, x = "Site1", fill="Genus", leg_size = 0.4, title="(A)")
p1

myco_bar = p1 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Acaulospora" = "#560d0d",
                               "Amanita" = "#a35151",
                               "Ambispora" = "#dba4a4", 
                               "Cenococcum" = "#cc1c1c",
                               "Ceratobasidium"="#111b77",
                               "Chloridium"="#283dff", 
                               "Chromelosporium"="#636bb7",
                               "Claroideoglomus"="#bfc5ff",
                               "Corymbiglomus"="#014443",
                               "Dentiscutata"="#195637",
                               "Diversispora"="#117744",
                               "Entoloma"="#60ffaf",
                               "Gigaspora"="#b7ffdb",
                               "Glomus"="#825121",
                               "Inocybe"="#ea7f17",
                               "Oidiodendron"="#fcb067",
                               "Paraglomus"="#ffe8d3",
                               "Peziza"="#d8d6d4",
                               "Racocetra"="#82807f",
                               "Ramaria" = "#3f3e3d",
                               "Rhizophagus" = "#000000",
                               "Rhizopogon" = "#5b5b19", 
                               "Scleroderma" = "#fcfc00",
                               "Scutellospora"="#ffff9e",
                               "Septoglomus"="#ffb7ef", 
                               "Serendipita"="#fa7efc",
                               "Stephanospora"="#ae09ea",
                               "Tomentella"="#521899",
                               "Wilcoxina"="#1e0047",
                               "unidentified1"="white",
                               "unidentified2"="white",
                               "unidentified2"="white")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Mycorrhizal", x="", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


#theme(legend.position="none")
  
  #theme(legend.text=element_text(size=10), legend.key.size = unit(1, "mm"))
  #theme(legend.position = "bottom") +
  #theme(legend.direction = "horizontal")

## now endophytic fungi

biom_ITS_endo = import_biom("endophyte_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_endo) <- map_ITS
colnames(tax_table(biom_ITS_endo)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_endo <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_endo <- merge_phyloseq(biom_ITS_endo, otus_rep_ITS_endo)
biom_ITS_endo
sample_data(biom_ITS_endo)
tax_table(biom_ITS_endo)
refseq(biom_ITS_endo)

biom_ITS_endo_gen = tax_glom(biom_ITS_endo, "Genus")
biom_ITS_endo_gen_merg = merge_samples(biom_ITS_endo_gen, "Site1")
sample_data(biom_ITS_endo_gen_merg)$Site1 <- levels(sample_data(biom_ITS_endo_gen)$Site1)
biom_ITS_endo_gen_merg_ab = transform_sample_counts(biom_ITS_endo_gen_merg, function(x) 100 * x/sum(x))

p2 <- plot_ordered_bar(biom_ITS_endo_gen_merg_ab, x = "Site1", fill="Genus", leg_size = 0.4, title="(B)")
p2

endo_bar = p2 + theme(legend.text = element_text(face="italic")) + 
  scale_fill_manual(values = c("Annulohypoxylon" = "#560d0d",
                               "Anthostomella" = "#a35151",
                               "Camarosporium" = "#dba4a4", 
                               "Capronia" = "#cc1c1c",
                               "Chalara"="#111b77",
                               "Cladosporium"="#283dff", 
                               "Colletotrichum"="#636bb7",
                               "Cytospora"="#bfc5ff",
                               "Diaporthe"="#014443",
                               "Kabatiella"="#195637",
                               "Lecythophora"="#117744",
                               "Microdochium"="#60ffaf",
                               "Paecilomyces"="#b7ffdb",
                               "Periconia"="#825121",
                               "Phialophora"="#ea7f17",
                               "Phomopsis"="#fcb067",
                               "Torula"="#ffe8d3",
                               "Trichoderma"="#d8d6d4",
                               "Williopsis"="#82807f",
                               "Xylaria" = "#3f3e3d")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Endophytic", x="Habitat", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")

## now pathogens

biom_ITS_pat = import_biom("pathogen_genus.biom")
map_ITS = import_qiime_sample_data("mapping_ITS.txt")
sample_data(biom_ITS_pat) <- map_ITS
colnames(tax_table(biom_ITS_pat)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
otus_rep_ITS_pat <- readDNAStringSet("otus_ITS.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
biom_ITS_pat <- merge_phyloseq(biom_ITS_pat, otus_rep_ITS_pat)
biom_ITS_pat
sample_data(biom_ITS_pat)
tax_table(biom_ITS_pat)
refseq(biom_ITS_pat)

biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
biom_ITS_pat_gen_merg = merge_samples(biom_ITS_pat_gen, "Site1")
biom_ITS_pat_gen_merg
tax_table(biom_ITS_pat_gen_merg)
otu_table(biom_ITS_pat_gen_merg)

sample_data(biom_ITS_pat_gen_merg)$Site1 <- levels(sample_data(biom_ITS_pat_gen)$Site1)
biom_ITS_pat_gen_merg_ab = transform_sample_counts(biom_ITS_pat_gen_merg, function(x) 100 * x/sum(x))


badTaxa = c("OTU_1494","OTU_2750","OTU_2559","OTU_3578","OTU_3405","OTU_3221","OTU_5205","OTU_3915","OTU_2624","OTU_3451"
            ,"OTU_3118","OTU_6904","OTU_4023","OTU_5200","OTU_5862","OTU_4387","OTU_5089","OTU_4523","OTU_4615","OTU_4945")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

biom_ITS_pat_gen_merg_ab_prun = pop_taxa(biom_ITS_pat_gen_merg_ab, badTaxa)
biom_ITS_pat_gen_merg_ab_prun

#top26 <- names(sort(taxa_sums(biom_ITS_pat_gen_merg_ab), TRUE)[1:26]) # sampling first 100 taxa
#biom_ITS_pat_gen_merg_ab_26 <- prune_taxa(top26, biom_ITS_pat_gen_merg_ab)

p3 <- plot_ordered_bar(biom_ITS_pat_gen_merg_ab_prun, x = "Site1", fill="Genus", leg_size = 0.4, title="(C)")
p3

pat_bar = p3 + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Aspergillus" = "#560d0d",
                               "Calonectria" = "#a35151",
                               "Ceratocystis" = "#dba4a4", 
                               "Clonostachys" = "#cc1c1c",
                               "Collophora"="#111b77",
                               "Coniothyrium"="#283dff", 
                               "Curvularia"="#636bb7",
                               "Cylindrocarpon"="#bfc5ff",
                               "Cylindrocladiella"="#014443",
                               "Endomelanconiopsis"="#195637",
                               "Fusarium"="#117744",
                               "Ganoderma"="#60ffaf",
                               "Gibberella"="#b7ffdb",
                               "Gibellulopsis"="#825121",
                               "Glomerella"="#ea7f17",
                               "Haematonectria"="#fcb067",
                               "Lasiodiplodia"="#ffe8d3",
                               "Magnaporthe"="#d8d6d4",
                               "Monographella"="#82807f",
                               "Mycoleptodiscus" = "#3f3e3d",
                               "Pestalotiopsis" = "#000000",
                               "Phaeoacremonium" = "#5b5b19", 
                               "Phoma" = "#c1c15d",
                               "Pyricularia"="#fcfc00",
                               "Rhizophydium"="#ffff9e", 
                               "Sclerotinia"="#ffb7ef",
                               "Spizellomyces"="#fa7efc",
                               "Thanatephorus"="#ae09ea",
                               "Trichothecium"="#521899",
                               "Volutella"="#1e0047")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Pathogenic", x="", y="") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")


library("ggpubr")
ggarrange(myco_bar,endo_bar,pat_bar,
          labels = c("A", "B","C"),
          widths = c(1,1,1.05),
          align = "h", ncol = 3, nrow = 1,
          legend="right",
          common.legend = FALSE)



biom_16S_gen = tax_glom(biom_16s_cv, "Genus")
biom_16S_gen
biom_16S_gen_merg = merge_samples(biom_16S_gen, "Site1")
biom_16S_gen_merg
sample_data(biom_16S_gen_merg)
otu_table(biom_16S_gen_merg)

sample_data(biom_16S_gen_merg)$Site1 <- levels(sample_data(biom_16S_gen)$Site1)
biom_16S_gen_merg_ab = transform_sample_counts(biom_16S_gen_merg, function(x) 100 * x/sum(x))
biom_16S_gen_merg_ab
otu_table(biom_16S_gen_merg_ab)
tax_table(biom_16S_gen_merg_ab)
sum(otu_table(biom_16S_gen_merg_ab))


badTaxa = c("OTU_30")

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}

biom_16s_gen_merg_ab_prun = pop_taxa(biom_16S_gen_merg_ab, badTaxa)
biom_16s_gen_merg_ab_prun


top_30_bact <- names(sort(taxa_sums(biom_16s_gen_merg_ab_prun), TRUE)[1:30]) # first 30 top OTUs 
top_30_bact
biom_16s_top30 <- prune_taxa(top_30_bact, biom_16s_gen_merg_ab_prun)
biom_16s_top30
tax_table(biom_16s_top30)

# get sequnces for BLAST against GenBank
write.csv(refseq(biom_16S_gen) , file = "refseq_16s_gen.csv")
# manual correction of the tax_table - OPTIONAL!
write.csv(tax_table(biom_16s_top30) , file = "tax_16s_top30.csv")
tax_biom_16s_top30 <- read.csv("tax_16s_top30.csv", header=T, row.names =1)
tax_table(biom_16s_top30) <- tax_table(as.matrix(tax_biom_16s_top30))
tax_table(biom_16s_top30)


source("../plot_ordered_bar_byGian.R")

p_bact <- plot_ordered_bar(biom_16s_top30, x = "Site1", fill="Genus", leg_size = 0.4, title="(D)")
p_bact

p_bact + theme(legend.text = element_text(face="italic")) +
  scale_fill_manual(values = c("Actinoallomurus" = "#560d0d",
                               "Aeribacillus" = "#a35151",
                               "Bacillus" = "#dba4a4", 
                               "Bradyrhizobium" = "#cc1c1c",
                               "Burkholderia"="#111b77",
                               "Candidatus Koribacter"="#283dff", 
                               "Candidatus Nitrososphaera"="#636bb7",
                               "Candidatus Solibacter"="#bfc5ff",
                               "Candidatus Xiphinematobacter"="#014443",
                               "Flavobacterium"="#195637",
                               "Gemmata"="#117744",
                               "Geobacillus"="#60ffaf",
                               "Hyphomicrobium"="#b7ffdb",
                               "Kaistobacter"="#825121",
                               "Ktedonobacter"="#ea7f17",
                               "Mesorhizobium"="#fcb067",
                               "Mycobacterium"="#ffe8d3",
                               "Nitrospira"="#d8d6d4",
                               "Nocardioides"="#82807f",
                               "Pedomicrobium" = "#3f3e3d",
                               "Phenylobacterium" = "#000000",
                               "Pirellula" = "#5b5b19", 
                               "Planctomyces" = "#c1c15d",
                               "Pseudonocardia"="#fcfc00",
                               "Rhizobium"="#ffff9e", 
                               "Rhodoplanes"="#ffb7ef",
                               "[Spartobacteria] DA101"="#fa7efc",
                               "Sphingomonas"="#ae09ea",
                               "Steroidobacter"="#521899",
                               "Streptomyces"="#1e0047")) + 
  theme(panel.background = element_blank()) + 
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") +
  guides(fill=guide_legend(ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  theme(legend.position="right") + guides(fill=guide_legend(reverse = TRUE, ncol=1,keywidth = 0.7, keyheight = 0.7)) +
  #facet_wrap(~Phylum, ncol=4) +
  #facet_grid(~Phylum, nrow=2) +
  ylim(0, 100) +
  labs(title="Bacteria and Archaea", x="Habitat", y="Relative abundance") +
  theme(panel.background = element_blank()) + 
  theme(axis.text.x = element_text(vjust=0.5, size=10)) +
  theme(axis.text.y = element_text(hjust=0.5, size=10)) +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  theme(axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  theme(legend.position="right")





#################################################################################################
#################################################################################################

## Calculate Relative Abundacnes

biom_ITS_cv
biom_ITS_phy = tax_glom(biom_ITS_cv, "Phylum")
otu_ITS_phy = taxa_sums(biom_ITS_phy)/sum(taxa_sums(biom_ITS_phy))*100
tax_ITS_phy <- data.frame(tax_table(biom_ITS_phy))
tax_ITS_phy <- tax_ITS_phy[c(1:2)]
tax_ITS_phy$abundance <- as.vector(otu_ITS_phy)
tax_ITS_phy <- tax_ITS_phy[order(tax_ITS_phy$abundance, decreasing = TRUE),] 
tax_ITS_phy

biom_LSU_cv
biom_LSU_phy = tax_glom(biom_LSU_cv, "Phylum")
otu_LSU_phy = taxa_sums(biom_LSU_phy)/sum(taxa_sums(biom_LSU_phy))*100
tax_LSU_phy <- data.frame(tax_table(biom_LSU_phy))
tax_LSU_phy <- tax_LSU_phy[c(1:2)]
tax_LSU_phy$abundance <- as.vector(otu_LSU_phy)
tax_LSU_phy <- tax_LSU_phy[order(tax_LSU_phy$abundance, decreasing = TRUE),] 
tax_LSU_phy

biom_16s_cv
biom_16s_phy = tax_glom(biom_16s_cv, "Phylum")
otu_16s_phy = taxa_sums(biom_16s_phy)/sum(taxa_sums(biom_16s_phy))*100
tax_16s_phy <- data.frame(tax_table(biom_16s_phy))
tax_16s_phy <- tax_16s_phy[c(1:2)]
tax_16s_phy$abundance <- as.vector(otu_16s_phy)
tax_16s_phy <- tax_16s_phy[order(tax_16s_phy$abundance, decreasing = TRUE),] 
tax_16s_phy




# different guilds

biom_ITS_symb
tax_table(biom_ITS_symb)
otu_table(biom_ITS_symb)

otu_ITS_symb_gen = taxa_sums(biom_ITS_symb)/sum(taxa_sums(biom_ITS_symb))*100
tax_ITS_symb_gen <- data.frame(tax_table(biom_ITS_symb))
tax_ITS_symb_gen <- tax_ITS_symb_gen[c(1:6)]
tax_ITS_symb_gen$abundance <- as.vector(otu_ITS_symb_gen)
tax_ITS_symb_gen <- tax_ITS_symb_gen[order(tax_ITS_symb_gen$abundance, decreasing = TRUE),] 
tax_ITS_symb_gen

otu_ITS_endo_gen = taxa_sums(biom_ITS_endo)/sum(taxa_sums(biom_ITS_endo))*100
tax_ITS_endo_gen <- data.frame(tax_table(biom_ITS_endo))
tax_ITS_endo_gen <- tax_ITS_endo_gen[c(1:6)]
tax_ITS_endo_gen$abundance <- as.vector(otu_ITS_endo_gen)
tax_ITS_endo_gen <- tax_ITS_endo_gen[order(tax_ITS_endo_gen$abundance, decreasing = TRUE),] 
tax_ITS_endo_gen

biom_ITS_pat_gen = tax_glom(biom_ITS_pat, "Genus")
otu_ITS_pat_gen = taxa_sums(biom_ITS_pat_gen)/sum(taxa_sums(biom_ITS_pat_gen))*100
tax_ITS_pat_gen <- data.frame(tax_table(biom_ITS_pat_gen))
tax_ITS_pat_gen <- tax_ITS_pat_gen[c(1:6)]
tax_ITS_pat_gen$abundance <- as.vector(otu_ITS_pat_gen)
tax_ITS_pat_gen <- tax_ITS_pat_gen[order(tax_ITS_pat_gen$abundance, decreasing = TRUE),] 
tax_ITS_pat_gen




## alpha diversity

plot_alpha_its <- plot_richness(biom_ITS_cv, x = "Site1", color = "Site1", measures = c("Observed","Simpson"))
plot_alpha_lsu <- plot_richness(biom_LSU_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))
plot_alpha_16s <- plot_richness(biom_16s_cv, x = "Site1", color = "Site1", measures = c("Observed","Shannon"))

rich_its = plot_alpha_its + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")

rich_lsu = plot_alpha_lsu + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")

rich_16s = plot_alpha_16s + geom_boxplot(outlier.colour="black", outlier.shape = 8) +
scale_colour_manual(values=palette_CB9) +
theme(legend.position="none")
rich_16s

# multiplot
source("../multiplot.R")
multiplot(rich_its, rich_lsu, rich_16s, cols=3)
dev.off()


## plotting heatmaps
library("superheat")
library("vegan")
data(dune)
dune

superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          pretty.order.cols = TRUE)

superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          pretty.order.rows = TRUE, 
          heat.pal = c("#b35806", "white", "#542788"))


superheat(t(dune), left.label.size = 0.3, bottom.label.size = 0.1,
          heat.pal = c("#b35806", "white", "#542788"), row.dendrogram = TRUE)




