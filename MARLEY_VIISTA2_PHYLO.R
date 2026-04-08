# Quick preliminary analysis of Marley data 

# SAMPLE INFO----
# dada2 pipeline run for 20260320_16SV4_PE250 (Seq3_VIISTA) and 20260313_16SV4_PE250
# This script assumes you have run "MARLEY_PRELIM_VIISTA.R" and merged runs
# This analysis is for VIISTA 2 only

# SETUP----
# Load libraries 
library(dada2); library(patchwork); library(ggplot2); library(phyloseq); library(data.table); library(dplyr)

# If you need to download phyloseq
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")

## Functions----
pretty.theme <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=14, color = "black", angle = 0, hjust = 1),
          axis.text.y=element_text(size=14, color = "black"),
          axis.title.x=element_text(size=20, color = "black"),             
          axis.title.y=element_text(size=20, color = "black"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), units = , "cm"),
          plot.title = element_text(size=20))
}

# Phyloseq setup----
## Import----
meta = read.table("PHYLOSEQ-INPUT/sample_metadata_all_REAL_redo.txt", header = TRUE, sep = "\t", row.names = 1)
taxa = readRDS("MARLEY_DATA/Intermediate/taxa-table-all.rds")
taxa_matrix <- as.matrix(taxa)
tax_phy <- tax_table(taxa_matrix)
seqtab.nochim = readRDS("MARLEY_DATA/Intermediate/seqtab_nochim_allDNA.rds")
phy <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                sample_data(meta), tax_table(taxa))

# This worked by reading in the metadata file as a text file so just DO THIS

# ## Troubleshooting-------
# # Issues with getting the names to match up
# # Check sample names in each component
# otu_sample_names <- sample_names(seqtab.nochim)
# meta_sample_names <- sample_names(meta)
# tax_sample_names <- sample_names(taxa)  # If you have a taxonomic table
# 
# # Print sample names for inspection
# print(otu_sample_names)
# print(meta_sample_names)
# print(tax_sample_names)
# 
# # Read in the old way
# # Loading into Phyloseq-----------------
# mat = read.table("PHYLOSEQ-INPUT/ASVs_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
# tax = read.table("PHYLOSEQ-INPUT/ASVs_taxonomy.tsv", header = TRUE, sep = "\t", row.names = 1)
# meta = read.table("PHYLOSEQ-INPUT/sample_metadata.txt", header = TRUE, sep = "\t", row.names = 1)
# 
# mat = as.matrix(mat)
# tax = as.matrix(tax)
# 
# OTU = otu_table(mat, taxa_are_rows = TRUE)
# TAX = tax_table(tax)
# META = sample_data(meta)
# sample_names(OTU)
# sample_names(TAX) # This should be null
# tail(META)
# head(TAX)
# head(OTU)
# 
# phy = phyloseq(OTU,TAX,META)
# sample_names(phy)

# Exploration----
# Look at statistics related to your data
num.reads = sum(sample_sums(phy))
lowest.sam = sort(sample_sums(phy)) 
mean.seq = mean(sample_sums(phy))  
std.seq = sd(sample_sums(phy))/sqrt(206) 
med.seq = median(sample_sums(phy))
phy.summary <- data.frame(num.reads, mean.seq, std.seq, med.seq)
phy.summary     

# Information on a per sample basis
seq.dt = data.table(as(sample_data(phy), "data.frame"),
                    TotalReads = sample_sums(phy), keep.rownames = TRUE)
seq.dt # information per sample after filtering out taxa 
write.csv(seq.dt, "MARLEY_DATA/Intermediate/reads-per-sample-VIISTA2.csv")

# Sample removal----
# Need to go back to run decontam
# Removing blank and mockDNA 
phy_fil <- prune_samples(!(sample_names(phy) %in% c("WBlank", "ZymoMockDNA")), phy)
phy_fil
sample_names(phy_fil)

## Filtering----
phy_fil %>%
  subset_taxa(Family != "Mitochondria",
              Order != "Chloroplast") -> phy_nochloro
phy_fil # 15660
phy_nochloro # 7263

phy_fil %>%
  subset_taxa(Family != "Mitochondria") -> phy_withchloro
phy_fil 
phy_withchloro 

# For now will proceeed with phy_nochloro

# Filter out low abundances
# lowabundnames = filter_taxa(phy_nochloro, function(x) sum(x) > 10) # which have <10 reads across all samples
# phy_nochloro_filt = prune_taxa(lowabundnames, phy_nochloro) # Remove taxa that meets these conditions
# ntaxa(phy_nochloro_filt) # 4260 remain

# Calculate relative abundance
per.f.final = transform_sample_counts(phy_nochloro, function (x) x/sum(x)*100)
saveRDS(per.f.final, "MARLEY_DATA/Intermediate/per_f_final.rds")

# Plotting----
# Ordination----
# Set depth to a factor
sample_data(per.f.final)$depth <- as.factor(sample_data(per.f.final)$depth)
sample_data(per.f.final)$station <- as.factor(sample_data(per.f.final)$station)

# NMDS by depth
BC_distance <- phyloseq::distance(per.f.final, "bray") 
bcOrd <- ordinate(per.f.final, "PCoA", BC_distance)

# Plot by depth
NMDS_plot_depth1 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = depth), size = 4, alpha = 0.9) +
  scale_color_viridis_d(option = "viridis", direction = -1) +
  pretty.theme() 
  # labs(fill = "depth")
NMDS_plot_depth1

# Plot by depth
NMDS_plot_depth1 <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = depth, shape = station), size = 4, alpha = 0.9) +
  scale_color_viridis_d(option = "viridis", direction = -1) +
  pretty.theme() 
# labs(fill = "depth")
NMDS_plot_depth1

# Plot by site
NMDS_plot_site <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = station), size = 4, alpha = 0.9) +
  scale_color_viridis_d(option = "rocket", direction = -1) +
  pretty.theme() 
NMDS_plot_site

# Plot by chl
NMDS_plot_chla <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color =chla), size = 4, alpha = 0.9) +
  scale_color_viridis_c(option = "mako", direction = -1) +
  pretty.theme() +
  labs(fill = "chla")
NMDS_plot_chla

# Plot by date
NMDS_plot_date <- plot_ordination(per.f.final, bcOrd, type = "samples") +
  geom_point(aes(color = collect_date), size = 4, alpha = 0.9) +
  scale_color_viridis_d(option = "magma", direction = -1) +
  pretty.theme() 
NMDS_plot_date

# Shannon Diversity----
# Estimate shannon diversity for the entire dataset 
# Then add to the metadata
shannon_df_all <- estimate_richness(phy_nochloro, measures = "Shannon")
shannon_df_all
write.csv(shannon_df_all, "MARLEY_DATA/Intermediate/shannon_all_NEW.csv")
# added metadata columns, do this better in the future

# read edited file back in
shannon <- read.csv("MARLEY_DATA/Intermediate/shannon_all_NEW_wMetadata.csv", header = T)

shannon$depth <- as.factor(shannon$depth)
shannon$station <- as.factor(shannon$station)

shannon_DNA_plot <- ggplot(shannon, aes(y = Shannon, x = depth, group = depth)) +
  geom_boxplot(aes(group = depth), outlier.shape = NA, alpha = 0.5) + 
  geom_jitter(aes(color = depth)) +
  scale_color_viridis_c(option = "viridis", direction = -1) +
  # geom_point(aes(fill = Site), shape = 21, color = "black", size = 2, stroke = 0.5) +  # Black outline for points
  #scale_fill_viridis_d(option = "viridis", direction = 1) +
  coord_flip() +
  scale_x_reverse() + 
  pretty.theme() +
  labs(x = "Depth (m)", y = "Shannon Diversity") +
  theme(legend.position = "none")
shannon_DNA_plot

shannon_DNA_plot_site <- ggplot(shannon, aes(y = Shannon, x = station)) +
  geom_boxplot(aes(group = station), outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = station)) + 
  scale_color_viridis_d(option = "rocket", direction = -1) +
  pretty.theme() 
 # labs(x = "Site", y = "Shannon Diversity")
shannon_DNA_plot_site

shannon_DNA_plot_chla <- ggplot(shannon, aes(y = Shannon, x = chla)) +
  geom_point(alpha = 0.5, size = 2) + 
 # geom_smooth(method = "lm") +
  # geom_point(aes(fill = Site), shape = 21, color = "black", size = 2, stroke = 0.5) +  # Black outline for points
  # coord_flip() +
  # scale_x_reverse() + 
  pretty.theme() +
  labs(x = "Chla (ug/L)", y = "Shannon Diversity")
shannon_DNA_plot_chla

mod <- lm(Shannon ~ chla, data = shannon)
summary(mod)

# Barplots----
## Combined plot with facet----
# Setup palette
pal_27 <- c(
  '#9e0142','#44AAAA','#fdae61','#808080','#117777','#4477AA','darkblue','#b15928',
  '#a6cee3','#771155','#CC99BB','#fee08b','#abdda4','#5e4fa2','#e7298a','#66c2a5',
  '#1f78b4',  # vibrant blue
  '#33a02c',  # bright green
  '#ff7f00',  # orange
  '#cab2d6',  # lavender
  '#a6d854',  # lime green
  '#ffd92f',  # bright yellow
  '#fb8072',  # coral pink
  '#8dd3c7',  # pale aqua
  '#bc80bd',  # muted purple
  '#b3de69',  # soft green
  '#fccde5'   # pink blush
)

phy.rel <- transform_sample_counts(phy_withchloro, function(x) x / sum(x) * 100)
phy.phylum <- tax_glom(phy.rel, taxrank = "Phylum")
df <- psmelt(phy.phylum)

# dna_barplot <- ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
#   geom_bar(stat = "identity") +
#   # coord_flip() +
#   # scale_x_reverse(breaks = seq(0,300, by = 50)) +
#   # scale_fill_manual(values = pal_shared) +
#   pretty.theme() +
#   theme(
#     axis.text = element_text(size = 12, color = "black"),
#     legend.position = "right",
#     legend.title = element_text(size = 13),
#     legend.text = element_text(size = 11)
#   ) +
#   labs(
#     x = "Depth (cm)",
#     y = "Relative Abundance (%)",
#     fill = "Phylum"
#   ) +
#   facet_wrap(.~depth)
# dna_barplot

# Ok try and average by depth
# Take average
df_avg_dna <- df %>%
  group_by(depth, Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  group_by(depth) %>%
  mutate(rel_abundance = mean_abundance / sum(mean_abundance) * 100) %>%
  filter(rel_abundance > 1) %>%
  ungroup()

# This is the plot averaged by depth
pal_12 <- c('#44AAAA','#fdae61','#808080','#4477AA','darkblue','#771155',
            '#ff7f00','#fee08b','#abdda4','#5e4fa2','#e7298a','#CC99BB')

pal_9 <- c('#44AAAA','#fdae61','#4477AA','#771155','#fee08b','#e7298a','#abdda4','#5e4fa2','#CC99BB')

dna_barplot <- ggplot(df_avg_dna, aes(x = depth, y = rel_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_x_reverse(breaks = seq(0,100, by = 10)) +
  scale_fill_manual(values = pal_9) +
  pretty.theme() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  ) +
  labs(
    x = "Depth (cm)",
    y = "Relative Abundance (%)",
    fill = "Phylum") 
dna_barplot

dna_barplot_bySite <- ggplot(df_avg_dna, aes(x = Site, y = rel_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = pal_8) +
  pretty.theme() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  ) +
  labs(
    x = "Depth (cm)",
    y = "Relative Abundance (%)",
    fill = "Phylum") 
dna_barplot_bySite

df_avg_dna_rare <- df %>%
  group_by(depth, Phylum) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  group_by(depth) %>%
  mutate(rel_abundance = mean_abundance / sum(mean_abundance) * 100) %>%
  filter(rel_abundance > 1) %>%
  ungroup()




