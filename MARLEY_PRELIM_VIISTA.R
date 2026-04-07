# Quick preliminary analysis of Marley data 
# dada2 pipeline

# SAMPLE INFO----
# Sequenced at MSU on an AVITI instrument **Need to account for this with the error models
  ## 20260306_16SV4_PE250
  ## 20260313_16SV4_PE250
  ## 20260320_16SV4_PE250 

# SETUP----
# Load libraries 
library(dada2); library(patchwork); library(ggplot2)

# If you need to install bioconductor and dada2 see below
# library(BiocManager)
# 
# # bioconductor
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.22")
# 
# # dada2
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2", version = "3.22")

# Set path (need to make this reproducible with here eventually)
path="/Users/ashley/Documents/r_local/draft-16S-workflow/MARLEY_DATA/Seq3_VIISTA/"
list.files(path)

# Read in your forward and reverse filenames
# Need to make sure this matches what your sequences look like

# 1st run----
# (From 20260320_16SV4_PE250) 

# Add the .gz if your files are still gunzipped
# make sure that it matches what comes after your "_"
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)
saveRDS(sample.names,"MARLEY_DATA/Intermediate/filenames.rds") # use this to save file at each step

# Inspect quality plots
P1 <- plotQualityProfile(fnFs[1:4])
P2 <- plotQualityProfile(fnRs[1:4])
P1 + P2

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trimming primers
# The primers used are: 16S V4 forward (515f): GTGCCAGCMGCCGCGGTAA 16S V4 reverse (806r): GGACTACHVGGGTWTCTAAT
# According to here: https://rtsf.natsci.msu.edu/genomics/technical-documents/amplicon-metagenomic-guide.aspx#16S-V4
# Will need to come back and use cutAdapt most likely
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-20

# Ensure that the lengths of fnFs and fnRs are equal
if(length(fnFs) != length(fnRs)) {
  stop("The number of forward and reverse reads must be equal.")
}

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN),
                     truncLen = c(240, 210),  # Example values
                     maxN = 0,  # Discard reads with any Ns
                     maxEE = c(2, 2),  # Allow max 2 expected errors
                     truncQ = 2,  # Truncate reads at the first quality score <= 2
                     compress = TRUE,  # Compress the output files
                     verbose = TRUE, 
                     multithread = TRUE) 
saveRDS(out,"MARLEY_DATA/Intermediate/out.rds")

# Learn error rates----
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"MARLEY_DATA/Intermediate/errF.rds")
saveRDS(errR,"MARLEY_DATA/Intermediate/errR.rds")

# Visualize the error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplication----
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepFs, "MARLEY_DATA/Intermediate/derepFs.rds")
saveRDS(derepRs, "MARLEY_DATA/Intermediate/derepRs.rds")
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Run dada2 algorithm----
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "MARLEY_DATA/Intermediate/dadaFs.rds")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs,"MARLEY_DATA/Intermediate/dadaRs.rds")
dadaFs[[1]]

# Output
# dada-class: object describing DADA2 denoising results
# 2122 sequence variants were inferred from 198566 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# DADA2 concluded that your dataset contains 2,122 ASVs that it believes represent 
# real biological sequences rather than sequencing/PCR error.

# Merge----
# Merge forward and reverse reads
# Name run1 so can track moving forward
mergers_run1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers_run1[[1]])
saveRDS(mergers_run1, "MARLEY_DATA/Intermediate/mergers_run1.rds")


# Construct a sequence table
seqtab_run1 <- makeSequenceTable(mergers_run1)
dim(seqtab_run1)
table(nchar(getSequences(seqtab_run1)))
saveRDS(seqtab_run1, "MARLEY_DATA/Intermediate/seqtab_dna-run1.rds")

# Continuing on with the protocol for now to get a complete dataset
# Eventually need to analyze seq run 2 and combine

# Remove chimeras----
# Chimeras are common but shouldn't be more than 3-4% 
# If high, then remove primers and redo analysis
seqtab.nochim <- removeBimeraDenovo(seqtab_run1, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab_run1) # This means that ~8.5% were chimeric
saveRDS(seqtab.nochim, "MARLEY_DATA/Intermediate/seqtab.nochim.rds")
# seqtab.nochim <- readRDS("seqtab-nochim.rds")

# Do a quick check in----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers_run1, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assigning taxonomy----
# Be sure to put the database in your working directory or another known location
# You can download proper database here: https://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/r_local/draft-16S-workflow/DB/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
saveRDS(taxa, "MARLEY_DATA/Intermediate/taxa-table1.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Now it is ready to import into Phyloseq

# Read in----
# Uncomment if you need to read the files back in
# sample.names <- readRDS("MARLEY_DATA/Intermediate/filenames.rds")
# out <- readRDS("MARLEY_DATA/Intermediate/out.rds")
# errF <- readRDS("MARLEY_DATA/Intermediate/errF.rds")
# errR <- readRDS("MARLEY_DATA/Intermediate/errR.rds")
# derepFs <- readRDS("MARLEY_DATA/Intermediate/derepFs.rds")
# derepRs <- readRDS("MARLEY_DATA/Intermediate/derepRs.rds")
# dadaFs <- readRDS("MARLEY_DATA/Intermediate/dadaFs.rds")
# dadaRs <- readRDS("MARLEY_DATA/Intermediate/dadaRs.rds")
# mergers_run1 <- readRDS("MARLEY_DATA/Intermediate/mergers_run1.rds")
# seqtab.nochim <- readRDS("MARLEY_DATA/Intermediate/seqtab.nochim.rds")
# taxa <- readRDS("MARLEY_DATA/Intermediate/taxa-table1.rds")
