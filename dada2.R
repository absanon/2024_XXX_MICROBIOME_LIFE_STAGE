# Author: Sanon Aboubakar
####### Install the required packages for your analysis
BiocManager::install("dada2")
BiocManager::install("DECIPHER")
BiocManager::install("microbiome")

install.packages("BiocManager")

install.packages("devtools")

library("devtools")

library(dada2)
packageVersion("dada2")

here::i_am("dada2.R")

############### Load your fastq data
path <- "C:/sanon"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq.gz", full.names = TRUE))
fnFs
fnRs

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

##################### Filter and trim
plotQualityProfile(fnFs[11:12])
plotQualityProfile(fnRs[11:12])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
  truncLen = c(290, 210),
  maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = FALSE
) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

##### Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim) / sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track
taxa <- assignTaxonomy(seqtab.nochim, "C:/sanon/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
taxa <- addSpecies(taxa, "C:/sanon/silva_species_assignment_v132.fa.gz")

# Removing sequence rownames for display only
taxa.print <- taxa

rownames(taxa.print) <- NULL

head(taxa.print)

taxa.print

tail(taxa.print)

View(taxa.print)

library(DECIPHER)
packageVersion("DECIPHER")

# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))

# creating Phyloseq object: OTU table, Tax table, Sample ID, Phylogeny, ref seq
# Extracting the standard goods from R

load("C:/sanon/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

ids <- IdTaxa(dna, trainingSet, strand = "top", processors = 10, verbose = TRUE) # use all processors

# ranks of interest
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxid <- x$taxon[m]
  # taxid[startsWith(taxid, "unclassified_")] <- NA
  taxid
}))

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
colnames(taxid)

# Removing sequence rownames for display only
taxid.print <- taxid

rownames(taxid.print) <- NULL

head(taxid.print)

#######
samples.out <- rownames(seqtab.nochim)

samples.out

sites <- sapply(strsplit(samples.out, "-"), `[`, 1)

sites

subject <- sapply(strsplit(samples.out, "_"), `[`, 1)

subject

breeding <- sapply(strsplit(samples.out, "-"), `[`, 2)

breeding

stage <- sapply(strsplit(samples.out, "-"), `[`, 3)

stage

samdf <- data.frame(Stage = stage, Breeding = breeding, Sites = sites)
rownames(samdf) <- samples.out

levels <- c("water", "larvae", "pupae", "adult")
samdf <- samdf %>%
  mutate(
    Stage = case_when(
      Stage == "wa" ~ "water",
      Stage == "lar" ~ "larvae",
      Stage == "pu" ~ "pupae",
      Stage == "ad" ~ "adult",
      TRUE ~ Stage
    ),
    Stage = factor(Stage, levels = levels),
    Breeding = case_when(
      Breeding == "plas" ~ "plastic",
      TRUE ~ Breeding
    ),
    Urbanisation = case_when(
      Sites == "LG" ~ "urban",
      Sites == "TD" ~ "peri-urban"
    ),
    breeding_urban = glue::glue("{Breeding}_{Urbanisation}")
  )

samdf


##########################################################################
##### Built the Phyloseq file
rownames(samdf)
ps <- phyloseq(
  otu_table(t(seqtab.nochim), taxa_are_rows = T),
  sample_data(samdf),
  tax_table(taxid)
)
ps
###################################################################################################################
### REMOVE CONTAMINANT########
theme_set(theme_bw())
library(decontam)

# Put sample_data into a ggplot-friendly data.frame
df <- as.data.frame(sample_data(ps))

df$LibrarySize <- sample_sums(ps)

df <- df[order(df$LibrarySize), ]

df$Index <- seq(nrow(df))

ggplot(data = df, aes(x = Index, y = LibrarySize, color = breeding)) +
  geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Breeding == "NC"

contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg")

table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)

table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples

ps.pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))

ps.pa.neg <- prune_samples(sample_data(ps.pa)$Breeding == "NC", ps.pa)

ps.pa.pos <- prune_samples(sample_data(ps.pa)$Breeding != "NC", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev05$contaminant
)

ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_jitter(alpha = .5, width = .3) +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)")

##### Remove negative controls and contaminants from phyloseq object

ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)

ps.noncontam

ps.noncontam <- prune_samples(sample_data(ps)$Breeding != "NC", ps.noncontam)

SANON <- ps.noncontam

tax_table(SANON)[1, ]

SANON

summary(SANON)

###### Summarize SANON

microbiome::summarize_phyloseq(SANON)

phyloseq::tax_table(SANON)[1:20, 1:6]

# Total number of individuals observed from each ASV

sum.check <- taxa_sums(SANON)

sum(sum.check) #

# Other descriptive statistics:

median(sample_sums(SANON))

SANON <- microbiome::add_refseq(SANON)

# Check if ref_seq slot is added to phyloseq object

print(SANON)

# now check taxa names are ASVids
taxa_names(SANON)[1:3]

phyloseq::tax_table(SANON)[1:6]

abou <- phyloseq::tax_table(SANON)
abou

# Add prevalence filter
## Calculate prevalence of each taxon
#prev <- apply(otu_table(SANON), 1, function(x) sum(x > 0) / length(x))
#
## Define prevalence threshold (5% prevalence)
#prevalence_threshold <- 0.1
#
## Identify taxa with prevalence greater than or equal to 5%
#taxa_to_keep <- names(prev[prev >= prevalence_threshold])
#
## Prune the phyloseq object to only keep those taxa
#physeq_filtered <- prune_taxa(taxa_to_keep, SANON)
#
#SANON <- physeq_filtered

# Save your workspace to a .RData file
save.image(file = "Sanon_16S_DADA2_data.RData")
