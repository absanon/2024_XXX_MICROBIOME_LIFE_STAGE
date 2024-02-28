# Author: Sanon Aboubakar
#######Install the required packages for your analysis
BiocManager:: install("dada2")
BiocManager::install("DECIPHER")
BiocManager::install("microbiome")

install.packages("BiocManager")

install.packages("devtools")

library("devtools")

library(dada2);

packageVersion("dada2")

library(here)

i_am("dada2.R")

###############Load your fastq data
path<-"C:/sanon"
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq.gz", full.names = TRUE))
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

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,truncLen=c(290,210), 
                     maxN=0,maxEE=c(2,2),truncQ=2,rm.phix =TRUE,
                    compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool="pseudo")
dadaFs[[1]]
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

##### Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

#### If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track
taxa <- assignTaxonomy(seqtab.nochim, "C:/sanon/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/sanon/silva_species_assignment_v132.fa.gz")

# Removing sequence rownames for display only
taxa.print <- taxa 

rownames(taxa.print) <- NULL

head(taxa.print)

taxa.print

tail(taxa.print)

View(taxa.print)

library(DECIPHER); packageVersion("DECIPHER")

# Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim)) 

# creating Phyloseq object: OTU table, Tax table, Sample ID, Phylogeny, ref seq
#Extracting the standard goods from R

load("C:/sanon/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

ids <- IdTaxa(dna, trainingSet, strand="top", processors=10, verbose=TRUE) # use all processors

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # ranks of interest
ranks
class(ranks)

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(taxid) <- ranks
rownames(taxid) <- getSequences(seqtab.nochim)
colnames(taxid)

# Removing sequence rownames for display only
taxa.print <- taxid

rownames(taxa.print) <- NULL

head(taxa.print)

####### 
samples.out <- rownames(seqtab.nochim)

samples.out

sites <-sapply(strsplit(samples.out, "-"), `[`, 1)

sites

subject <- sapply(strsplit(samples.out, "_"), `[`, 1)

subject

breeding <- sapply(strsplit(samples.out, "-"), `[`, 2)

breeding

stage <- sapply(strsplit(samples.out, "-"), `[`, 3)

stage

samdf <- data.frame(Stage=stage, Breeding=breeding, Sites=sites)

samdf

rownames(samdf) <- samples.out
########################################################################## 
##### Built the Phylosequ file
rownames(samdf)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps
###################################################################################################################
###REMOVE CONTAMINANT########
theme_set(theme_bw())
library(decontam)

# Put sample_data into a ggplot-friendly data.frame
df <- as.data.frame(sample_data(ps)) 

df$LibrarySize <- sample_sums(ps)

df <- df[order(df$LibrarySize),]

df$Index <- seq(nrow(df))

ggplot(data=df, aes(x=Index, y=LibrarySize, color=breeding)) + geom_point()

sample_data(ps)$is.neg <- sample_data(ps)$Breeding == "NC"

contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")

table(contamdf.prev$contaminant)

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)

table(contamdf.prev05$contaminant)

# Make phyloseq object of presence-absence in negative controls and true samples

ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))

ps.pa.neg <- prune_samples(sample_data(ps.pa)$Breeding == "NC", ps.pa)

ps.pa.pos <- prune_samples(sample_data(ps.pa)$Breeding != "NC", ps.pa)

# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

##### Remove negative controls and contaminants from phyloseq object

ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps)

ps.noncontam

ps.noncontam <- prune_samples(sample_data(ps)$Breeding!='NC', ps.noncontam)

Sanon <- ps.noncontam

Sanon

Sanon_DECIPHER <- Sanon

tax_table(Sanon_DECIPHER) <- taxid

Sanon_DECIPHER

# Save your workspace to a .RData file
save.image(file = "C:/sanon/Sanon_16S_DADA2_data.RData")
