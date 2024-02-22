#######Install the required packages for your analysis
BiocManager:: install.packages("dada2", version = "3.11")

install.packages("BiocManager")

install.packages("devtools")

library("devtools")

library(dada2);

packageVersion("dada2")

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

if (!requireNamespace("BiocManager", quietly=TRUE))

install.packages("BiocManager")

BiocManager::install("DECIPHER")

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

colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
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

# Load your workspace from a .RData file
load("C:/sanon/Sanon_16S_DADA2_data.RData")
##################################################################################################################
##### clean up your sequence
library(microbiome)

ps.fam <- subset_taxa(Sanon, !is.na(Family) & !Family %in% c("", "uncharacterized"))

SANON<-ps.fam

SANON

tax_table(SANON)[1,]

SANON

summary(SANON)

###### Summarize SANON

microbiome::summarize_phyloseq(SANON)

phyloseq::tax_table(Sanon)[1:20,1:6]

# Total number of individuals observed from each ASV

sum.check <- taxa_sums(SANON)

sum(sum.check)#

# Other descriptive statistics:

median(sample_sums(SANON))

SANON <- microbiome::add_refseq(SANON)

# Check if ref_seq slot is added to phyloseq object

print(SANON)

# now check taxa names are ASVids
taxa_names(SANON)[1:3]

# Check 
taxa_names(SANON)[1:10]

phyloseq::tax_table(SANON)[1:6]

abou<-phyloseq::tax_table(SANON)
abou

#### Sort samples as you want they appear in your graph
levels <- c("wa", "lar", "pu", "ad" ) 

## Alpha diversity measures per developmental stages

plot_richness(SANON, x="Stage", measures=c("Simpson", "Shannon"), color = "Stage") +
  geom_boxplot() +
  theme_classic() +
 # levels(c("wa", "lar", "pu", "ad"))+
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
  #                   method = "wilcox.test", hide.ns = F, tip.length = 0)+
  stat_pwc(method = "wilcox.test", p.adjust.method = "BH", label="p.adj.format", tip.length = 0)+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90 ))

# Pick relative abundances (compositional) and sample metadata 

SANON.rel <- microbiome::transform(SANON, "compositional")

otu <- microbiome::abundances(SANON.rel)

meta <- meta(SANON.rel)


#Permanova for Beta diversity

#library(vegan)
permM<- adonis2(t(otu)~Sites*Stage,data= meta, permutations=999, 
                        method="bray", by= "terms",na.rm=T)
permM
permM1<- adonis2(t(otu)~Sites+Stage,data= meta, permutations=999, 
                method="bray", by= "terms",na.rm=T)
permM1
permM2<- adonis2(t(otu)~Breeding*Stage,data= meta, permutations=999, 
                method="bray", by= "terms",na.rm=T)
permM2
permM3<- adonis2(t(otu)~Breeding+Stage,data= meta, permutations=999, 
                 method="bray", by= "terms",na.rm=T)
permM3
permM4<- adonis2(t(otu)~Sites*Stage*Breeding,data= meta, permutations=999, 
                 method="bray", by= "terms",na.rm=T)
permM4

#### Plots for Beta diversity 
# Transform data to proportions as appropriate for Bray-Curtis distances

ps.prop <- transform_sample_counts(SANON, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ord.pca.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

plot_ordination(ps.prop, ord.pca.bray, color="Breeding",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Stage",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Sites",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Stage",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Stage", shape="Sites", title="Bray NMDS")+geom_point(size=3)

##### Community composition per sites
#### Phylum
ps.rel <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

ps.rel
# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Phylum', NArm = FALSE)

ps.melt <- psmelt(glom)

# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(Stage, Phylum) %>%
  mutate(median=median(Abundance))

# select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])

ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "<1%"

#to get the same rows together
ps.melt_sum <- ps.melt %>%
  group_by(Sample, Stage,Phylum) %>%
  summarise(Abundance=sum(Abundance))

ps.melt_sum <- ps.melt %>%
  group_by(Sample, Stage,Phylum) %>%
  summarise(Abundance=sum(Abundance))

# Setting factor levels
ps.melt_sum$Stage <- factor(ps.melt_sum$Stage, levels = c("wa", "lar", "pu", "ad"))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum)) + 
  labs(x="Samples", y="Relative abundance (%)") +
  facet_wrap(~Stage, scales= "free_x", nrow=1) +
  theme_classic() +
  scale_fill_brewer("Phylum", palette = "Paired")+
      theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))
#####################################################################################################################"
####### Class

ps.rel1 = transform_sample_counts(SANON, function(x) x/sum(x)*100)

# agglomerate taxa
glom <- tax_glom(ps.rel1, taxrank = 'Class', NArm = FALSE)

ps.melt1 <- psmelt(glom)

# change to character for easy-adjusted level
ps.melt1$Class <- as.character(ps.melt1$Class)

ps.melt1 <- ps.melt1 %>%
  group_by(Stage, Class) %>%
  mutate(median=median(Abundance))

# select group mean > 1
keep <- unique(ps.melt1$Class[ps.melt1$median > 1])

ps.melt1$Class[!(ps.melt1$Class %in% keep)] <- "< 1%"

#to get the same rows together
ps.melt_sum <- ps.melt1 %>%
  group_by(Sample,Stage,Class) %>%
  summarise(Abundance=sum(Abundance))

# Setting factor levels
ps.melt_sum$Stage <- factor(ps.melt_sum$Stage, levels = c("wa", "lar", "pu", "ad"))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill=Class)) + 
  labs(x="Samples", y="Relative abundance (%)") +
  facet_wrap(~Stage, scales= "free_x", nrow=1) +
  theme_classic() + 
  scale_fill_brewer("Class", palette = "Paired")+
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90,size = 10, face = "bold"))
##########################################################################################################################

####Community composition per stages and sites
#### Phylum
 ps.relA <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

 ps.relA

 ## agglomerate taxa
glom <- tax_glom(ps.relA, taxrank = 'Phylum', NArm = FALSE)

 ps.melt <- psmelt(glom)

  # change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

 ps.melt <- ps.melt %>%

      group_by(Stage, Breeding, Phylum) %>%
   
   mutate(median=median(Abundance))

 # select group median > 1
keep <- unique(ps.melt$Phylum[ps.melt$median >1 ])

 ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "<1%"

  #to get the same rows together
 
ps.melt_sum <- ps.melt %>%

    group_by(Sample,Stage,Breeding,Phylum) %>%
  
   summarise(Abundance=sum(Abundance))

ps.melt_sum

 ps.melt_sum$Stage <- factor(ps.melt_sum$Stage, levels = c("wa", "lar", "pu", "ad"))
 
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
   geom_bar(stat = "identity", aes(fill=Phylum)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Stage+Breeding, scales= "free_x", nrow=1) +
   theme_classic() + 
   scale_fill_brewer("Phylum", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))

##################################################################################################################  #### Class
 ps.relB <- transform_sample_counts(SANON, function(x) x/sum(x)*100)
 
 ps.relB
 
 ## agglomerate taxa
 glom <- tax_glom(ps.relB, taxrank = 'Class', NArm = FALSE)
 
 ps.melt <- psmelt(glom)
 
 # change to character for easy-adjusted level
 ps.melt$Class <- as.character(ps.melt$Class)
 
 ps.melt <- ps.melt %>%
   group_by(Stage, Breeding, Class) %>%
   mutate(median=median(Abundance))
 
 # select group median > 1
 keep <- unique(ps.melt$Class[ps.melt$median > 1])
 ps.melt$Class[!(ps.melt$Class %in% keep)] <- "<1%"
 
 #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Stage,Breeding,Class) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum

 #setting level 
 ps.melt_sum$Stage <- factor(ps.melt_sum$Stage, levels = c("wa", "lar", "pu", "ad"))
 
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
   geom_bar(stat = "identity", aes(fill=Class)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Stage+Breeding, scales= "free_x", nrow=1) +
  # scale_fill_manual(values = mycolors)+
   theme_classic() + 
   scale_fill_brewer("Class", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold")) 

 ###########################################################################################################################
####Community composition per sites

 #### Phylum
 
 ps.relC <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

 ps.relC
 
 ## agglomerate taxa
 glom <- tax_glom(ps.relC, taxrank = 'Phylum', NArm = FALSE)
 ps.melt <- psmelt(glom)

# change to character for easy-adjusted level
ps.melt$Phylum <- as.character(ps.melt$Phylum)

 ps.melt <- ps.melt %>%

      group_by(Sites, Phylum) %>%
   mutate(median=median(Abundance))

# select group median > 1
 keep <- unique(ps.melt$Phylum[ps.melt$median > 1])

  ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "<1%"
 
  #to get the same rows together
ps.melt_sum <- ps.melt %>%

     group_by(Sample,Sites,Phylum) %>%

  summarise(Abundance=sum(Abundance))

ps.melt_sum

 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
   geom_bar(stat = "identity", aes(fill=Phylum)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Sites, scales= "free_x", nrow=1) +
   theme_classic() + 
   scale_fill_brewer("Phylum", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))
#######################################################################
 ###### Composition per breeding material 
 ps.relK <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

  ps.relK
 
## agglomerate taxa
 glom <- tax_glom(ps.relK, taxrank = 'Phylum', NArm = FALSE)

  ps.melt <- psmelt(glom)
 
  # change to character for easy-adjusted level
 ps.melt$ Phylum <- as.character(ps.melt$Phylum)
ps.melt <- ps.melt %>%
   group_by(Breeding, Phylum) %>%
   mutate(median=median(Abundance))

 # select group median > 1
 keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
 ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "<1%"

  #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Breeding,Phylum) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
   geom_bar(stat = "identity", aes(fill=Phylum)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Breeding, scales= "free_x", nrow=1) +
   theme_classic() + 
   scale_fill_brewer("Phylum", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))
 ######################################################################################################"
 ps.relP <- transform_sample_counts(SANON, function(x) x/sum(x)*100)
 
 ps.relP
 
 ## agglomerate taxa
 glom <- tax_glom(ps.relK, taxrank = 'Family', NArm = FALSE)
 ps.melt <- psmelt(glom)
 
 # change to character for easy-adjusted level
 ps.melt$Family<- as.character(ps.melt$Family)
 ps.melt <- ps.melt %>%
   group_by(Breeding, Family) %>%
   mutate(median=median(Abundance))
 
 # select group median > 1
 keep <- unique(ps.melt$Family[ps.melt$median > 0.5])
 ps.melt$Family[!(ps.melt$Family %in% keep)] <- "<0.5%"

 #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Breeding,Family) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
   geom_bar(stat = "identity", aes(fill=Family)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Breeding, scales= "free_x", nrow=1) +
   theme_classic() + 
   scale_fill_brewer("Family", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))
 #######################################################################################
 
 ####Community composition per site and breeding material
 ps.relC1 <- transform_sample_counts(SANON, function(x) x/sum(x)*100)
 
 ps.relC1

## agglomerate taxa
glom <- tax_glom(ps.relC1, taxrank = 'Phylum', NArm = FALSE)
 ps.melt <- psmelt(glom)

# change to character for easy-adjusted level
 ps.melt$Phylum <- as.character(ps.melt$Phylum)
 ps.melt <- ps.melt %>%
   group_by( Stage,Sites, Phylum) %>%
   mutate(median=median(Abundance))

  # select group median > 1
 keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
 ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "<1%"
 
 #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Stage, Sites, Phylum) %>%
   summarise(Abundance=sum(Abundance))

ps.melt_sum

 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
   geom_bar(stat = "identity", aes(fill=Phylum)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Stage+Sites, scales= "free_x", nrow=1) +
   theme_classic() + 
   scale_fill_brewer("Phylum", palette = "Paired")+
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))
 
 ###############################################################################################################
 #### Class

  ps.relD <- transform_sample_counts(Sanon, function(x) x/sum(x)*100)
 
 ps.relD

 ## agglomerate taxa
 glom <- tax_glom(ps.relB, taxrank = 'Class', NArm = FALSE)
 ps.melt <- psmelt(glom)
 
# change to character for easy-adjusted level
 ps.melt$Class <- as.character(ps.melt$Class)
 ps.melt <- ps.melt %>%
   group_by(Sites, Genus) %>%
   mutate(median=median(Abundance))

  # select group median > 1
 keep <- unique(ps.melt$Class[ps.melt$median > 1])
 ps.melt$Class[!(ps.melt$Class %in% keep)] <- "<1%"
 
#to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Sites,Genus) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
   geom_bar(stat = "identity", aes(fill=Genus)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Sites, scales= "free_x", nrow=1) +
   theme_classic() + 
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold")) 
######################################################################################################################### 
 #### Class

 ps.relF <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

  ps.relF
 
  ## agglomerate taxa
 glom <- tax_glom(ps.relF, taxrank = 'Class', NArm = FALSE)
ps.melt <- psmelt(glom)

# change to character for easy-adjusted level
 ps.melt$Class <- as.character(ps.melt$Class)
 ps.melt <- ps.melt %>%
   group_by(Breeding, Class) %>%
   mutate(median=median(Abundance))

  # select group median > 1
 keep <- unique(ps.melt$Class[ps.melt$median > 1])
 ps.melt$Class[!(ps.melt$Class %in% keep)] <- "<1%"

 #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Breeding,Class) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
   geom_bar(stat = "identity", aes(fill=Class)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Breeding, scales= "free_x", nrow=1) +
   theme_classic() + 
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))  
########################################################################################################################

  #### Class
 ps.relM <- transform_sample_counts(SANON, function(x) x/sum(x)*100)
 
 ps.relM
 
 ## agglomerate taxa
 glomM <- tax_glom(ps.relM, taxrank = 'Class', NArm = FALSE)
 ps.melt <- psmelt(glomM)
 
 # change to character for easy-adjusted level
 ps.melt$Class <- as.character(ps.melt$Class)
 ps.melt <- ps.melt %>%
   group_by(Breeding+Sites, Class) %>%
   mutate(median=median(Abundance))
 
 # select group median > 1
 keep <- unique(ps.melt$Class[ps.melt$median > 1])
 ps.melt$Class[!(ps.melt$Class %in% keep)] <- "<1%"
 
 #to get the same rows together
 ps.melt_sum <- ps.melt %>%
   group_by(Sample,Breeding+Sites,Class) %>%
   summarise(Abundance=sum(Abundance))
 ps.melt_sum
 ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
   geom_bar(stat = "identity", aes(fill=Class)) + 
   labs(x="Samples", y="Relative abundance (%)") +
   facet_wrap(~Breeding+sites, scales= "free_x", nrow=1) +
   theme_classic() + 
   theme(strip.background = element_blank(), 
         axis.text.x.bottom = element_text(angle = -90, size = 10, face = "bold"))  
 #########################################################################################################################