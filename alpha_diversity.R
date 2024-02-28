here::i_am("alpha_diversity.R")

# Load your workspace from a .RData file
load("Sanon_16S_DADA2_data.RData")
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