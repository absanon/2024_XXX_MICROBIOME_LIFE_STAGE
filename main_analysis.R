library(here)
i_am("main_analysis.R")

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