library(tidyverse)
library(vegan)
library(glue)
library(RLdbRDA)
here::i_am("beta_diversity.R")
load("Sanon_16S_DADA2_data.RData")

df <- as_tibble(sample_data(SANON))
df$LibrarySize <- sample_sums(SANON)
df<- df[order(df$LibrarySize),]

min_depth <- df %>% 
  select(LibrarySize) %>% 
  min()

vegan_avgdist <- avgdist(as.data.frame(t(otu_table(SANON))), 
                         sample=min_depth, iterations = 1000)

# PCoA
pcoa.ord <- ape::pcoa(vegan_avgdist)

pcoa_df <- data.frame(pcoa.ord$vectors) %>% 
  rownames_to_column("Sample") %>% 
  left_join(samdf %>% rownames_to_column("Sample"))

pcoa_rel_eig <- round(pcoa.ord$values$Relative_eig*100, digits = 2)

ggplot(pcoa_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(color=Stage, shape=breeding_urban), size=2.5)+
  scale_color_brewer(palette = "Oranges")+
  scale_shape_manual(values=c(15, 16, 17, 8))+
  labs(x=glue("PCoA1 ({pcoa_rel_eig[1]}%)"), y=glue("PCoA2 ({pcoa_rel_eig[2]}%)"),
       shape="Breeding material/Urbanisation")+
  theme_bw()

ggsave("figures/PCoA.pdf", dpi=300)

# NMDS
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress),
       xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress",
       main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

NMDS.scree(vegan_avgdist)

# Scree plot shows that from 4 dimensions the stress is significantly reduced
nmds.ord <- metaMDS(vegan_avgdist, k=4)

nmds_df <- data.frame(nmds.ord$points) %>% 
  rownames_to_column("Sample") %>% 
  left_join(samdf %>% rownames_to_column("Sample"))

ggplot(nmds_df, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=Stage, shape=breeding_urban), size=2.5)+
  scale_color_brewer(palette = "Oranges")+
  scale_shape_manual(values=c(15, 16, 17, 8))+
  labs(x="NMDS1", y="NMDS2", shape="Breeding material/Urbanisation")+
  theme_bw()

ggsave("figures/NMDS.pdf", dpi=300)


# PERMANOVA

meta <- data.frame(sample_data(SANON))

adonis2(vegan_avgdist~Sites+Stage+Breeding, data=meta, permutations=999)

adonis2(vegan_avgdist ~ Sites + Sites:Breeding + Sites:Breeding:Stage, data = meta, permutations = 999)

CTRL <- how(plots=Plots(strata = meta$Stage, type="free"),
            within=Within(type="free"),
            nperm=999)
set.seed(4)
CTRL <- how(blocks=meta$Stage,
            nperm=999)

adonis2(vegan_avgdist~Sites+Breeding, data=meta, permutations=CTRL)

# db-RDA
source("custom_rldbrda.R")

meta_overall_rda  <- samdf %>% 
  select(-Urbanisation, -breeding_urban) %>% 
  filter(Breeding != "NC")

rda <- custom_rldbrda(
  vegan_avgdist, 
  meta_overall_rda
  )

plot_data <- prepare_plot_data(rda)
plot_data

g <- plot_dbrda(plot_data)
g+
  labs(y="")+
  scale_y_discrete(labels=c("Stage"="Mosquito\nlife stage", 
                            "Breeding"="Breeding\nmaterial", 
                            "Sites"="Location"))+
  scale_fill_grey(start = .8, end=0.2,
                  labels=c(bquote(R^2), bquote('Cumulative' ~ R^2)))+
  theme_bw()

ggsave("figures/dbRDA.pdf", dpi=300)


#' # Per stage
set.seed(1234)
final_rda <- data.frame(matrix(ncol = 14, nrow = 0))
for (i in c("water", "larvae", "pupae", "adult")){
  names <- samdf %>% 
    filter(Breeding!="NC", 
           Stage == i) %>% 
    rownames_to_column("Sample") %>% 
    pull(Sample)
  
  distm <- as.data.frame(as.matrix(vegan_avgdist)) %>% 
    select(all_of(names)) %>% 
    filter(rownames(.) %in% names) %>% 
    as.dist(diag = T, upper = T)
  
  meta_rda <- meta_overall_rda %>% 
    select(-Stage) %>% 
    filter(rownames(.) %in% names)
  
  rda <- custom_rldbrda(
    distm,
    meta_rda
  )
  
  if (is.null(rda)){
    rda <- data.frame(matrix(ncol = 14, nrow = 2))
    colnames(rda) <- colnames(final_rda)
    rownames(rda) <- c("Breeding", "Sites")
    rda["Stage"] <- i
  } else {
    rda <- rda %>% 
      mutate(Stage = i)
  }
    
  final_rda <- rbind(final_rda, rda)
}

plot_data <- prepare_plot_data(final_rda)
plot_data

plot_data %>% 
  filter(variable != "RDAcumul_R2.adj") %>% 
  plot_dbrda()+
  facet_wrap(~factor(Stage, levels = c("water", "larvae", "pupae", "adult")))+
  scale_fill_grey(start = 0.4,
                  labels=c(bquote(R^2), bquote('Cumulative' ~ R^2)))+
  scale_y_discrete(labels=c("Breeding"="Breeding\nmaterial", 
                            "Sites"="Location"))+
  labs(x=bquote("Effect size (adjusted "*R^2*")"))+
  guides(fill="none")+  
  theme_bw()
ggsave("figures/dbRDA_per_stage.pdf", dpi=300)

# permM<- adonis2(t(otu)~Sites*Stage,data= meta, permutations=999, 
#                         method="bray", by= "terms",na.rm=T)
# permM
# permM1<- adonis2(t(otu)~Sites+Stage,data= meta, permutations=999, 
#                 method="bray", by= "terms",na.rm=T)
# permM1
# permM2<- adonis2(t(otu)~Breeding*Stage,data= meta, permutations=999, 
#                 method="bray", by= "terms",na.rm=T)
# permM2
# permM3<- adonis2(t(otu)~Breeding+Stage,data= meta, permutations=999, 
#                  method="bray", by= "terms",na.rm=T)
# permM3
# permM4<- adonis2(t(otu)~Sites*Stage*Breeding,data= meta, permutations=999, 
#                  method="bray", by= "terms",na.rm=T)
# permM4
# 
# #### Plots for Beta diversity 
# # Transform data to proportions as appropriate for Bray-Curtis distances
# 
# #ps.prop <- transform_sample_counts(SANON, function(otu) otu/sum(otu))
# 
# ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
# 
# ord.pca.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Breeding",title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Stage",title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Sites",title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Stage",title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Sites", title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)
# 
# plot_ordination(ps.prop, ord.pca.bray, color="Stage", shape="Sites", title="Bray NMDS")+geom_point(size=3)