library(tidyverse)
library(vegan)
library(glue)
library(RLdbRDA)
library(patchwork)
library(phyloseq)
here::i_am("beta_diversity.R")
load("Sanon_16S_DADA2_data.RData")

df <- as_tibble(sample_data(SANON))
df$LibrarySize <- sample_sums(SANON)
df <- df[order(df$LibrarySize),]

min_depth <- df %>% 
  select(LibrarySize) %>% 
  min()

set.seed(1234)
bray_avgdist <- vegan::avgdist(as.data.frame(t(otu_table(SANON))), 
                         sample=min_depth, iterations = 1000)

# PCoA
pcoa.ord <- ape::pcoa(bray_avgdist)

## Plot loadings on PCoA plot
#vf <- vegan::envfit(pcoa.ord$vectors, meta[1:3])
#
## Select factors you want to plot
#spp.scrs <- as.data.frame(scores(vf, display = "factors"))
#spp.scrs$Group <- rownames(spp.scrs)
#
#spp.scrs <- spp.scrs %>% 
#  mutate(Group=gsub("Stage|Breeding|Sites", "", Group))

pcoa_df <- data.frame(pcoa.ord$vectors) %>% 
  rownames_to_column("Sample") %>% 
  left_join(samdf %>% rownames_to_column("Sample"))

pcoa_rel_eig <- round(pcoa.ord$values$Relative_eig*100, digits = 2)

p_pcoa <- ggplot(pcoa_df, aes(x=Axis.1, y=Axis.2))+
  geom_point(aes(color=Stage, shape=breeding_urban), size=2.5)+
  #scale_color_brewer(palette = "Oranges")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"))+
  scale_shape_manual(values=c(0, 15, 1, 16))+
  #geom_segment(data = spp.scrs,
  #             aes(x = 0, xend = Axis.1, y = 0, yend = Axis.2), linetype=2,
  #             arrow = arrow(length = unit(0.25, "cm")), colour = "grey") + 
  #geom_text(data = spp.scrs, aes(x = Axis.1, y = Axis.2, label = Group), colour="black",
  #          size = 3, nudge_y = .01)+
  labs(x=glue("PCoA1 ({pcoa_rel_eig[1]}%)"), y=glue("PCoA2 ({pcoa_rel_eig[2]}%)"),
       shape="Breeding material/\nUrbanisation")+
  theme_bw()
p_pcoa

ggsave("figures/PCoA.pdf", dpi=300)

# NMDS
source("custom_rldbrda.R")
NMDS_scree(bray_avgdist)

# Scree plot shows that from 4 dimensions the stress is significantly reduced
nmds.ord <- vegan::metaMDS(bray_avgdist, k=4)

nmds_df <- data.frame(nmds.ord$points) %>% 
  rownames_to_column("Sample") %>% 
  left_join(samdf %>% rownames_to_column("Sample"))

p_nmds <- ggplot(nmds_df, aes(x=MDS1, y=MDS2))+
  geom_point(aes(color=Stage, shape=breeding_urban), size=2.5)+
  #scale_color_brewer(palette = "Oranges")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"))+
  scale_shape_manual(values=c(0, 15, 1, 16))+
  labs(x="NMDS1", y="NMDS2", shape="Breeding material/\nUrbanisation")+
  theme_bw()
p_nmds

ggsave("figures/NMDS.pdf", dpi=300)


## PERMANOVA
#
#meta <- data.frame(sample_data(SANON))
#
#adonis2(bray_avgdist~Sites+Stage+Breeding, data=meta, permutations=999)
#
#adonis2(bray_avgdist ~ Sites + Sites:Breeding + Sites:Breeding:Stage, data = meta, permutations = 999)
#
#CTRL <- how(plots=Plots(strata = meta$Stage, type="free"),
#            within=Within(type="free"),
#            nperm=999)
#set.seed(4)
#CTRL <- how(blocks=meta$Stage,
#            nperm=999)
#
#adonis2(bray_avgdist~Sites+Breeding, data=meta, permutations=CTRL)

# db-RDA on all explanatory variables
set.seed(1234)
meta_overall_rda  <- samdf %>% 
  select(-Urbanisation, -breeding_urban) %>% 
  filter(Breeding != "NC")

rda <- custom_rldbrda(
  bray_avgdist, 
  meta_overall_rda
  )

plot_data <- RLdbRDA::prepare_plot_data(rda)
plot_data

g <- RLdbRDA::plot_dbrda(plot_data)
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

rda <- custom_rldbrda(
  bray_avgdist, 
  meta_overall_rda
)

final_rda <- data.frame(matrix(ncol = 14, nrow = 0))
for (i in c("water", "larvae", "pupae", "adult")){
  names <- samdf %>% 
    filter(Breeding!="NC", 
           Stage == i) %>% 
    rownames_to_column("Sample") %>% 
    pull(Sample)
  
  distm <- as.data.frame(as.matrix(bray_avgdist)) %>% 
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

p_dbrda <- plot_data %>% 
  #filter(variable != "RDAcumul_R2.adj") %>% 
  plot_dbrda()+
  facet_wrap(~factor(Stage, levels = c("water", "larvae", "pupae", "adult")), nrow=1)+
  #scale_fill_grey(start = 0.1, end=.85, 
  #                labels=c(bquote(R^2), bquote('Cumulative' ~ R^2)))+
  scale_y_discrete(labels=c("Breeding"="Breeding\nmaterial", 
                            "Sites"="Location"))+
  labs(x=bquote("Effect size (adjusted "*R^2*")"))+
  #guides(fill="none")+  
  theme_bw()+
  theme(axis.title.y = element_blank())
p_dbrda
ggsave("figures/dbRDA_per_stage.pdf", dpi=300)

# Combine pcoa + dbRDA plots
p_pcoa / p_dbrda +
  plot_layout(heights = c(4, 1))+ 
  plot_annotation(tag_levels = 'A')&
  theme(plot.tag = element_text(face="bold"))
ggsave("figures/beta_diversity.pdf", dpi=300)

p1 <- alpha +
  #labs(x="")+
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) 

p2 <- (p_pcoa / p_dbrda) +
  plot_layout(heights = c(4, 1))


free(p1) + p2 +
  plot_layout(widths = c(1, 1))+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
ggsave("figures/combined_diversity.pdf", dpi=300, width = 10, height = 6)
