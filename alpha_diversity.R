here::i_am("alpha_diversity.R")
library(tidyverse)
library(microbiome)
library(vegan)
library(ggpubr)

# Load your workspace from a .RData file
load("Sanon_16S_DADA2_data.RData")
##################################################################################################################
##### clean up your sequence

#### Sort samples as you want they appear in your graph
## Alpha diversity measures per developmental stages

#plot_richness(SANON, x="Stage", measures=c("Simpson", "Shannon"), color = "Stage") +
#  geom_boxplot() +
#  theme_classic() +
#  ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", label="p.adj.format", tip.length = 0)+
#  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90 ))


alpha_rarefied <- function(ab_table, sequencing_depth) {
  df <- ab_table %>%
    t() %>% 
    rrarefy(., sample=sequencing_depth) %>% # rrafefy samples from rows, not from columns!
    as_tibble(rownames="sample") %>%
    group_by(sample) %>%
    pivot_longer(-sample) %>%
    summarize(Observed = specnumber(value),
              Shannon = diversity(value, index="shannon"),
              Simpson = diversity(value, index="simpson")
    ) %>%
    as.data.frame() %>%
    column_to_rownames("sample")
  df
}

#' ## Alpha diversity
df <- as_tibble(sample_data(SANON))
df$LibrarySize <- sample_sums(SANON)
df<- df[order(df$LibrarySize),]

#threshold <- quantile(df$LibrarySize, 0.05)

min_depth <- df %>% 
  select(LibrarySize) %>% 
  min()

ggplot(data=df, aes(x=LibrarySize, y="samples")) + 
  geom_jitter(height = 0.02)+
  ggtitle("Library sizes")+
  geom_vline(xintercept = min_depth, linetype="dashed", color="grey")+
  labs(y="")+
  #scale_x_log10()+
  theme_bw()

set.seed(1234)
alpha_df_list <- purrr::map(1:1000, ~alpha_rarefied(ab_table=as_tibble(otu_table(SANON)), sequencing_depth=min_depth))
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>% 
  mutate(Observed=round(Observed))

alpha <- alpha_average_df %>% 
  rownames_to_column("Sample") %>% 
  select(-Observed) %>% 
  left_join(samdf %>% rownames_to_column("Sample")) %>% 
  pivot_longer(c(-Sample, -Stage, -Breeding, -Sites, -Urbanisation), names_to = "Metric", values_to = "Diversity") %>%
  ggplot(aes(x=Stage, y=Diversity, color=Stage))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)+
  facet_wrap(~Metric, nrow=1, scales="free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.text.align = 0)+
  scale_color_brewer(palette = "Set1")+
  stat_pwc(method = "wilcox.test", p.adjust.method="BH",
           label = "p.adj.format", hide.ns="p.adj", show.legend = F, tip.length = 0.01)
alpha

ggsave("figures/alpha_diversity.pdf", dpi=300)
