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

# plot_richness(SANON, x="Stage", measures=c("Simpson", "Shannon"), color = "Stage") +
#  geom_boxplot() +
#  theme_classic() +
#  ggpubr::stat_pwc(method = "wilcox.test", p.adjust.method = "BH", label="p.adj.format", tip.length = 0)+
#  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90 ))


alpha_rarefied <- function(ab_table, sequencing_depth) {
  df <- ab_table %>%
    t() %>%
    vegan::rrarefy(., sample = sequencing_depth) %>% # rrafefy samples from rows, not from columns!
    as_tibble(rownames = "sample") %>%
    group_by(sample) %>%
    pivot_longer(-sample) %>%
    summarize(
      Observed = vegan::specnumber(value),
      Shannon = vegan::diversity(value, index = "shannon"),
      Simpson = vegan::diversity(value, index = "simpson")
    ) %>%
    as.data.frame() %>%
    column_to_rownames("sample")
  df
}

#' ## Alpha diversity
df <- as_tibble(sample_data(SANON))
df$LibrarySize <- sample_sums(SANON)
df <- df[order(df$LibrarySize), ]

# threshold <- quantile(df$LibrarySize, 0.05)

min_depth <- df %>%
  select(LibrarySize) %>%
  min()

ggplot(data = df, aes(x = LibrarySize, y = "samples")) +
  geom_jitter(height = 0.02) +
  ggtitle("Library sizes") +
  geom_vline(xintercept = min_depth, linetype = "dashed", color = "grey") +
  labs(y = "") +
  # scale_x_log10()+
  theme_bw()

set.seed(1234)
alpha_df_list <- purrr::map(1:1000, ~ alpha_rarefied(ab_table = as_tibble(otu_table(SANON)), sequencing_depth = min_depth))
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(samdf %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(
      -Sample, -Stage, -Breeding, -Sites,
      -Urbanisation, -breeding_urban
    ),
    names_to = "Metric", values_to = "Diversity"
  ) %>%
  ggplot(aes(x = Stage, y = Diversity, color = Stage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.text.align = 0
  ) +
  # scale_color_brewer(palette="Oranges")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF")) +
  ggpubr::stat_pwc(
    method = "wilcox.test", p.adjust.method = "BH",
    label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )
alpha

ggsave("figures/alpha_diversity.pdf", dpi = 300)

# Plot shape based on location/breeding material
alpha_shape <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(samdf %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(
      -Sample, -Stage, -Breeding, -Sites,
      -Urbanisation, -breeding_urban
    ),
    names_to = "Metric", values_to = "Diversity"
  ) %>%
  ggplot(aes(x = Stage, y = Diversity, color = Stage)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=gsub(pattern = "_", replacement = "/", x = breeding_urban)), width = 0.3) +
  facet_wrap(~Metric, nrow = 1, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.text.align = 0
  ) +
  # scale_color_brewer(palette="Oranges")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF")) +
  scale_shape_manual(values = c(0, 15, 1, 16)) +
  labs(shape = "Breeding material/\nLocation")+
  ggpubr::stat_pwc(
    method = "wilcox.test", p.adjust.method = "BH",
    label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )
alpha_shape
ggsave("figures/alpha_diversity_shape.pdf", dpi = 300, height=7, width=9)

# Alpha diversity per breeding/urbanisation

alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(samdf %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(
      -Sample, -Stage, -Breeding, -Sites,
      -Urbanisation, -breeding_urban
    ),
    names_to = "Metric", values_to = "Diversity"
  ) %>%
  filter(Metric == "Simpson") %>%
  ggplot(aes(x = Urbanisation, y = Diversity, color = Urbanisation)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=gsub(pattern = "_", replacement = "/", x = breeding_urban)), width = 0.3) +
  facet_nested_wrap(~Stage, nrow=2, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.text.align = 0
  ) +
  scale_shape_manual(values = c(0, 15, 1, 16)) +
  labs(shape = "Breeding material/\nLocation")+
  ggpubr::stat_pwc(
    method = "wilcox.test", p.adjust.method = "BH",
    label = "p.adj.format", hide.ns = "p.adj", show.legend = F, tip.length = 0.01
  )
