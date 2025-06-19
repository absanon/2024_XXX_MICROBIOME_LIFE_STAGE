library(tidyverse)
library(ggalluvial)
library(ggtext)
library(ggh4x)
library(phyloseq)
library(patchwork)

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

cols <- c(
  "water" = "#3B9AB2",
  "eggs" = "#78B7C5",
  "larvae" = "#EBCC2A",
  "pupae" = "#F21A00",
  "adult" = "#7A0403FF"
)

# Bennett
load("bennett_dada2.RData")

# alpha diversity
min_depth <- min(sample_sums(ps)[sample_sums(ps) > 1000])
set.seed(1234)
alpha_df_list <- purrr::map(
  1:1000,
  ~ alpha_rarefied(
    ab_table = as_tibble(otu_table(ps)),
    sequencing_depth = min_depth
  )
)
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha_bennett <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(bennett_meta %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Diversity"
  ) %>%
  mutate(
    life_stage = case_when(
      life_stage == "Adult" ~ "adult",
      life_stage == "Larvae" ~ "larvae",
      T ~ life_stage
    ),
    life_stage = forcats::fct_relevel(life_stage, "larvae", "adult")
  ) %>%
  ggplot(aes(x = life_stage, y = Diversity, color = life_stage)) +
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
  scale_color_manual(values = cols) +
  ggpubr::stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )
alpha_bennett
ggsave("bennett_alpha.pdf", dpi = 300, width = 7, height = 5)


# Hernandez
load("hernandez_dada2.RData")

# alpha diversity
min_depth <- min(sample_sums(ps)[sample_sums(ps) > 1000])
set.seed(1234)
alpha_df_list <- purrr::map(
  1:1000,
  ~ alpha_rarefied(
    ab_table = as_tibble(otu_table(ps)),
    sequencing_depth = min_depth
  )
)
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha_hernandez <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(hernandez_meta %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Diversity"
  ) %>%
  mutate(
    life_stage = case_when(
      life_stage == "larves" ~ "larvae",
      T ~ life_stage
    ),
    life_stage = forcats::fct_relevel(life_stage, "eggs", "larvae", "adult")
  ) %>%
  ggplot(aes(x = life_stage, y = Diversity, color = life_stage)) +
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
  scale_color_manual(values = cols) +
  ggpubr::stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )
alpha_hernandez
ggsave("hernandez_alpha.pdf", dpi = 300, width = 7, height = 5)

# Rodpai
load("rodpai_dada2.RData")

# alpha diversity
min_depth <- min(sample_sums(ps)[sample_sums(ps) > 1000])
set.seed(1234)
alpha_df_list <- purrr::map(
  1:1000,
  ~ alpha_rarefied(
    ab_table = as_tibble(otu_table(ps)),
    sequencing_depth = min_depth
  )
)
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha_rodpai <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(rodpai_meta %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Diversity"
  ) %>%
  mutate(
    life_stage = case_when(
      isolation_source == "Mosquito abdomen" ~ "adult",
      isolation_source == "Mosquito larvae" ~ "larvae",
      isolation_source == "Mosquito larvae habitat water" ~ "water",
      T ~ isolation_source
    )
  ) %>%
  mutate(
    life_stage = forcats::fct_relevel(life_stage, "water", "larvae", "adult")
  ) %>%
  ggplot(aes(x = life_stage, y = Diversity, color = life_stage)) +
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
  scale_color_manual(values = cols) +
  ggpubr::stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )
alpha_rodpai
ggsave("rodpai_alpha.pdf", dpi = 300, width = 7, height = 10)

# Hery
load("hery_dada2.RData")

# alpha diversity
min_depth <- min(sample_sums(ps)[sample_sums(ps) > 1000])
set.seed(1234)
alpha_df_list <- purrr::map(
  1:1000,
  ~ alpha_rarefied(
    ab_table = as_tibble(otu_table(ps)),
    sequencing_depth = min_depth
  )
)
alpha_average_df <- Reduce(`+`, alpha_df_list) / length(alpha_df_list)
alpha_average_df <- alpha_average_df %>%
  mutate(Observed = round(Observed))

alpha_hery <- alpha_average_df %>%
  rownames_to_column("Sample") %>%
  select(-Observed) %>%
  left_join(hery_meta %>% rownames_to_column("Sample")) %>%
  pivot_longer(
    c(Shannon, Simpson),
    names_to = "Metric",
    values_to = "Diversity"
  ) %>%
  mutate(life_stage = forcats::fct_relevel(life_stage, "water", "larvae")) %>%
  ggplot(aes(x = life_stage, y = Diversity, color = life_stage)) +
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
  scale_color_manual(values = cols) +
  ggpubr::stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )
alpha_hery
ggsave("hery_alpha.pdf", dpi = 300, width = 7, height = 5)

alpha_plot <- ((alpha_bennett +
  ggtitle("Bennett et al. (2019)") |
  alpha_hery + ggtitle("Hery et al. (2021)")) /
  (alpha_rodpai +
    ggtitle("Rodpai et al. (2023)") |
    alpha_hernandez + ggtitle("Hernandez et al. (2024)"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold"),
    plot.title = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.spacing.y = unit(0, "mm"),
    legend.spacing.x = unit(0, "mm"),
    legend.position = "none"
  )

ggsave("alpha_plot.pdf", alpha_plot, dpi = 300, width = 7, height = 8)
