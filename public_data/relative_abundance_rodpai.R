library(tidyverse)
library(ggpubr)
library(phyloseq)
library(glue)
library(ggnested)
library(ggh4x)
library(RColorBrewer)
library(ggtext)
library(ggalluvial)
library(patchwork)
library(microbiome)

setwd("/lustre1/scratch/337/vsc33750/")
load("16S_Aeg_public/rodpai/rodpai_dada2.RData")

rodpai_meta <- read_tsv("16S_Aeg_public/rodpai.tsv")
rodpai_meta <- as.data.frame(rodpai_meta) %>%
  filter(
    Experiment %in% sample_names(otu_table(t(seqtab.nochim), taxa_are_rows = T))
  )
rownames(rodpai_meta) <- rodpai_meta$Experiment

ps <- phyloseq(
  otu_table(t(seqtab.nochim), taxa_are_rows = T),
  sample_data(rodpai_meta),
  tax_table(taxid)
)
ps
ps <- microbiome::add_refseq(ps)

# Custom legend plot function
addSmallLegend <- function(
  myPlot,
  pointSize = 0.75,
  textSize = 6,
  spaceLegend = 0.1
) {
  myPlot +
    guides(
      shape = guide_legend(override.aes = list(size = pointSize)),
      color = guide_legend(override.aes = list(size = pointSize))
    ) +
    theme(
      legend.title = element_markdown(size = textSize),
      legend.text = element_markdown(size = textSize),
      legend.key.size = unit(spaceLegend, "lines")
    )
}

# Relative abundance plot
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x) * 100)
ps.melt <- psmelt(ps.rel)

# Define a function to replace "unclassified" values with NA
replace_unclassified <- function(x) {
  ifelse(grepl(".*unclassified.*", x), NA, x)
}

# Define a function to replace "unclassified NA" with a specified string
replace_unclassified_na <- function(x, replacement) {
  str_replace(x, "unclassified NA", replacement)
}

ps.melt_clean_tax <- ps.melt |>
  # Replace "unclassified" values with NA
  mutate(
    Family = replace_unclassified(Family),
    Order = replace_unclassified(Order),
    Class = replace_unclassified(Class),
    Phylum = replace_unclassified(Phylum),
    Domain = replace_unclassified(Domain)
  ) |>
  # Create a new column "Family" with appropriate values
  mutate(
    Family = ifelse(
      is.na(Family),
      glue("unclassified {coalesce(Order, Class, Phylum)}"),
      glue("<i>{Family}</i>")
    )
  ) |>
  # Replace NA values in Domain, Phylum, Class, and Order with values from Family
  mutate(across(
    c(Domain, Phylum, Class, Order),
    ~ if_else(is.na(.), Family, .)
  )) |>
  # Replace "unclassified NA" with "unclassified Bacteria ({OTU})"
  mutate(across(
    c(Domain, Phylum, Class, Order, Family),
    ~ replace_unclassified_na(., glue("unclassified Bacteria ({OTU})"))
  )) |>
  # Replace Phylum values if Family starts with "unclassified Bacteria"
  mutate(
    Phylum = if_else(
      startsWith(Family, "unclassified Bacteria"),
      "Others",
      Phylum
    )
  )

ps.melt_sum <- ps.melt_clean_tax |>
  group_by(Sample, Family, Genus) |>
  mutate(F_Abundance = sum(Abundance)) |>
  ungroup()

ps.melt_clean <- ps.melt_sum

# Calculate max relative abundance
rel_abundance_df <- ps.melt_clean |>
  filter(!is.na(Abundance)) %>%
  group_by(Sample, Phylum) |>
  mutate(phylum_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Phylum) |>
  mutate(max_phylum_abundance = max(phylum_abundance)) |>
  ungroup() |>
  group_by(Sample, Phylum, Family) |>
  mutate(family_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Family) |>
  mutate(max_family_abundance = max(family_abundance)) |>
  ungroup() |>
  group_by(Sample, Phylum, Order) |>
  mutate(order_abundance = sum(Abundance)) |>
  ungroup() |>
  group_by(Order) |>
  mutate(max_order_abundance = max(order_abundance)) |>
  ungroup()

# Clean up Phylum, Order and Family names based on relative abundance
rel_abundance_clean <- rel_abundance_df |>
  group_by(Phylum) |>
  mutate(
    clean_Phylum = case_when(
      max_phylum_abundance < 1 ~ "Others",
      T ~ Phylum
    ),
    clean_Family = case_when(
      max_family_abundance < 5 &
        max_family_abundance < max(max_family_abundance) &
        !startsWith(Family, "unclassified Bacteria") ~
        glue::glue("Other {Phylum}"),
      clean_Phylum == "Others" & startsWith(Family, "unclassified Bacteria") ~
        "unclassified Bacteria",
      T ~ Family
    )
  ) |>
  mutate(
    clean_Family = if_else(
      clean_Phylum == "Others" & clean_Family != "unclassified Bacteria",
      "Others",
      clean_Family
    )
  ) |>
  ungroup() |>
  group_by(Order) |>
  mutate(
    clean_Order = case_when(
      max_order_abundance < 1 ~ "Others",
      T ~ Order
    )
  ) |>
  mutate(
    clean_Order = if_else(
      clean_Phylum == "Others" &
        (!startsWith(Order, "unclassified Bacteria") |
          max_family_abundance < 5),
      "Others",
      clean_Order
    )
  ) |>
  ungroup()

# Setting factor levels
rel_abundance_clean <- rel_abundance_clean %>%
  mutate(
    isolation_source = case_when(
      isolation_source == "Mosquito abdomen" ~ "Adult",
      isolation_source == "Mosquito larvae" ~ "Larvae",
      isolation_source == "Mosquito larvae habitat water" ~ "Water",
      T ~ isolation_source
    )
  )

rel_abundance_clean$isolation_source <- factor(
  rel_abundance_clean$isolation_source,
  levels = c("Water", "Larvae", "Adult")
)
rel_abundance_clean$clean_Phylum <- factor(
  rel_abundance_clean$clean_Phylum,
  levels = c(
    sort(unique(rel_abundance_clean$clean_Phylum[
      rel_abundance_clean$clean_Phylum != "Others"
    ])),
    "Others"
  )
)
rel_abundance_clean$clean_Order <- factor(
  rel_abundance_clean$clean_Order,
  levels = c(
    sort(unique(rel_abundance_clean$clean_Order[
      rel_abundance_clean$clean_Order != "Others"
    ])),
    "Others"
  )
)
rel_abundance_clean$clean_Family <- factor(
  rel_abundance_clean$clean_Family,
  levels = c(
    sort(unique(rel_abundance_clean$clean_Family[
      !startsWith(rel_abundance_clean$clean_Family, "Other")
    ])),
    sort(unique(rel_abundance_clean$clean_Family[startsWith(
      rel_abundance_clean$clean_Family,
      "Other"
    )]))
  )
)

# Create relative abundance plot
pal <- c(
  viridisLite::viridis(
    length(unique(rel_abundance_clean$clean_Phylum)) - 1,
    direction = -1
  ),
  "grey70"
)

names(pal) <- levels(rel_abundance_clean$clean_Phylum)

p <- rel_abundance_clean |>
  select(-OTU, -Abundance) |>
  distinct() %>%
  #group_by(container_description, province) %>%
  #filter(all(c("Larvae", "Adult") %in% life_stage)) |>
  #ungroup() %>%
  ggnested(
    aes(
      x = Sample,
      y = F_Abundance,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_nested(
    ~ collection_date + lat_lon + isolation_source,
    scales = "free_x",
    space = "free",
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  ggtitle("Rodpai et al. (2023)") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text = element_text(size = 6),
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown(),
    plot.title = element_text(size = 8)
  )
addSmallLegend(p, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )
ggsave(
  "16S_Aeg_public/rodpai.pdf",
  dpi = 300,
  width = 169,
  height = 140,
  units = "mm"
)

# Alluvial plot
ps.melt_clean <- ps.melt_clean |>
  mutate(
    life_stage = case_when(
      isolation_source == "Mosquito abdomen" ~ "Adult",
      isolation_source == "Mosquito larvae" ~ "Larvae",
      isolation_source == "Mosquito larvae habitat water" ~ "Water",
      T ~ isolation_source
    )
  )

rel_abundance_clean <- rel_abundance_clean |>
  mutate(
    life_stage = case_when(
      isolation_source == "Mosquito abdomen" ~ "Adult",
      isolation_source == "Mosquito larvae" ~ "Larvae",
      isolation_source == "Mosquito larvae habitat water" ~ "Water",
      T ~ isolation_source
    )
  )

top_ASV <- ps.melt_clean |>
  filter(life_stage == "Adult") |>
  select(OTU, life_stage, Abundance) |>
  drop_na() |>
  group_by(OTU) |>
  summarise(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  slice_max(order_by = mean, n = 15) |>
  pull(OTU)

# Alluvial plot of top ASVs per metadata variables
grouped_relabund <- rel_abundance_clean |>
  group_by(OTU, life_stage, lat_lon, collection_date) |>
  mutate(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  select(
    OTU,
    clean_Phylum,
    life_stage,
    lat_lon,
    collection_date,
    mean,
    median
  ) |>
  distinct()

grouped_relabund2 <- grouped_relabund |>
  group_by(life_stage, lat_lon, collection_date) |>
  slice_max(mean, n = 15) |>
  pull(OTU) |>
  unique()

get_top_ASV_stage <- function(df, stage) {
  df |>
    group_by(life_stage, lat_lon) |>
    slice_max(mean, n = 15) |>
    filter(life_stage == stage) |>
    pull(OTU) |>
    unique()
}

# stage_list <- lapply(c("eggs", "larves", "adult"), function(stage) get_top_ASV_stage(grouped_relabund, stage))

df_alluvial <- grouped_relabund |>
  filter(OTU %in% grouped_relabund2) |>
  mutate(mean_larvae = if_else(life_stage == "Larvae", mean, NA_real_)) |>
  mutate(
    life_stage = factor(life_stage, levels = c("Water", "Larvae", "Adult")),
    OTU = fct_reorder(OTU, mean_larvae, .na_rm = TRUE)
  ) |>
  select(-mean_larvae)

alluvial_grid <- df_alluvial |>
  ggplot(aes(
    x = life_stage,
    y = mean,
    alluvium = OTU,
    stratum = OTU,
    fill = clean_Phylum
    # fill=as.integer(OTU)
  )) +
  scale_fill_manual(values = pal) +
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE, linewidth = .1) +
  facet_grid2(
    vars(lat_lon),
    vars(collection_date),
    axes = "all",
    scales = "free_y"
  ) +
  labs(y = "Mean relative abundance (%)") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 1)),
    guide = guide_axis_color(color = "grey50"),
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "lightgrey", linetype = 0),
    strip.placement = "outside",
    strip.text = element_text(size = 8),
    strip.switch.pad.grid = unit(.5, "cm"),
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 8, color = "grey50"),
    axis.title.y.right = element_text(color = "blue", vjust = -17),
    legend.text = element_markdown()
  )
addSmallLegend(alluvial_grid, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )
ggsave("16S_Aeg_public/rodpai_alluvial.pdf", dpi = 300, width = 7, height = 10)

save.image("16S_Aeg_public/rodpai/rodpai_dada2.RData")
