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
load("16S_Aeg_public/hery/hery_dada2.RData")

hery_meta <- read_tsv("16S_Aeg_public/hery.tsv")
hery_meta <- as.data.frame(hery_meta) %>%
  filter(
    Experiment %in% sample_names(otu_table(t(seqtab.nochim), taxa_are_rows = T))
  )
hery_meta <- hery_meta %>%
  separate(
    isolation_source,
    into = c("life_stage", "container"),
    sep = "_",
    remove = FALSE
  ) %>%
  mutate(
    life_stage = str_remove_all(life_stage, "[0-9[:space:]b]"),
    container = str_remove_all(container, "[0-9[:space:]]")
  )
rownames(hery_meta) <- hery_meta$Experiment

ps <- phyloseq(
  otu_table(t(seqtab.nochim), taxa_are_rows = T),
  sample_data(hery_meta),
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
rel_abundance_clean$life_stage <- factor(
  rel_abundance_clean$life_stage,
  levels = c("water", "larvae")
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
  #  group_by(container, geographic_location, collection_date) %>%
  #  filter(all(c("water", "larvae") %in% life_stage)) |>
  #  ungroup() %>%
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
    ~ geographic_location + collection_date + container + life_stage,
    scales = "free_x",
    space = "free",
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  ggtitle("Hery et al. (2021)") +
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
  "16S_Aeg_public/hery.pdf",
  dpi = 300,
  width = 169,
  height = 200,
  units = "mm"
)

# Alluvial plot

top_ASV <- ps.melt_clean |>
  filter(life_stage == "larvae") |>
  select(OTU, life_stage, Abundance) |>
  drop_na() |>
  group_by(OTU) |>
  summarise(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  slice_max(order_by = mean, n = 15) |>
  pull(OTU)

# Alluvial plot of top ASVs per metadata variables
grouped_relabund <- rel_abundance_clean |>
  group_by(OTU, life_stage, container, geographic_location, collection_date) |>
  mutate(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  select(
    OTU,
    clean_Phylum,
    life_stage,
    container,
    geographic_location,
    collection_date,
    mean,
    median
  ) |>
  distinct()

grouped_relabund2 <- grouped_relabund |>
  group_by(life_stage, container, geographic_location, collection_date) |>
  slice_max(mean, n = 15) |>
  pull(OTU) |>
  unique()

# get_top_ASV_stage <- function(df, stage) {
#  df |>
#    group_by(life_stage, lat_lon) |>
#    slice_max(mean, n = 15) |>
#    filter(life_stage == stage) |>
#    pull(OTU) |>
#    unique()
# }

# stage_list <- lapply(c("eggs", "larves", "adult"), function(stage) get_top_ASV_stage(grouped_relabund, stage))

df_alluvial <- grouped_relabund |>
  filter(OTU %in% grouped_relabund2) |>
  mutate(mean_water = if_else(life_stage == "water", mean, NA_real_)) |>
  mutate(
    life_stage = factor(life_stage, levels = c("water", "larvae")),
    OTU = fct_reorder(OTU, mean_water, .na_rm = TRUE)
  ) |>
  select(-mean_water)

alluvial_grid_g <- df_alluvial |>
  filter(geographic_location == "Guadeloupe") |>
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
    vars(container),
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
addSmallLegend(alluvial_grid_g, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )
ggsave(
  "16S_Aeg_public/hery_alluvial_guadeloupe.pdf",
  dpi = 300,
  width = 7,
  height = 10
)

alluvial_grid_fg <- df_alluvial |>
  filter(geographic_location == "French Guiana") |>
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
    vars(container),
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
addSmallLegend(alluvial_grid_fg, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )
ggsave(
  "16S_Aeg_public/hery_alluvial_french_guiana.pdf",
  dpi = 300,
  width = 7,
  height = 10
)


save.image("16S_Aeg_public/hery/hery_dada2.RData")

group_vars <- c("life_stage", "geographic_location")

# Ensure variables exist
if (!all(group_vars %in% colnames(sample_data(ps)))) {
  stop("One or more grouping variables not found in sample_data.")
}

# Create a combined group column in sample_data
sample_data(ps)$GroupID <- apply(
  sample_data(ps)[, group_vars],
  1,
  paste,
  collapse = "_"
)

ps_merged <- merge_samples(ps, group = "GroupID")

ps_group_rel <- transform_sample_counts(ps_merged, function(x) x / sum(x) * 100)

ps.melt_group <- psmelt(ps_group_rel)

ps.melt_clean_tax <- ps.melt_group |>
  select(
    -geographic_location,
    -life_stage,
    -sample_name,
    -isolation_source,
    -container
  ) |>
  tidyr::separate_wider_delim(
    cols = Sample,
    delim = "_",
    names = c("life_stage", "geographic_location"),
    cols_remove = FALSE
  ) |>
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
rel_abundance_clean$life_stage <- factor(
  rel_abundance_clean$life_stage,
  levels = c("water", "larvae")
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
  #  group_by(container, geographic_location, collection_date) %>%
  #  filter(all(c("water", "larvae") %in% life_stage)) |>
  #  ungroup() %>%
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
    ~ geographic_location + life_stage,
    scales = "free_x",
    space = "free",
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  ggtitle("Hery et al. (2021)") +
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
  "16S_Aeg_public/hery_grouped.pdf",
  dpi = 300,
  width = 169,
  height = 140,
  units = "mm"
)
