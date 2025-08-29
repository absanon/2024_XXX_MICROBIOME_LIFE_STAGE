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
library(ComplexUpset)
here::i_am("relative_abundance.R")
load("Sanon_16S_DADA2_data.RData")

#' Custom legend plot function
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

#' ## Relative abundance plot
ps.rel <- transform_sample_counts(SANON, function(x) x / sum(x) * 100)

# agglomerate taxa
# glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)

# ps.melt <- psmelt(glom)
ps.melt <- psmelt(ps.rel)

# Clean up taxonomy

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
  ) |>
  # Replace "plas" with "plastic" in the Breeding column
  mutate(
    Breeding = case_when(Breeding == "plas" ~ "plastic", TRUE ~ Breeding),
    Urbanisation = case_when(
      Urbanisation == "urban" ~ "U",
      Urbanisation == "peri-urban" ~ "PU"
    )
  )

ps.melt_sum <- ps.melt_clean_tax |>
  group_by(Sample, Family, Genus) |>
  mutate(F_Abundance = sum(Abundance)) |>
  ungroup()

ps.melt_clean <- ps.melt_sum


# Calculate max relative abundance
rel_abundance_df <- ps.melt_clean |>
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
rel_abundance_clean$Stage <- factor(
  rel_abundance_clean$Stage,
  levels = c("water", "larvae", "pupae", "adult")
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

# Create table with percentages per phylum
phylum_table <- rel_abundance_clean |>
  group_by(clean_Phylum, Sample) |>
  summarize(clean_phylum_sum = round(sum(Abundance), 2)) |>
  pivot_wider(names_from = Sample, values_from = clean_phylum_sum)
write_delim(
  phylum_table,
  "data/phylum_relative_abundance.tsv",
  delim = "\t",
  col_names = T
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

strip <- strip_nested(
  text_x = elem_list_text(
    face = c(rep("bold", 4), rep("plain", 24)),
    size = rep(6, 28)
  )
)

p <- rel_abundance_clean |>
  select(-OTU, -Abundance) |>
  distinct() |>
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
    ~ Stage + Breeding + Urbanisation,
    scales = "free_x",
    strip = strip,
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown()
  )
p

rel_ab_plot <- addSmallLegend(p, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )

ggsave(
  "figures/relative_abundance.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)

# Look at 'human activity' hypothesis in plastic urban (Others classified ASVs)
rel_abundance_clean |>
  filter(
    Breeding == "plastic" & Urbanisation == "U",
    clean_Phylum == "Others"
  ) |>
  select(-OTU, -Abundance) |>
  distinct() |>
  mutate(
    Domain = ifelse(
      startsWith(Domain, "unclassified Bacteria"),
      "unclassified Bacteria",
      Domain
    )
  ) |>
  ggplot(aes(
    x = Sample,
    y = F_Abundance,
    fill = Domain
  )) +
  geom_bar(stat = "identity")
ggsave("test.pdf", dpi = 300, width = 169, height = 120, units = "mm")

# Create Proteobacteria plot

pal2 <- c(
  viridisLite::viridis(
    length(unique(rel_abundance_clean$clean_Order[
      rel_abundance_clean$Phylum == "Proteobacteria"
    ])) -
      1,
    direction = -1,
    option = "plasma"
  ),
  "grey90"
)

rel_abundance_protebacteria <- rel_abundance_clean |>
  filter(Phylum == "Proteobacteria") |>
  mutate(
    clean_Family = case_when(
      family_abundance > 1 ~ Family,
      str_detect(clean_Family, "Other") & clean_Order != "Others" ~
        paste("Other", clean_Order),
      T ~ clean_Family
    ),
    clean_Family = as.character(clean_Family),
  ) |>
  select(-OTU, -Abundance) |>
  distinct()

rel_abundance_protebacteria$clean_Family <- factor(
  rel_abundance_protebacteria$clean_Family,
  levels = c(
    sort(unique(rel_abundance_protebacteria$clean_Family[
      !startsWith(rel_abundance_protebacteria$clean_Family, "Other")
    ])),
    sort(unique(rel_abundance_protebacteria$clean_Family[startsWith(
      rel_abundance_protebacteria$clean_Family,
      "Other"
    )]))
  )
)


p2 <- rel_abundance_protebacteria |>
  ggnested(
    aes(
      x = Sample,
      y = F_Abundance,
      main_group = clean_Order,
      sub_group = clean_Family
    ),
    main_palette = pal2,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_nested(
    ~ Stage + Breeding + Urbanisation,
    scales = "free_x",
    strip = strip,
    switch = "x"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.ticks.x = element_blank(),
    legend.text = element_markdown()
  )

rel_ab_proteo_plot <- addSmallLegend(
  p2,
  spaceLegend = .5,
  pointSize = .5,
  textSize = 4
) +
  theme(
    legend.title = element_blank(),
    # legend.position = "right",
    axis.title.y = element_text(size = 6),
    legend.box.spacing = unit(0, "pt")
  )

ggsave(
  "figures/relative_abundance_proteobacteria.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)

#' ### Longitudinal top ASVs

top_ASV <- ps.melt_clean |>
  filter(Stage == "adult") |>
  group_by(OTU, Breeding, Sites) |>
  summarise(mean = mean(Abundance), median = median(Abundance)) |>
  # group_by(Breeding, Sites) |>
  slice_max(order_by = mean, n = 15, with_ties = F) |>
  pull(OTU) |>
  unique()

ps.melt_clean |>
  filter(OTU %in% top_ASV) |>
  group_by(Stage, Breeding, Sites) |>
  mutate(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  ggplot(aes(x = Stage, y = mean)) +
  geom_smooth(aes(group = -1), method = "lm", color = "grey40") +
  geom_point(aes(shape = Sites, color = Breeding)) +
  scale_color_viridis_d() +
  theme_bw()

#' Relative abundance per stage
top_ASV_water <- ps.melt_clean |>
  filter(Stage == "water") |>
  group_by(OTU) |>
  summarise(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  slice_max(order_by = median, n = 15) |>
  pull(OTU)

top_ASV_adult <- ps.melt_clean |>
  filter(Stage == "adult") |>
  group_by(OTU) |>
  summarise(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  slice_max(order_by = median, n = 15) |>
  pull(OTU)

top_ASV <- c(top_ASV_water, top_ASV_adult) |> unique()

#' Venn diagram
venn_result <- VennDiagram::venn.diagram(
  x = list(top_ASV_water, top_ASV_adult),
  category.names = c("water", "adult"),
  scaled = F,
  filename = NULL, # "figures/venn_diagram.tiff",
  # imagetype = "tiff",
  disable.logging = T,
  height = 480,
  width = 480,
  resolution = 300,
  # compression = "lzw",
  main = "Top 15 ASVs",
  main.cex = .4,
  main.pos = c(.5, .95),
  main.fontfamily = "sans",
  lwd = 1,
  col = c("steelblue", "#D94701"),
  fill = scales::alpha(c("steelblue", "#D94701"), .2),
  # Numbers
  cex = .4,
  cat.default.pos = "text",
  cat.cex = .4,
  cat.pos = c(15, 345),
  cat.dist = c(0.1, 0.1),
  cat.fontfamily = "sans",
  cat.col = c("steelblue", "#D94701"),
  label.col = c("steelblue", "#976559", "#D94701"),
  fontfamily = "sans"
)

# Alluvial plot
alluvial_data <- ps.melt_clean |>
  filter(OTU %in% top_ASV) |>
  group_by(OTU, Stage) |>
  mutate(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  select(OTU, Phylum, Stage, mean, median) |>
  distinct() |>
  mutate(mean_water = if_else(Stage == "water", mean, NA_real_)) |>
  mutate(OTU = fct_reorder(OTU, mean_water, .na_rm = TRUE)) |>
  select(-mean_water)

alluvial_plot <- alluvial_data |>
  # mutate(top50_ASV=case_when(OTU %in% intersect(top_ASV_adult, top_ASV_water) ~ "both",
  #                         !OTU %in% top_ASV_adult ~"only water",
  #                         !OTU %in% top_ASV_water ~"only adult")) |>
  ggplot(aes(x = Stage, y = mean, alluvium = OTU, stratum = OTU, fill = OTU)) +
  # scale_fill_brewer(type = "qual", palette = "Paired")+
  scale_fill_viridis_d(option = "turbo") +
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE, linewidth = .1) +
  geom_vline(aes(xintercept = 2.25), linetype = "dashed") +
  annotate(
    "text",
    x = 2.25,
    y = 55,
    label = "Larvae stop eating",
    angle = 90,
    vjust = -1,
    hjust = 1,
    size = 2
  ) +
  labs(y = "Mean relative abundance (%)") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 0.1))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    legend.text = element_markdown()
  )
alluvial_plot
ggsave(
  "figures/alluvial_plot_top15asv.pdf",
  dpi = 300,
  width = 169,
  height = 100,
  units = "mm"
)

# Alluvial plot of top ASVs per metadata variables
grouped_relabund <- rel_abundance_clean |>
  group_by(OTU, Stage, Breeding, Urbanisation) |>
  mutate(mean = mean(Abundance), median = median(Abundance)) |>
  ungroup() |>
  select(OTU, clean_Phylum, Stage, Breeding, Urbanisation, mean, median) |>
  distinct()

grouped_relabund2 <- grouped_relabund |>
  group_by(Stage, Breeding, Urbanisation) |>
  slice_max(median, n = 15) |>
  pull(OTU) |>
  unique()

get_top_ASV_stage <- function(df, stage) {
  df |>
    group_by(Stage, Breeding, Urbanisation) |>
    slice_max(median, n = 15) |> # top 15 ASVs
    filter(Stage == stage) |>
    pull(OTU) |>
    unique()
}

stage_list <- lapply(c("water", "adult", "larvae", "pupae"), function(stage) {
  get_top_ASV_stage(grouped_relabund, stage)
})

# Venn diagram
venn_stage <- VennDiagram::venn.diagram(
  x = stage_list,
  category.names = c("water", "adult", "larvae", "pupae"),
  scaled = F,
  filename = NULL, # "figures/venn_diagram.tiff",
  # imagetype = "tiff",
  disable.logging = T,
  height = 500,
  width = 500,
  resolution = 300,
  # compression = "lzw",
  main = "Top ASV distribution",
  main.cex = .4,
  main.pos = c(.5, .95),
  main.fontfamily = "sans",
  lwd = 1,
  col = c("#3B9AB2", "#7A0403FF", "#EBCC2A", "#F21A00"),
  fill = scales::alpha(c("#3B9AB2", "#7A0403FF", "#EBCC2A", "#F21A00"), .2),
  # Numbers
  cex = .4,
  cat.default.pos = "text",
  cat.cex = .4,
  # cat.pos=c(15, 345),
  cat.dist = c(0.2, 0.2, 0.1, 0.1),
  cat.fontfamily = "sans",
  cat.col = c("#3B9AB2", "#7A0403FF", "#EBCC2A", "#F21A00"),
  # label.col=c("steelblue", "#976559", "#D94701"),
  fontfamily = "sans"
)
grid::grid.newpage()
grid::grid.draw(venn_stage)

# Venn top 15 ASVs adult per location/material
get_top_adult_ASV <- function(df, breeding, urbanisation) {
  df |>
    group_by(Stage, Breeding, Urbanisation) |>
    slice_max(median, n = 50) |>
    filter(
      Stage == "adult",
      Breeding == breeding,
      Urbanisation == urbanisation
    ) |>
    pull(OTU) |>
    unique()
}

## Use expand.grid to generate all combinations of breeding and urbanisation
combinations <- expand.grid(
  breeding = c("plastic", "tire"),
  urbanisation = c("U", "PU")
)

## Apply the function to each combination of breeding and urbanisation
adult_list <- lapply(1:nrow(combinations), function(i) {
  breeding <- combinations$breeding[i]
  urbanisation <- combinations$urbanisation[i]
  get_top_adult_ASV(
    grouped_relabund,
    #rel_abundance_clean,
    breeding = breeding,
    urbanisation = urbanisation
  )
})
names(adult_list) <- apply(combinations, 1, function(row) {
  paste(row, collapse = "/")
})

## Venn diagram
venn_adult <- VennDiagram::venn.diagram(
  x = adult_list,
  category.names = c(
    "plastic/urban",
    "tire/urban",
    "plastic/peri-urban",
    "tire/peri-urban"
  ),
  scaled = F,
  filename = NULL, # "figures/venn_diagram.tiff",
  # imagetype = "tiff",
  disable.logging = T,
  height = 500,
  width = 500,
  resolution = 300,
  # compression = "lzw",
  main = "Top adult ASV distribution",
  main.cex = .4,
  main.pos = c(.5, .95),
  main.fontfamily = "sans",
  lwd = 1,
  col = c("#DFC27D", "#A6611A", "#80CDC1", "#018571"),
  fill = scales::alpha(c("#DFC27D", "#A6611A", "#80CDC1", "#018571"), .2),
  # Numbers
  cex = .4,
  cat.default.pos = "text",
  cat.cex = .4,
  # cat.pos=c(15, 345),
  cat.dist = c(0.2, 0.2, 0.1, 0.1),
  cat.fontfamily = "sans",
  cat.col = c("#DFC27D", "#A6611A", "#80CDC1", "#018571"),
  fontfamily = "sans"
)
grid::grid.newpage()
grid::grid.draw(venn_adult)

# Alluvial plot of top ASVs per metadata variables
df_alluvial <- grouped_relabund |>
  filter(OTU %in% grouped_relabund2) |>
  mutate(mean_water = if_else(Stage == "water", mean, NA_real_)) |>
  mutate(
    OTU = fct_reorder(OTU, mean_water, .na_rm = TRUE),
    Urbanisation = case_when(
      Urbanisation == "U" ~ "urban",
      Urbanisation == "PU" ~ "peri-urban",
      T ~ Urbanisation
    )
  ) |>
  select(-mean_water)

weird_top_asv <- ps.melt_clean |>
  filter(
    OTU %in% stage_list[[1]],
    Phylum %in% c("Actinobacteriota", "Bacteroidota")
  )

alluvial_grid <- df_alluvial |>
  ggplot(aes(
    x = Stage,
    y = mean,
    alluvium = OTU,
    stratum = OTU,
    fill = clean_Phylum
    # fill=as.integer(OTU)
  )) +
  # scale_fill_viridis_c(option="turbo")+
  scale_fill_manual(values = pal) +
  # scale_fill_gradientn(colors=c("#858b8cff", "#F2E191", "#F21A00"))+
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE, linewidth = .1) +
  geom_vline(aes(xintercept = 2.25), linetype = "dashed") +
  annotate(
    "text",
    x = 2.25,
    y = 75,
    label = "Larvae stop eating",
    angle = 90,
    vjust = -1,
    hjust = 1,
    size = 1.5
  ) +
  geom_point(
    data = filter(df_alluvial, OTU %in% top_ASV_adult),
    aes(y = mean * 10),
    color = "black",
    show.legend = F,
    size = 1
  ) +
  geom_smooth(
    data = filter(df_alluvial, OTU %in% top_ASV_adult),
    method = "lm",
    se = T,
    aes(x = Stage, y = mean * 10, group = 1),
    color = "blue",
    fill = "grey",
    show.legend = F,
    inherit.aes = F
  ) +
  facet_grid2(
    vars(Urbanisation),
    vars(Breeding),
    axes = "all",
    scales = "free_y"
  ) +
  labs(y = "Mean relative abundance (%)") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 1)),
    guide = guide_axis_color(color = "grey50"),
    sec.axis = sec_axis(
      ~ . / 10,
      name = "Mean relative abundance of top adult ASVs (%)",
      guide = guide_axis_color(color = "blue")
    )
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

ag <- addSmallLegend(alluvial_grid, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )

#ag +
#  inset_element(venn,
#    0.55,
#    .75,
#    .75,
#    1,
#    ignore_tag = T
#  )

# Convert the Venn diagram to a grob
venn_grobs <- grid::grobTree(venn_stage)
venn_groba <- grid::grobTree(venn_adult)

# Wrap the Venn diagram into a ggplot object using annotation_custom()
venn_stage_plot <- ggplot() +
  annotation_custom(
    venn_grobs,
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  ) +
  theme_void()
venn_adult_plot <- ggplot() +
  annotation_custom(
    venn_groba,
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  ) +
  theme_void()

((venn_adult_plot /
  venn_stage_plot) |
  ag) +
  plot_layout(widths = c(.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(
  "figures/alluvial_grid.pdf",
  dpi = 300,
  width = 169,
  units = "mm",
  height = 120
)

# Upset plot of the top 50 ASVs in adults
# create "set" identifiers = Stage/Breeding combination
df <- grouped_relabund |>
  mutate(set = paste(Breeding, Urbanisation, sep = "/")) |>
  filter(
    Stage == "adult",
    (OTU %in%
      adult_list$`plastic/U` &
      Breeding == "plastic" &
      Urbanisation == "U") |
      (OTU %in%
        adult_list$`tire/U` &
        Breeding == "tire" &
        Urbanisation == "U") |
      (OTU %in%
        adult_list$`plastic/PU` &
        Breeding == "plastic" &
        Urbanisation == "PU") |
      (OTU %in%
        adult_list$`tire/PU` &
        Breeding == "tire" &
        Urbanisation == "PU")
  )

# 2) build membership matrix: OTU x set (logical)
# xtabs will count occurrences; > 0 turns it into logical membership
membership_mat <- as.data.frame.matrix(xtabs(~ OTU + set, data = df) > 0)

# membership_mat currently has rownames = OTU; move OTU into a column
membership_mat$OTU <- rownames(membership_mat)
rownames(membership_mat) <- NULL

# 3) create OTU-level metadata (one row per OTU) -- choose how to collapse if multiple rows:
# here we take the (mean of mean) and take the first Phylum if consistent
otu_meta <- df %>%
  group_by(OTU) %>%
  summarize(
    clean_Phylum = first(clean_Phylum),
    mean_abund = mean(mean, na.rm = TRUE),
    .groups = "drop"
  )

# 4) join membership + metadata
membership_df <- membership_mat %>% left_join(otu_meta, by = "OTU")

# 5) determine set column names for ComplexUpset
set_cols <- setdiff(
  colnames(membership_df),
  #c("OTU", "clean_Phylum", "mean_abund")
  c(
    "mean_abund",
    "clean_Phylum",
  )
)

# ComplexUpset plot
p <- upset(
  membership_df,
  intersect = set_cols,
  min_size = 1,
  width_ratio = 0.2,
  set_sizes = FALSE, # hide left set-size plot per your request
  base_annotations = list(
    #  # intersection_size uses geom="bar" with a weight aesthetic to sum fractional contributions
    'ASV count' = intersection_size(
      counts = T,
      mapping = aes(fill = clean_Phylum),
    ) +
      scale_fill_manual(values = pal) +
      theme(legend.title = element_blank(), axis.text.x = element_blank())
  ),
  annotations = list(
    # extra panel: quasirandom points for mean, with mean point
    'Mean relative\n abundance (%)' = ggplot(mapping = aes(y = mean_abund)) +
      ggbeeswarm::geom_quasirandom()
  ),
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(axis.title.x = element_blank())
    )
  )
)
print(p)

ggsave(
  "figures/upset_top50_adult.pdf",
  dpi = 300,
  width = 169,
  height = 100,
  units = "mm"
)

#' Compare ASVs only in adult and water per location/breeding material
asv_stage_mean <- ps.melt_clean |>
  select(-Sample) |>
  pivot_wider(names_from = Stage, values_from = Abundance, values_fn = mean) |>
  group_by(OTU) |>
  summarise(across(c(pupae, adult, water, larvae), \(x) mean(x, na.rm = TRUE)))

water <- asv_stage_mean |>
  filter(water > 0) |>
  pull(OTU) |>
  unique()

larvae <- asv_stage_mean |>
  filter(water == 0 & larvae > 0) |>
  pull(OTU) |>
  unique()

pupae <- asv_stage_mean |>
  filter(water == 0 & larvae == 0 & pupae > 0) |>
  pull(OTU) |>
  unique()

adult <- asv_stage_mean |>
  filter(adult > 0 & water == 0 & larvae == 0 & pupae == 0) |>
  pull(OTU) |>
  unique()

stage_relative_ab <- ps.melt_clean |>
  mutate(
    fill_color = case_when(
      OTU %in% water ~ "water",
      OTU %in% larvae ~ "larvae",
      OTU %in% pupae ~ "pupae",
      OTU %in% adult ~ "adult"
    )
  ) |>
  group_by(fill_color, Sample, Stage, Breeding, Urbanisation) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  mutate(
    fill_color = factor(
      fill_color,
      levels = c("water", "larvae", "pupae", "adult")
    )
  )

ASV_ra_per_stage_p <- stage_relative_ab |>
  # filter(Stage != "water") |>
  ggplot(aes(x = Sample, y = Abundance, fill = fill_color)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Relative abundance (%)") +
  facet_nested(
    ~ Stage + Breeding + Urbanisation,
    scales = "free_x",
    strip = strip,
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  # scale_fill_brewer(palette = "Oranges", name="ASVs present starting from:")+
  scale_fill_manual(
    values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"),
    breaks = c("water", "larvae", "pupae", "adult"),
    name = "ASVs present starting from:"
  ) +
  theme_classic() +
  guides(
    shape = guide_legend(override.aes = list(size = .75)),
    color = guide_legend(override.aes = list(size = .75))
  ) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside",
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.title = element_markdown(size = 6),
    legend.text = element_markdown(size = 6),
    legend.key.size = unit(.5, "lines"),
    legend.box.spacing = unit(0, "pt")
  )
ASV_ra_per_stage_p
ggsave(
  "figures/ASV_relative_abundance_per_stage.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)

# Which taxonomy for per stage ASVs?
p3 <- rel_abundance_clean |>
  filter(OTU %in% c(larvae, pupae, adult)) |>
  distinct() |>
  mutate(
    stage_present = case_when(
      OTU %in% larvae ~ "larvae",
      OTU %in% pupae ~ "pupae",
      OTU %in% adult ~ "adult"
    ),
    stage_present = factor(
      stage_present,
      levels = c("larvae", "pupae", "adult")
    )
  ) |>
  group_by(
    Sample,
    Stage,
    Breeding,
    clean_Phylum,
    Urbanisation,
    clean_Family,
    stage_present
  ) |>
  summarise(ab = sum(Abundance)) |>
  ungroup() |>
  select(
    Sample,
    Stage,
    Breeding,
    clean_Phylum,
    Urbanisation,
    clean_Family,
    ab,
    stage_present
  ) |>
  distinct() |>
  ggnested(
    aes(
      x = Sample,
      y = ab,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal,
    gradient_type = "tints",
    max_l = 1
  ) +
  ggnewscale::new_scale_color() +
  geom_bar(aes(color = stage_present), stat = "identity", lwd = .3, lty = 2) +
  labs(x = "", y = "Relative abundance (%)") +
  facet_nested(
    ~ Stage + Breeding + Urbanisation,
    scales = "free_x",
    strip = strip,
    switch = "x",
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  # scale_color_brewer(palette="Oranges")+
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown()
  )

addSmallLegend(p3, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  ) +
  guides(color = "none")

ggsave(
  "figures/relative_abundance_of_stage_ASVs.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)
#' combine plots

rel_ab_plot /
  ((free(alluvial_plot) +
    inset_element(venn_result, 0.02, .5, .32, .9, ignore_tag = T) +
    ASV_ra_per_stage_p) +
    plot_layout(widths = c(.7, 1))) +
  plot_layout(heights = c(1, .5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(
  "figures/combined_relative_abundance.pdf",
  dpi = 300,
  width = 7,
  height = 8
)

# Absolute quantification with 16S qPCR
#qpcr <- read_tsv("data/16S_qPCR_results.tsv")
qpcr <- read_csv("data/16S_qPCR_results.csv") |>
  mutate(Sample = gsub(" ", "-", Sample)) |>
  filter(Task != "Standard")

copies_per_stage <- qpcr |>
  left_join(samdf |> rownames_to_column("Sample")) |>
  filter(!is.na(Stage)) |>
  ggplot(aes(x = Stage, y = Quantity, color = Stage)) +
  geom_boxplot(outlier.shape = NA, show.legend = F) +
  geom_jitter(
    aes(shape = gsub(pattern = "_", replacement = "/", x = breeding_urban)),
    width = 0.3
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.text.align = 0
  ) +
  scale_y_log10() +
  scale_color_manual(values = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF")) +
  scale_shape_manual(values = c(0, 15, 1, 16)) +
  labs(shape = "Breeding material/\nLocation", y = "copies/µl") +
  stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )
ggsave("figures/copies_per_stage.pdf", dpi = 300, height = 7, width = 5)

qpcr |>
  left_join(samdf |> rownames_to_column("Sample")) |>
  filter(!is.na(Stage)) |>
  ggplot(aes(x = Breeding, y = Quantity, color = Breeding)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    aes(shape = gsub(pattern = "_", replacement = "/", x = breeding_urban)),
    width = 0.3
  ) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    axis.title.x = element_blank(),
    legend.text.align = 0
  ) +
  scale_y_log10() +
  labs(shape = "Breeding material/\nLocation") +
  facet_wrap(~Stage, ncol = 2, nrow = 2) +
  scale_shape_manual(values = c(0, 15, 1, 16)) +
  stat_pwc(
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.adj.format",
    hide.ns = "p.adj",
    show.legend = F,
    tip.length = 0.01
  )

absolute_quantity_df <- rel_abundance_clean |>
  left_join(qpcr) |>
  mutate(abs_quantity = (Abundance / 100) * Quantity) |>
  group_by(OTU, Stage, Breeding, Urbanisation) |>
  mutate(mean = mean(abs_quantity), median = median(abs_quantity)) |>
  ungroup() |>
  select(
    OTU,
    clean_Phylum,
    Stage,
    Breeding,
    Urbanisation,
    mean,
    median,
    abs_quantity
  ) |>
  #distinct() |>
  filter(OTU %in% stage_list[[2]]) |>
  mutate(mean_water = if_else(Stage == "water", mean, NA_real_)) |>
  mutate(
    OTU = fct_reorder(OTU, mean_water, .na_rm = TRUE),
    Urbanisation = case_when(
      Urbanisation == "U" ~ "urban",
      Urbanisation == "PU" ~ "peri-urban",
      T ~ Urbanisation
    )
  ) |>
  select(-mean_water)

# Absolute abundance barplot
p_abs <- rel_abundance_clean |>
  left_join(qpcr) |>
  mutate(F_Abundance_abs = ((F_Abundance / 100) * Quantity)) |>
  select(-OTU, -Abundance) |>
  distinct() |>
  ggnested(
    aes(
      x = Sample,
      y = F_Abundance_abs,
      main_group = clean_Phylum,
      sub_group = clean_Family
    ),
    main_palette = pal,
    gradient_type = "tints",
    max_l = 1
  ) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Absolute abundance (copies/µl") +
  facet_nested(
    Stage ~ Breeding + Urbanisation,
    scales = "free",
    #strip = strip,
    switch = "x",
    independent = T,
    nest_line = element_line(color = "black", linewidth = .2)
  ) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0.1)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    # ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
    panel.spacing.x = unit(.15, "lines"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    legend.text = element_markdown()
  )
p_abs

abs_ab_plot <- addSmallLegend(p_abs, spaceLegend = .5) +
  theme(
    legend.title = element_blank(),
    axis.title.y = element_text(size = 8),
    legend.box.spacing = unit(0, "pt")
  )

ggsave(
  "figures/absolute_abundance.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)

# Top 15 ASV absolute abundance
top_ASV_abs_quant <- absolute_quantity_df |>
  select(-abs_quantity) |>
  distinct() |>
  ggplot(aes(
    x = Stage,
    y = mean,
    alluvium = OTU,
    stratum = OTU,
    fill = clean_Phylum
    # fill=as.integer(OTU)
  )) +
  # scale_fill_viridis_c(option="turbo")+
  scale_fill_manual(values = pal) +
  # scale_fill_gradientn(colors=c("#3B9AB2", "#F2E191", "#F21A00"))+
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE, linewidth = .1) +
  geom_vline(aes(xintercept = 2.25), linetype = "dashed") +
  #annotate("text",
  #  x = 2.25, y = 61, label = "Larvae stop eating",
  #  angle = 90, vjust = -1, hjust = 1, size = 2
  #) +
  facet_grid2(
    vars(Urbanisation),
    vars(Breeding),
    axes = "all",
    scales = "free",
    independent = "all"
  ) +
  labs(y = "Mean absolute abundance (copies/µl)") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(add = c(0, 1))
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
    legend.text = element_markdown()
  )

ggsave(
  "figures/absolute_quantity_top_adult_asv.pdf",
  dpi = 300,
  width = 10,
  height = 7
)

(free(
  addSmallLegend(copies_per_stage) +
    theme(
      legend.position = "bottom",
      # legend.title = element_text(vjust = 1),
      legend.box.just = "left",
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.box = "vertical",
      legend.box.spacing = unit(1, "mm"),
      legend.spacing = unit(1, "mm"),
      legend.margin = margin(t = 0.1, unit = 'cm')
    )
) +
  guides(
    color = guide_legend(nrow = 2, title.position = "top"), # Number of rows for Group 1 legend
    shape = guide_legend(nrow = 2, title.position = "top") # Number of rows for Group 2 legend
  ) |
  addSmallLegend(top_ASV_abs_quant) +
    theme(legend.key.size = unit(2, "mm"), legend.margin = margin(0, 0, 0, 0)) +
    theme(legend.title = element_blank())) +
  plot_layout(widths = c(.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("figures/absolute_combined.pdf", dpi = 300, width = 7.4, height = 5)

# Venn all ASVs

get_ASV_stage <- function(df, stage, abundance_filter) {
  df |>
    filter(Abundance > abundance_filter, Stage == stage) |>
    pull(OTU) |>
    unique()
}

venn_list <- lapply(c("water", "larvae", "pupae", "adult"), function(stage) {
  get_ASV_stage(ps.melt, stage, 0)
})

venn_0 <- VennDiagram::venn.diagram(
  x = venn_list,
  category.names = c("water", "larvae", "pupae", "adult"),
  scaled = F,
  filename = NULL, # "figures/venn_diagram.tiff",
  # print.mode = c("raw", "percent"),
  # imagetype = "tiff",
  disable.logging = T,
  height = 480,
  width = 480,
  resolution = 300,
  # compression = "lzw",
  main = "ASV distribution",
  main.cex = .4,
  main.pos = c(.5, .95),
  main.fontfamily = "sans",
  lwd = 1,
  col = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"),
  fill = scales::alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"), .2),
  # Numbers
  cex = .4,
  cat.default.pos = "text",
  cat.cex = .4,
  # cat.pos=c(15, 345),
  cat.dist = c(0.2, 0.2, 0.1, 0.1),
  cat.fontfamily = "sans",
  cat.col = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"),
  # label.col=c("steelblue", "#976559", "#D94701"),
  fontfamily = "sans"
)

venn_list <- lapply(c("water", "larvae", "pupae", "adult"), function(stage) {
  get_ASV_stage(ps.melt, stage, 1)
})

venn_1 <- VennDiagram::venn.diagram(
  x = venn_list,
  category.names = c("water", "larvae", "pupae", "adult"),
  scaled = F,
  filename = NULL, # "figures/venn_diagram.tiff",
  # print.mode = c("raw", "percent"),
  # imagetype = "tiff",
  disable.logging = T,
  height = 480,
  width = 480,
  resolution = 300,
  # compression = "lzw",
  main = "ASV distribution",
  main.cex = .4,
  main.pos = c(.5, .95),
  main.fontfamily = "sans",
  lwd = 1,
  col = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"),
  fill = scales::alpha(c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"), .2),
  # Numbers
  cex = .4,
  cat.default.pos = "text",
  cat.cex = .4,
  # cat.pos=c(15, 345),
  cat.dist = c(0.2, 0.2, 0.1, 0.1),
  cat.fontfamily = "sans",
  cat.col = c("#3B9AB2", "#EBCC2A", "#F21A00", "#7A0403FF"),
  # label.col=c("steelblue", "#976559", "#D94701"),
  fontfamily = "sans"
)

for (i in 9:23) {
  venn_0[[i]]$label <- paste0(venn_0[[i]]$label, "\n(", venn_1[[i]]$label, ")")
}

# Convert the Venn diagram to a grob
venn_grob <- grid::grobTree(venn_0)

# Wrap the Venn diagram into a ggplot object using annotation_custom()
venn_plot <- ggplot() +
  annotation_custom(
    venn_grob,
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  ) +
  theme_void()

rel_ab_plot /
  ((ASV_ra_per_stage_p +
    free(venn_plot)) +
    plot_layout(widths = c(1, .7))) +
  plot_layout(heights = c(1, .5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(
  "figures/combined_relative_abundance.pdf",
  dpi = 300,
  width = 7,
  height = 8
)

(venn_plot | rel_ab_proteo_plot) +
  plot_layout(widths = c(.5, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))
ggsave(
  "figures/combined_relative_abundance_proteo.pdf",
  dpi = 300,
  width = 169,
  height = 120,
  units = "mm"
)
