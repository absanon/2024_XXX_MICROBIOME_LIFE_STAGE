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
here::i_am("relative_abundance.R")
load("Sanon_16S_DADA2_data.RData")

#' Custom legend plot function
addSmallLegend <- function(myPlot, pointSize = 0.75, textSize = 6, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_markdown(size = textSize), 
          legend.text  = element_markdown(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

#' ## Relative abundance plot
ps.rel <- transform_sample_counts(SANON, function(x) x/sum(x)*100)

# agglomerate taxa
glom <- tax_glom(ps.rel, taxrank = 'Family', NArm = FALSE)

ps.melt <- psmelt(glom)
#ps.melt <- psmelt(ps.rel)

# change to character for easy-adjusted level
ps.melt$Family <- as.character(ps.melt$Family)

# Clean up taxonomy
# TODO: unitalicize non-family names for the unclassified instances

ps.melt_clean <- ps.melt %>% 
  mutate(Family=gsub(".*unclassified.*", NA, Family),
         Order=gsub(".*unclassified.*", NA, Order),
         Class=gsub(".*unclassified.*", NA, Class),
         Phylum=gsub(".*unclassified.*", NA, Phylum),
         Domain=gsub(".*unclassified.*", NA, Domain),
         Family=ifelse(is.na(Family),
                       glue("unclassified <i>{coalesce(Order, Class, Phylum)}</i>"),
                       glue("<i>{Family}</i>"))
         ) %>% 
  mutate(across(c(Domain, Phylum, Class, Order), ~ if_else(is.na(.), Family, .))) %>%
  mutate(across(c(Domain, Phylum, Class, Order, Family), ~str_replace(., "unclassified <i>NA</i>", glue("unclassified Bacteria ({OTU})")))) %>% 
  mutate(Phylum=if_else(startsWith(Family, "unclassified Bacteria"), "Others", Phylum)) %>% 
  mutate(Breeding=case_when(Breeding == "plas" ~ "plastic",
                            T ~ Breeding))

# Calculate max relative abundance
rel_abundance_df <- ps.melt_clean %>% 
  group_by(Sample, Phylum) %>% 
  mutate(phylum_abundance=sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Phylum) %>% 
  mutate(max_phylum_abundance=max(phylum_abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, Phylum, Family) %>% 
  mutate(family_abundance=sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Family) %>% 
  mutate(max_family_abundance=max(family_abundance)) %>% 
  ungroup() %>% 
  group_by(Sample, Phylum, Order) %>% 
  mutate(order_abundance=sum(Abundance)) %>% 
  ungroup() %>% 
  group_by(Order) %>% 
  mutate(max_order_abundance=max(order_abundance)) %>% 
  ungroup()

# Clean up Phylum, Order and Family names based on relative abundance
rel_abundance_clean <- rel_abundance_df %>% 
  group_by(Phylum) %>% 
  mutate(clean_Phylum = case_when(max_phylum_abundance < 1 ~ "Others",
                                T ~ Phylum),
         clean_Family = case_when(max_family_abundance < 5 & 
                                    max_family_abundance < max(max_family_abundance) ~ glue("Other {Phylum}"),
                                T ~ Family)) %>% 
  mutate(clean_Family = if_else(clean_Phylum == "Others",
                                "Others", clean_Family)) %>% 
  ungroup() %>% 
  group_by(Order) %>% 
  mutate(clean_Order = case_when(max_order_abundance < 1 ~ "Others",
                                 T ~ Order)) %>% 
  mutate(clean_Order = if_else(clean_Phylum == "Others" & 
                                 (!startsWith(Order, "unclassified Bacteria") | max_family_abundance < 5), 
                               "Others", clean_Order)) %>% 
  ungroup()

# Setting factor levels
rel_abundance_clean$Stage <- factor(rel_abundance_clean$Stage, levels = c("water", "larvae", "pupae", "adult"))
rel_abundance_clean$clean_Phylum <- factor(rel_abundance_clean$clean_Phylum, 
                                     levels = c(sort(unique(rel_abundance_clean$clean_Phylum[rel_abundance_clean$clean_Phylum!="Others"])), 
                                                                            "Others"))
rel_abundance_clean$clean_Order <- factor(rel_abundance_clean$clean_Order, 
                                           levels = c(sort(unique(rel_abundance_clean$clean_Order[rel_abundance_clean$clean_Order!="Others"])), 
                                                      "Others"))
rel_abundance_clean$clean_Family <- factor(rel_abundance_clean$clean_Family, 
                                           levels = c(sort(unique(rel_abundance_clean$clean_Family[!startsWith(rel_abundance_clean$clean_Family, "Other")])), 
                                             sort(unique(rel_abundance_clean$clean_Family[startsWith(rel_abundance_clean$clean_Family, "Other")]))))

# Create relative abundance plot
# TODO: color facet lines on life stage
pal <- c(viridisLite::viridis(
  length(unique(rel_abundance_clean$clean_Phylum))-1, 
  direction = -1), "grey90")

strip <- strip_nested(text_x = elem_list_text(face=c(rep("bold", 4), rep("plain", 24)), size=rep(6, 28)))

p <- ggnested(rel_abundance_clean, 
              aes(x = Sample, 
                  y = Abundance, 
                  main_group=clean_Phylum, 
                  sub_group = clean_Family),
         main_palette = pal,
         gradient_type = "tints",
         max_l = 1) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Relative abundance (%)") +
  facet_nested(~Stage+Breeding+Sites, scales= "free_x", 
               strip = strip, switch="x",
               nest_line = element_line(color = "black", linewidth = .2)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.1))+
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
        panel.spacing.x = unit(.15, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        legend.text = element_markdown())
p

rel_ab_plot <- addSmallLegend(p, spaceLegend = .1)+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size=8))

ggsave("figures/relative_abundance.pdf", dpi=300, width = 7, height = 5)

#' Create Proteobacteria plot

# TODO: Replace "Other Proteobacteria" with Order name
# TODO: add max family to order

pal <- c(viridisLite::viridis(
  length(unique(rel_abundance_clean$clean_Order[rel_abundance_clean$Phylum=="Proteobacteria"]))-1, 
  direction = -1, option = "plasma"), "grey80")

p2 <- rel_abundance_clean %>% 
  filter(Phylum=="Proteobacteria") %>% 
  ggnested(aes(x = Sample, 
               y = Abundance, 
               main_group=clean_Order, 
               sub_group = clean_Family),
           main_palette = pal,
           gradient_type = "tints",
           max_l = 1) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Relative abundance (%)") +
  facet_nested(~Stage+Breeding+Sites, scales= "free_x", 
               strip=strip, switch="x") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.1))+
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
        panel.spacing.x = unit(.15, "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown())

rel_ab_proteo_plot <- addSmallLegend(p2, spaceLegend = .1)+
  theme(legend.title = element_blank(),
        #legend.position = "right",
        axis.title.y = element_text(size=8))

ggsave("figures/relative_abundance_proteobacteria.pdf", dpi=300, width = 7, height = 5)

#' ### Longitudinal top ASVs
# TODO: do we see top ASV's from adult in other stages?
asv_rel <- psmelt(ps.rel)

top_ASV <- asv_rel %>% 
  filter(Stage == "adult") %>% 
  group_by(OTU, Breeding, Sites) %>% 
  summarise(mean=mean(Abundance), median=median(Abundance)) %>% 
  #group_by(Breeding, Sites) %>%
  slice_max(order_by=mean, n=15, with_ties = F) %>% 
  pull(OTU) %>% 
  unique()
  
asv_rel %>% 
  filter(OTU %in% top_ASV) %>% 
  ggplot(aes(x=Sample, y=Abundance, fill=OTU))+
  geom_bar(stat = "identity") + 
  labs(x="Samples", y="Relative abundance (%)") +
  facet_nested(~Stage+Breeding+Sites, scales= "free_x", 
               strip=strip, switch="x") +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0.1))+
  scale_fill_viridis_d()+
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
        strip.placement = "outside",
        panel.spacing.x = unit(.15, "lines"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.y = element_text(size=8),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown())

asv_rel %>% 
  filter(OTU %in% top_ASV) %>% 
  group_by(Stage, Breeding, Sites) %>% 
  mutate(mean=mean(Abundance), median=median(Abundance)) %>% 
  ungroup() %>% 
  ggplot(aes(x=Stage, y=mean))+
  geom_smooth(aes(group=-1), method = "lm", color="grey40")+
  geom_point(aes(shape=Sites, color=Breeding))+
  scale_color_viridis_d()+
  theme_bw()



# Alluvial plot
top_ASV_water <- asv_rel %>% 
  filter(Stage == "water") %>% 
  group_by(OTU) %>% 
  summarise(mean = mean(Abundance), median = median(Abundance)) %>% 
  ungroup() %>% 
  slice_max(order_by = median, n = 15) %>%
  pull(OTU)

top_ASV_adult <- asv_rel %>% 
  filter(Stage == "adult") %>% 
  group_by(OTU) %>% 
  summarise(mean = mean(Abundance), median = median(Abundance)) %>% 
  ungroup() %>%
  slice_max(order_by = median, n = 15) %>%
  pull(OTU)

top_ASV <- c(top_ASV_water, top_ASV_adult) %>% unique()


alluvial_data <- asv_rel %>% 
  filter(OTU %in% top_ASV) %>% 
  group_by(OTU, Stage) %>% 
  mutate(mean=mean(Abundance), median=median(Abundance)) %>% 
  ungroup() %>% 
  select(OTU, Order, Stage, mean, median) %>%
  distinct() %>% 
  mutate(mean_water = if_else(Stage == "water", mean, NA_real_)) %>%
  mutate(OTU = fct_reorder(OTU, mean_water, .na_rm = TRUE)) %>% 
  select(-mean_water)

alluvial_plot <- alluvial_data %>% 
  #mutate(top50_ASV=case_when(OTU %in% intersect(top_ASV_adult, top_ASV_water) ~ "both",
  #                         !OTU %in% top_ASV_adult ~"only water",
  #                         !OTU %in% top_ASV_water ~"only adult")) %>% 
  ggplot(aes(x = Stage, y = mean, alluvium = OTU, stratum=OTU, fill=OTU)) +
  #scale_fill_brewer(type = "qual", palette = "Paired")+
  scale_fill_viridis_d(option="turbo")+
  geom_flow(decreasing = TRUE) +
  geom_stratum(decreasing = TRUE, linewidth=.1) +
  geom_vline(aes(xintercept = 2.25),
             linetype="dashed")+
  annotate("text", x=2.25, y=55, label="Larvae stop eating", 
           angle=90, vjust = -1,hjust=1, size=2)+
  labs(y="Median relative abundance (%)")+
  scale_y_continuous(limits=c(0,NA), 
                     expand = expansion(add=c(0, 0.1)))+
  theme_classic()+
  theme(legend.position = "none",
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=6),
        axis.text.y = element_text(size=6),
        axis.title.y = element_text(size=8),
        legend.text = element_markdown())
alluvial_plot
ggsave("figures/alluvial_plot_top15asv.pdf", dpi=300, width=7, height = 4)

#' combine plots

rel_ab_plot / (((rel_ab_proteo_plot + free(alluvial_plot))+plot_layout(widths = c(1, .5))) / guide_area() + 
                  plot_layout(guides = 'collect'))+
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face="bold"))
ggsave("figures/combined_relative_abundance.pdf", dpi=300, width=7, height = 8)


#' Compare ASVs only in adult and water per location/breeding material
asv_stage_mean <- asv_rel %>% 
  select(-Sample) %>% 
  pivot_wider(names_from = Stage, values_from = Abundance, values_fn=mean) %>% 
  group_by(OTU) %>% 
  summarise(across(c(pupae, adult, water, larvae), \(x) mean(x, na.rm = TRUE)))

water <- asv_stage_mean %>% 
  filter(water > 0) %>% 
  pull(OTU) %>% 
  unique()

larvae <- asv_stage_mean %>% 
  filter(water == 0 & larvae > 0) %>% 
  pull(OTU) %>% 
  unique()

pupae <- asv_stage_mean %>% 
  filter(water == 0 & larvae == 0 & pupae > 0) %>% 
  pull(OTU) %>% 
  unique()

adult <- asv_stage_mean %>% 
  filter(adult > 0 & water == 0 & larvae == 0 & pupae == 0) %>% 
  pull(OTU) %>% 
  unique()

ASV_subset <- c(only_adult, only_water)

stage_relative_ab <- asv_rel %>% 
  mutate(fill_color=case_when(OTU %in% water ~ "water", 
                              OTU %in% larvae ~ "larvae", 
                              OTU %in% pupae ~ "pupae",
                              OTU %in% adult ~ "adult")) %>% 
  group_by(fill_color, Sample, Stage, Breeding, Sites) %>% 
  summarise(Abundance=sum(Abundance), .groups = "drop") %>% 
  mutate(fill_color=factor(fill_color, levels = c("water", "larvae", "pupae", "adult")))

stage_relative_ab %>% 
  #filter(Stage != "water") %>% 
  ggplot(aes(x = Sample, y = Abundance, fill=fill_color)) + 
  geom_bar(stat = "identity") + 
  labs(x="", y="Relative abundance (%)") +
  facet_nested(~Stage+Breeding+Sites, scales= "free_x", 
               strip = strip, switch="x",
               nest_line = element_line(color = "black", linewidth = .2)) +
  scale_y_continuous(limits = c(0,101), expand = c(0, 0))+
  scale_fill_brewer(palette = "Oranges", name="ASVs present from following stage:")+
  theme_classic() + 
  theme(legend.position = "bottom", 
        strip.background = element_blank(),
        strip.placement = "outside",
        #ggh4x.facet.nestline = element_line(color = list("orange", rep("black", 11)), linewidth = .2),
        panel.spacing.x = unit(.15, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        legend.text = element_markdown())
ggsave("figures/ASV_relative_abundance_per_stage.pdf", dpi=300, width=7, height=5)
