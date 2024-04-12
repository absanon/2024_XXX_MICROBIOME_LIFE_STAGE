library(tidyverse)
library(ggpubr)
library(phyloseq)
library(glue)
library(ggnested)
library(ggh4x)
library(RColorBrewer)
library(ggtext)
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

# change to character for easy-adjusted level
ps.melt$Family <- as.character(ps.melt$Family)

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

## select group mean > 1
#keep <- unique(ps.melt_clean$Family[ps.melt_clean$Abundance > 1])
#
#ps.melt_clean$Family[!(ps.melt_clean$Family %in% keep)] <- "< 1%"
#ps.melt_clean$Phylum[!(ps.melt_clean$Family %in% keep)] <- "< 1%"
#ps.melt_clean$Order[!(ps.melt_clean$Family %in% keep)] <- "< 1%"
#
###to get the same rows together
##ps.melt_sum <- ps.melt_clean %>%
##  group_by(Sample,Stage,Family) %>%
##  summarise(Abundance=sum(Abundance))
#
#unique(ps.melt_clean$Family)
#unique(ps.melt_clean$Phylum)
#unique(ps.melt_clean$Order)


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

pal <- c(viridisLite::viridis(9, direction = -1), "grey90")

strip <- strip_nested(text_x = elem_list_text(face=c(rep("bold", 4), rep("plain", 24)), size=rep(6, 28)))

p <- ggnested(rel_abundance_clean, aes(x = Sample, y = Abundance, main_group=clean_Phylum, sub_group = clean_Family),
         main_palette = pal,
         gradient_type = "tints",
         max_l = 1) + 
  geom_bar(stat = "identity") + 
  labs(x="Samples", y="Relative abundance (%)") +
  facet_nested_wrap(~Stage+Breeding+Sites, scales= "free_x", nrow=1, strip = strip) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
        panel.spacing.x = unit(.15, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown())
p

addSmallLegend(p, spaceLegend = .1)+
  theme(legend.title = element_blank(),
        axis.title = element_text(size=8))

ggsave("figures/relative_abundance.pdf", dpi=300, width = 7, height = 5)

#' ### Proteobacteria
pal <- c(viridisLite::viridis(10, direction = -1, option = "plasma"), "grey80")

p2 <- rel_abundance_clean %>% 
  filter(Phylum=="Proteobacteria") %>% 
  ggnested(aes(x = Sample, y = Abundance, main_group=clean_Order, sub_group = clean_Family),
           main_palette = pal,
           gradient_type = "tints",
           max_l = 1) + 
  geom_bar(stat = "identity") + 
  labs(x="Samples", y="Relative abundance (%)") +
  facet_nested_wrap(~Stage+Breeding+Sites, scales= "free_x", nrow=1, strip=strip) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0))+
  theme_classic() + 
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(color = "black", linewidth = .2),
        panel.spacing.x = unit(.15, "lines"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_markdown())

addSmallLegend(p2, spaceLegend = .1)+
  theme(legend.title = element_blank(),
        #legend.position = "right",
        axis.title = element_text(size=8))

ggsave("figures/relative_abundance_proteobacteria.pdf", dpi=300, width = 7, height = 5)
