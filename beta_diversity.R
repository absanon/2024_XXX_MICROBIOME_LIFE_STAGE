library(vegan)
here::i_am("beta_diversity.R")

source("alpha_diversity.R")

permM<- adonis2(t(otu)~Sites*Stage,data= meta, permutations=999, 
                        method="bray", by= "terms",na.rm=T)
permM
permM1<- adonis2(t(otu)~Sites+Stage,data= meta, permutations=999, 
                method="bray", by= "terms",na.rm=T)
permM1
permM2<- adonis2(t(otu)~Breeding*Stage,data= meta, permutations=999, 
                method="bray", by= "terms",na.rm=T)
permM2
permM3<- adonis2(t(otu)~Breeding+Stage,data= meta, permutations=999, 
                 method="bray", by= "terms",na.rm=T)
permM3
permM4<- adonis2(t(otu)~Sites*Stage*Breeding,data= meta, permutations=999, 
                 method="bray", by= "terms",na.rm=T)
permM4

#### Plots for Beta diversity 
# Transform data to proportions as appropriate for Bray-Curtis distances

ps.prop <- transform_sample_counts(SANON, function(otu) otu/sum(otu))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ord.pca.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

plot_ordination(ps.prop, ord.pca.bray, color="Breeding",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Stage",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Sites",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Stage",title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Breeding", shape="Sites", title="Bray NMDS")+geom_point(size=3)

plot_ordination(ps.prop, ord.pca.bray, color="Stage", shape="Sites", title="Bray NMDS")+geom_point(size=3)