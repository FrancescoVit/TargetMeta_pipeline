#' ---
#' title: "Script 3 of pipeline: read in workspace from script1, run beta diversity analysis"
#' author: "Francesco Vitali"
#' output: html_document
#' ---

#' ```{r setup library load, include=FALSE}


library(phyloseq)
library(tidyverse)
library(data.table)
library(ggsci)
library(ggpubr)
library(microbiome)
library(patchwork)
library(metagenomeSeq)
library(vegan)
library(reshape2)
library(psych)
library(corrplot)
library(rstatix)
library(pairwiseAdonis)

# red in workspace from script 1

load ("./WP3WP4_workspace")

palette_custom <- c("#1F78B4", "#E31A1C", "#FF7F00", "#33A02C") 

#' knitr::opts_chunk$set(fig.width=unit(15,"cm"), fig.height=unit(11,"cm"))

#' ```

#' ```

#' *BACTERIA*

#############

#' **WP3 BACTERIA**

#############
#' ```{r betadiv bact WP3, include=TRUE, echo = FALSE}

# how to transform data for beta diversity analysis? 


# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP3_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP3_initial_bact), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)
# variation in read depth is 7X 

# plot reads to find outliers

sdt %>%
  filter(AnimalID != "E149-40") %>%
  ggboxplot(.,y = "TotalReads")

sdt %>%
  filter(AnimalID != "E149-40") %>%
  summary()


# NORMALIZATION: open problematic of which to choose

# see https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y
# see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003531
# see https://www.nature.com/articles/s41522-020-00160-w

# new kid on the block! : avgdist() supported by Schloss 
# (see https://www.youtube.com/watch?v=xyufizOpc5I) also (https://www.youtube.com/watch?v=ht3AX5uZTTQ) also (https://www.youtube.com/watch?v=t5qXPIS-ECU)
# So, the whole thing is about rarefaction or not. He support rarefaction as is the only method that is removing some relation between BC distance and 
# depth of sampling (n reads). To note, however, is that he is not doing the "classical" rarefaction (i.e. as in script2.r for alpha-div) using subsampling 
# and replacement, as it is not that good, as it really all depends on the seed and on which sequences get subsampled.
# Instead, in the videos above, he use a nice combination of rarefaction and vegdist: he use rarefaction analysis followed by BC calculation multiple time, 
# so that the final distance matrix is an average od the distance matrices obtained over the different replicates and the relation with n° seqs disappear, but 
# without the risk of normal (single) rarefaction of loosing OTUs only by chance of the random procedure.
# So, he is using un-rarefied and un-transformed counts, but rather performs this normalization step at the level of distance calculation. 
# From the video, and also conceptually, seems very convincing.

# CSS remain one of the most interesting and supported, but rarefaction is still an option. 
# Evaluate both Bray and Jaccard eventually, in one case to evaluate changes in abundance, in the other
# changes as presence/absence

# Data driven evaluation: transform with rarefaction and with CSS, calculate bray and jaccard, compare with Mantel

# generate a 10% prevalece filter dataset (see https://www.nature.com/articles/s41467-022-28034-z)

# calculate prevalence for each ASV
WP3.bact.ASVprev <- prevalence(subset_samples(WP3_initial_bact, AnimalID != "E149-40"), detection=0, sort=TRUE, count=F)
summary(names(WP3.bact.ASVprev[WP3.bact.ASVprev > 0.1]))
summary(names(WP3.bact.ASVprev[WP3.bact.ASVprev < 0.1]))
# create an object with only the ASV in more than 10% of samples
WP3_initial_bact_prevalence <- prune_taxa(names(WP3.bact.ASVprev[WP3.bact.ASVprev > 0.1]), subset_samples(WP3_initial_bact, AnimalID != "E149-40")) 
#check
min(prevalence(WP3_initial_bact_prevalence, detection=0, sort=TRUE, count=F))

# rarefaction on the full dataset
WP3_raref_bact<- rarefy_even_depth(subset_samples(WP3_initial_bact, AnimalID != "E149-40"),
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP3_initial_bact, AnimalID != "E149-40")))),
                                   replace=F)

# rarefaction on the prevalence dataset
WP3_raref_bact.prevalence<- rarefy_even_depth(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40"),
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40")))),
                                   replace=F)

# CSS scaling (not removing singletons) on the full dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(subset_samples(WP3_initial_bact, AnimalID != "E149-40"))
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP3_CSS_bact <- subset_samples(WP3_initial_bact, AnimalID != "E149-40")
otu_table(WP3_CSS_bact) <- otu_table(data.CSS, taxa_are_rows = T)

# CSS scaling (not removing singletons) on the prevalence dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40"))
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP3_CSS_bact.prevalence <- subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40")
otu_table(WP3_CSS_bact.prevalence) <- otu_table(data.CSS, taxa_are_rows = T)


# calculate distance matrices
bray_WP3_raref <- phyloseq::distance(WP3_raref_bact, "bray")
bray_WP3_CSS <- phyloseq::distance(WP3_CSS_bact, "bray")

bray_WP3_raref.prevalence <- phyloseq::distance(WP3_raref_bact.prevalence, "bray")
bray_WP3_CSS.prevalence <- phyloseq::distance(WP3_CSS_bact.prevalence, "bray")

# check if order of the matrix is the same
summary(melt(as.matrix(bray_WP3_raref))[1] == melt(as.matrix(bray_WP3_CSS))[1])

# perform mantel test
set.seed(35264)
mantel(xdis = bray_WP3_CSS, ydis = bray_WP3_raref, permutations = 999, method = "pearson")

set.seed(35264)
mantel(xdis = bray_WP3_CSS, ydis = bray_WP3_CSS.prevalence, permutations = 999, method = "pearson")

set.seed(35264)
mantel(xdis = bray_WP3_raref, ydis = bray_WP3_raref.prevalence, permutations = 999, method = "pearson")

mds.css <- monoMDS(bray_WP3_CSS)
mds.raref <- monoMDS(bray_WP3_raref)
WP3.bact.proc <- procrustes(mds.raref, mds.css)
summary(WP3.bact.proc)
plot(WP3.bact.proc)

# see if biological interpretation is the same or not
df <- as(sample_data(WP3_CSS_bact), "data.frame")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df, AnimalID)

set.seed(12387)
adonis2(bray_WP3_CSS ~ Time * Diet * Sacrif_dfference, data = df, permutations = perm)
set.seed(12387)
adonis2(bray_WP3_CSS.prevalence ~ Time * Diet * Sacrif_dfference, data = df, permutations = perm)

set.seed(12387)
adonis2(bray_WP3_raref ~ Time * Diet * Sacrif_dfference, data = df, permutations = perm)
set.seed(12387)
adonis2(bray_WP3_raref.prevalence ~ Time * Diet * Sacrif_dfference, data = df, permutations = perm)

# significance and R2 only show minor changes, but the Mantel correl between the two dist matrices is not high
# overall, higher R2 in all factors (time diet and sacrifice difference) with CSS.prevalence dataset

# This set the basis of choosing the CSS transformation on the prevalence (10%) filtered dataset as this seems
# to show higher power to discriminate over the exp factors

# Perform ordination now

# NMDS
set.seed(34624)
WP3_bact_ord <- ordinate(WP3_CSS_bact.prevalence, "NMDS", distance = "bray")
p.WP3_bact_ord <- plot_ordination(physeq = WP3_CSS_bact.prevalence, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8, aes(fill = Diet, shape = Time)) #+ geom_text(mapping = aes(label = Gender), size = 3, nudge_x = 0.015, nudge_y = -0.015)
p.WP3_bact_ord <- p.WP3_bact_ord + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(21,22,23,24))+
  xlab("Dim1") + ylab("Dim2") + labs( fill = "Diet") + theme(legend.position="right")

p.WP3_bact_ord

df <- as(sample_data(WP3_CSS_bact.prevalence), "data.frame")
perm <- how(nperm = 999)
setBlocks(perm) <- with(df, AnimalID)

bray_WP3_CSS.prevalence <- phyloseq::distance(WP3_CSS_bact.prevalence, "bray")

set.seed(12387)
adonis2(bray_WP3_CSS.prevalence ~ Time + Diet + Sacrif_dfference, data = df, permutations = perm)

# NMDS with Schloss approach

# starting data is the untransformed DF, will still remove the lower outlier and remove low prevalence OTUs
phylo_schloss <- subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40")
otu_table_schloss <- otu_table(phylo_schloss)@.Data

otu_table_schloss %>%
  t() %>%
  as.data.frame() -> otu_table_schloss

# adding some metadata
row.names(otu_table_schloss) == row.names(sample_data(phylo_schloss))
otu_table_schloss <- cbind(sample_data(phylo_schloss)[,3:4], otu_table_schloss)

# calculate depth
nSeqs <- min(sample_sums(phylo_schloss))

# calculate rarefaction of BC distance and plot MDS
set.seed(43655826)
raref_BC <- avgdist(x = otu_table_schloss[,-c(1,2)], sample = nSeqs, dmethod = "bray", iterations = 100)
set.seed(62853)
schloss_MDS <- metaMDS(raref_BC)

scores(schloss_MDS) %>%
  as_tibble(rownames = "Group") %>%
  ggplot(aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 5, stroke = 0.8, aes(fill = otu_table_schloss$Diet, shape = otu_table_schloss$Time))+
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(21,22,23,24))+
  xlab("Dim1") + ylab("Dim2") + labs( fill = "Diet") + theme(legend.position="right")

# PCOA

WP3_CSS_bact.prevalenceT3 <- subset_samples(WP3_CSS_bact.prevalence, Time == "T3")

WP3_bact_ord <- ordinate(WP3_CSS_bact.prevalenceT3, "PCoA", distance = "bray")
p.WP3_bact_ord <- plot_ordination(physeq = WP3_CSS_bact.prevalenceT3, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + 
  geom_point(size = 3, stroke = 0.8, aes(fill = Diet), shape = 22) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  xlab("PCoA1 (17.4%)") + 
  ylab("PCoA2 (11.6%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP3_bact_ord_1

p.WP3_bact_ord -> p.WP3_bact_ord_final

# pairwise adonis
df <- as(sample_data(WP3_CSS_bact.prevalenceT3), "data.frame")
d.bacteria = phyloseq::distance(WP3_CSS_bact.prevalenceT3, "bray")
set.seed(12387)
pw_ad_bactT3 <- pairwise.adonis(d.bacteria, factors = df$Diet, p.adjust.m = "BH", perm = 999)

pw_ad_bactT3$R2

# PCoA using the Schloss avgdist method
WP3_bact_ord <- ordinate(phylo_schloss, "PCoA",distance = raref_BC)
p.WP3_bact_ord <- plot_ordination(physeq = phylo_schloss, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8, aes(fill = Diet, shape = Time)) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(21,22,23,24))+
  xlab("PCoA1 (21.2%)") + 
  ylab("PCoA2 (8.6%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP3_bact_ord_1



# plot on distance from sacrifice
p.WP3_bact_ord <- plot_ordination(physeq = WP3_CSS_bact.prevalence, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8, aes(fill = Sacrif_dfference, shape = Time)) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_viridis_c(option = "C") +
  scale_shape_manual(values = c(21,22,23,24))+
  xlab("PCoA1 (17.3%)") + 
  ylab("PCoA2 (11.5%)") + 
  labs(fill = "Days from \nfirst sacrif.", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP3_bact_ord_2

p.WP3_bact_ord_1 / p.WP3_bact_ord_2 -> p_suppl_ordi

ggsave('./Supplementary figure 2.png', p_suppl_ordi)

# Which component capture best experimental factors?

plot_scree(WP3_bact_ord, title = NULL)

# extract coordinates on the first 20 axis, mount a df with exp factors
df <- as(sample_data(WP3_CSS_bact.prevalence), "data.frame")
ordi.coordinates <- as.data.frame(WP3_bact_ord$vectors[,1:10])
row.names(df) == row.names(ordi.coordinates) # row name is in same order
ordi.coordinates <- cbind(ordi.coordinates, df[,c(3,4,9:14)]) # attach some variables

# BEST AXIS FOR TIME
ordi.coordinates %>%
  wilcox_effsize(Axis.1 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.2 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.3 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.4 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.5 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.6 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.7 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.8 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.9 ~ Time)
ordi.coordinates %>%
  wilcox_effsize(Axis.10 ~ Time)


ordi.coordinates %>%
  wilcox_test(Axis.1 ~ Time)
ordi.coordinates %>%
  wilcox_test(Axis.2 ~ Time)

# BEST AXIS FOR Diet
ordi.coordinates %>%
  wilcox_effsize(Axis.1 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.2 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.3 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.4 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.5 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.6 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.7 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.8 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.9 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.10 ~ Diet)


## Other numeric variables with correlation

ordi.coordinates %>%
  select(-c(Diet, start_body_weight, Time)) %>%
  corr.test(scale(.),method = "spearman", adjust = "none") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r),  order="original", method = "number",tl.cex = 0.8,number.cex = 1, type = "lower",
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)



# Maybe it has more sense to see beta diversity at T3, where the diet should have had an effect

WP3_CSS_bact.prevalence.T3 <- subset_samples(WP3_CSS_bact.prevalence, Time == "T3") %>% prune_taxa(taxa_sums(.) > 0, .)

# PCOA
WP3_bact_ord <- ordinate(WP3_CSS_bact.prevalence.T3, "PCoA", distance = "bray")
p.WP3_bact_ord <- plot_ordination(physeq = WP3_CSS_bact.prevalence.T3, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8, shape = 21, aes(fill = Diet)) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  xlab("PCoA1 (17.4%)") + 
  ylab("PCoA2 (11.6%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23)))  -> p.WP3_bact_ord.T3

p.WP3_bact_ord

## With Schloss
# starting data is the untransformed DF, will still remove the lower outlier and remove low prevalence OTUs
phylo_schloss <- subset_samples(subset_samples(WP3_initial_bact_prevalence, Time == "T3"), AnimalID != "E149-40")
otu_table_schloss <- otu_table(phylo_schloss)@.Data

otu_table_schloss %>%
  t() %>%
  as.data.frame() -> otu_table_schloss

# adding some metadata
row.names(otu_table_schloss) == row.names(sample_data(phylo_schloss))
otu_table_schloss <- cbind(sample_data(phylo_schloss)[,3:4], otu_table_schloss)

# calculate depth
nSeqs <- min(sample_sums(phylo_schloss))

# calculate rarefaction of BC distance and plot MDS
set.seed(43655826)
raref_BC <- avgdist(x = otu_table_schloss[,-c(1,2)], sample = nSeqs, dmethod = "bray", iterations = 100)

# PCoA using the Schloss avgdist method
WP3_bact_ord <- ordinate(phylo_schloss, "PCoA",distance = raref_BC)
  p.WP3_bact_ord <- plot_ordination(physeq = phylo_schloss, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8, aes(fill = Diet, shape = Time)) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(21,22,23,24))+
  xlab("PCoA1 (24.1%)") + 
  ylab("PCoA2 (12.3%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP3_bact_ord_1


# see if biological interpretation is the same or not
bray_WP3_CSS.prevalence.T3 <- phyloseq::distance(WP3_CSS_bact.prevalence.T3, "bray")
df <- as(sample_data(WP3_CSS_bact.prevalence.T3), "data.frame")
set.seed(12387)
adonis2(bray_WP3_CSS.prevalence.T3 ~ Diet + Sacrif_dfference + end_body_weight, data = df, permutations = 999)


# plot on distance from sacrifice
p.WP3_bact_ord <- plot_ordination(physeq = WP3_CSS_bact.prevalence.T3, WP3_bact_ord, axes = c(1,2))
p.WP3_bact_ord <- p.WP3_bact_ord + geom_point(size = 5, stroke = 0.8,shape = 21, aes(fill = Sacrif_dfference)) 
p.WP3_bact_ord <- p.WP3_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_viridis_c(option = "C") +
  xlab("PCoA1 (17.3%)") + 
  ylab("PCoA2 (11.5%)") + 
  labs(fill = "Days from \nfirst sacrif.", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP3_bact_ord.T3

# Which component capture best experimental factors?

plot_scree(WP3_bact_ord, title = NULL)

# extract coordinates on the first 20 axis, mount a df with exp factors
df <- as(sample_data(WP3_CSS_bact.prevalence.T3), "data.frame")
ordi.coordinates <- as.data.frame(WP3_bact_ord$vectors[,1:10])
row.names(df) == row.names(ordi.coordinates) # row name is in same order
ordi.coordinates <- cbind(ordi.coordinates, df[,c(4,9:14)]) # attach some variables

# for sacrifice difference and other numeric var, I could do a scatterplot matrix

ordi.coordinates %>%
  select(-c(Diet, start_body_weight)) %>%
  corr.test(scale(.),method = "spearman", adjust = "none") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r),  order="original", method = "number",tl.cex = 0.8,number.cex = 1,
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)

# sacrifice difference strongly correlate with axis.2 

# BEST AXIS FOR Diet


ordi.coordinates %>%
  wilcox_test(Axis.1 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.2 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.3 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.4 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.5 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.6 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.7 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.8 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.9 ~ Diet)
ordi.coordinates %>%
  wilcox_test(Axis.10 ~ Diet)

# axis 1:5 have a significant anova respect to diet, should calculate which has the higher effect size

ordi.coordinates %>%
  wilcox_effsize(Axis.1 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.2 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.3 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.4 ~ Diet)
ordi.coordinates %>%
  wilcox_effsize(Axis.5 ~ Diet)


#############

#' **WP4 BACTERIA**

###########

#' ```{r betadiv bact WP4, include=TRUE, echo = FALSE}

# Optimization of transformation and of methods was performed only for WP3, here take the same 

# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP4_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP4_initial_bact), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)

# plot reads to find outliers

sdt %>%
  filter(AnimalID != "WP4-928") %>%
  ggboxplot(.,y = "TotalReads")

sdt %>%
  filter(AnimalID != "WP4-928") %>%
  summary()


# generate a 10% prevalece filter dataset (see https://www.nature.com/articles/s41467-022-28034-z)

# calculate prevalence for each ASV
WP4.bact.ASVprev <- prevalence(subset_samples(WP4_initial_bact, AnimalID != "WP4-928"), detection=0, sort=TRUE, count=F)
summary(names(WP4.bact.ASVprev[WP4.bact.ASVprev > 0.1]))
summary(names(WP4.bact.ASVprev[WP4.bact.ASVprev < 0.1]))
# create an object with only the ASV in more than 10% of samples
WP4_initial_bact_prevalence <- prune_taxa(names(WP4.bact.ASVprev[WP4.bact.ASVprev > 0.1]), subset_samples(WP4_initial_bact, AnimalID != "WP4-928")) 
#check
min(prevalence(WP4_initial_bact_prevalence, detection=0, sort=TRUE, count=F))

# CSS scaling (not removing singletons) on the prevalence dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(subset_samples(WP4_initial_bact_prevalence, AnimalID != "WP4-928"))
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP4_CSS_bact.prevalence <- subset_samples(WP4_initial_bact_prevalence, AnimalID != "WP4-928")
otu_table(WP4_CSS_bact.prevalence) <- otu_table(data.CSS, taxa_are_rows = T)

# Perform ordination now
# PCOA

WP4_bact_ord <- ordinate(WP4_CSS_bact.prevalence, "PCoA", distance = "bray")
p.WP4_bact_ord <- plot_ordination(physeq = WP4_CSS_bact.prevalence, WP4_bact_ord, axes = c(1,2))
p.WP4_bact_ord <- p.WP4_bact_ord + geom_point(size = 3, stroke = 0.8, aes(fill = Inoculum_diet, shape = Time)) 
p.WP4_bact_ord <- p.WP4_bact_ord + 
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(24,21,22,23))+
  xlab("PCoA1 (24%)") + 
  ylab("PCoA2 (14.5%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) -> p.WP4_bact_ord_1

## Plot with WP3
figure2 <- p.WP3_bact_ord_final / p.WP4_bact_ord + plot_layout(guides = "collect") & theme(legend.position = 'left') 

#ggsave('./Figure 2.png', figure2)

# pairwise adonis
df <- as(sample_data(subset_samples(WP4_CSS_bact.prevalence, Time == "T3" & Type != "Inoculum")), "data.frame")
d.bacteria = phyloseq::distance(subset_samples(WP4_CSS_bact.prevalence, Time == "T3"  & Type != "Inoculum"), "bray")
set.seed(12387)
pw_ad_bactWP4 <- pairwise.adonis(d.bacteria, factors = df$Diets, p.adjust.m = "BH", perm = 999)

mat1 <- data.frame(matrix(nrow=4,ncol=4,byrow=TRUE))

mat1[1,1] <- 0
mat1[1,2] <- pw_ad_bactT3$R2[1]
mat1[1,3] <-pw_ad_bactT3$R2[2]
mat1[1,4] <-pw_ad_bactT3$R2[4]
mat1[2,1] <- pw_ad_bactT3$R2[1]
mat1[2,2] <- 0
mat1[2,3] <-pw_ad_bactT3$R2[4]
mat1[2,4] <-pw_ad_bactT3$R2[5]
mat1[3,1] <-pw_ad_bactT3$R2[2] 
mat1[3,2] <- pw_ad_bactT3$R2[4]
mat1[3,3] <- 0
mat1[3,4] <-pw_ad_bactT3$R2[6]
mat1[4,1] <- pw_ad_bactT3$R2[3]
mat1[4,2] <- pw_ad_bactT3$R2[5]
mat1[4,3] <-pw_ad_bactT3$R2[6]
mat1[4,4] <- 0

colnames(mat1) <- c("MBDT", "PVD", "MBD", "CTR")
row.names(mat1) <- c("MBDT", "PVD", "MBD", "CTR")

corrplot(as.matrix(mat1), type = "lower", diag = F, method = 'number', col = "black")

###########

#' **WP3 WP4 BACTERIA**

###########

# here analysis of WP3 at T3 and WP4 all timepoint, to highlight the "passage" in the FTM
# optionally, one could select only the ID that actually were used for the pool

WP3WP4_initial_bact <- subset_samples(phyloseq_obj_initial_bact,  
                                      AnimalID != "E149-40" &  AnimalID != "WP4-928") %>% prune_taxa(taxa_sums(.) > 0, .)



# select sample subset of interest (i.e. only T3 of WP3)

WP3WP4_initial_bact <- subset_samples(WP3WP4_initial_bact, WP != "WP3" | Time != "T0")

# obtain rarefied object for alpha diversity

WP3WP4_raref_bact<- rarefy_even_depth(WP3WP4_initial_bact,
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(WP3WP4_initial_bact))),
                                   replace=F)

# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP3WP4_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP3WP4_initial_bact), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)

# plot reads to find outliers

sdt %>%
  ggboxplot(.,y = "TotalReads")

sdt %>%
  summary()

# generate a 10% prevalece filter dataset (see https://www.nature.com/articles/s41467-022-28034-z)

# calculate prevalence for each ASV
WP34.bact.ASVprev <- prevalence(WP3WP4_initial_bact, detection=0, sort=TRUE, count=F)
summary(names(WP34.bact.ASVprev[WP34.bact.ASVprev > 0.1]))
summary(names(WP34.bact.ASVprev[WP34.bact.ASVprev < 0.1]))
# create an object with only the ASV in more than 10% of samples
WP34_initial_bact_prevalence <- prune_taxa(names(WP34.bact.ASVprev[WP34.bact.ASVprev > 0.1]), WP3WP4_initial_bact) 
#check
min(prevalence(WP34_initial_bact_prevalence, detection=0, sort=TRUE, count=F))

# CSS scaling (not removing singletons) on the prevalence dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(WP34_initial_bact_prevalence)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP34_CSS_bact.prevalence <- WP34_initial_bact_prevalence
otu_table(WP34_CSS_bact.prevalence) <- otu_table(data.CSS, taxa_are_rows = T)

sample_data(WP34_CSS_bact.prevalence)$Time

WP34_bact_ord <- ordinate(WP34_CSS_bact.prevalence, "PCoA", distance = "bray")
p.WP34_bact_ord <- plot_ordination(physeq = WP34_CSS_bact.prevalence, WP34_bact_ord, axes = c(1,2))
p.WP34_bact_ord <- p.WP34_bact_ord + 
  geom_point(size = 3, stroke = 0.8, aes(fill = Diets, shape = Time))+
  #geom_line(aes(group = AnimalID)) + 
  geom_vline(xintercept= 0, linetype="dotted") +
  geom_hline(yintercept= 0, linetype="dotted") +
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(24,21,22,25))+
  xlab("PCoA1 (24.1%)") + 
  ylab("PCoA2 (11.9%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) + 
  theme(legend.position="none")

p.WP34_bact_ord

# Calculate alpha-div
# richness
alphadiv_WP34_bact <- microbiome::richness(WP3WP4_raref_bact)
df <- as(sample_data(WP3WP4_raref_bact), "data.frame")
df <- cbind(df, alphadiv_WP34_bact)

row.names(alphadiv_WP34_bact$vectors) == row.names(df)

df <- cbind(df, WP34_bact_ord$vectors[,1:2])

ggscatter(data = df, x = "Axis.1", y = "observed", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA1") + 
  ylab("Richness (n° of observed ASVs)")

ggscatter(data = df, x = "Axis.2", y = "observed", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA2") + 
  ylab("Richness (n° of observed ASVs)")

ggboxplot(data = df, x = "WP", y = "observed", add = "jitter", 
          fill = "WP", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(label.y = 100) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  xlab("") + 
  ylab("Richness (n° of observed ASVs)") + 
  coord_flip() + 
  theme(legend.position="none") -> sottoplot_b_ord


# evenness
alphadiv_WP34_bact <- microbiome::evenness(WP3WP4_raref_bact)
df <- as(sample_data(WP3WP4_raref_bact), "data.frame")
df <- cbind(df, alphadiv_WP34_bact)

row.names(alphadiv_WP34_bact$vectors) == row.names(df)

df <- cbind(df, WP34_bact_ord$vectors[,1:2])

ggscatter(data = df, x = "Axis.1", y = "pielou", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA1") + 
  ylab("Richness (n° of observed ASVs)")

ggscatter(data = df, x = "Axis.2", y = "pielou", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA2") + 
  ylab("Richness (n° of observed ASVs)")

ggboxplot(data = df, x = "WP", y = "pielou", add = "jitter", 
          fill = "WP", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(label.y = .3, label.x = 2) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  xlab("") + 
  ylab("Pielou's Evenness") + 
  coord_flip() + 
  theme(legend.position="none") -> sottoplot_b_ord_eve


## Plot with other single 

figure2 <- p.WP34_bact_ord + p.WP3_bact_ord_final/p.WP4_bact_ord + plot_layout(guides = "collect") & theme(legend.position = 'bottom')

#ggsave('./Figure 2.png', figure2)

layout <- "
A
A
A
B
C
"

p.WP34_bact_ord / sottoplot_b_ord / sottoplot_b_ord_eve+ 
  plot_layout(design = layout) -> panel_ord_bact

#####

#' *FUNGI*

#' **WP3 FUNGI**

#' ```{r betadiv fungi WP3, include=TRUE, echo = FALSE}


#' **WP4 FUNGI**

#' ```{r betadiv fungi WP4, include=TRUE, echo = FALSE}

# how to transform data for alpha diversity analysis? 

#########################
#' **WP3 WP4 FUNGI**
######################

# here analysis of WP3 at T3 and WP4 all timepoint, to highlight the "passage" in the FTM
# optionally, one could select only the ID that actually were used for the pool

WP3WP4_initial_fung <- subset_samples(phyloseq_obj_initial_fungi,  
                                      AnimalID != "E149-40" &  AnimalID != "WP4-928") %>% prune_taxa(taxa_sums(.) > 0, .)

# select sample subset of interest (i.e. only T3 of WP3)

WP3WP4_initial_fung <- subset_samples(WP3WP4_initial_fung, WP != "WP3" | Time != "T0")

# obtain rarefied object for alpha diversity

WP3WP4_raref_fung<- rarefy_even_depth(WP3WP4_initial_fung,
                                      rngseed=1234,
                                      sample.size=round(0.99*min(sample_sums(WP3WP4_initial_fung))),
                                      replace=F)


# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP3WP4_initial_fung), "data.frame"),
                 TotalReads = sample_sums(WP3WP4_initial_fung), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)

# plot reads to find outliers

sdt %>%
  ggboxplot(.,y = "TotalReads")

sdt %>%
  summary()

# generate a 10% prevalece filter dataset (see https://www.nature.com/articles/s41467-022-28034-z)

# calculate prevalence for each ASV
WP34.fung.ASVprev <- prevalence(WP3WP4_initial_fung, detection=0, sort=TRUE, count=F)
summary(names(WP34.fung.ASVprev[WP34.fung.ASVprev > 0.01]))
summary(names(WP34.fung.ASVprev[WP34.fung.ASVprev < 0.01]))
# create an object with only the ASV in more than 10% of samples
WP34_initial_fung_prevalence <- prune_taxa(names(WP34.fung.ASVprev[WP34.fung.ASVprev > 0.1]), WP3WP4_initial_fung) 
#check
min(prevalence(WP34_initial_fung_prevalence, detection=0, sort=TRUE, count=F))

# CSS scaling (not removing singletons) on the prevalence dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(WP34_initial_fung_prevalence)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP34_CSS_fung.prevalence <- WP34_initial_fung_prevalence
otu_table(WP34_CSS_fung.prevalence) <- otu_table(data.CSS, taxa_are_rows = T)

sample_data(WP34_CSS_fung.prevalence)$Time

  
# Ordination and plot
  
WP34_fung_ord <- ordinate(WP34_CSS_fung.prevalence, "PCoA", distance = "bray")
p.WP34_fung_ord <- plot_ordination(physeq = WP34_CSS_fung.prevalence, WP34_fung_ord, axes = c(1,2))
p.WP34_fung_ord <- p.WP34_fung_ord + 
  geom_point(size = 3, stroke = 0.8, aes(fill = Diets, shape = Time))+
  #geom_line(aes(group = AnimalID)) + 
  geom_vline(xintercept= 0, linetype="dotted") +
  geom_hline(yintercept= 0, linetype="dotted") +
  theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank()) + 
  scale_fill_manual(values = palette_custom) +
  scale_shape_manual(values = c(24,21,22,25))+
  xlab("PCoA1 (40.9%)") + 
  ylab("PCoA2 (5.4%)") + 
  labs(fill = "Diet", shape = "Time") + 
  theme(legend.position="right") +
  guides(fill = guide_legend(override.aes = list(shape = 23))) + 
  theme(legend.position="none")

p.WP34_fung_ord


# Calculate alpha-div
alphadiv_WP34_fung <- microbiome::richness(WP3WP4_raref_fung)
df <- as(sample_data(WP3WP4_raref_fung), "data.frame")
df <- cbind(df, alphadiv_WP34_fung)

row.names(WP34_fung_ord$vectors) == row.names(df)

df <- cbind(df, WP34_fung_ord$vectors[,1:2])

ggscatter(data = df, x = "Axis.1", y = "observed", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA1") + 
  ylab("Richness (n° of observed ASVs)")

ggscatter(data = df, x = "Axis.2", y = "observed", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA2") + 
  ylab("Richness (n° of observed ASVs)")

ggboxplot(data = df, x = "WP", y = "observed", add = "jitter", 
          fill = "WP", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(label.y = 90) +
  #ggtitle(label = "Fungi", subtitle =  "Richness (n° of ASVs)") +
  xlab("") + 
  ylab("Richness (n° of observed ASVs)") + 
  coord_flip() + 
  theme(legend.position="none") -> sottoplot_f_ord


# evenness
alphadiv_WP34_fung <- microbiome::evenness(WP3WP4_raref_fung)
df <- as(sample_data(WP3WP4_raref_fung), "data.frame")
df <- cbind(df, alphadiv_WP34_fung)

row.names(alphadiv_WP34_fung$vectors) == row.names(df)

df <- cbind(df, WP34_fung_ord$vectors[,1:2])

ggscatter(data = df, x = "Axis.1", y = "pielou", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA1") + 
  ylab("Pielou's Evenness")

ggscatter(data = df, x = "Axis.2", y = "pielou", add = "reg.line", add.params = list(linetype = "dotted")) + 
  stat_cor(method = "pearson") + 
  xlab("PCoA2") + 
  ylab("Pielou's Evenness")

ggboxplot(data = df, x = "WP", y = "pielou", add = "jitter", 
          fill = "WP", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(label.y = .2, label.x = 2) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  xlab("") + 
  ylab("Pielou's Evenness") + 
  coord_flip() + 
  theme(legend.position="none") -> sottoplot_f_ord_eve


  
layout <- "
A
A
A
B
C
"

p.WP34_fung_ord / sottoplot_f_ord / sottoplot_f_ord_eve+ 
  plot_layout(design = layout) -> panel_ord_fung


ggarrange(panel_ord_bact, panel_ord_fung, 
          labels = c("A", "B"),
          ncol = 2) -> figure2


ggsave('./Figure 2.png', figure2)

#' ```

#' ```{r write images, include=FALSE}
#' knitr::opts_chunk$set(echo = TRUE)

# use ggsave for images

#' ```

