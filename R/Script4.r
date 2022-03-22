#' ---
#' title: "Script 4 of pipeline: differential abundance anaysis"
#' author: "Francesco Vitali"
#' output: html_document
#' ---

#' ```{r setup library load, include=FALSE}


library(phyloseq)
library(tidyverse)
library(data.table)
library(ggsci)
library(ggpubr)
library(Maaslin2)
library(microbiome)
library(genefu)
library(metagenomeSeq)

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

# generate a 10% prevalece filter dataset (see https://www.nature.com/articles/s41467-022-28034-z)
# calculate prevalence for each ASV
WP3.bact.ASVprev <- prevalence(subset_samples(WP3_initial_bact, AnimalID != "E149-40"), detection=0, sort=TRUE, count=F)
summary(names(WP3.bact.ASVprev[WP3.bact.ASVprev > 0.1]))
summary(names(WP3.bact.ASVprev[WP3.bact.ASVprev < 0.1]))
# create an object with only the ASV in more than 10% of samples
WP3_initial_bact_prevalence <- prune_taxa(names(WP3.bact.ASVprev[WP3.bact.ASVprev > 0.1]), subset_samples(WP3_initial_bact, AnimalID != "E149-40"))
#check
min(prevalence(WP3_initial_bact_prevalence, detection=0, sort=TRUE, count=F))

# NORMALIZATION: CSS

# CSS scaling (not removing singletons) on the full dataset
data.metagenomeSeq = phyloseq_to_metagenomeSeq(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40"))
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
WP3_CSS_bact <- subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40")
otu_table(WP3_CSS_bact) <- otu_table(data.CSS, taxa_are_rows = T)


# NORMALIZATION: TSS

WP3_TSS_bact <- transform(WP3_initial_bact_prevalence, "compositional")

# NORMALIZATION: rarefaction
WP3_raref_bact<- rarefy_even_depth(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40"),
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP3_initial_bact_prevalence, AnimalID != "E149-40")))),
                                   replace=F)

WP3_raref_TSS_bact <- transform(WP3_raref_bact, "compositional")

# RUN MAASLIN2 WITH CSS

# build input data: df with ASV
df_input_data <- as.data.frame(otu_table(tax_glom(WP3_CSS_bact, taxrank = "Genus")))
# build input data: change row.names from OTU to phylum_genus
tax_input_data <- as.data.frame(tax_table(tax_glom(WP3_CSS_bact, taxrank = "Genus"))) 
tax_input_data$label <- paste(tax_input_data$Phylum, tax_input_data$Genus, sep = ".")
# remove spaces
tax_input_data$label <- gsub(" ", replacement = ".", x = tax_input_data$label)
# solve duplicates
duplicated <- genefu::rename.duplicate(x=tax_input_data$label, verbose=TRUE)
tax_input_data$label <- duplicated$new.x
# check
row.names(tax_input_data) == row.names(df_input_data)
# rename
row.names(df_input_data) <- tax_input_data$label
#transpose
df_input_data <- t(df_input_data)

# build input metadata:

df_input_metadata <- as(sample_data(WP3_CSS_bact), "data.frame")

# build new variable pasting the diet and time variables, to use and test for interaction between the two;
# maaslin has no proper interaction term like there would be for example in a ANOVA

df_input_metadata$TimeDiet <- paste(df_input_metadata$Time, df_input_metadata$Diet, sep = "_")

# run maaslin on all data, this allows to not loose statistical power

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "fit_data_WP3_bact_CSS_LM",
  fixed_effects = c("Time", "Diet", "TimeDiet"),
  random_effects = c("AnimalID", "end_body_weight"),
  reference = "Diet,CTR;Time,T0;TimeDiet,T0_CTR",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  cores = 7)


# Select only T3 then run with numerical vars

df_input_metadata %>%
  filter(Time == "T3") -> df_input_metadata_T3

df_input_data[row.names(df_input_data) %in% row.names(df_input_metadata_T3),] -> df_input_data_T3

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data_T3,
  input_metadata = df_input_metadata_T3,
  output = "fit_data_WP3_bact_CSS_LM_colonTum",
  fixed_effects = c("Colon_tumor"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  cores = 7)

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data_T3,
  input_metadata = df_input_metadata_T3,
  output = "fit_data_WP3_bact_CSS_CPLM_colonTum",
  fixed_effects = c("Colon_tumor"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "CPLM",
  cores = 7)



# fit_data_WP3_bact <- Maaslin2(
#   input_data = df_input_data,
#   input_metadata = df_input_metadata,
#   output = "fit_data_WP3_bact_CSS_ZINB",
#   fixed_effects = c("Time", "Diet", "TimeDiet"),
#   random_effects = c("AnimalID", "end_body_weight"),
#   reference = "Diet,CTR;Time,T0;TimeDiet,T0_CTR",
#   normalization = "NONE",
#   transform = "NONE",
#   analysis_method = "ZINB",
#   cores = 7)
# 
# fit_data_WP3_bact <- Maaslin2(
#   input_data = df_input_data,
#   input_metadata = df_input_metadata,
#   output = "fit_data_WP3_bact_CSS_NEGBIN",
#   fixed_effects = c("Time", "Diet", "TimeDiet"),
#   random_effects = c("AnimalID", "end_body_weight"),
#   reference = "Diet,CTR;Time,T0;TimeDiet,T0_CTR",
#   normalization = "NONE",
#   transform = "NONE",
#   analysis_method = "NEGBIN",
#   cores = 7)

# RUN MAASLIN2 WITH TSS

# build input data: df with ASV
df_input_data <- as.data.frame(otu_table(tax_glom(WP3_TSS_bact, taxrank = "Genus")))
# build input data: change row.names from OTU to phylum_genus
tax_input_data <- as.data.frame(tax_table(tax_glom(WP3_TSS_bact, taxrank = "Genus"))) 
tax_input_data$label <- paste(tax_input_data$Phylum, tax_input_data$Genus, sep = ".")
# remove spaces
tax_input_data$label <- gsub(" ", replacement = ".", x = tax_input_data$label)
# solve duplicates
duplicated <- genefu::rename.duplicate(x=tax_input_data$label, verbose=TRUE)
tax_input_data$label <- duplicated$new.x
# check
row.names(tax_input_data) == row.names(df_input_data)
# rename
row.names(df_input_data) <- tax_input_data$label
#transpose
df_input_data <- t(df_input_data)

# build input metadata:

df_input_metadata <- as(sample_data(WP3_TSS_bact), "data.frame")

# build new variable pasting the diet and time variables, to use and test for interaction between the two;
# maaslin has no proper interaction term like there would be for example in a ANOVA

df_input_metadata$TimeDiet <- paste(df_input_metadata$Time, df_input_metadata$Diet, sep = "_")

# run maaslin on all data, this allows to not loose statistical power

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "fit_data_WP3_bact_TSS_LM",
  fixed_effects = c("Time", "Diet", "TimeDiet"),
  random_effects = c("AnimalID", "end_body_weight"),
  reference = "Diet,CTR;Time,T0;TimeDiet,T0_CTR",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  cores = 7)

# Select only T3 then run with numerical vars

df_input_metadata %>%
  filter(Time == "T3") -> df_input_metadata_T3

df_input_data[row.names(df_input_data) %in% row.names(df_input_metadata_T3),] -> df_input_data_T3

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data_T3,
  input_metadata = df_input_metadata_T3,
  output = "fit_data_WP3_bact_TSS_LM_colonTum_abb",
  fixed_effects = c("Colon_tumor"),
  normalization = "NONE",
  min_abundance = 0.001,
  transform = "NONE",
  analysis_method = "LM",
  cores = 7)

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data_T3,
  input_metadata = df_input_metadata_T3,
  output = "fit_data_WP3_bact_TSS_CPLM_colonTum",
  fixed_effects = c("Colon_tumor"),
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "CPLM",
  cores = 7)


# RUN MAASLIN2 WITH RAREFACTION

# build input data: df with ASV
df_input_data <- as.data.frame(otu_table(tax_glom(WP3_raref_TSS_bact, taxrank = "Genus")))
# build input data: change row.names from OTU to phylum_genus
tax_input_data <- as.data.frame(tax_table(tax_glom(WP3_raref_TSS_bact, taxrank = "Genus"))) 
tax_input_data$label <- paste(tax_input_data$Phylum, tax_input_data$Genus, sep = ".")
# remove spaces
tax_input_data$label <- gsub(" ", replacement = ".", x = tax_input_data$label)
# solve duplicates
duplicated <- genefu::rename.duplicate(x=tax_input_data$label, verbose=TRUE)
tax_input_data$label <- duplicated$new.x
# check
row.names(tax_input_data) == row.names(df_input_data)
# rename
row.names(df_input_data) <- tax_input_data$label
#transpose
df_input_data <- t(df_input_data)

# build input metadata:

df_input_metadata <- as(sample_data(WP3_raref_TSS_bact), "data.frame")

# build new variable pasting the diet and time variables, to use and test for interaction between the two;
# maaslin has no proper interaction term like there would be for example in a ANOVA

df_input_metadata$TimeDiet <- paste(df_input_metadata$Time, df_input_metadata$Diet, sep = "_")

# run maaslin on all data, this allows to not loose statistical power

fit_data_WP3_bact <- Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = "fit_data_WP3_bact_raref_LM",
  fixed_effects = c("Time", "Diet", "TimeDiet"),
  random_effects = c("AnimalID", "end_body_weight"),
  reference = "Diet,CTR;Time,T0;TimeDiet,T0_CTR",
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  cores = 7)


# 
# 
# # analysis have to be done inside each diet, as binary variable is better dealt with.
# 
# # MBD
# 
# df_input_metadata %>%
#   filter(Diet == "MBD") -> df_input_metadata_MBD
# 
# df_input_data[row.names(df_input_data) %in% row.names(df_input_metadata_MBD),] -> df_input_data_MBD
# 
# fit_data_WP3_bact <- Maaslin2(
#   input_data = df_input_data_MBD,
#   input_metadata = df_input_metadata_MBD,
#   output = "fit_data_WP3_bact_MBD",
#   fixed_effects = c("Time"),
#   random_effects = c("AnimalID"),
#   normalization = "TSS",
#   transform = "LOG",
#   analysis_method = "LM",
#   cores = 7)
# 
# 
# #

























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
