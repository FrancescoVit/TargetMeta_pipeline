#' ---
#' title: "Script 2 of pipeline: read in workspace from script1, run alpha diversity analysis"
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
library(vegan)
library(patchwork)
library(psych)
library(corrplot)

# red in workspace from script 1

load ("./WP3WP4_workspace")

palette_custom <- c("#1F78B4", "#E31A1C", "#FF7F00", "#33A02C") 

#' knitr::opts_chunk$set(fig.width=unit(15,"cm"), fig.height=unit(11,"cm"))

#' ```

#' ```

#' **Here we calculate alpha diversity measure, make some comparisons, and plot**

#' *BACTERIA*

#############

#' **WP3 BACTERIA**

#############
#' ```{r alphadiv bact WP3, include=TRUE, echo = FALSE}

# how to transform data for alpha diversity analysis? 

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

# see rarefaction curves
rarecurve(t(otu_table(WP3_initial_bact)), step=50, cex=0.5)

# most of the sample are in the 40k-70k reads, the min sample is fairly different from the other
# create rarefied phyloseq object, remove sample "E149.40.T0", both time point as is lower outlier for total reads
WP3_raref_bact<- rarefy_even_depth(subset_samples(WP3_initial_bact, AnimalID != "E149-40"),
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP3_initial_bact, AnimalID != "E149-40")))),
                                   replace=F)

# now, following P. Schloss (https://www.biorxiv.org/content/10.1101/2020.12.11.422279v1.abstract) alpha
# diversity would be impacted by removing of rare ASVs to a certain treshold, and to minimize different 
# library size effect, he perform rarefaction. Sound fair conceptually, at least for alpha diversity

# calculate all indices. 
alphadiv_WP3_bact <- microbiome::alpha(WP3_raref_bact)
df <- as(sample_data(WP3_raref_bact), "data.frame")
df <- cbind(df, alphadiv_WP3_bact)

df1 <- as(sample_data(WP3_raref_bact), "data.frame")
df1 <- cbind(df1[,c(3,4,10,11,12,13)], alphadiv_WP3_bact)

write.csv2(x =  df1, file = "alphadiv_WP3.csv")


df1 %>%
  filter(Time == "T3", Diet == "PVD" | Diet == "CTR")  %>%
  select(-c(Time, Diet, start_body_weight)) %>%
  corr.test(method = "pearson", adjust = "BH") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r), type="lower", order="hclust", 
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)


ggboxplot(data = df, x = "Time", y = "observed", add = "jitter", 
          fill = "Time", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test", comparisons = list(c("T0","T3")), paired = T) + 
  facet_grid(.~Diet) +
  ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness") + guides(fill="none") #+ ylim(c(100,1000))

ggboxplot(data = df, x = "Time", y = "chao1", add = "jitter",
          fill = "Time", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test", comparisons = list(c("T0","T3")), paired = T) + 
  facet_grid(.~Diet) +
  ggtitle(label = "Bacteria", subtitle =  "Richness (ChaoI index)") +
  ylab("Richness") + guides(fill="none") #+ ylim(c(100,1000))

ggboxplot(data = df, x = "Time", y = "evenness_pielou", add = "jitter", 
          fill = "Time", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test", comparisons = list(c("T0","T3")), paired = T) + 
  facet_grid(.~Diet) +
  ggtitle(label = "Bacteria", subtitle =  "Pielou's Evenness") +
  ylab("Evenness") + guides(fill="none") #+ ylim(c(100,1000))

ggboxplot(data = df, x = "Time", y = "diversity_shannon", add = "jitter", 
          fill = "Time", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test", comparisons = list(c("T0","T3")), paired = T) + 
  facet_grid(.~Diet) +
  ggtitle(label = "Bacteria", subtitle =  "Shannon's index") +
  ylab("Shannon") + guides(fill="none") #+ ylim(c(100,1000))

## comparing at T3 respect to diet

df$Diet <- factor(df$Diet, levels = c("CTR", "MBD", "MBDT", "PVD"))

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "observed", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,1150) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) -> p_rich_2

df %>%
  filter(Diet == "PVD") %>%
  ggplot(aes(x = Time, y = observed)) +
  geom_boxplot(outlier.shape = NA, col = "grey50")+
  geom_jitter(width = 0.02, size = 2.5) +
  coord_flip() +
  theme_bw() +
  stat_compare_means(method = "t.test",comparisons = list(c("T0", "T3"))) +
  xlab("PVD subset") +
  ylab("") + theme(
    panel.background = element_rect(fill = "grey80", colour = "lightblue",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()) -> p_rich_1

p_rich_2 + inset_element(p_rich_1, left = 0.05, bottom = 0.1, right = 0.95, top = 0.05) -> p_rich_final

ggsave('./PVD_bact_rich_WP3.png', p_rich_final)

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "observed", add = "jitter", 
            fill = "Diet", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test",ref.group = "CTR") +
  ggtitle(label = "Bacteria", subtitle =  "Richness") +
  ylab("Richness") + guides(fill="none") 


df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "chao1", add = "jitter", 
            fill = "Diet", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test",ref.group = "CTR") +
  ggtitle(label = "Bacteria", subtitle =  "Richness (ChaoI index)") +
  ylab("Richness") + guides(fill="none") 

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "evenness_pielou", add = "jitter", 
            fill = "Diet", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test",ref.group = "CTR") +
  ggtitle(label = "Bacteria", subtitle =  "Pielou's Evenness") +
  ylab("Evenness") + guides(fill="none") 

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "diversity_shannon", add = "jitter", 
            fill = "Diet", palette ="grey", rotate = F, dot.size = 4) + theme_bw() + 
  stat_compare_means(method = "t.test",ref.group = "CTR") +
  ggtitle(label = "Bacteria", subtitle =  "Shannon's index") +
  ylab("Shannon") + guides(fill="none")

# inspect correlation between diet and alphadiv

df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "observed", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(12, 960), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  ylab("Richness (n° of observed ASVs)") +
  xlab("Tumours in the Colon") +
  stat_cor(aes(color = Diet), label.x = 12, method = "pearson") + 
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))  -> p_correl_rich

df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "diversity_shannon", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(12, 5), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  ylab("Shannon's Index") +
  xlab("Tumours in the Colon") +
  stat_cor(aes(color = Diet), label.x = 12, method = "pearson") + 
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) -> p_correl_shan


p_correl_rich/p_correl_shan + 
  plot_layout(guides = "collect") & theme(legend.position = 'top') -> panel_correlation

p_rich_final + panel_correlation -> figure1

ggsave('./Figure1.png', figure1)

df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "observed", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(15, 960), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  ylab("Richness (n° of observed ASVs)") +
  xlab("Colon tumours") +
  stat_cor(aes(color = Diet), label.x = 15, method = "pearson") + 
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) -> p_correl_rich_2

df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "diversity_shannon", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(15, 5), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  ylab("Shannon's Index") +
  xlab("Colon tumours") +
  stat_cor(aes(color = Diet), label.x = 15, method = "pearson") + 
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) -> p_correl_shan_2


p_correl_rich_2/p_correl_shan_2 +
  plot_layout(guides = "collect") & theme(legend.position = 'top') -> panel_correlation_2


df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "dominance_simpson", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(15, 0.70), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  #ylab("Shannon's Index") +
  xlab("Colon tumours") +
  stat_cor(aes(color = Diet), method = "pearson") + 
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"))

df %>%
  filter(Time == "T3") %>%
  ggscatter(y = "rarity_log_modulo_skewness", x = "SI_tumour", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(8, 2.061), cor.coef.size = 4,
            add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  ylim(c(2.061, 2.062)) +
  #ylab("Shannon's Index") +
  xlab("SI tumour") +
  stat_cor(aes(color = Diet), method = "pearson") + 
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"))


df1 %>%
  filter(Time == "T3", Diet == "PVD" | Diet == "CTR") %>% 
  ggscatter(y = "dominance_simpson", x = "Colon_tumor", shape = 21, fill = "Diet",col = "Diet", size = 3, cor.coef.coord = c(15, 0.7), cor.coef.size = 4,
                                                                 add = "reg.line", cor.coef = T, cor.method = "pearson") +  
  scale_fill_manual(values = palette_custom) +
  scale_color_manual(values = palette_custom) +
  #ylab("Shannon's Index") +
  xlab("Colon tumours") +
  stat_cor(aes(color = Diet), method = "pearson") + 
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"))




#############

#' **WP4 BACTERIA**

###########

#' ```{r alphadiv bact WP4, include=TRUE, echo = FALSE}

# how to transform data for alpha diversity analysis? 

# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP4_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP4_initial_bact), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)
# variation in read depth is 8 

# plot reads to find outliers; WP4-928 is lower outlier

sdt %>%
  filter(Type == "Feces", AnimalID != "WP4-928") %>%
  ggboxplot(.,y = "TotalReads", label = "AnimalID")

# see rarefaction curves
rarecurve(t(otu_table(WP4_initial_bact)), step=50, cex=0.5)

# most of the sample are in the 40k-70k reads, the min sample is fairly different from the other
# create rarefied phyloseq object, remove sample "E149.40.T0", both time point as is lower outlier for total reads
WP4_raref_bact<- rarefy_even_depth(subset_samples(WP4_initial_bact, AnimalID != "WP4-928"),
                                   rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP4_initial_bact, AnimalID != "WP4-928")))),
                                   replace=F)

# now, following P. Schloss (https://www.biorxiv.org/content/10.1101/2020.12.11.422279v1.abstract) alpha
# diversity would be impacted by removing of rare ASVs to a certain treshold, and to minimize different 
# library size effect, he perform rarefaction. Sound fair conceptually, at least for alpha diversity

# calculate all indices. 
alphadiv_WP4_bact <- microbiome::alpha(WP4_raref_bact)
df <- as(sample_data(WP4_raref_bact), "data.frame")
df <- cbind(df, alphadiv_WP4_bact)


df$Time <- factor(df$Time, levels = c("Inoculum", "Tftm", "T0", "T3"))

df$Inoculum_diet <- factor(df$Inoculum_diet, levels = c("CTR", "MBD", "MBDT", "PVD"))

df1 <- as(sample_data(WP4_raref_bact), "data.frame")
df1 <- cbind(df1[,c(3,5,6,15)], alphadiv_WP4_bact)

write.csv2(x =  df1, file = "alphadiv_WP4.csv")

# calculate correlation matrix

df1 %>%
  filter(Time == "T3")  %>%
  select(-c(Time, Inoculum_diet,Type)) %>%
  corr.test(method = "pearson", adjust = "BH") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r), type="lower", order="hclust", 
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)

# In WP4 no correlations between MDF and diversity, maybe MDF are too early 

# So, exploration here is more in the sense of time and diet of inoculum


df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Inoculum_diet", y = "observed", add = "jitter", 
            fill = "Inoculum_diet", palette ="grey", rotate = F, dot.size = 6) + 
  theme_bw() +
  #geom_hline(data = df[df$Type == "Inoculum",], aes(yintercept = observed), linetype = 2) +
  facet_wrap(. ~ Time) + 
  stat_compare_means(method = "t.test",comparisons = list(c("Tftm", "T0"),c("T0","T3"), c("Tftm","T3"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") + 
  guides(fill="none") + 
  xlab("") +
  ylim(350,1250) +
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) #-> p_observed_WP4

df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Time", y = "observed", add = "jitter", 
            fill = "Time", palette ="grey", rotate = F, dot.size = 6) + 
  theme_bw() +
  geom_hline(data = df[df$Type == "Inoculum",], aes(yintercept = observed), linetype = 2) +
  facet_wrap(. ~ Inoculum_diet) + 
  stat_compare_means(method = "t.test",comparisons = list(c("Tftm", "T0"),c("T0","T3"), c("Tftm","T3"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") + 
  guides(fill="none") + 
  xlab("") +
  ylim(350,1250) +
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) #-> p_observed_WP4

df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Time", y = "evenness_pielou", add = "jitter", 
            fill = "Time", palette ="grey", rotate = F, dot.size = 6) + 
  theme_bw() +
  geom_hline(data = df[df$Type == "Inoculum",], aes(yintercept = evenness_pielou), linetype = 2) +
  facet_wrap(. ~ Inoculum_diet) + 
  stat_compare_means(method = "t.test",comparisons = list(c("Tftm", "T0"),c("T0","T3"), c("Tftm","T3"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Pielou's Evenness") + 
  guides(fill="none") + 
  xlab("") +
  ylim(0.55,0.75) +
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) #-> p_eve_WP4

df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Time", y = "diversity_shannon", add = "jitter", 
            fill = "Time", palette ="grey", rotate = F, dot.size = 6) + 
  theme_bw() +
  geom_hline(data = df[df$Type == "Inoculum",], aes(yintercept = diversity_shannon), linetype = 2) +
  facet_wrap(. ~ Inoculum_diet) + 
  stat_compare_means(method = "t.test",comparisons = list(c("Tftm", "T0"),c("T0","T3"), c("Tftm","T3"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Shannon's Index") + 
  guides(fill="none") + 
  xlab("Diet of donor rat")  +
  ylim(3,6) +
  theme(axis.text.x = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 14,  hjust = .5, vjust = .5, face = "plain")) #-> p_eve_WP4


### compare overall the diets of donor

df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Inoculum_diet", y = "observed", shape = "Time", 
            fill = "Inoculum_diet", palette ="grey", rotate = F, dot.size = 4) + 
  theme_bw() +
  #stat_compare_means(method = "t.test",comparisons = list(c("CTR", "MBD"),c("CTR", "MBDT"), c("CTR", "PVD")), label = "p.signif") +
  #stat_compare_means(method = "t.test",comparisons = list(c("MBD", "MBDT"),c("MBD", "PVD"), c("MBDT", "PVD"))) +
    #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") + 
  guides(fill="none", shape = "none") + 
  xlab("Diet of donor rat") +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) -> p_WP4_1


df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Inoculum_diet", y = "evenness_pielou", shape = "Time", 
            fill = "Inoculum_diet", palette ="grey", rotate = F, dot.size = 4) + 
  theme_bw() +
  #stat_compare_means(method = "t.test",comparisons = list(c("CTR", "MBD"),c("CTR", "MBDT"), c("CTR", "PVD")), label = "p.signif") +
  #stat_compare_means(method = "t.test",comparisons = list(c("MBD", "MBDT"),c("MBD", "PVD"), c("MBDT", "PVD"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Pielou's Evenness index") + 
  guides(fill="none", shape = "none") + 
  xlab("Diet of donor rat") +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))  -> p_WP4_2

df %>%
  filter(Time != "Inoculum") %>%
  ggboxplot(data = ., x = "Inoculum_diet", y = "diversity_shannon", shape = "Time", 
            fill = "Inoculum_diet", palette ="grey", rotate = F, dot.size = 4) + 
  theme_bw() +
  #stat_compare_means(method = "t.test",comparisons = list(c("CTR", "MBD"),c("CTR", "MBDT"), c("CTR", "PVD")), label = "p.signif") +
  #stat_compare_means(method = "t.test",comparisons = list(c("MBD", "MBDT"),c("MBD", "PVD"), c("MBDT", "PVD"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Shannon's Diversity index") + 
  guides(fill="none", shape = "none") + 
  xlab("Diet of donor rat")  +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))  -> p_WP4_3

p_WP4_1 + p_WP4_2 + p_WP4_3 -> p_alpha_WP4

# plot with WP3

(p_rich_final + panel_correlation) / p_alpha_WP4 + plot_layout(heights = c(2, 1)) -> figure1

ggsave('./Figure1.png', figure1)


# Use ordination analysis to evaluate if overall time or diet has the highest contribution to the alpha diversity

df %>%
  filter(Time != "Inoculum") -> df_ord
  

df %>%
  filter(Time != "Inoculum") %>%
  select(-c(1,2,4:15)) %>%
  mutate(across(Time, as.numeric)) %>%
  scale(.) %>%
  PCA(na.omit(.), graph = FALSE) %>%
  fviz_pca_biplot(.,pointsize = 4, invisible = "quali", habillage = df_ord$Inoculum_diet, 
                  col.ind = df_ord$Inoculum_diet,
                  col.var = "gray", alpha.var =0.8, repel = T,
                  label = "var") +
  geom_text_repel(label = df_ord$Time, nudge_x = 0.2, nudge_y = 0.2)

#####

#' *FUNGI*

#' **WP3 FUNGI**

#' ```{r alphadiv fungi WP3, include=TRUE, echo = FALSE}

# how to transform data for alpha diversity analysis? 

# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP3_initial_fung), "data.frame"),
                 TotalReads = sample_sums(WP3_initial_fung), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)
# variation in read depth is 4 

# plot reads to find outliers

sdt %>%
  filter(AnimalID != "E149-60" & AnimalID != "E149-15") %>%
  #filter(AnimalID != "E149-15") %>%
  ggboxplot(.,y = "TotalReads", label = "AnimalID")

# see rarefaction curves
rarecurve(t(otu_table(WP3_initial_fung)), step=50, cex=0.5)
# most of the sample are in the 40k-70k reads, the min sample is fairly different from the other
# create rarefied phyloseq object, remove sample "E149.40.T0", both time point 
WP3_raref_fung<- rarefy_even_depth(subset_samples(WP3_initial_fung, AnimalID != "E149-60" & AnimalID != "E149-15"), rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(subset_samples(WP3_initial_fung, AnimalID != "E149-60" & AnimalID != "E149-15")))),
                                   replace=F)

# now, following P. Schloss (https://www.biorxiv.org/content/10.1101/2020.12.11.422279v1.abstract) alpha
# diversity would be impacted by removing of rare ASVs to a certain treshold, and to minimize different 
# library size effect, he perform rarefaction. Sound fair conceptually, at least for alpha diversity

# calculate all indices. 
alphadiv_WP3_fung <- microbiome::alpha(WP3_raref_fung)
df <- as(sample_data(WP3_raref_fung), "data.frame")
df <- cbind(df, alphadiv_WP3_fung)

df1 <- as(sample_data(subset_samples(WP3_initial_fung, AnimalID != "E149-60" & AnimalID != "E149-15")), "data.frame")
df1 <- cbind(df1[,c(3,4,10,11,12,13)], alphadiv_WP3_fung)


df1 %>%
  filter(Time == "T3")  %>%
  select(-c(Time, Diet, start_body_weight)) %>%
  corr.test(method = "pearson", adjust = "BH") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r), type="lower", order="hclust", 
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)


df$Diet <- factor(df$Diet, levels = c("CTR", "MBD", "MBDT", "PVD"))

## no correlation found between fungi and tumour counts, so exploration to see difference in diet

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "observed", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,70) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "observed", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Richness (n° of observed ASVs)") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,70) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


###

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "diversity_shannon", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Shannon's Index") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,3.5) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "diversity_shannon", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Shannon's Index") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,3.5) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))

###

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "evenness_pielou", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Pielou's Evenness Index") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,1) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "evenness_pielou", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Pielou's Evenness Index") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,1) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))

#' **WP4 FUNGI**

#' ```{r alphadiv fungi WP4, include=TRUE, echo = FALSE}

# how to transform data for alpha diversity analysis? 

# calculate the variation in read obtained per sample
sdt = data.table(as(sample_data(WP4_initial_fung), "data.frame"),
                 TotalReads = sample_sums(WP4_initial_fung), keep.rownames = TRUE)
summary(sdt$TotalReads)
max(sdt$TotalReads)/min(sdt$TotalReads)
# variation in read depth is 5 

# plot reads to find outliers

sdt %>%
  #filter(AnimalID != "E149-60" & AnimalID != "E149-15") %>%
  #filter(AnimalID != "E149-15") %>%
  ggboxplot(.,y = "TotalReads", label = "AnimalID")

# see rarefaction curves
rarecurve(t(otu_table(WP3_initial_fung)), step=50, cex=0.5)

# most of the sample are in the 40k-70k reads, the min sample is fairly different from the other
# create rarefied phyloseq object, remove sample "E149.40.T0", both time point 
WP4_raref_fung<- rarefy_even_depth(WP4_initial_fung, rngseed=1234,
                                   sample.size=round(0.99*min(sample_sums(WP4_initial_fung))),
                                   replace=F)

# now, following P. Schloss (https://www.biorxiv.org/content/10.1101/2020.12.11.422279v1.abstract) alpha
# diversity would be impacted by removing of rare ASVs to a certain treshold, and to minimize different 
# library size effect, he perform rarefaction. Sound fair conceptually, at least for alpha diversity

# calculate all indices. 
alphadiv_WP3_fung <- microbiome::alpha(WP3_raref_fung)
df <- as(sample_data(WP3_raref_fung), "data.frame")
df <- cbind(df, alphadiv_WP3_fung)

df1 <- as(sample_data(subset_samples(WP3_initial_fung, AnimalID != "E149-60" & AnimalID != "E149-15")), "data.frame")
df1 <- cbind(df1[,c(3,4,10,11,12,13)], alphadiv_WP3_fung)


df1 %>%
  filter(Time == "T3")  %>%
  select(-c(Time, Diet, start_body_weight)) %>%
  corr.test(method = "pearson", adjust = "BH") ->  cor_test_mat  # Apply corr.test function

corrplot(as.matrix(cor_test_mat$r), type="lower", order="hclust", 
         p.mat = as.matrix(cor_test_mat$p), sig.level = 0.05)


df$Diet <- factor(df$Diet, levels = c("CTR", "MBD", "MBDT", "PVD"))

## no correlation found between fungi and tumour counts, so exploration to see difference in diet

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "observed", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Richness (n° of observed ASVs)") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,70) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "observed", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Richness (n° of observed ASVs)") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,70) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


###

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "diversity_shannon", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Shannon's Index") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,3.5) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "diversity_shannon", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Shannon's Index") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,3.5) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))

###

df %>%
  filter(Time == "T3") %>%
  ggboxplot(data = ., x = "Diet", y = "evenness_pielou", add = "jitter", 
            fill = "Diet", rotate = F, dot.size = 6) + theme_bw() + 
  stat_compare_means(method = "anova") +
  #stat_compare_means(method = "t.test",comparisons = list(c("PVD", "CTR"), c("PVD","MBD"), c("PVD", "MBDT"),
  #                                                        c("MBDT", "MBD"), c("MBDT","CTR"),
  #                                                        c("MBD", "CTR"))) +
  #ggtitle(label = "Bacteria", subtitle =  "Richness (n° of ASVs)") +
  ylab("Pielou's Evenness Index") +
  scale_fill_manual(values = palette_custom) +
  guides(fill="none") + 
  xlab("") +
  ylim(0,1) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain")) 


df %>%
  ggboxplot(data = ., x = "Time", y = "evenness_pielou", add = "jitter", 
            fill = "Time", rotate = F, dot.size = 6) + theme_bw() + 
  facet_wrap(. ~ Diet) +
  stat_compare_means(method = "t.test") +
  ylab("Pielou's Evenness Index") +
  scale_fill_grey() +
  guides(fill="none") + 
  xlab("") +
  ylim(0,1) +
  theme(axis.text.x = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"),
        axis.title.y = element_text(size = 12,  hjust = .5, vjust = .5, face = "plain"))



#' ```

#' ```{r write images, include=FALSE}
#' knitr::opts_chunk$set(echo = TRUE)

# use ggsave for images

#' ```

