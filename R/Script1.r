#' ---
#' title: "Script 1 of pipeline: read in data, prepare phyloseq obj, visualize reads, save environment"
#' author: "Francesco Vitali"
#' output: html_document
#' ---

#' ```{r setup library load, include=FALSE}


library(phyloseq)
library(ape)
library(tidyverse)
library(data.table)
library(ggsci)
library(ggpubr)

#' knitr::opts_chunk$set(fig.width=unit(15,"cm"), fig.height=unit(11,"cm"))

#' ```

#' ```{r read file, include=FALSE}
## BACTERIA
#Uploading data tables
metadata_bact <- read.table(file = "./WP3WP4_metadata.csv", header = T, sep = "\t", row.names = 1)
otutable_bact <- read.table(file = "./otutable_16S.txt", header = T, sep = "", row.names = 1)
taxonomy_bact <- read.table(file = "./taxa_SILVA_R.csv", sep = "\t", row.names = 1)
tree_bact <- read.tree("./WP3WP4_16S_tree_rooted.tree")
colnames(taxonomy_bact) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# R converts the "-" in the otu table as ".", change
row.names(metadata_bact) <- gsub("-", ".", row.names(metadata_bact)) 

# mounting phyloseq object
tax_table_bact <- tax_table(as.matrix(taxonomy_bact))
otu_table_bact <- otu_table(as.matrix(otutable_bact), taxa_are_rows = T)
meta_table_bact <- sample_data(metadata_bact)
# mount the phyloseq object
phyloseq_obj_initial_bact <- merge_phyloseq(tax_table_bact,otu_table_bact,meta_table_bact, tree_bact)

# divide WP3 and WP4 objects
WP3_initial_bact <- subset_samples(phyloseq_obj_initial_bact, WP == "WP3") %>% prune_taxa(taxa_sums(.) > 0, .)
WP4_initial_bact <- subset_samples(phyloseq_obj_initial_bact, WP == "WP4") %>% prune_taxa(taxa_sums(.) > 0, .)
# 
# # correct time point name
# sample_data(WP3_initial_bact)$Time <- gsub("T1", "T3",sample_data(WP3_initial_bact)$Time)
# # for WP4 is a bit more diffcult, as the time is by weeks can change in
# # T0 <- Tftm
# # T2 <- T0
# # T14 <- T3
# sample_data(WP4_initial_bact)$Time <- gsub("T0", "Tftm",sample_data(WP4_initial_bact)$Time)
# sample_data(WP4_initial_bact)$Time <- gsub("T2", "T0",sample_data(WP4_initial_bact)$Time)
# sample_data(WP4_initial_bact)$Time <- gsub("T14", "T3",sample_data(WP4_initial_bact)$Time)


## FUNGI
#Uploading data tables
metadata_fungi <- read.table(file = "./WP3WP4_metadata_fungi.csv", header = T, sep = "\t", row.names = 1)
otutable_fungi <- read.table(file = "./otutable_ITS.txt", header = T, sep = "", row.names = 1)
taxonomy_fungi <- read.table(file = "./taxa_R_ITS.csv", sep = "\t", row.names = 1)
tree_fungi <- read.tree("./WP3WP4_ITS_tree.tree")
colnames(taxonomy_fungi) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# R converts the "-" in the otu table as "."
row.names(metadata_fungi) <- gsub("-", ".", row.names(metadata_fungi))
# mounting phyloseq object
tax_table_fungi <- tax_table(as.matrix(taxonomy_fungi))
otu_table_fungi <- otu_table(as.matrix(otutable_fungi), taxa_are_rows = T)
meta_table_fungi <- sample_data(metadata_fungi)
# mount the phyloseq object
phyloseq_obj_initial_fungi <- merge_phyloseq(tax_table_fungi,otu_table_fungi,meta_table_fungi, tree_fungi)

# divide WP3 and WP4 objects
WP3_initial_fung <- subset_samples(phyloseq_obj_initial_fungi, WP == "WP3") %>% prune_taxa(taxa_sums(.) > 0, .)
WP4_initial_fung <- subset_samples(phyloseq_obj_initial_fungi, WP == "WP4") %>% prune_taxa(taxa_sums(.) > 0, .)
# 
# # correct time point name
# sample_data(WP3_initial_fung)$Time <- gsub("T1", "T3",sample_data(WP3_initial_fung)$Time)
# # for WP4 is a bit more diffcult, as the time is by weeks can change in
# # T0 <- Tftm
# # T2 <- T0
# # T14 <- T3
# sample_data(WP4_initial_fung)$Time <- gsub("T0", "Tftm",sample_data(WP4_initial_fung)$Time)
# sample_data(WP4_initial_fung)$Time <- gsub("T2", "T0",sample_data(WP4_initial_fung)$Time)
# sample_data(WP4_initial_fung)$Time <- gsub("T14", "T3",sample_data(WP4_initial_fung)$Time)
# 

#' ```

#' **Here we verify correct data import, and plot or test some initial things:**

#' *BACTERIA*

#' **WP3 BACTERIA**

#' ```{r read visualization bact WP3, include=TRUE, echo = FALSE}
WP3_initial_bact

sdt = data.table(as(sample_data(WP3_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP3_initial_bact), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

cat("Figure shows distribution of read among diets at the two time points. Anova test was used for
overall significance, t-test was used for each diet in comparison to CTR diet")

means <- aggregate(TotalReads ~  Diet, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = "Diet", y = "TotalReads", add = "jitter", 
                     color = "Diet", palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  ylim(c(0,150000)) + 
  stat_compare_means(method = "anova", label.y = 3000, label.x = 4.3) + 
  stat_compare_means(ref.group = "CTR", method = "t.test", label.y = 130000, label = "p.signif") + 
  facet_wrap(.~ Time)
p_reads1

#' ```

#' **WP4 BACTERIA**
 
#' ```{r read visualization bact WP4, include=TRUE, echo = FALSE}
WP4_initial_bact

sdt = data.table(as(sample_data(WP4_initial_bact), "data.frame"),
                 TotalReads = sample_sums(WP4_initial_bact), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

cat("Figure shows distribution of read among diets at the two time points. Anova test was used for
overall significance, t-test was used for each diet in comparison to CTR diet")

means <- aggregate(TotalReads ~  Diet, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = "Inoculum_diet", y = "TotalReads", add = "jitter", 
                     color = "Inoculum_diet", palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  ylim(c(0,70000)) + 
  stat_compare_means(method = "anova", label.y = 3000, label.x = 4.3) + 
  facet_wrap(.~ Time)
p_reads1

#' ```

#' *FUNGI*

#' **WP3 FUNGI**

#' ```{r read visualization fungi WP3, include=TRUE, echo = FALSE}
WP3_initial_fung

sdt = data.table(as(sample_data(WP3_initial_fung), "data.frame"),
                 TotalReads = sample_sums(WP3_initial_fung), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

cat("Figure shows distribution of read among diets at the two time points. Anova test was used for
overall significance, t-test was used for each diet in comparison to CTR diet")

means <- aggregate(TotalReads ~  Diet, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = "Diet", y = "TotalReads", add = "jitter", 
                      color = "Diet", palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  ylim(c(0,150000)) + 
  stat_compare_means(method = "anova", label.y = 3000, label.x = 4.3) + 
  stat_compare_means(ref.group = "CTR", method = "t.test", label.y = 130000, label = "p.signif") + 
  facet_wrap(.~ Time)
p_reads1

#' ```

#' **WP4 FUNGI**

#' ```{r read visualization fungi WP4, include=TRUE, echo = FALSE}
WP4_initial_fung

sdt = data.table(as(sample_data(WP4_initial_fung), "data.frame"),
                 TotalReads = sample_sums(WP4_initial_fung), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

cat("Figure shows distribution of read among diets at the two time points. Anova test was used for
overall significance, t-test was used for each diet in comparison to CTR diet")

means <- aggregate(TotalReads ~  Diet, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = "Inoculum_diet", y = "TotalReads", add = "jitter", 
                      color = "Inoculum_diet", palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  ylim(c(0,110000)) + 
  stat_compare_means(method = "anova", label.y = 3000, label.x = 4.3) + 
  facet_wrap(.~ Time)
p_reads1

#' ```

#' ```{r save workspace, include=FALSE}
#' knitr::opts_chunk$set(echo = TRUE)

save.image("./WP3WP4_workspace")

#' ```

