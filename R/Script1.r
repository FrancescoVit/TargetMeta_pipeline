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

#' ```{create indices variables and similar, include=FALSE}

# here we set some variables that can be used to index the scripts. 
# All variables will be stored in environment, so will be available also
# for the other scripts

workdir <- setwd("/...") # change accordingly to working directory
projname <- "" # set an overall short project name
type <- # set if 16S, ITS, o both (value 1,2,3)
readvar <- "" # set a variable to use for read distribution vizualization

# ideally the folder in which rad data are collected
rawdatadir <- "/..." 

# CONTENT OF RAW DATA DIR

# If csv (i.e. MICCA pipeline)
# '$projname'_16S_metadata.csv
# '$projname'_otutable_16S.txt
# '$projname'_taxa_SILVA_R.csv
# '$projname'_16S_tree_rooted.tree

#' ```

#' ```{r read file, include=FALSE}

# usually, one could start either with txt/csv files or by parsing data from .qza 
# at some point, here an if statement with one or the other way would be good

# moreover, the first if statement would be on the value of type variable. If equal to 
# 1 only upload 16S, if 2 only ITS, if 3 upload both

## BACTERIA
#Uploading data tables
metadata_bact <- read.table(file = "./WP3WP4_metadata.csv", header = T, sep = "\t", row.names = 1)
otutable_bact <- read.table(file = "./otutable_16S.txt", header = T, sep = "", row.names = 1)
taxonomy_bact <- read.table(file = "./taxa_SILVA_R.csv", sep = "\t", row.names = 1)
tree_bact <- read.tree("./WP3WP4_16S_tree_rooted.tree")
colnames(taxonomy_bact) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# mounting phyloseq object
tax_table_bact <- tax_table(as.matrix(taxonomy_bact))
otu_table_bact <- otu_table(as.matrix(otutable_bact), taxa_are_rows = T)
meta_table_bact <- sample_data(metadata_bact)
# mount the phyloseq object
phyloseq_obj_initial_bact <- merge_phyloseq(tax_table_bact,otu_table_bact,meta_table_bact, tree_bact)

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

#' ```

#' **Here we verify correct data import, and plot or test some initial things:**

#' *BACTERIA*

#' ```{r read visualization, include=TRUE, echo = FALSE}

sdt = data.table(as(sample_data(phyloseq_obj_initial_bact), "data.frame"),
                 TotalReads = sample_sums(phyloseq_obj_initial_bact), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

means <- aggregate(TotalReads ~  readvar, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = readvar, y = "TotalReads", add = "jitter", 
                     color = readvar, palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  stat_compare_means(method = "anova", label.y = 3000, label.x = 4.3)
p_reads1

#' ```

#' *FUNGI*
#' ```{r read visualization fungi WP3, include=TRUE, echo = FALSE}
sdt = data.table(as(sample_data(phyloseq_obj_initial_fungi), "data.frame"),
                 TotalReads = sample_sums(phyloseq_obj_initial_fungi), keep.rownames = TRUE)
setnames(sdt, "rn", "Sample")

cat("A total of", sum(sdt$TotalReads), "were obtained")

means <- aggregate(TotalReads ~  readvar, sdt, mean)
p_reads1 <- ggboxplot(data = sdt, x = readvar, y = "TotalReads", add = "jitter", 
                      color = readvar, palette ="npg", rotate = T, dot.size = 3) + 
  theme_cleveland() + 
  stat_compare_means(method = "anova")
p_reads1

#' ```

#' ```{r save workspace, include=FALSE}
#' knitr::opts_chunk$set(echo = TRUE)

save.image("./$projname_workspace")

#' ```

