#' ---
#' title: "Script 0 of pipeline: check if needed packages are installed, else install them"
#' author: "Francesco Vitali"
#' output: html_document
#' ---

#' ```{r setup library load, include=FALSE}

#' For a relly fresh system, some dependencies would be on the linux install

sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev cmake

#' here a list of all the packages that are needed by other scripts
#' When opened in RStudio, this script should prompt a message abut pks
#' required but not installed, and ask if should install.
#' On a fresh system, or on a R update, o is situations with missing pks, this
#' be enough for setting up the environment

library(ape)
library(data.table)
library(ggsci)
library(tidyverse)
library(vegan)
library(patchwork)
library(psych)
library(corrplot)
library(ggpubr)
library(reshape2)
library(rstatix)
library(BiocManager)

install.packages(c("ape","data.table","ggsci","tidyverse","vegan","patchwork","psych","corrplot","ggpubr","reshape2","rstatix"))

#' Actions are still required for installing bioconductor pks or pks compiled from
#' github

BiocManager::install("Biostrings")
BiocManager::install("phyloseq")
BiocManager::install("microbiome")
BiocManager::install("metagenomeSeq")

library(phyloseq)
library(microbiome)
library(metagenomeSeq)

#' ```
