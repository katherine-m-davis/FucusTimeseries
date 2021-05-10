############################################################################################################# 
#### INDVAL by habitat type (F. distichus, Rock substrate, Seawater) at the genus- or family-level ##########
# Indicator species analyses from the indicspecies package. 
# https://cran.r-project.org/web/packages/indicspecies/vignettes/indicspeciesTutorial.pdf
############################################################################################################# 

# Indicator Statistics
# ‘A’ is sample estimate of the probability that the surveyed site belongs to the target site 
# group given the fact that the species has been found. This conditional probability 
# is called the specificity or positive predictive value of the species as indicator of the site group. 
# ‘B’ is sample estimate of the probability of finding the species
# in sites belonging to the site group. This second conditional probability is
# called the fidelity or sensitivity of the species as indicator of the target site group.

# Load packages
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(phyloseq)
library(dplyr)
library(vegan)
library(ggplot2)
library(indicspecies)
library(sjmisc)

# Set up environment
theme_set(theme_bw())
setwd("~/Desktop/Desktop2020/FQ-TS/Data/")
############################################################################################################# 


#### 16S #####################################################################################################
# Load data
ps.16.df <- readRDS(file = '~/Desktop/Desktop2020/FQ-TS/Data/fts.df.16.filtered.rarefied_10K.May2021.RDS')

# Select and summarise by sample number, Genus, Abundance
ps.16.g <- ps.16.df %>% group_by(Sample, Genus) %>% summarise(g_abund = sum(Abundance))

#### Make table where samples are rows and genera are columns ####
ps.16.g.w <- ps.16.g %>% spread(Genus, g_abund) %>% replace(is.na(.), 0)
fts.16.samples <- ps.16.g.w$Sample
rownames(ps.16.g.w) <- fts.16.samples
ps.16.g.w <- ps.16.g.w[,-1]  
ps.16.mat <- as.matrix(ps.16.g.w)

#### 16S Indicator analysis by sample type (habitat) ####
ps.16.meta <- ps.16.df %>% select(Sample, Year:QC.ed) %>% distinct()
rownames(ps.16.meta) <- ps.16.meta$Sample
ps.16.meta <-ps.16.meta[ps.16.g.w$Sample,]

sample_type <- as.character(ps.16.meta$Sample_type) # vector of months matching sample order in otu table
indval.16 <- multipatt(ps.16.mat, sample_type, duleg = TRUE, control = how(nperm=999))
capture.output(summary(indval.16, indvalcomp = T), file="indval_16S_genus_sample_type_fts.txt")
summary(indval.16, indvalcomp = T) #output results
############################################################################################################# 

#### 18S ####################################################################################################
# Load data
ps.18.df <- readRDS(file="~/Desktop/Desktop2020/FQ-TS/Data/fts.df.18.filtered.rarefied_750.Apr2021.RDS")

# Select and summarise by sample number, Genus, Abundance
ps.18.g <- ps.18.df %>% group_by(Sample, Genus) %>% summarise(g_abund = sum(Abundance))

# Make table where samples are rows and genera are columns 
ps.18.g.w <- ps.18.g %>% spread(Genus, g_abund) %>% replace(is.na(.), 0)
fts.18.samples <- ps.18.g.w$Sample
rownames(ps.18.g.w) <- fts.18.samples
ps.18.g.w <- ps.18.g.w[,-1]  
ps.18.mat <- as.matrix(ps.18.g.w)

#### 18S Indicator analysis by sample type (habitat) by species ####
ps.18.meta <- ps.18.df %>% select(Sample, Year:QC.ed) %>% distinct()
rownames(ps.18.meta) <- ps.18.meta$Sample
ps.18.meta <-ps.18.meta[fts.18.samples,]

sample_type <- as.character(ps.18.meta$Sample_type) # vector of months matching sample order in otu table
indval.18 <- multipatt(ps.18.mat, sample_type, duleg = TRUE, control = how(nperm=999))
capture.output(summary(indval.18, indvalcomp = T), file="indval_18S_genus_sample_type_fts.txt")
summary(indval.18, indvalcomp = T) #output results


#### Analysis for 18S by genus####
ps.18.fam <- ps.18.df %>% ungroup() %>% group_by(Sample, Family) %>% summarise(f_abund = sum(Abundance))

ps.18.f.w <- ps.18.fam %>% spread(Family, f_abund) %>% replace(is.na(.), 0)
fts.18.samples <- ps.18.f.w$Sample
rownames(ps.18.f.w) <- fts.18.samples
ps.18.f.w <- ps.18.f.w[,-1]  
ps.18.fam.mat <- as.matrix(ps.18.f.w)

ps.18.meta <- ps.18.df %>% select(Sample, Year:QC.ed) %>% distinct()
rownames(ps.18.meta) <- ps.18.meta$Sample
ps.18.meta <-ps.18.meta[fts.18.samples,]

sample_type <- as.character(ps.18.meta$Sample_type) # vector of months matching sample order in otu table
indval.18.fam <- multipatt(ps.18.fam.mat, sample_type, duleg = TRUE, control = how(nperm=999))
capture.output(summary(indval.18.fam, indvalcomp = T), file="indval_18S_family_sample_type_fts.txt")
summary(indval.18.fam, indvalcomp = T) #output results

