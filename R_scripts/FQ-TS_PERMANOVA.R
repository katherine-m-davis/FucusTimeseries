##################################################################################################################
#### PERMANOVA to test community variation by explained by habitat type, temporal, spatial, and host factors #####
##################################################################################################################

# Load packages
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(phyloseq)
library(dplyr)
library(vegan)
library(ggplot2)
library(RVAideMemoire)
library(BiodiversityR)
library(purrr)
library(ecodist)

#setup environment
theme_set(theme_bw())
setwd("~/Desktop/Desktop2020/FQ-TS/Data/")

# Load metadata
fts.meta <- readRDS(file = "~/Desktop/Desktop2020/FQ-TS_new/Data/fts_metadata_v6.rds")
fts.meta <-fts.meta %>% ungroup() %>% mutate_if(is.factor, as.character)

#### Read distance matrices ####
dist.bc.16 <- readRDS("distance_matrix_bc_16s.rds")
dist.bc.18 <-readRDS("distance_matrix_bc_18s.rds")
dist.bc.mg <- readRDS("distance_matrix_bc_metag.rds")

##################################################################################################################
#### Combine all distance matrices for habitat analysis ##########################################################
dist.list <- list(dist.bc.16, dist.bc.18, dist.bc.mg)

# Get rownames
dist.samples <- map(dist.list, row.names)

#### Combine and filter metadata ####
dist.meta <- replicate(3, fts.meta, simplify = FALSE) # repeat metadata as 3 lists
dist.meta.new <- map2(dist.meta, dist.samples, ~ .x %>% 
                         filter(Swab_number %in% .y))  # filter to appropriate samples for 16S, 18S, and shotgun metagenomics

# Make distance matrices ordered by metadata sample number order
dist.list.ordered <- map2(dist.list, dist.meta.new, ~ .x[.y$Swab_number, .y$Swab_number])

#### Run PERMANOVA by Habitat Type ####
permanova.habitat <- map2(dist.list.ordered, dist.meta.new, ~adonis(.x ~  Sample_type, data = .y))
permanova.habitat
capture.output(permanova.habitat, file="~/Desktop/Desktop2020/FQ-TS/Data/FQ-TS_PERMANOVA_habitat_type.txt")

#### Run pairwise PERMANOVA by Habitat Type ####

##################################################################################################################

##################################################################################################################
#### Combine taxa distance matrices for temporal,spatial, and host factor analysis ##############################
dist.list.taxa <- list(dist.bc.16, dist.bc.18)

# Get rownames
dist.taxa.samples <- map(dist.list.taxa, row.names)

#### Combine and filter metadata ####
fts.meta.f <- fts.meta %>% filter(Sample_type == "Fucus") %>% filter(! is.na(Total_dichotomies_whole_plant_binned))
dist.meta.taxa <- replicate(2, fts.meta.f, simplify = FALSE) # repeat metadata as 3 lists
dist.meta.taxa <- map2(dist.meta.taxa, dist.taxa.samples, ~ .x %>% 
                        filter(Swab_number %in% .y))  # filter to appropriate samples for 16S, 18S, and shotgun metagenomics

# Make distance matrices ordered by metadata sample number order
dist.list.taxa.ordered <- map2(dist.list.taxa, dist.meta.taxa, ~ .x[.y$Swab_number, .y$Swab_number])

#### Run PERMANOVA by for time, space, and host factors ####
permanova.tsh <- map2(dist.list.taxa.ordered, dist.meta.taxa, ~adonis(.x ~ Month*quadrat*is.repro*Total_dichotomies_whole_plant_binned, data = .y))
permanova.tsh
capture.output(permanova.tsh, file="~/Desktop/Desktop2020/FQ-TS/Data/FQ-TS_PERMANOVA_month_quadrat_host-factors.txt")

