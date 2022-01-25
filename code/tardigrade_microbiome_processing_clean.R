# Laura Tibbs-Cortes
# Sept 7, 2021

# Tardigrade microbiome analysis.

# Using CornCob, DivNet, phyloseq packages

# NOTE: originally run in R 4.0.3

# # First time only, install packages: 

# # For phyloseq and decontam, use:
# if (!requireNamespace("BiocManager"))
# install.packages("BiocManager")
# BiocManager::install("phyloseq")
# BiocManager::install("decontam")

# For all other required packages, use install.packages("packagename")

# load packages
library(phyloseq)
# library(ggplot2)
library(vegan)
# library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(data.table)
library("xlsx")
library(stringr)
library(tibble)
library(splitstackshape)
library(viridis)
library(remotes)
library(DescTools)

library(tidyverse)
library(decontam)
library(gridExtra)
library(writexl)

set.seed(734095) # seed from random number generator for reproducibility

# FULL DATA -----------------------------------------------------------
# working with all data together here, from running Mothur on all sequences at once.

# NOTE: need to first download and unzip the taxonomy file (mothur_output.tar.gz) 
# found at https://github.com/LTibbs/tardigrade_microbiome/blob/main/taxonomy/mothur_output.tar.gz
# After downloading, unzip the file by typing "tar -xzvf mothur_output.tar.gz" in a terminal

# creating a taxonomy object for later use ----------------------------

#reading in the cons.taxonomy file produced by mothur:
tax <- read.table("../taxonomy/stability.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.01.cons.taxonomy", header = TRUE)
#making a list of characters representing the bootstrap values produced in the cons.taxonomy file
bootstraps <- c("(100)", "(99)", "(98)", "(97)", "(96)", "(95)", "(94)", "(93)", "(92)", "(91)", "(90)", "(89)", "(88)", "(87)", "(86)", "(85)", "(84)", "(85)", "(84)", "(83)", "(82)", "(81)", "(80)", "(79)", "(78)", "(77)", "(76)", "(75)", "(74)", "(73)", "(72)", "(71)", "(70)", "(69)", "(68)", "(67)", "(66)", "(65)", "(64)", "(63)", "(62)", "(61)", "(60)", "(59)", "(58)", "(57)", "(56)", "(55)", "(54)", "(53)", "(52)", "(51)")
#creating a function that will allow a grep command to escape and important characters.
#we will want to exscape the "()" in the bootstrap values so they will be considered a character and deleted.
regex.escape <- function(string) {
  gsub("([][{}()+*^${|\\\\?])", "\\\\\\1", string)
}
#making a new list including the escapes and a pipe character between each value
bootstrap_input <- paste0("\\b", paste0(regex.escape(bootstraps), collapse="|"), "\\b")
#using gsub to find and replace all values in the input and replace them with nothing
tax$Taxonomy <- gsub(bootstrap_input, "", tax$Taxonomy)
#spliting the file apart into columns by the semicolon
tax <- cSplit(data.frame(tax), "Taxonomy", sep=";", drop=TRUE)
#renaming the columns
names(tax)[3:8] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

##### importing mothur output into a physeq object  -----------------

data <- import_mothur(mothur_shared_file = "../taxonomy/complete.stability.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.shared", 
                      mothur_constaxonomy_file = "../taxonomy/stability.trim.contigs.good.unique.good.filter.precluster.pick.opti_mcc.0.01.cons.taxonomy")

# read in key between group names and map names:
key <- read.xlsx("../taxonomy/mothur_group_vs_sample_code.xlsx",1)
stopifnot(all.equal(key$mothur_ID, sample_names(data)))
sample_names(data) <- key$map_ID

#read your meta data in as well
map <- read.csv("../design/Map_final.csv")
save.rep <- map$rep.ID
# add columns to the map:
map_df <- map %>%
  mutate(Sample.ID=str_sub(rep.ID, 1, -2)) %>%
  separate(col = rep.ID, sep = "_", into=c("Location", "Year", "Tree", "Subsample")) %>%
  mutate(Year=paste0(20,Year))
map_df$rep.ID <- save.rep

#convert the metadata into a phyloseq object
map <- sample_data(map_df)

#naming the rows to correspond with the SampleID column within the metadatafile.
rownames(map) <- map$old.sample.ID

#merging the shared-taxonomy phyloseq object with the metadata file (map) and viewing it
data_merge <- merge_phyloseq(data, map)

# NOTE: removing the A1 and A2 samples here; they aren't part of this experiment
data_merge <- prune_samples(!(data_merge@sam_data$Subsample %in% c("liA1", "liA2")), data_merge)

data_merge

#assigning taxonomic rankings
colnames(tax_table(data_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")

# # output the metadata and taxonomy:
# fwrite(file = "../design/metadata.csv", x = data_merge@sam_data)
# write.csv(data_merge@tax_table, "taxonomy.csv", row.names = T)


# make logical sample order for figures etc
sample.order <- c("L1_19_Tr1_li", "L1_19_Tr2_li", "L1_19_Tr3_li", 
                  "L1_20_Tr1_li", "L1_20_Tr2_li", "L1_20_Tr3_li",  "L1_20_Tr4_li",
                  "L1_20_Tr1_sub",  "L1_20_Tr2_sub", "L1_20_Tr3_sub", "L1_20_Tr4_sub",
                  "L2_19_Tr1_li", "L2_19_Tr2_li", "L2_19_Tr3_li", "L2_19_Tr4_li", 
                  "L2_19_Tr1_mo", "L2_19_Tr2_mo", "L2_19_Tr3_mo", 
                  "L3_19_Tr1_li",  "L3_19_Tr2_li", "L3_19_Tr3_li", 
                  "L4_19_Tr1_li", "L4_19_Tr2_li", "L4_19_Tr3_li", "L4_19_Tr4_li",
                  "L5_19_Tr1_li", "L5_19_Tr2_li", "L5_19_Tr3_li", "L5_19_Tr4_li",
                  "L6_19_Tr1_li", "L6_19_Tr2_li", "L6_19_Tr3_li", "L6_19_Tr4_li", "L6_19_Tr5_li")

# Decontam ----------------------------------------------------------------
# Use decontam package (https://github.com/benjjneb/decontam) to remove contaminant sequences

# # first time only: install decontam
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("decontam")

# follow decontam vignette from https://benjjneb.github.io/decontam/vignettes/decontam_intro.html:
df <- as.data.frame(sample_data(data_merge)) # pull sample data into df
df$LibrarySize <- sample_sums(data_merge) # get library size for all samples
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

# Because tardigrades are low biomass, need to use prevalence only method:
sample_data(data_merge)$is.neg <- sample_data(data_merge)$Sample_or_Control=="Control" # make into logical vector

# run isContaminant, prevalence method (on ALL experimental data):
# use 0.25 threshold instead of default
contamdf.prev <- isContaminant(data_merge, method="prevalence", neg="is.neg",  
                               batch=data_merge@sam_data$plate,
                               threshold=0.25)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))
fwrite(contamdf.prev, "../output/decontam_output.csv", row.names = T)

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(data_merge, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

# remove the contaminants (and controls, now we're done with them) from the data:
decon <- prune_taxa(!(df.pa$contaminant), data_merge)
decon <- prune_samples(data_merge@sam_data$Sample_or_Control == "Sample", decon)

# # output the removed contaminants:
# contam <- prune_taxa((df.pa$contaminant), data_merge)
# fwrite(as_tibble(rownames(contam@otu_table)),"removed_contaminants.csv")

data_merge <- decon # make the decontaminated data set our main data set

# Removing OTUs that represent less than 10 reads  ---------------

#checking the number of OTUs before pruning
taxa_sum_df <- data.frame(sum = taxa_sums(data_merge))

#pruning the taxa so only OTUs with more than 10 sequences per read remain
data_sub <- data_merge %>%
  prune_taxa(taxa_sums(.) > 9, .)

#checking the number of OTUs left after pruning
pruned_taxa_sum_df <- data.frame(sum = taxa_sums(data_sub))

data_merge <- data_sub
data_merge

# Alpha diversity measures: DivNet ------------------------------------------------

# install packages:
# remotes::install_github("adw96/breakaway")
# remotes::install_github("adw96/DivNet")
# devtools::install_github("bryandmartin/corncob")
library(DivNet)
library(magrittr)

# This section takes more time and computing power than available on my laptop.
# Used HPC instead; time and memory required are shown in comments for each contrast.

# first, prune the samples to only include relevant ones for a given contrast
data_C1 <- prune_samples(data_merge@sam_data$E1==T, data_merge)
data_C2 <- prune_samples(data_merge@sam_data$E2==T, data_merge)
data_C3 <- prune_samples(data_merge@sam_data$E3==T, data_merge)
data_C4 <- prune_samples(data_merge@sam_data$E4==T, data_merge)

# Now, run test for the differences:
# contrast 1--location; took ~40GB and 12hr 11 min
divnet_C1 <- data_C1 %>%
  divnet(formula = ~Location,
         base = "Otu000001")
save(divnet_C1, file="output/divnet_C1")


# contrast 2 -- moss vs lichen
# took 4 GB, 18 minutes
# So, can also run on laptop
divnet_C2 <- data_C2 %>%
  divnet(formula = ~moss.lichen,
         base = "Otu000001")
save(divnet_C2, file="output/divnet_C2")

# # contrast 3 -- tardigrades vs substrate
# NOT possible in regular DivNet, need divnet-rs
# divnet_C3 <- data_C3 %>%
#   divnet(formula = ~Sample_Type,
#          network="diagonal",
#          base = "Otu000001")
# save(divnet_C3, file="output/divnet_C3")

# contrast 4 -- 2019 vs 2020
# took 13 GB, 1 hr 38 min
divnet_C4 <- data_C4 %>%
  divnet(formula = ~Year)
save(divnet_C4, file="output/divnet_C4")

# divnet-rs ---------------------------------------------------------------

# some data too large to run in regular divnet
# On HPC, run:
# module load openblas
# git clone https://github.com/mooreryan/divnet-rs.git
# cd divnet-rs/
# cargo build --release

# make files for divnet-rs:

# full data count table as csv:
otu.rs.full <- as_tibble(data_merge@otu_table)
otu.rs.full$OTU <- rownames(data_merge@otu_table)
otu.rs.full <- otu.rs.full %>%
  select(OTU, everything())
fwrite(otu.rs.full, "full_otu_divnetrs.csv")

# contrast3 count table as csv:
otu.rs3 <- as_tibble(data_C3@otu_table)
otu.rs3$OTU <- rownames(data_C3@otu_table)
otu.rs3 <- otu.rs3 %>%
  select(OTU, everything())
fwrite(otu.rs3, "C3_otu_divnetrs.csv")

# contrast3 sample data:
sam.data3 <- model.matrix(object=~Sample_Type, data=as_tibble(data_C3@sam_data)) # use model.matrix
fwrite(cbind(tibble(sample=data_C3@sam_data$old.sample.ID), 
             Sample_TypeTardigrade=as.tibble(sam.data3)[,2]),
       "C3_samdata_divnetrs.csv")

# full sample data:
sam.data.full <- model.matrix(object=~Sample.ID, data=as_tibble(data_merge@sam_data))
fwrite(cbind(tibble(sample=data_merge@sam_data$old.sample.ID), 
             as.tibble(sam.data.full)[,-1]),
       "full_samdata_divnetrs.csv")

# On server: run config.toml 
# full data took: 14 GB and 5 hr 
# contrast 3 took: 12 GB and 23 hr

# back on laptop: ---------------------------------------------------------
# first, prune the samples to only include relevant ones for a given contrast
data_C1 <- prune_samples(data_merge@sam_data$E1==T, data_merge)
data_C2 <- prune_samples(data_merge@sam_data$E2==T, data_merge)
data_C3 <- prune_samples(data_merge@sam_data$E3==T, data_merge)
data_C4 <- prune_samples(data_merge@sam_data$E4==T, data_merge)


# test for differences in alpha diversity across contrasts:
# Shannon:
attach("output/divnet_C1")
shannon.C1 <- betta(divnet_C1$shannon %>% summary %$% estimate,
                    sqrt(divnet_C1$`shannon-variance`),
                    breakaway::make_design_matrix(data_C1, "Location"))
shannon.C1$table
shannon.C1$global[2]

attach("output/divnet_C2")
shannon.C2 <- betta(divnet_C2$shannon %>% summary %$% estimate,
                    sqrt(divnet_C2$`shannon-variance`),
                    breakaway::make_design_matrix(data_C2, "moss.lichen"))
shannon.C2$table

attach("output/divnet_C4")
shannon.C4 <- betta(divnet_C4$shannon %>% summary %$% estimate,
                    sqrt(divnet_C4$`shannon-variance`),
                    breakaway::make_design_matrix(data_C4, "Year"))
shannon.C4$table

# Simpson:
simpson.C1 <- betta(divnet_C1$simpson %>% summary %$% estimate,
                    sqrt(divnet_C1$`simpson-variance`),
                    breakaway::make_design_matrix(data_C1, "Location"))
simpson.C1$table
simpson.C1$global[2]

simpson.C2 <- betta(divnet_C2$simpson %>% summary %$% estimate,
                    sqrt(divnet_C2$`simpson-variance`),
                    breakaway::make_design_matrix(data_C2, "moss.lichen"))
simpson.C2$table

simpson.C4 <- betta(divnet_C4$simpson %>% summary %$% estimate,
                    sqrt(divnet_C4$`simpson-variance`),
                    breakaway::make_design_matrix(data_C4, "Year"))
simpson.C4$table

# get pairwise comparisons of locations' alpha diversity:
source("https://raw.githubusercontent.com/adw96/breakaway/master/R/betta_lincom.R")
# based on 
# https://github.com/adw96/breakaway/issues/111 and
# https://github.com/adw96/DivNet/issues/37

# Shannon:
# pairwise comparisons: 95% CI for difference between pairs of levels 
# NOTE: any contrast using L1 isn't going to work correctly because L1 is the REFERENCE level;
# So for L1 just use the results of the original shannon.C1 betta
(L2.L3 <- betta_lincom(shannon.C1, c(0,1,-1,0,0,0))) # 95% CI for L2-L3
(L2.L4 <- betta_lincom(shannon.C1, c(0,1,0,-1,0,0))) # 95% CI for L2-L4
(L2.L5 <- betta_lincom(shannon.C1, c(0,1,0,0,-1,0))) # 95% CI for L2-L5
(L2.L6 <- betta_lincom(shannon.C1, c(0,1,0,0,0,-1))) # 95% CI for L2-L6
(L3.L4 <- betta_lincom(shannon.C1, c(0,0,1,-1,0,0))) # 95% CI for L3-L4
(L3.L5 <- betta_lincom(shannon.C1, c(0,0,1,0,-1,0))) # 95% CI for L3-L5
(L3.L6 <- betta_lincom(shannon.C1, c(0,0,1,0,0,-1))) # 95% CI for L3-L6
(L4.L5 <- betta_lincom(shannon.C1, c(0,0,0,1,-1,0))) # 95% CI for L4-L5
(L4.L6 <- betta_lincom(shannon.C1, c(0,0,0,1,0,-1))) # 95% CI for L4-L6
(L5.L6 <- betta_lincom(shannon.C1, c(0,0,0,0,1,-1))) # 95% CI for L5-L6
# multiple testing correction:
p.adjust(c(shannon.C1$table[2,3], shannon.C1$table[3,3], shannon.C1$table[4,3], shannon.C1$table[5,3], shannon.C1$table[6,3],
           as.numeric(L2.L3$`p-values`),as.numeric(L2.L4$`p-values`),as.numeric(L2.L5$`p-values`),as.numeric(L2.L6$`p-values`),
           as.numeric(L3.L4$`p-values`),as.numeric(L3.L5$`p-values`),as.numeric(L3.L6$`p-values`),
           as.numeric(L4.L5$`p-values`),as.numeric(L4.L6$`p-values`),
           as.numeric(L5.L6$`p-values`)), method="BH")

# Simpson:
# pairwise comparisons: 95% CI for difference between pairs of levels 
# NOTE: any contrast using L1 isn't going to work correctly because L1 is the REFERENCE level;
# So for L1 just use the results of the original simpson.C1 betta
(L2.L3 <- betta_lincom(simpson.C1, c(0,1,-1,0,0,0))) # 95% CI for L2-L3
(L2.L4 <- betta_lincom(simpson.C1, c(0,1,0,-1,0,0))) # 95% CI for L2-L4
(L2.L5 <- betta_lincom(simpson.C1, c(0,1,0,0,-1,0))) # 95% CI for L2-L5
(L2.L6 <- betta_lincom(simpson.C1, c(0,1,0,0,0,-1))) # 95% CI for L2-L6
(L3.L4 <- betta_lincom(simpson.C1, c(0,0,1,-1,0,0))) # 95% CI for L3-L4
(L3.L5 <- betta_lincom(simpson.C1, c(0,0,1,0,-1,0))) # 95% CI for L3-L5
(L3.L6 <- betta_lincom(simpson.C1, c(0,0,1,0,0,-1))) # 95% CI for L3-L6
(L4.L5 <- betta_lincom(simpson.C1, c(0,0,0,1,-1,0))) # 95% CI for L4-L5
(L4.L6 <- betta_lincom(simpson.C1, c(0,0,0,1,0,-1))) # 95% CI for L4-L6
(L5.L6 <- betta_lincom(simpson.C1, c(0,0,0,0,1,-1))) # 95% CI for L5-L6
# multiple testing correction:
p.adjust(c(simpson.C1$table[2,3], simpson.C1$table[3,3], simpson.C1$table[4,3], simpson.C1$table[5,3], simpson.C1$table[6,3],
           as.numeric(L2.L3$`p-values`),as.numeric(L2.L4$`p-values`),as.numeric(L2.L5$`p-values`),as.numeric(L2.L6$`p-values`),
           as.numeric(L3.L4$`p-values`),as.numeric(L3.L5$`p-values`),as.numeric(L3.L6$`p-values`),
           as.numeric(L4.L5$`p-values`),as.numeric(L4.L6$`p-values`),
           as.numeric(L5.L6$`p-values`)), method="BH")

# Now, divnet-rs result processing:
# based on https://github.com/mooreryan/divnet-rs/blob/master/test_files/lee_data_walkthrough/3_import_divnet_rs_data.R:
library(tictoc)
setwd("../output")

# set this manually!
nreplicates <- 5

# For contrast 3: # # # #
# Read in data file.
divnet_rs <- fread( cmd = 'findstr "^[^#]" C3_divnet_output.csv', 
                    sep = ",", header = TRUE, data.table = F )

# Replicate 0 is actually the estimates for the real data.
rep0 <- divnet_rs[divnet_rs$replicate == 0, -1]
rownames(rep0) <- rep0$sample
rep0$sample <- NULL

rep0_shannon <- apply(rep0, 1, DivNet::shannon_true)

# Now calculate the shannon index for the replicates.
reps_shannon <- sapply(1:nreplicates, function (i) {
  d <- divnet_rs[divnet_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL
  
  apply(d, 1, DivNet::shannon_true)
})

# What we want is the variance in the diversity estimates for the replicates.
reps_shannon_error <- t(apply(reps_shannon, 1, function (x) {
  c(var(x), sd(x))
}))
shannon.C3 <- betta(chats=unname(rep0_shannon),
                    sqrt(unname(reps_shannon_error[,1])),
                    breakaway::make_design_matrix(data_C3, "Sample_Type"))
shannon.C3$table

# Get the Simpson index for the actual data.
rep0_simpson <- apply(rep0, 1, DivNet::simpson_true)

# Now calculate the simpson index for the replicates.
reps_simpson <- sapply(1:nreplicates, function (i) {
  d <- divnet_rs[divnet_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL
  
  apply(d, 1, DivNet::simpson_true)
})

# What we want is the variance in the diversity estimates for the replicates.
reps_simpson_error <- t(apply(reps_simpson, 1, function (x) {
  c(var(x), sd(x))
}))
colnames(reps_simpson_error) <- c("variance", "sd")

simpson.C3 <- betta(chats=unname(rep0_simpson),
                    sqrt(unname(reps_simpson_error[,1])),
                    breakaway::make_design_matrix(data_C3, "Sample_Type"))
simpson.C3$table

# output estimates:
shannon.tibble=tibble(old.sample.ID=names(rep0_shannon),
                      shannon.est=rep0_shannon,
                      shannon.var=unname(reps_shannon_error[,1]))
write.xlsx(left_join(shannon.tibble, 
                     as_tibble(data_merge@sam_data) %>% select(old.sample.ID, rep.ID)),
           file="alpha_div_by_contrast.xlsx", sheetName="shannon_C3", append=T)
simpson.tibble=tibble(old.sample.ID=names(rep0_simpson),
                      simpson.est=rep0_simpson,
                      simpson.var=unname(reps_simpson_error[,1]))
write.xlsx(left_join(simpson.tibble, 
                     as_tibble(data_merge@sam_data) %>% select(old.sample.ID, rep.ID)),
           file="alpha_div_by_contrast.xlsx", sheetName="simpson_C3", append=T)

# For the full data: # # # #
nreplicates <- 5

divnet_rs <- fread("grep -v '^#' full_divnet_output.csv", header=T, data.table=F)

# Replicate 0 is actually the estimates for the real data.
rep0 <- divnet_rs[divnet_rs$replicate == 0, -1]
rownames(rep0) <- rep0$sample
rep0$sample <- NULL

rep0_shannon <- apply(rep0, 1, DivNet::shannon_true)

# Now calculate the shannon index for the replicates.
reps_shannon <- sapply(1:nreplicates, function (i) {
  d <- divnet_rs[divnet_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL
  
  apply(d, 1, DivNet::shannon_true)
})

# What we want is the variance in the diversity estimates for the replicates.
reps_shannon_error <- t(apply(reps_shannon, 1, function (x) {
  c(var(x), sd(x))
}))
shannon.full <- betta(chats=unname(rep0_shannon),
                      sqrt(unname(reps_shannon_error[,1])),
                      breakaway::make_design_matrix(data_merge, "rep.ID"))
shannon.full$table
shannon.full$global[2]

all.equal(names(rep0_shannon), rownames(reps_shannon_error))

shannon.tibble=tibble(old.sample.ID=names(rep0_shannon),
                      shannon.est=rep0_shannon,
                      shannon.var=unname(reps_shannon_error[,1]))
fwrite(left_join(shannon.tibble, 
                 as_tibble(data_merge@sam_data) %>% select(old.sample.ID, rep.ID)),
       "full_shannon.csv")

# Get the Simpson index for the actual data.
rep0_simpson <- apply(rep0, 1, DivNet::simpson_true)

# Now calculate the simpson index for the replicates.
reps_simpson <- sapply(1:nreplicates, function (i) {
  d <- divnet_rs[divnet_rs$replicate == i, -1]
  rownames(d) <- d$sample
  d$sample <- NULL
  
  apply(d, 1, DivNet::simpson_true)
})

# What we want is the variance in the diversity estimates for the replicates.
reps_simpson_error <- t(apply(reps_simpson, 1, function (x) {
  c(var(x), sd(x))
}))
colnames(reps_simpson_error) <- c("variance", "sd")

simpson.full <- betta(chats=unname(rep0_simpson),
                      sqrt(unname(reps_simpson_error[,1])),
                      breakaway::make_design_matrix(data_merge, "rep.ID"))
simpson.full$table

all.equal(names(rep0_simpson), rownames(reps_simpson_error))

simpson.tibble=tibble(old.sample.ID=names(rep0_simpson),
                      simpson.est=rep0_simpson,
                      simpson.var=unname(reps_simpson_error[,1]))
fwrite(left_join(simpson.tibble, 
                 as_tibble(data_merge@sam_data) %>% select(old.sample.ID, rep.ID)),
       "full_simpson.csv")

# Corncob: modeling and testing relative abundance -------------------------------------

library(corncob)
library(magrittr)

# Make function to pull estimated relative abundances:

get_relative_abundance <- function(da_output, phyloseq_data_full, taxa_level) {
  my.ra <- as_tibble(phyloseq_data_full@sam_data) %>% select(Sample.ID, old.sample.ID)
  for (i in 1:(nrow(da_output$data@tax_table))) {
    if((i-1)%%100==0) {print(i)}
    # for (i in 1:3) {
    # get current formula
    current.OTU <- rownames(da_output$data@tax_table)[i]
    if(taxa_level=="Phylum") {
      current.full <- paste0(as_tibble(da_output$data@tax_table@.Data)$Phylum[i])
    } else if (taxa_level=="Genus") {
      current.full <- paste0(as_tibble(phyloseq_data_full@tax_table@.Data)$Phylum[i], "_",
                             as_tibble(phyloseq_data_full@tax_table@.Data)$Class[i], "_",
                             as_tibble(phyloseq_data_full@tax_table@.Data)$Order[i], "_",
                             as_tibble(phyloseq_data_full@tax_table@.Data)$Family[i], "_",
                             as_tibble(phyloseq_data_full@tax_table@.Data)$Genus[i])
      
      
    } else if (taxa_level=="OTU") {
      current.full <- current.OTU
    } else {warning("check taxa level input")}
    
    tryCatch({
      temp <- tibble(ra=unname(da_output$all_models[[i]]$mu.resp), 
                     old.sample.ID=rownames(da_output$all_models[[i]]$mu.resp)) %>%
        rename_with(.fn=~paste0(., "_", current.full), .cols=-old.sample.ID)
      # stopifnot(all.equal(temp$old.sample.ID, my.ra$old.sample.ID))
      
      my.ra <- left_join(my.ra,
                         temp, by="old.sample.ID")
      # my.ra <- cbind(my.ra, 
      #                    temp %>% select(-old.sample.ID))
    },
    error=function(e){# if warning, skip this phylum
      # print(e)
    })
    
    rm(current.OTU, current.full)
  }
  return(my.ra)
}

# run this on laptop: takes several hours to do all, but doable
data_phylum <- data_merge %>% # make small version with only phyla
  tax_glom("Phylum")
# first, prune the samples to only include relevant ones for a given contrast
phylum_C1 <- prune_samples(data_phylum@sam_data$E1==T, data_phylum)
phylum_C2 <- prune_samples(data_phylum@sam_data$E2==T, data_phylum)
phylum_C3 <- prune_samples(data_phylum@sam_data$E3==T, data_phylum)
phylum_C4 <- prune_samples(data_phylum@sam_data$E4==T, data_phylum)

# test for differential mean and variability in phyla across my contrasts:
# Contrast1: location
c1_phylum_da <- differentialTest(formula = ~ Location,
                                 phi.formula = ~ Location,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Location,
                                 test = "Wald", boot = FALSE,
                                 data = phylum_C1,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_phylum_da$significant_taxa, data = phylum_C1)
plot(c1_phylum_da)

c1_phylum_dv <- differentialTest(formula = ~ Location,
                                 phi.formula = ~ Location,
                                 formula_null = ~ Location,
                                 phi.formula_null = ~ 1,
                                 data = phylum_C1,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_phylum_dv$significant_taxa, data = phylum_C1)
plot(c1_phylum_dv)

# Contrast2: moss.lichen
c2_phylum_da <- differentialTest(formula = ~ moss.lichen,
                                 phi.formula = ~ moss.lichen,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ moss.lichen,
                                 test = "Wald", boot = FALSE,
                                 data = phylum_C2,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_phylum_da$significant_taxa, data = phylum_C2)
plot(c2_phylum_da)

c2_phylum_dv <- differentialTest(formula = ~ moss.lichen,
                                 phi.formula = ~ moss.lichen,
                                 formula_null = ~ moss.lichen,
                                 phi.formula_null = ~ 1,
                                 data = phylum_C2,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_phylum_dv$significant_taxa, data = phylum_C2)
plot(c2_phylum_dv)

# Contrast3: Sample_Type
c3_phylum_da <- differentialTest(formula = ~ Sample_Type,
                                 phi.formula = ~ Sample_Type,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Sample_Type,
                                 test = "Wald", boot = FALSE,
                                 data = phylum_C3,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_phylum_da$significant_taxa, data = phylum_C3)
plot(c3_phylum_da)

c3_phylum_dv <- differentialTest(formula = ~ Sample_Type,
                                 phi.formula = ~ Sample_Type,
                                 formula_null = ~ Sample_Type,
                                 phi.formula_null = ~ 1,
                                 data = phylum_C3,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_phylum_dv$significant_taxa, data = phylum_C3)
plot(c3_phylum_dv)

# Contrast4: Year
c4_phylum_da <- differentialTest(formula = ~ Year,
                                 phi.formula = ~ Year,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ Year,
                                 test = "Wald", boot = FALSE,
                                 data = phylum_C4,
                                 fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c4_phylum_da$significant_taxa, data = phylum_C4)
plot(c4_phylum_da)

c4_phylum_dv <- differentialTest(formula = ~ Year,
                                 phi.formula = ~ Year,
                                 formula_null = ~ Year,
                                 phi.formula_null = ~ 1,
                                 data = phylum_C4,
                                 test = "LRT", boot = FALSE,
                                 fdr_cutoff = 0.05)
c4_phylum_dv$significant_taxa
# otu_to_taxonomy(OTU = c4_phylum_dv$significant_taxa, data = phylum_C4)
plot(c4_phylum_dv)


data_genus <- data_merge %>% # make small version with only genus
  tax_glom("Genus")
# first, prune the samples to only include relevant ones for a given contrast
genus_C1 <- prune_samples(data_genus@sam_data$E1==T, data_genus)
genus_C2 <- prune_samples(data_genus@sam_data$E2==T, data_genus)
genus_C3 <- prune_samples(data_genus@sam_data$E3==T, data_genus)
genus_C4 <- prune_samples(data_genus@sam_data$E4==T, data_genus)

# test for differential mean and variability in genera across my contrasts:
# Contrast1: location
c1_genus_da <- differentialTest(formula = ~ Location,
                                phi.formula = ~ Location,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Location,
                                test = "Wald", boot = FALSE,
                                data = genus_C1,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_genus_da$significant_taxa, data = genus_C1)
plot(c1_genus_da)

c1_genus_dv <- differentialTest(formula = ~ Location,
                                phi.formula = ~ Location,
                                formula_null = ~ Location,
                                phi.formula_null = ~ 1,
                                data = genus_C1,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_genus_dv$significant_taxa, data = genus_C1)
plot(c1_genus_dv)

# Contrast2: moss.lichen
c2_genus_da <- differentialTest(formula = ~ moss.lichen,
                                phi.formula = ~ moss.lichen,
                                formula_null = ~ 1,
                                phi.formula_null = ~ moss.lichen,
                                test = "Wald", boot = FALSE,
                                data = genus_C2,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_genus_da$significant_taxa, data = genus_C2)
plot(c2_genus_da)

c2_genus_dv <- differentialTest(formula = ~ moss.lichen,
                                phi.formula = ~ moss.lichen,
                                formula_null = ~ moss.lichen,
                                phi.formula_null = ~ 1,
                                data = genus_C2,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_genus_dv$significant_taxa, data = genus_C2)
plot(c2_genus_dv)

# Contrast3: Sample_Type
c3_genus_da <- differentialTest(formula = ~ Sample_Type,
                                phi.formula = ~ Sample_Type,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Sample_Type,
                                test = "Wald", boot = FALSE,
                                data = genus_C3,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_genus_da$significant_taxa, data = genus_C3)
plot(c3_genus_da)

c3_genus_dv <- differentialTest(formula = ~ Sample_Type,
                                phi.formula = ~ Sample_Type,
                                formula_null = ~ Sample_Type,
                                phi.formula_null = ~ 1,
                                data = genus_C3,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_genus_dv$significant_taxa, data = genus_C3)
plot(c3_genus_dv)

# Contrast4: Year
c4_genus_da <- differentialTest(formula = ~ Year,
                                phi.formula = ~ Year,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Year,
                                test = "Wald", boot = FALSE,
                                data = genus_C4,
                                fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c4_genus_da$significant_taxa, data = genus_C4)
plot(c4_genus_da)

c4_genus_dv <- differentialTest(formula = ~ Year,
                                phi.formula = ~ Year,
                                formula_null = ~ Year,
                                phi.formula_null = ~ 1,
                                data = genus_C4,
                                test = "LRT", boot = FALSE,
                                fdr_cutoff = 0.05)
c4_genus_dv$significant_taxa
otu_to_taxonomy(OTU = c4_genus_dv$significant_taxa, data = genus_C4)
plot(c4_genus_dv)

data_otu <- data_merge
# first, prune the samples to only include relevant ones for a given contrast
otu_C1 <- prune_samples(data_otu@sam_data$E1==T, data_otu)
otu_C2 <- prune_samples(data_otu@sam_data$E2==T, data_otu)
otu_C3 <- prune_samples(data_otu@sam_data$E3==T, data_otu)
otu_C4 <- prune_samples(data_otu@sam_data$E4==T, data_otu)

# test for differential mean and variability in OTUs across my contrasts:
# Contrast1: location
c1_otu_da <- differentialTest(formula = ~ Location,
                              phi.formula = ~ Location,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Location,
                              test = "Wald", boot = FALSE,
                              data = otu_C1,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_otu_da$significant_taxa, data = otu_C1)
plot(c1_otu_da)

c1_otu_dv <- differentialTest(formula = ~ Location,
                              phi.formula = ~ Location,
                              formula_null = ~ Location,
                              phi.formula_null = ~ 1,
                              data = otu_C1,
                              test = "LRT", boot = FALSE,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c1_otu_dv$significant_taxa, data = otu_C1)
plot(c1_otu_dv)

# Contrast2: moss.lichen
c2_otu_da <- differentialTest(formula = ~ moss.lichen,
                              phi.formula = ~ moss.lichen,
                              formula_null = ~ 1,
                              phi.formula_null = ~ moss.lichen,
                              test = "Wald", boot = FALSE,
                              data = otu_C2,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_otu_da$significant_taxa, data = otu_C2)
plot(c2_otu_da)

c2_otu_dv <- differentialTest(formula = ~ moss.lichen,
                              phi.formula = ~ moss.lichen,
                              formula_null = ~ moss.lichen,
                              phi.formula_null = ~ 1,
                              data = otu_C2,
                              test = "LRT", boot = FALSE,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c2_otu_dv$significant_taxa, data = otu_C2)
plot(c2_otu_dv)

# Contrast3: Sample_Type
c3_otu_da <- differentialTest(formula = ~ Sample_Type,
                              phi.formula = ~ Sample_Type,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Sample_Type,
                              test = "Wald", boot = FALSE,
                              data = otu_C3,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_otu_da$significant_taxa, data = otu_C3)
plot(c3_otu_da)

c3_otu_dv <- differentialTest(formula = ~ Sample_Type,
                              phi.formula = ~ Sample_Type,
                              formula_null = ~ Sample_Type,
                              phi.formula_null = ~ 1,
                              data = otu_C3,
                              test = "LRT", boot = FALSE,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c3_otu_dv$significant_taxa, data = otu_C3)
plot(c3_otu_dv)

# Contrast4: Year
c4_otu_da <- differentialTest(formula = ~ Year,
                              phi.formula = ~ Year,
                              formula_null = ~ 1,
                              phi.formula_null = ~ Year,
                              test = "Wald", boot = FALSE,
                              data = otu_C4,
                              fdr_cutoff = 0.05)
otu_to_taxonomy(OTU = c4_otu_da$significant_taxa, data = otu_C4)
plot(c4_otu_da)

c4_otu_dv <- differentialTest(formula = ~ Year,
                              phi.formula = ~ Year,
                              formula_null = ~ Year,
                              phi.formula_null = ~ 1,
                              data = otu_C4,
                              test = "LRT", boot = FALSE,
                              fdr_cutoff = 0.05)
c4_otu_dv$significant_taxa
otu_to_taxonomy(OTU = c4_otu_dv$significant_taxa, data = otu_C4)
plot(c4_otu_dv)

# output results:
diff.abund.out <- list(contrast=rep(c(1:4), 6),
                       taxa_level=c(rep("phylum", 4),
                                    rep("genus", 4),
                                    rep("otu", 4)),
                       differentially_abundant=rep(NA, 12),
                       differentially_variable=rep(NA, 12))
for(i in 1:12) {
  if(diff.abund.out$taxa_level[i] != "otu") {
    diff.abund.out$differentially_abundant[i] <- ifelse(length(get(paste0("c", diff.abund.out$contrast[i],
                                                                          "_", diff.abund.out$taxa_level[i], "_da"))$significant_taxa)>0,
                                                        list(otu_to_taxonomy(OTU=get(paste0("c", diff.abund.out$contrast[i],
                                                                                            "_", diff.abund.out$taxa_level[i], "_da"))$significant_taxa,
                                                                             data=get(paste0(diff.abund.out$taxa_level[i], "_C", diff.abund.out$contrast[i])))),
                                                        0)
    diff.abund.out$differentially_variable[i] <- ifelse(length(get(paste0("c", diff.abund.out$contrast[i],
                                                                          "_", diff.abund.out$taxa_level[i], "_dv"))$significant_taxa)>0,
                                                        list(otu_to_taxonomy(OTU=get(paste0("c", diff.abund.out$contrast[i],
                                                                                            "_", diff.abund.out$taxa_level[i], "_da"))$significant_taxa,
                                                                             data=get(paste0(diff.abund.out$taxa_level[i], "_C", diff.abund.out$contrast[i])))),
                                                        0)                                                      
    
  } else {
    diff.abund.out$differentially_abundant[i] <- ifelse(length(get(paste0("c", diff.abund.out$contrast[i],
                                                                          "_", diff.abund.out$taxa_level[i], "_da"))$significant_taxa)>0,
                                                        list(get(paste0("c", diff.abund.out$contrast[i],
                                                                        "_", diff.abund.out$taxa_level[i], "_da"))$significant_taxa),
                                                        0)
    diff.abund.out$differentially_variable[i] <- ifelse(length(get(paste0("c", diff.abund.out$contrast[i],
                                                                          "_", diff.abund.out$taxa_level[i], "_dv"))$significant_taxa)>0,
                                                        list(get(paste0("c", diff.abund.out$contrast[i],
                                                                        "_", diff.abund.out$taxa_level[i], "_dv"))$significant_taxa),
                                                        0)
  }
}

fwrite(diff.abund.out, "output/corncob_output.csv")

# get relative abundance values by contrast
write.xlsx(get_relative_abundance(da_output = c1_phylum_da, phyloseq_data_full = data_phylum, taxa_level = "Phylum"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="phylum_C1")
write.xlsx(get_relative_abundance(da_output = c2_phylum_da, phyloseq_data_full = data_phylum, taxa_level = "Phylum"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="phylum_C2", append = T)
write.xlsx(get_relative_abundance(da_output = c3_phylum_da, phyloseq_data_full = data_phylum, taxa_level = "Phylum"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="phylum_C3", append = T)
write.xlsx(get_relative_abundance(da_output = c4_phylum_da, phyloseq_data_full = data_phylum, taxa_level = "Phylum"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="phylum_C4", append = T)
write.xlsx(get_relative_abundance(da_output = c1_genus_da, phyloseq_data_full = data_genus, taxa_level = "Genus"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="genus_C1", append=T)
write.xlsx(get_relative_abundance(da_output = c2_genus_da, phyloseq_data_full = data_genus, taxa_level = "Genus"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="genus_C2", append = T)
write.xlsx(get_relative_abundance(da_output = c3_genus_da, phyloseq_data_full = data_genus, taxa_level = "Genus"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="genus_C3", append = T)
write.xlsx(get_relative_abundance(da_output = c4_genus_da, phyloseq_data_full = data_genus, taxa_level = "Genus"),
           file="output/rel_abund_by_contrast.xlsx", sheetName="genus_C4", append = T)
otu.c1.out <- get_relative_abundance(da_output = c1_otu_da, phyloseq_data_full = data_otu, taxa_level = "OTU")
otu.c2.out <- get_relative_abundance(da_output = c2_otu_da, phyloseq_data_full = data_otu, taxa_level = "OTU")
otu.c3.out <- get_relative_abundance(da_output = c3_otu_da, phyloseq_data_full = data_otu, taxa_level = "OTU")
otu.c4.out <- get_relative_abundance(da_output = c4_otu_da, phyloseq_data_full = data_otu, taxa_level = "OTU")
fwrite(otu.c1.out, "output/otu_rel_abund_C1.csv")
fwrite(otu.c2.out, "output/otu_rel_abund_C2.csv")
fwrite(otu.c3.out, "output/otu_rel_abund_C3.csv")
fwrite(otu.c4.out, "output/otu_rel_abund_C4.csv")

# OTU:
data_OTU <- data_merge
full_OTU_da <- differentialTest(formula = ~ Sample.ID,
                                phi.formula = ~ Sample.ID,
                                formula_null = ~ 1,
                                phi.formula_null = ~ Sample.ID,
                                test = "Wald", boot = FALSE,
                                data = data_OTU,
                                fdr_cutoff = 0.05)
# output significantly different OTUs:
writeLines(c("Differentially abundant OTUs across Sample.ID are: ",
             full_OTU_da$significant_taxa),
           "output/diff_abund_OTUs_Sample.ID.txt") 
# Now output the relative abundance
OTU_ra <- as_tibble(data_merge@sam_data) %>% select(Sample.ID, old.sample.ID)
for (i in 1:(nrow(data_OTU@tax_table))) {
  # for (i in c(1:2)) {
  # get current formula
  print(i)
  current.OTU <- rownames(data_OTU@tax_table)[i]
  current.OTU_full <-otu_to_taxonomy(OTU = current.OTU, data = data_OTU)
  
  tryCatch({
    temp <- tibble(ra=full_OTU_da$all_models[[i]]$mu.resp, 
                   old.sample.ID=rownames(full_OTU_da$all_models[[i]]$mu.resp)) %>%
      rename_with(.fn=~paste0(., "_", current.OTU), .cols=-old.sample.ID)
    
    OTU_ra <- left_join(OTU_ra, 
                        temp, by="old.sample.ID")
  },
  error=function(e){# if warning, skip this OTU
    print(e)})
  
  rm(current.OTU, current.OTU)
}
fwrite(OTU_ra %>% select(-old.sample.ID) %>% distinct,
       "output/OTU_relative_abundance_only_SampleID.csv")

