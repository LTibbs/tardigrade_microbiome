# Calculating Bray-Curtis and Euclidean distance from Divnet-rs output
# NOTE: this was NOT used in the published analysis but was requested by a divnet-rs user
# Based on https://github.com/mooreryan/divnet-rs/issues/3 

# load libraries
library(tidyverse)
library(data.table)
library(DivNet)
library(magrittr)

# you will need to set your working directory to the location of the divnet-rs output

nreplicates <- 5 # set to number of replicates used in divnet-rs

# this reads in the output file from divnet-rs. Change file name to match your own divnet-rs output
divnet_rs <- fread( cmd = 'findstr "^[^#]" full_divnet_output.csv', 
                    sep = ",", header = TRUE, data.table = F )

# Replicate 0 is actually the estimates for the real data.
rep0 <- divnet_rs[divnet_rs$replicate == 0, -1]
rownames(rep0) <- rep0$sample
rep0$sample <- NULL

# pull Bray-Curtis distance:
BC <- matrix(nrow=nrow(rep0), ncol=nrow(rep0))
colnames(BC) <- rownames(rep0)
rownames(BC) <- rownames(rep0)

# loop through and calculate the BC for each pair of samples:
# Sys.time()
for (i in 1:nrow(rep0)) {
  print(paste0("i=",i))
  for(j in 1:nrow(rep0)) {
    if(j>i) {break}
    BC[i,j] <- DivNet::bc_fast(rep0[i,], rep0[j,])
  }
}
# Sys.time()
fwrite(BC, "Bray_Curtis.csv")

# pull Euclidean distance:
Euc <- matrix(nrow=nrow(rep0), ncol=nrow(rep0))
colnames(Euc) <- rownames(rep0)
rownames(Euc) <- rownames(rep0)

for (i in 1:nrow(rep0)) {
  print(paste0("i=",i))
  for(j in 1:nrow(rep0)) {
    if(j>i) {break}
    Euc[i,j] <- DivNet::euc_fast(rep0[i,], rep0[j,])
  }
}
fwrite(Euc, "Euclidean_distance.csv")