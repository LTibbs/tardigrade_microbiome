# Tardigrade Microbiome
Contains code used for analysis in the paper "Tardigrade community microbiomes in North American orchards includes putative endosymbionts and plant pathogens" (Laura E. Tibbs-Cortes, Bienvenido W. Tibbs-Cortes, and Stephan Schmitz-Esser 2022).

## Taxonomy folder
This folder contains the mothur output that will be read by the R code and a key (`mothur_group_vs_sample_code.xlsx`) to convert between the sample names used in Mothur and the Map file (in `design` subfolder).

To unzip the tar.gz file, type tar -xzvf mothur_output.tar.gz in a terminal.

## Design
This folder contains the experimental design data for each sample.

## Code
This folder contains code used. Use the `tardigrade_microbiome_processing_clean.R` script to process the mothur output (found in `taxonomy` folder). The `config.toml` file is used partway through the R file.
