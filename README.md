# Description

Code to accompany the paper "Inducing Novel Endosymbioses by Bacteria Implantation in Fungi".

# Details

## bacteria_comparison.py

The snakemake pipeline used to compare bacterial samples.

## braker.sh

The job submission script for running BRAKER to annotate the fungus assembly.

## fungus_comparison.py

The snakemake pipeline used to compare fungal samples.

## illumina_assembly.py

The snakemake pipeline used to independently assemble all bacterial samples.

## pacbio_assembly.sh

The job submission script for assembling the fungus Pacbio data.

## snp_analysis.r

The R script used to extract information about SNPs occuring over time along the different lineages.

## summarise_pile.r

The R script used to check for SNP frequencies over time along different lineages.
