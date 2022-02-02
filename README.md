# Vaccinia intron coordinates

*Michael Ly, Hannah M Burgess, Sahil B. Shah, Ian Mohr, Britt A. Glaunsinger*

# Script overview

This script extracts intron coordinates based on exon data from gtf files and outputs a new gtf file containing all genomic coordinates including introns. The RNA-seq data to replicate our work in the paper, "Vaccinia virus D10 has broad decapping activity that is regulated by mRNA splicing", published on ... can be found using the GEO accession number: GSE185520. However, this script can be used generally on any gtf file of interest.

## Requirements

Most packages that this script relies upon are from the standard scientific/numeric python distributions. However, the AGEpy package developed at the Bioinformatics Core Facility of the Max Planck Institute for Biology of Ageing was also used to output python data frames as a gtf file. The instructions to install this package can be found here: <https://github.com/mpg-age-bioinformatics/AGEpy>. 

## Input files

The only required input file needed is a gtf file. The original gtf file was created using RNA-seq data from the first replicate of the "WT_NoDox" and "WT_PlusDox" conditions.
