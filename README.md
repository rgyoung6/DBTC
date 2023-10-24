# Dada-BLAST-Taxon Assign-Condense (DBTC) package
DBTC is a R package for metabarcode analysis. This package is a non-shiny version of the DBTCShiny package. For instructions and specifics on using this command line version of the DBTCShiny please see the DBTCShiny [README](https://github.com/rgyoung6/DBTCShiny/tree/main#readme)

# Description

This repository contains the DBTC package located at rgyoung6/DBTC. The Dada-BLAST-Taxon Assign-Condense  package contains DBTC functions to process metabarcode data. These functions are executed in a R terminal using the command line. The DBTC functions have four main outcomes...

  - Fastq file processing using Dada in R
  - Using the Basic Local Alignment Search Tool (BLAST) amplicon sequence variants (ASV) can be searched against local NCBI or custom libraries
  - Assign taxa to the unique reads using NCBI taxon database through taxa names and/or taxaID's
  - Condense the resulting ASV taxonomic assignment tables to unique taxa with the ability to combine datasets (using different sequence libraries for the same reads, or results from the same samples for different molecular regions) into a combined results table

# Installation via GitHub

Run the following commands in your R terminal...<br/>

```
install.packages("devtools")
library(devtools)
devtools::install_github("rgyoung6/DBTC")
library(DBTC)
```

**Note:** the first command to install the "devtools" may not be necessary if already installed.<br/>
