# README - HZAR_pipeline
Instructions on how to run the R program [HZAR](https://cran.r-project.org/web/packages/hzar/hzar.pdf) (Hybrid Zone Analysis in R) on RADseq genomic data with reference scripts

Latest update: 5/20/21

This doc is notes, instructions, reference scripts, and a few quality of life changes and additonal graphs to the original HZAR software by [Graham Derryberry et. al](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12209). HZAR is an R package designed to do 1 dimensinal geographic cline analyses in hybrid zones. This document serves as an additional resource to the manual provided by Graham Derryberry on [CRAN](https://cran.rstudio.com/web/packages/hzar/hzar.pdf) and the source code on [GitHub](https://github.com/cran/hzar). The reference scipts here and the pipeline was written for use with Restriction site associated DNA as the genetic markers for running the geographic clines on a hybrid zone between the white-collared manakin and golden-collared manakin. 

## Overall Pipeline Summary

raw RAD reads (fastq file) -> see [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) -> hzar input file from `populations` -> `HZAR` -> additional graphs

## Pipeline Steps

### Step 1: Process raw RAD reads
See [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) for instructions on how to process raw RAD reads after you get your giant fastq.gz file from the sequencing facility. At the step for running `populations` in `stacks`, you will need to make sure you add the flag `--hzar` to have the `populations` module output a hzar file, your input file for HZAR.

### Step 2: Format starting input file for HZAR
`stacks` will output a file called `populations.hzar.csv` that you will need to make some minor changes to before running `HZAR`. You first need to remove the first row, as this is an added header from `stacks`. Then you will need to fill in the geographic distances between your populations into the second column. 

### Step 3: Running HZAR


### Step 4: Making some additional graphs
