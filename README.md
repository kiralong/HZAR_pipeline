# README - HZAR_pipeline
Instructions on how to run the R program [HZAR](https://cran.r-project.org/web/packages/hzar/hzar.pdf) (Hybrid Zone Analysis in R) on RADseq genomic data with reference scripts

Latest update: 5/20/21

This doc is notes, instructions, reference scripts, and a few quality of life changes and additonal graphs to the original HZAR software by [Graham Derryberry et. al](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12209). HZAR is an R package designed to do 1 dimensinal geographic cline analyses in hybrid zones. This document serves as an additional resource to the manual provided by Graham Derryberry on [CRAN](https://cran.rstudio.com/web/packages/hzar/hzar.pdf) and the source code on [GitHub](https://github.com/cran/hzar). The reference scipts here and the pipeline was written for use with Restriction site associated DNA as the genetic markers for running the geographic clines on a hybrid zone between the white-collared manakin and golden-collared manakin. 

## Overall Pipeline Summary

raw RAD reads (fastq file) -> see [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) -> hzar input file from `populations` -> `HZAR` -> additional graphs

## Pipeline Steps

### Step 1: Process raw RAD reads
See [RADseq_pipeline](https://github.com/kiralong/RADseq_pipeline) for instructions on how to process raw RAD reads after you get your giant fastq.gz file from the sequencing facility. At the step for running `populations` in `stacks`, you will need to make sure you add the flag `--hzar` to have the `populations` module output a hzar file, your input file for HZAR.

Note that if you are running HZAR on an entire chromosome or genome, you may want to only select SNPs that are fixed differences between the parental species. These clines will be at a frequency of 0 in one parental and 1 in the other. These clines will be much cleaner and easier to interpret and you should still have hundreds if not thousands of them in your RAD dataset. Use the script [find_fixed_parental_snps.py](find_fixed_parental_snps.py) on the `populations.sumstats.tsv` file from `populations`. Note that you will need to edit thhis python script to use your population names instead of my manakin parental populations. You can then use the output of [find_fixed_parental_snps.py](find_fixed_parental_snps.py) to create a whitelist of only fixed parental snps to then rerun `populations` with the `--hzar` flag to get a smaller, cleaner dataset fro HZAR to run on an entire chromosome or your whole genome if you'd like. 

### Step 2: Format starting input file for HZAR
`stacks` will output a file called `populations.hzar.csv` that you will need to make some minor changes to before running `HZAR`. You first need to remove the first row, as this is an added header from `stacks`. Then you will need to fill in the geographic distances between your populations into the second column. 

### Step 3: Running HZAR
Use the script [hzarscript_loop.R](hzarscript_loop.R) to run `HZAR`. If you are running 1,000 snps or more, I suggest running this in a computing cluster. Running about 1,000 SNPS will take `HZAR` about 2 days. This script is written to take in 2 data files: a list of your SNP IDs and the `HZAR` input file you get out of `populations` in `stacks`. The rest of the scrit is written with variables so you only need to change the inputs and population ID names and distances at the top of the script. After running [hzarscript_loop.R](hzarscript_loop.R) you should have 3 output files: individual_snp_clines.pdf, snp_clines_parameters.tsv, and all_clines.tsv.

Individual_snp_clines.pdf is a pdf graphing each individual geographic cline for every snp you had in your snp ID's file. snp_cline_parameters.tsv has a bunch of parameter values such as your cline center location and the confidence intervals. Lastly all_clines.tsv has the graphing shape information for all the clines you ran so you can graph them all on the same graphing space. 

### Step 4: Making some additional graphs
If you want to graph all the clines you ran on the same graphing space, instead of having them all on their own individual pdf page, use the script [plot_all_clines.R](plot_all_clines.R). You can also highlight clines of special interest in different colors to show if their cline centers are vastly shifted, for example. 

