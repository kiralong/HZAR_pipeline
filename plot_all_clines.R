#Plotting all HZAR Geographic clines together on one plot
library(scales)

#Set working directory
setwd("~/your/working/directory/path")

all_clines <- read.delim("filtered_clines.tsv")

# Make the basic plot
plot(all_clines$kms, all_clines[,2],
     type='n',
     xlab='Distance from M. vitellinus (Km)',
     ylab='Allele Frequency',
     main='Geographic Clines Near BCO2, filtered for CI 10km',
     ylim=c(0,1),
     xlim=c(0,100),
     las=1)

# Add all background plots
for (i in 2:ncol(all_clines)){
  kms <- all_clines$kms
  afq <- all_clines[,i]
  if (afq[1] < afq[length(afq)]){
    afq <- 1-afq
  }
  lines(kms, afq,
        col=alpha('black',0.25),
        lwd=2)
}
# Add the a SNP of interest you want to highlight on the graph
snp  <- 'X122394_27'
kms <- all_clines$kms
afq <- all_clines[,snp]
lines(kms, afq, lwd=4, col=alpha('firebrick',0.8))


