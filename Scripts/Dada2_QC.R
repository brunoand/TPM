#!/usr/local/bin/Rscript

#Importing libraries
library(dada2)
library(ggplot2)

oldw <- getOption("warn")
options(warn = -1)
args = commandArgs(trailingOnly=TRUE)

#Input is the path to the Filtered reads (R1 and R2 in the same folder)
path <- args[1]
Output <- args[2]
#LocR1 <- args[3]
#LocR2 <- args[4]
#Trim <- args[5]

#These variables will hold the list with the full paths for each R of the paired-end reads
fnFs <- sort(list.files(path, pattern="_R1_001.(fastq.gz|fastq|fq|fq.gz)", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.(fastq.gz|fastq|fq|fq.gz)", full.names = TRUE))

print(fnRs)

#Generating quality plots
#plot_qualF <- plotQualityProfile(fnFs, aggregate = TRUE)
#ggsave(paste0(basename(Output), '_R1.pdf'), plot_qualF, device="pdf")


#plot_qualR <- plotQualityProfile(fnRs, aggregate = TRUE)
#ggsave(paste0(basename(Output),'_R2.pdf'), plot_qualR, device="pdf")
