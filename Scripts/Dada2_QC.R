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


#Generating quality plots
plot_qualF <- plotQualityProfile(fnFs[1])
ggsave(paste0(dirname(Output), '/', basename(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(Output))),'_R1.png'), plot_qualF, device="png")


plot_qualR <- plotQualityProfile(fnRs[1])
ggsave(paste0(dirname(Output), '/', basename(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(Output))),'_R2.png'), plot_qualF, device="png")
