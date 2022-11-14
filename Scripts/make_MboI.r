#!/usr/bin/env Rscript

library(stringr)

args <- commandArgs(TRUE) 
workdir <- args[1]
genome <- args[2]
print(workdir)

filepath <- paste(workdir, '/SourceData/MboI/', genome, '_MboI.txt', sep='') 
print(filepath)
MboI <- read.delim(filepath, sep=' ', header=FALSE)
workdir1 <- paste(workdir, '/SourceData/results_', str_replace_all(genome, '-', '_'), '/gmapped_parsed_sorted_chunks', sep='')
print(workdir1)
setwd(workdir1)

scaff1 <- t(MboI[1, ])
scaff1 <- scaff1[-1, ]
scaff1 <- data.frame(scaff1)
scaff1 <- apply(scaff1, 1, as.numeric)
scaff1 <- data.frame(scaff1)
 
scaff1$start <- scaff1$scaff1 - 25
scaff1$end <- scaff1$scaff1 + 25
scaff1$chrom <- 'scaffold_1'
scaff <- scaff1[, c(4, 2, 3, 1)]
colnames(scaff) <- c('chrom', 'start', 'end', 'position')
for (i in 2:25){
     scaffi <- t(MboI[i, ])
     scaffi <- scaffi[-1, ]
     scaffi <- data.frame(scaffi)
     scaffi <- apply(scaffi, 1, as.numeric)
     scaffi <- data.frame(scaffi)
     
     scaffi$start <- scaffi$scaffi - 25
     scaffi$end <- scaffi$scaffi + 25
     scaffi$chrom <- paste('scaffold_', i , sep = '')
     scaffi <- scaffi[, c(4, 2, 3, 1)]
     colnames(scaffi) <- c('chrom', 'start', 'end', 'position')
     scaff <- rbind(scaff, scaffi)
}
scaff <- scaff[!is.na(scaff$start), ]
scaff[scaff$start < 1, ] <- 1
filepath_out <- paste('../MBoI25.bed', sep = '')

write.table(scaff[,c(-4)], file = filepath_out, sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

