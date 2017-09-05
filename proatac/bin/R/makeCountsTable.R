#!/usr/bin/env Rscript

# input: path to .bam files and peak file
# author: Caleb Lareau 5 September 2017

suppressMessages(suppressWarnings(require(chromVAR)))
suppressMessages(suppressWarnings(require(SummarizedExperiment)))
suppressMessages(suppressWarnings(require(tools)))
"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

bamdir <- args[i+1]
peak_file <- args[i+2]
by_rg_in <- args[i+3]
outdir <- args[i+4]
name <- args[i+5]

if(by_rg_in == "True"){
  by_rg <- TRUE
} else {
  by_rg <- FALSE
}

peaks <- getPeaks(peak_file)

bamfiles <- list.files(bamdir, full.names = TRUE)
sampleNames <- file_path_sans_ext(basename(bamfiles))
countsSE <- getCounts(bamfiles, peaks, paired =  TRUE,  by_rg = by_rg, format = "bam", 
                              colData = DataFrame(sampleNames = c(sampleNames)))

counts <- data.frame(data.matrix(assays(countsSE)[["counts"]]))
if(!by_rg) colnames(counts) <- sampleNames
write.table(counts,
            file = paste0(outdir, "/", name, ".counts.tsv"),
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

