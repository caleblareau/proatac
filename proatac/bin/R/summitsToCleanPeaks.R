#!/usr/bin/env Rscript

# input: macs2 summit calls, blacklist regions, integer for bp padding and top n
# output: top n refined peak calls with fixed width, no blacklist
# author: Caleb Lareau 9 June 2017, Broad Institute
# inspirted by function in chromVAR by A. Schepp

suppressMessages(suppressWarnings(require(GenomicRanges)))
suppressMessages(suppressWarnings(require(data.table)))

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly = TRUE)
bedfile <- args[1]
blacklist <- args[2]
pad <- args[3]
n <- args[4]

# For dev purposes
#bedfile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/proatac/dirOut/04_qc/proatacProject_summits2.bed"
#blacklist <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/proatac/proatac/anno/blacklist/hg19.full.blacklist.bed"
#pad <- 250
#n <- 1240

# Make GRanges of peaks and blacklist
peakdf <- data.frame(fread(paste0(input = bedfile), header = FALSE))
peaks <- makeGRangesFromDataFrame(setNames(data.frame(peakdf[, 1], peakdf[, 2]-pad, peakdf[, 3]+pad, peakdf[, 5]),
  c("seqnames", "start", "end", "score")), keep.extra.columns = TRUE)

bdf <- data.frame(fread(paste0(input = blacklist), header = FALSE))
bg <- makeGRangesFromDataFrame(setNames(data.frame(
  bdf[, 1], bdf[, 2], bdf[, 3]), c("seqnames", "start", "end")))

# Remove blacklist
peaks <- peaks[!(1:length(peaks) %in% data.frame(findOverlaps(peaks, bg))$queryHits)]
peaks <- sort(peaks)

# Filter peaks based on summit score
keep_peaks <- 1:length(peaks)
while (!(isDisjoint(peaks[keep_peaks]))) {
  
  # Fast variable access
  chr_names <- as.character(seqnames(peaks[keep_peaks]))
  starts <- start(peaks[keep_peaks])
  ends <- end(peaks[keep_peaks])
  scores <- mcols(peaks)$score
  
  # See if consecutive peaks are overlapping
  overlap_next <- intersect(
    which(chr_names[1:(length(keep_peaks) - 1)] == chr_names[2:(length(keep_peaks))]), 
    which(ends[1:(length(keep_peaks) - 1)] >= starts[2:(length(keep_peaks))] )
  )
  
  # Compare consectuive peaks
  overlap_previous <- overlap_next + 1
  overlap_comparison <- scores[keep_peaks[overlap_previous]] > scores[keep_peaks[overlap_next]]
  discard <- keep_peaks[c(overlap_previous[!overlap_comparison], overlap_next[overlap_comparison])]
  keep_peaks <- keep_peaks[keep_peaks %ni% discard]
}

# Export the final result by making a data frame; getting the top (or as many) n peaks
# based on the score and then resort based on genomic position.
fP <- data.frame(peaks[keep_peaks], rank = 1:length(keep_peaks))
nout <- min(n, dim(fP)[1])
odf <- head(fP[order(fP$score, decreasing = TRUE),], nout)

# Write to file
write.table(odf[sort(odf$rank, decreasing = FALSE, index.return = TRUE)$ix, c(1,2,3)],
            file = paste0(dirname(bedfile), "/cleanedRefinedPeaks.bed"),
            col.names = FALSE,row.names = FALSE, sep = "\t", quote = FALSE)

