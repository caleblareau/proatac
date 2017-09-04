#!/usr/bin/env Rscript

# input: macs2 summit calls, blacklist regions, integer for bp padding, chromosome sizes and top n
# output: top n refined peak calls with fixed width, no blacklist
# author: Caleb Lareau 4 September 2017

suppressMessages(suppressWarnings(require(GenomicRanges)))
suppressMessages(suppressWarnings(require(data.table)))

"%ni%" <- Negate("%in%")

args <- commandArgs(trailingOnly = TRUE)
summit_files <- args[1]
peak_width <- args[2]
blacklist <- args[3]
chrom_sizes <- args[4]
n <- args[5]
fdr_threshold <- args[6]
outdir <- args[7]
name <- args[8]

makePeaksDF <- function(summit_files, peak_width, blacklist, chrom_sizes, n = 999999999, fdr_threshold = 0.01){
  
  # Make GRanges of peaks and blacklist
  pad <- round(as.numeric(peak_width)/2)
  l <- lapply(summit_files, function(file){
    if(tools::file_ext(file) == "gz"){
      dt <- fread(paste0("zcat < ", file), stringsAsFactors = TRUE)
    } else if(tools::file_ext(file) %in% c("txt", "csv", "tsv")){
      dt <- fread(paste0(file), stringsAsFactors = TRUE)
    }
    dt
  })
  peakdf <- data.frame(rbindlist( l ))
  peakdf <- peakdf[peakdf$V5 > -1*log10(fdr_threshold), ]
  peakdf <- peakdf[peakdf[,1] %in% chrs,]
  peaks <- makeGRangesFromDataFrame(setNames(data.frame(peakdf[, 1], peakdf[, 2]-pad+1, peakdf[, 3]+pad, peakdf[, 5]),
                                             c("seqnames", "start", "end", "score")), keep.extra.columns = TRUE)
  
  bdf <- data.frame(fread(paste0(input = blacklist), header = FALSE))
  bg <- makeGRangesFromDataFrame(setNames(data.frame(
    bdf[, 1], bdf[, 2], bdf[, 3]), c("seqnames", "start", "end")))
  
  # Remove blacklist and peaks out of bounds
  peaks <- peaks[!(1:length(peaks) %in% data.frame(findOverlaps(peaks, bg))$queryHits)]
  sizedf <- read.table(chrom_sizes, stringsAsFactors = FALSE); names(sizedf) <- c("seqnames", "end"); sizedf$start <- 0
  chrranges <- makeGRangesFromDataFrame(sizedf)
  peaks <- subsetByOverlaps(peaks, chrranges, type = "within")
  peaks <- sortSeqlevels(peaks); peaks <- sort(peaks)
  
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
  peaks <- sortSeqlevels(peaks); peaks <- sort(peaks)
  
  # Export the final result by making a data frame; getting the top (or as many) n peaks
  # based on the score and then resort based on genomic position.
  fP <- data.frame(peaks[keep_peaks], rank = 1:length(keep_peaks))
  nout <- min(as.numeric(n), dim(fP)[1])
  odf <- head(fP[order(fP$score, decreasing = TRUE),], nout)
  return(odf)
}

odf <- makePeaksDF(summit_files, peak_width, blacklist, chrom_sizes, n = n, fdr_threshold = fdr_threshold)

# Write to file
write.table(odf[sort(odf$rank, decreasing = FALSE, index.return = TRUE)$ix, c(1,2,3)],
            file = paste0(outdir, "/", name, ".fixedpeaks.bed"),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

