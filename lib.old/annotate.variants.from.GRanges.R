# Annotate a list of variants with aritrary ranges ------------------------


## Generic function that takes a list of variants and a list of arbitrary ranges and annotates the variants based on the ovelap with the ranges
## Returns a vector of annotations
## If there is a n.tissue column in the region, also return the corresponding number

annotate.variants.from.GRanges <- function(variants, regions) {
  require(GenomicRanges)
  if(length(grep("chr", seqlevels(variants))) == 0) {
    seqlevels(regions) <- gsub("chr", "", seqlevels(regions))
  }
  annotations <- rep(0, length(variants))
  olaps <- findOverlaps(variants, regions)
  annotations[queryHits(olaps)] <- 1
  if (!is.null(regions$n.tissues)) {
    n.tissues <- rep(0, length(variants))
    n.tissues[queryHits(olaps)] <- regions$n.tissues[subjectHits(olaps)]
    return(as.data.frame(cbind(annotations, n.tissues)))
  } else {
    return(annotations)
  }
}
