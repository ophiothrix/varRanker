annotate.GRanges <- function(variants, tissue.id) {
	require(GenomicRanges)
	
	## Annotate variants with genomic features, poly-tissue GEP and ChromHMM sets, TF motif scores
	## From all pre-computed scores at the moment it is only feasible to extract a score for ANY variant for CADD. The others are limited to dbSNP variants, more or less...
	## For the moment, leave out the pre-computed scores
	source("./lib/feature.motif.annotation.R")
	variants <- feature.motif.annotation(variants)
	
	## Annotate variants with conservation scores
	source("./lib/add.conservation.scores.R")
	variants <- annotate.GERP(variants)
	variants <- annotate.phastCons(variants)
	variants <- annotate.phyloP(variants)
	
	## Annotate variants with tissue-specific epigenetic marks and ChromHMM and GEP predicted states
	source("./lib/epigenome.annotation.R")
	variants <- epigenome.annotation(variants, tissue.id)
	# colnames(mcols(variants))
	
	## Save the variant file
	return(variants)
}
