annotate.AS.DHS <- function(variant.file.id, tissue.id) {
	require(GenomicRanges)
	variants <- GRanges(read.table(paste0("./cache/snps.multicell.bySample/snps.multicell.", variant.file.id, ".bed"), col.names = c("chr", "start", "end", "AS_DHS"), stringsAsFactors = F))
	## Convert to 1-based coordinates
	start(variants) <- end(variants)
	summary(variants$AS_DHS)
	
	## The tissue-specific variant file doesn't have allele information. Load it from the full variant file
	if (file.exists("./cache/allele.biased.DHS.variants.all.rds")) {
		all.as.dhs <- readRDS("./cache/allele.biased.DHS.variants.all.rds")
	} else {
		source("./lib/prepare.AS_DHS.variants.R")
	}
	
	if (any(is.na(match(variants, all.as.dhs)))) {
		warning("Some of the tisue-specific variants are not found in the global file and can't be assigned REF and ALT alleles.")
	}
	variants$REF <- all.as.dhs[match(variants, all.as.dhs)]$REF
	variants$ALT <- all.as.dhs[match(variants, all.as.dhs)]$ALT
	
	## Subset to the variants that are found to be imbalanced in the given tissue
	variants <- variants[variants$AS_DHS]
	
	## Order variant object by genomic coordinates
	variants <- variants[order(variants)]
	
	## Annotate variants with genomic features, poly-tissue GEP and ChromHMM sets, TF motif scores
	## From all pre-computed scores at the moment it is only feasible to extract a score for ANY variant for CADD. The others are limited to dbSNP variants, more or less...
	## For the moment, leave out the pre-computed scores
	source("./lib/feature.motif.annotation.R")
	variants <- feature.motif.annotation(variants, tissue.id)
	
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
	saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", variant.file.id, ".variants.rds"))
	return(variants)
}
