#### Annotate variants from a VCF file ####

### Function to annotate variants from a vcf file with a given tissue-specific annotation and convert them into a dataframe compatible with prediction.
path.to.vcf <- "./cache/BRCA-US.small.vcf.gz"
tissue.id <- "E119"
annotate.vcf <- function(path.to.vcf, tissue.id) {
	require(GenomicRanges)
	require(data.table)
	if (length(grep("gz", path.to.vcf)) == 0) {
		variants <- fread(path.to.vcf, stringsAsFactors = F)
	} else {
		variants <- fread(paste0("gunzip -c ", path.to.vcf), stringsAsFactors = F)
	}
	head(variants, 12)
	## Check if there is a SNP id column
	if (length(grep("rs", variants$V3)) == 0) {
		colnames(variants)[1:4] <- c("chr", "start", "REF", "ALT")
	} else {
		colnames(variants)[1:5] <- c("chr", "start", "rsID", "REF", "ALT")
	}
	
	## Extract allele frequencies
	colnames(variants)[grep("AF", variants[1,])] <- "info"
	variants[, AF:=gsub(".*AF=([0-9\\.,]+);.*", "\\1", variants[["info"]])]
	
	## Remove extra columns
	variants[, grep("^V[0-9]+", colnames(variants)):=NULL]
	variants[, c("info", "rsID"):=NULL]
	
	## Convert to UCSC chromosome names, if not already
	if (length(grep("chr", variants$chr)) == 0) {
		variants[, chr:=paste0("chr", chr)]
	}
	
	## Check if it's a multi-allele vcf. If so, split into single alleles
	if (length(grep(",", variants$alt)) != 0) {
		print("The vcf file is multi-allelic. Splitting into single alleles.")
		split.alt <- strsplit(variants$alt, ",")
		split.af <- strsplit(variants$AF, ",")
		n.alleles <- unlist(lapply(split.alt, length))
		variants <- variants[rep(1:nrow(variants), n.alleles)]
		variants$alt <- unlist(split.alt)
		variants$AF <- unlist(split.af)
	}
	
	## Add end coordinate
	variants[,end:=start]
	head(variants, 12)
	
	## Convert to GRanges object
	variants <- GRanges(variants)
	
	## Order variant object by genomic coordinates
	variants <- variants[order(variants)]
	
	## Remove indels (for now)
	variants <- variants[nchar(variants$REF) == 1 & nchar(variants$ALT) == 1]
	
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
	saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))
	return(variants)
}

