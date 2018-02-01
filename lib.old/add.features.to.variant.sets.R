# CADD scores were designed to assess both coding and regulatory variants, GWAVA was designed specifically for regulatory variants. eQTLs don't distinguish between coding and non-coding variants. The problem may arise when generating negative sets: if there are a lot of exonic variants in those sets, the CADD score is likely to be biased towards high compared to the regulatory variants. To see if there's any realtionship between the scores and the corresponding feature of the variant, let's assign features to all the negative sets.


# Assign variants to the features -----------------------------------------

## Function to assign variants to the corresponding feature. Taked a GRanges object of variants and assigns each variant to intergenic, exon, or intron regions. As we're mostly interested if something is an exon or not the precedence is as following
# exon > promoter > intron > UTRs3 > UTRs5 > intergenic
# That is to say, if a variant is in the exon and promoter (+/- 1kb of TSS), it is assigned to the exon, if it is in the promoter region, then it gets assigned to promoter even if it is also in the intergenic or intron or UTR space.

# load("./cache/Whole.Blood.eQTL.damage.scores.rda")
# table(eQTL.scores$feature)
# variants <- eQTL.scores

assign.features <- function(variants) {
	require(GenomicRanges)
	require(GenomicFeatures)
	require(TxDb.Hsapiens.UCSC.hg19.knownGene)
	
	#### Assign the eQTLs to genomic features
	print("Assigning variants to genomic features")
	variants$feature <- NA
	# Get intergenic vs genic regions for hg19
	genic.regions <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
	# If the genes overlap (including genes on the opposite strands) consider them a single genic region - change strand to * and reduce.
	strand(genic.regions) <- "*"
	genic.regions <- reduce(genic.regions)
	# Get intergenic regions by calculating the gaps between the genes
	intergenic.regions <- gaps(genic.regions)
	# Gaps forces strand information, i.e. whole chr1 positive and negative strands are considered gaps relative to the strandless GRanges object. Remove the stranded gaps.
	intergenic.regions <- intergenic.regions[strand(intergenic.regions) == "*"]
	variants$feature[queryHits(findOverlaps(query = variants, subject = intergenic.regions))] <- "intergenic"
	
	# Get intron coordinates
	intronic.regions <- unlist(intronsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
	# Remove strand information and reduce
	strand(intronic.regions) <- "*"
	intronic.regions <- reduce(intronic.regions)
	variants$feature[queryHits(findOverlaps(query = variants, subject = intronic.regions))] <- "intron"
	
	# Get exon coordinates
	exonic.regions <- unlist(exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "tx"))
	# Remove strand information and reduce
	strand(exonic.regions) <- "*"
	exonic.regions <- reduce(exonic.regions)
	variants$feature[queryHits(findOverlaps(query = variants, subject = exonic.regions))] <- "exon"
	
	# Get 5' UTRs
	UTRs5 <- unlist(fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
	# Remove strand information and reduce
	strand(UTRs5) <- "*"
	UTRs5 <- reduce(UTRs5)
	variants$feature[queryHits(findOverlaps(query = variants, subject = UTRs5))] <- "UTRs5"
	
	# Get 3' UTRs
	UTRs3 <- unlist(threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg19.knownGene))
	# Remove strand information and reduce
	strand(UTRs3) <- "*"
	UTRs3 <- reduce(UTRs3)
	variants$feature[queryHits(findOverlaps(query = variants, subject = UTRs3))] <- "UTRs3"
	
	# Get Promoter regions (+/- 3kb of TSS)
	promoter.regions <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene, upstream = 1000, downstream = 1000)
	# Remove strand information and reduce
	strand(promoter.regions) <- "*"
	promoter.regions <- reduce(promoter.regions)
	variants$feature[queryHits(findOverlaps(query = variants, subject = promoter.regions))] <- "promoter"
	
	# Get coding sequence coordinates
	coding.regions <- unlist(cdsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "tx"))
	# Remove strand information and reduce
	strand(coding.regions) <- "*"
	coding.regions <- reduce(coding.regions)
	variants$feature[queryHits(findOverlaps(query = variants, subject = coding.regions))] <- "coding"
	
	
	return(variants)    
}

