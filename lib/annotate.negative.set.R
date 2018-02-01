#### Function that extracts 1KG variants with MAF > 25% matched to a given positive variant set and annotates them in the same manner as the positive set. ####
## Takes as arguments positive variant file id and ROADMAP tissue id used to annotate positive variants.
## Extracts common 1KG variants within 2 kb of the positive variants, as well as the same number of randomly sampled variants across genome.
## Returns annotated set of matched and random negative variants.
## Saves GRanges object containing annotated variants

annotate.negative.set <- function(variant.file.id, tissue.id, window.to.match = 1000, maf.cutoff = 0.05) {
	require(GenomicRanges)
	print("Putting together a set of negative variants...")	
	path.to.positive.set <- paste0("./cache/", tissue.id, ".annotated.", variant.file.id, ".variants.rds")
	
	### We'll need to exclude any variants that has been found to show imbalanced DHS or been found to act as an eQTL. ###
	
	## Load allele imbalanced DHS variants
	if (file.exists("./cache/allele.biased.DHS.variants.all.rds")) {
		all.as.dhs <- readRDS("./cache/allele.biased.DHS.variants.all.rds")
	} else {
		source("./lib/prepare.AS_DHS.variants.R")
	}
	## save not imbalances DHS variants into a separate object
	negative.as.dhs <- all.as.dhs[all.as.dhs$significance.level == "not_imbalanced"]
	all.as.dhs <- all.as.dhs[all.as.dhs$significance.level != "not_imbalanced"]
	
	## Load eQTLs
	# for fls in `ls *snpgenes`
	# do 
	#   echo $fls
	#   cut -f14,15,20,21 $fls | tail -n +2 >> all.eQTLS.txt
	# done
	# sort -u all.eQTLS.txt | gzip -c all.eQTLS.txt.gz
	# rm all.eQTLS.txt
	if (file.exists("./data/GTEx/all.eQTLs.rds")) {
		eQTLs <- readRDS("./data/GTEx/all.eQTLs.rds")
	} else {
		source("./lib/process.GTEx.data.R")
		eQTLs <- readRDS("./data/GTEx/all.eQTLs.rds")
	}
	
	## Load variants from positive training set
	variants <- readRDS(path.to.positive.set)
	
	## Load 1KG dataset
	variants.1KG <- readRDS("./data/control.variants/variants.1KG.full.rds")
	## Subset to variants with MAF > 0.25
	summary(variants.1KG$MAF > maf.cutoff)
	common.variants <- variants.1KG[variants.1KG$MAF > maf.cutoff]
	## Let's not limit the 1000 genomes variants to high frequencies. Probably doesn't matter much anyway, as we only include very few variants from this set.
	# common.variants <- variants.1KG
	
	## get 1KG variants within 2 kb of the positive set variants
	search.space <- GenomicRanges::shift(variants, -(window.to.match/2))
	width(search.space) <- window.to.match
	
	olaps <- findOverlaps(common.variants, search.space)
	olaps
	matched.1KG.variants <- common.variants[queryHits(olaps)]
	
	## Sample random 1KG variants, the same number as the matched variants
	## Remove 1KG variants that overlap with reported eQTLs and AS DHS
	# summary(matched.1KG.variants %in% all.as.dhs)
	# summary(matched.1KG.variants %in% eQTLs)
	olaps <- findOverlaps(matched.1KG.variants, eQTLs)
	matched.1KG.variants <- matched.1KG.variants[-queryHits(olaps)]
	olaps <- findOverlaps(matched.1KG.variants, all.as.dhs)
	matched.1KG.variants <- matched.1KG.variants[-queryHits(olaps)]
	matched.1KG.variants$source <- "matched"
	set.seed(15052017)
	ids <- sample(1:length(common.variants), size = length(matched.1KG.variants)*1.3, replace = F)
	random.1KG.variants <- common.variants[ids]
	# summary(random.1KG.variants %in% all.as.dhs | random.1KG.variants %in% eQTLs)
	# summary(random.1KG.variants %in% eQTLs)
	olaps <- findOverlaps(random.1KG.variants, eQTLs)
	random.1KG.variants <- random.1KG.variants[-queryHits(olaps)]
	olaps <- findOverlaps(random.1KG.variants, all.as.dhs)
	random.1KG.variants <- random.1KG.variants[-queryHits(olaps)]
	random.1KG.variants$source <- "random"
	## Remove MAF column
	mcols(matched.1KG.variants) <- mcols(matched.1KG.variants)[,c(1,2,4)]
	mcols(random.1KG.variants) <- mcols(random.1KG.variants)[,c(1,2,4)]
	
	
	## Finally, sample randomly from the DHS variants with no evidence for allele bias (q value > 0.99)
	## Additionally, the variant could be not detected as imbalanced, if there are too few het individuals.
	## Add additional filter to require at least 10 hets
	summary(negative.as.dhs$q.value > 0.90 & negative.as.dhs$numhets >= 10)
	dnase.variants <- negative.as.dhs[negative.as.dhs$q.value > 0.90 & negative.as.dhs$numhets >= 10]
	## Sample the same number of DNase variants as for random variants. Or take all of them, if random variant size is greater
	if (length(random.1KG.variants) < length(dnase.variants)) {
		dnase.variants <- dnase.variants[sample(1:length(dnase.variants), size = length(random.1KG.variants), replace = F)]
	}
	dnase.variants$source <- "dnase"
	mcols(dnase.variants) <- mcols(dnase.variants)[,11:13]
	
	## Combine matched and random variants
	rm(variants)
	variants <- c(matched.1KG.variants, random.1KG.variants, dnase.variants)
	## Reorder 1KG variants
	variants <- variants[order(variants)]
	
	### Annotate negative variant set with the same tissue as the positive set ###
	
	## Subset to SNVs only
	require(stringr)
	all(nchar(variants$REF) == 1)
	all(nchar(variants$ALT) == 1)
	
	## Annotate variants with genomic features, poly-tissue GEP and ChromHMM sets, TF motif scores, pre-computed damage scores
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
	
	## Save the variant file
	saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", variant.file.id, ".negative.set.rds"))
	return(variants)
}
