### Function to Prepare the input data ###
## Takes paths to GRanges objects for positive and negative sets.
## Removes coding variants, converts features into appropriate formats
## Unifies and combines into a single training set

# path.to.positive.set <- "./cache/E120.annotated.fMuscle.variants.rds"
# path.to.negative.set <- "./cache/E120.annotated.fMuscle.negative.set.rds"
require(GenomicRanges)
prepare.training.set <- function(path.to.positive.set, path.to.negative.set, prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05, negative.positive.ratio = 3) {
	## load the datasets
	positive.set <- readRDS(path.to.positive.set)
	negative.set <- readRDS(path.to.negative.set)
	
	## We'd like to be able to identify the variants, especially, when combining them from multiple tissues. So let's create a variant id from coordinates and alternate allele (reference allele should be the same for a given coordinate)
	positive.set$varID <- paste0(as.character(seqnames(positive.set)), ":", start(positive.set), ":", positive.set$ALT)
	negative.set$varID <- paste0(as.character(seqnames(negative.set)), ":", start(negative.set), ":", negative.set$ALT)
	
	
	## Convert to data.frames
	positive.set <- as.data.frame(mcols(positive.set))
	negative.set <- as.data.frame(mcols(negative.set))
	dim(positive.set)
	dim(negative.set)
	
	## Resample the negative set so that:
	# Total number of negative variants == X times the number of positive
	# The number of negative variants from different sources is deternined by the proportions set
	table(negative.set$source)
	n.matched <- round(nrow(positive.set) * prop.matched * negative.positive.ratio)
	n.dnase <- round(nrow(positive.set) * prop.dnase * negative.positive.ratio)
	n.random <- round(nrow(positive.set) * prop.random * negative.positive.ratio)
	
	## Subsample the negative variants
	set.seed(08062017)
	negative.set <- negative.set[c(sample(x = which(negative.set$source == "matched"), size = n.matched, replace = F),
	  sample(x = which(negative.set$source == "dnase"), size = n.dnase, replace = F),
	  sample(x = which(negative.set$source == "random"), size = n.random, replace = F)),]
	table(negative.set$source)
	
	## Subset to a common set of columns
	## We'd like to keep variant source for debugging, so add source column to the postitive variants
	positive.set$source <- "positive"
	common.cols <- intersect(colnames(positive.set), colnames(negative.set))
	positive.set <- positive.set[,colnames(positive.set) %in% common.cols]
	negative.set <- negative.set[,colnames(negative.set) %in% common.cols]
	
	## Mark the variants with the outcome
	positive.set$regulatory <- 1
	negative.set$regulatory <- 0
	
	## Combine variants into a single data frame
	training.set <- rbind(positive.set, negative.set)
	head(training.set)
	table(training.set$regulatory)
	
	## Remove non-predictive columns
	training.set <- training.set[,!(colnames(training.set) %in% c("id", "REF", "ALT"))]
	
	## Remove coding variants
	table(training.set$feature)
	training.set <- training.set[training.set$feature != "coding",]
	
	## Check if any columns contain NAs
	training.set <- training.set[,!(apply(training.set, 2, function(x) any(is.na(x))))]
	
	## explicitly convert to a data frame
	training.set <- as.data.frame(training.set)
	
	## Convert character features into factors and numeric as necessary
	training.set$feature <- as.factor(training.set$feature)
	training.set$chromHMM.state <- as.factor(training.set$chromHMM.state)
	training.set$GERP.score <- as.numeric(training.set$GERP.score)
	training.set$regulatory <- as.factor(training.set$regulatory)
	if (length(grep("IDEAS.state", colnames(training.set))) != 0) {
		training.set$IDEAS.state <- as.factor(training.set$IDEAS.state)
	}
	if (length(grep("jaspar.motif.hit", colnames(training.set))) != 0) {
		training.set$jaspar.motif.hit <- as.factor(training.set$jaspar.motif.hit)
	}
	if (length(grep("hocomoco.motif.hit", colnames(training.set))) != 0) {
		training.set$hocomoco.motif.hit <- as.factor(training.set$hocomoco.motif.hit)
	}
	training.set$prom.Core <- as.factor(training.set$prom.Core)
	training.set$prom.Domain <- as.factor(training.set$prom.Domain)
	
	## Do the same for all the binary epigenome values
	for (ftre in grep("bin", colnames(training.set), value = T)) {
		training.set[,ftre] <- as.factor(training.set[,ftre])
	}

	## Save prepared training set
	return(training.set)
}
