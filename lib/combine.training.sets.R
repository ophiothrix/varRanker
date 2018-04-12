### Combine multiple sets into a single training set ###
combine.training.sets <- function(positive.sets, dnase = 0.5, matched = 0.5, random = 0.05) {
	### Function that combines several positive (and corresponding negative) sets into a single mega training set
	#' @param positive.sets vector of paths to the positive sets to be included in the mega-set
	#' @param dnase proportion of negative DNase variants to include
	#' @param matched proportion of matched negative variants to include
	#' @param random proportion of randomly sampled genomic variants to include
	#' @return Training set data.frame containing positive and negative variants from the target sets in pre-defined proportions.
	source("./lib/prepare.training.set.R")
	
	# ## List the positive training sets
	# positive.sets <- list.files("./cache", "variants.rds", full.names = T)
	# positive.sets <- positive.sets[-grep("cato", positive.sets)]
	# positive.sets <- positive.sets[-grep("E119.annotated.fMuscle.variants.rds", positive.sets)]
	# length(positive.sets)

	training.sets <- list()
	for (positive.set in positive.sets) {
		print(positive.set)
		training.sets[[positive.set]] <- prepare.training.set(path.to.positive.set = positive.set, path.to.negative.set = gsub("variants", "negative.set", positive.set), prop.matched = matched, prop.dnase = dnase, prop.random = random)
	}
	col.intersect <- colnames(training.sets[[which.min(unlist(lapply(training.sets, ncol)))]])
	
	combined.set <- data.frame()
	for (s in 1:length(training.sets)) {
		print(names(training.sets)[s])
		combined.set <- rbind(combined.set, training.sets[[s]][,col.intersect])
	}
	
	## Remove duplicated entries
	combined.set <- combined.set[!duplicated(combined.set),]
	
	return(combined.set)
}
