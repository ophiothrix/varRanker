#### Prepare the training set for model development ####
### Script to prepare the training set by combining variants from multiple tissues, removing variants that are being used for validation and testing, removing duplicated negative variants and re-balancing the them.

## The script returns training set and validation set objects as well as the list containing test sets

### Note that this script generates the biggest possible set of variants. I.e. there are no limitations on the closeness of annotation match (both perfectly and imperfectly matched annotations are considered) or on the available annotations (H3K9ac is not required for the variant set to be considered).

# rm(list = ls())
# gc()
prep.sets <- function(path.to.full.set, path.to.test.sets, test.set.tissues) {
	source("./lib/prepare.training.set.R")
	## Select the variant set(s) that are to be used for testing
	# test.set.tissues <- c("HSMM")
	# test.set.tissues <- c("fMuscle")
	# test.set.tissues <- c("fMuscle", "HSMM")
	
	## Load the training set (the combination of all variant sets)
	# train.set <- readRDS("./cache/all.variants.partial.annotation.rds")
	train.set <- readRDS(path.to.full.set)
	dim(train.set)
	table(train.set$regulatory)
	
	### Remove test set variants from the training set ###
	## We want to avoid using the same variants as in the test sets, so remove all the variants with matching variants IDs from the training set
	
	all.sets <- character()
	## Get all annotated variant files to be used as tests (positive sets only)
	for (t.set in test.set.tissues) {
		all.sets <- c(all.sets, list.files(path.to.test.sets, paste0(t.set, ".*variants.rds"), full.names = T))
	}
	
	## For each variant set - annotation combo, make test set and remove those variants from the training set
	sets.list <- list()
	for (t.set in all.sets) {
		print(t.set)
		pos.path <- t.set
		neg.path <- gsub("variants.rds", "negative.set.rds", t.set)
		## Make the test set
		test.set <- prepare.training.set(path.to.positive.set = pos.path, path.to.negative.set = neg.path, prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05, negative.positive.ratio = 2)
		## Subsample negative set to make it the same size as positive
		ids <- sample(which(test.set$regulatory == 0), -diff(table(test.set$regulatory)), replace = F)
		test.set <- test.set[-ids,]
		## Remove the variants in the test set from the training set
		train.set <- train.set[!(train.set$varID %in% test.set$varID),]
		## Save the test set for later
		sets.list[[basename(t.set)]] <- test.set
	}
	names(sets.list)
	table(train.set$regulatory)
	
	
	## Get validation set
	pos.ids <- sample(which(train.set$regulatory == 1), 1000, replace = F)
	neg.ids <- sample(which(train.set$regulatory == 0), 1000, replace = F)
	
	valid.set <- train.set[c(pos.ids, neg.ids),]
	train.set <- train.set[-c(pos.ids, neg.ids),]
	
	## Remove duplicated variants from the validation set
	valid.set <- valid.set[!duplicated(valid.set$varID),]
	table(valid.set$regulatory)
	table(valid.set$source)
	
	## Remove variants from the training set that match variants IDs in the validation set
	train.set <- train.set[!(train.set$varID %in% valid.set$varID),]
	table(train.set$source)
	
	
	## We've deliberately over-sampled the negative examples, so we can have some freedom to remove duplicated negative variants. Let's use that freedom
	summary(duplicated(train.set$varID[train.set$regulatory == 0]))
	train.set <- train.set[-which(duplicated(train.set$varID) & train.set$regulatory == 0),]
	
	## Finally, let's subsample the training set to go back to the original ratios of positive and negative variants
	dnase.ids <- sample(x = which(train.set$source == "dnase"), size = round(sum(train.set$regulatory == 1) * 0.5, 0), replace = F)
	
	matched.ids <- sample(x = which(train.set$source == "matched"), size = round(sum(train.set$regulatory == 1) * 0.45, 0), replace = F)
	
	random.ids <- sample(x = which(train.set$source == "random"), size = round(sum(train.set$regulatory == 1) * 0.05, 0), replace = F)
	
	pos.ids <- which(train.set$regulatory == 1)
	
	## Combine the negative variants into a new training set
	train.set <- train.set[c(pos.ids, dnase.ids, matched.ids, random.ids),]
	table(train.set$source)
	table(train.set$regulatory)
	
	## Add the training and validation sets to the list
	sets.list[["training"]] <- train.set
	sets.list[["validation"]] <- valid.set
	names(sets.list)
	return(sets.list)
}
