#### The minimal script to test model performance on the pre-specified training set

##!! To be run after every change to the script !!##

## Setting up
# rm(list = ls())
# gc()

# path.to.positive <- "./cache/E126.annotated.fSkin_fibro.variants.rds"
# path.to.negative <- "./cache/E126.annotated.fSkin_fibro.negative.set.rds"

model.check.from.path <- function(path.to.positive, path.to.negative) {
	require(pROC)
	require(h2o)
	require(RColorBrewer)
	## Initiate H2O instance
	localH2O <- h2o.init(nthreads=-1)
	nfolds <- 5
	clrs <- brewer.pal(5, "Dark2")
	### Prepare the input data ###
	source("./lib/prepare.training.set.R")
	training.set <- prepare.training.set(path.to.positive.set = path.to.positive, path.to.negative.set = path.to.negative, prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
	
	head(training.set)
	table(training.set$regulatory)
	dim(training.set)
	colnames(training.set)
	
	## Remove qval features, as they are redundant with pval
	if( length(grep("qval", colnames(training.set))) != 0){
		training.set <- training.set[-grep("qval", colnames(training.set))]
	}
	
	
	## Try h2o random forest with DNase
	train.set.h2o <- as.h2o(training.set)
	# train.set.h2o <- as.h2o(no.dnase.training.set)
	target <- "regulatory"
	features <- setdiff(colnames(train.set.h2o), target)
	## Remove cadd score from training features
	features <- features[features != "cadd.score"]
	# features <- features[-grep("global", features)]
	# features <- features[-grep("local", features)]
	
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	# summary(model_drf)
	# h2o.performance(model_drf, xval = T)
	# h2o.performance(model_drf, train = T)
	h2o.auc(model_drf, xval = T)
	# 0.7752843 with cadd; 0.7729757 without
	# pdf("./graphs/feature.importance.15.12.17.pdf", 12, 8)
	h2o.varimp_plot(model_drf, num_of_features = 40)
	as.data.frame(h2o.varimp(model_drf))
	# dev.off()
	return(h2o.auc(model_drf, xval = T))
}

model.check.from.df <- function(training.set) {
	require(pROC)
	require(h2o)
	require(RColorBrewer)
	## Initiate H2O instance
	localH2O <- h2o.init(nthreads=-1)
	nfolds <- 5
	clrs <- brewer.pal(5, "Dark2")
	### Prepare the input data ###
	
	head(training.set)
	table(training.set$regulatory)
	dim(training.set)
	colnames(training.set)
	
	## Remove qval features, as they are redundant with pval
	if( length(grep("qval", colnames(training.set))) != 0){
		training.set <- training.set[-grep("qval", colnames(training.set))]
	}
	
	
	## Try h2o random forest with DNase
	train.set.h2o <- as.h2o(training.set)
	# train.set.h2o <- as.h2o(no.dnase.training.set)
	target <- "regulatory"
	features <- setdiff(colnames(train.set.h2o), target)
	## Remove cadd score from training features
	features <- features[features != "cadd.score"]
	# features <- features[-grep("global", features)]
	# features <- features[-grep("local", features)]
	
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	# summary(model_drf)
	# h2o.performance(model_drf, xval = T)
	# h2o.performance(model_drf, train = T)
	h2o.auc(model_drf, xval = T)
	# 0.7752843 with cadd; 0.7729757 without
	# pdf("./graphs/feature.importance.15.12.17.pdf", 12, 8)
	h2o.varimp_plot(model_drf, num_of_features = 40)
	as.data.frame(h2o.varimp(model_drf))
	# dev.off()
	return(h2o.auc(model_drf, xval = T))
}

