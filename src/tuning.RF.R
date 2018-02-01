## Load training set
rm(list = ls())
train.set <- readRDS("./cache/train.full.rds")
dim(train.set)
table(train.set$regulatory)

## Get validation set
set.seed(17012018)
pos.ids <- sample(which(train.set$regulatory == 1), 1000, replace = F)
neg.ids <- sample(which(train.set$regulatory == 0), 1000, replace = F)

valid.set <- train.set[c(pos.ids, neg.ids),]
train.set <- train.set[-c(pos.ids, neg.ids),]

## Train a random forest model
require(pROC)
require(h2o)
require(RColorBrewer)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
nfolds <- 5
clrs <- brewer.pal(5, "Dark2")

## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(train.set))) != 0){
	train.set <- train.set[-grep("qval", colnames(train.set))]
}


train.set.h2o <- as.h2o(train.set)
valid.set.h2o <- as.h2o(valid.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

## With cross-validation
# model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
## Specifying a validation frame
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = 200, max_depth = 40, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)


summary(model_drf)

# h2o.performance(model_drf, xval = T)
# h2o.performance(model_drf, train = T)
h2o.auc(model_drf, xval = T)
as.data.frame(h2o.varimp(model_drf))
# 0.7752843 with cadd; 0.7729757 without
# pdf("./graphs/feature.importance.15.12.17.pdf", 12, 8)
h2o.varimp_plot(model_drf, num_of_features = 40)
as.data.frame(h2o.varimp(model_drf))
# dev.off()
return(h2o.auc(model_drf, xval = T))
}

model.check.from.df <- function(train.set) {
	require(pROC)
	require(h2o)
	require(RColorBrewer)
	## Initiate H2O instance
	localH2O <- h2o.init(nthreads=-1)
	nfolds <- 5
	clrs <- brewer.pal(5, "Dark2")
	### Prepare the input data ###
	
	head(train.set)
	table(train.set$regulatory)
	dim(train.set)
	colnames(train.set)
	
	## Remove qval features, as they are redundant with pval
	if( length(grep("qval", colnames(train.set))) != 0){
		train.set <- train.set[-grep("qval", colnames(train.set))]
	}
	
	
	## Try h2o random forest with DNase
	train.set.h2o <- as.h2o(train.set)
	# train.set.h2o <- as.h2o(no.dnase.train.set)
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
