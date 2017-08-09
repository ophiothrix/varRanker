rm(list = ls())
gc()
require(pROC)
require(h2o)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
nfolds <- 5

require(RColorBrewer)
clrs <- brewer.pal(5, "Dark2")

### Prepare the input data ###
source("./lib/prepare.training.set.R")
# training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds")

head(training.set)
dim(training.set)
if( length(grep("IDEAS.state", colnames(training.set))) == 1){
	training.set$IDEAS.state <- as.factor(training.set$IDEAS.state)
	table(training.set$IDEAS.state, training.set$regulatory)
}

## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}


## Build h2o random forest model with signed motif damage values
train.set.h2o <- as.h2o(training.set)
# train.set.h2o <- as.h2o(no.dnase.training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
h2o.performance(model_drf, xval = T)
h2o.varimp_plot(model_drf, num_of_features = 30)


train.set.h2o.abs.MDS <- train.set.h2o

## convert damage scores to absolute values
train.set.h2o.abs.MDS$hocomoco.damage.score <- abs(train.set.h2o.abs.MDS$hocomoco.damage.score)
train.set.h2o.abs.MDS$jaspar.damage.score <- abs(train.set.h2o.abs.MDS$jaspar.damage.score)

model_drf_abs_MDS <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o.abs.MDS, model_id = "h2o_drf_abs_MDS", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
h2o.performance(model_drf_abs_MDS, xval = T)
h2o.varimp_plot(model_drf_abs_MDS, num_of_features = 30)




## vs old set full value = 0.837168
## vs new set full value = 0.8079437
## vs new set absolute values = 0.7990063
## vs old set absolute values = 0.8014541
