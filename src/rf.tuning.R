## Setting up
rm(list = ls())
gc()
require(pROC)
require(h2o)
require(RColorBrewer)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
nfolds <- 5
clrs <- brewer.pal(5, "Dark2")

### Prepare the training set ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)

## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}

train.set.h2o <- as.h2o(training.set)
# train.set.h2o <- as.h2o(no.dnase.training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)
## Remove cadd score from training features
features <- features[features != "cadd.score"]

## Optimise number of trees
auc.xval <- auc.train <- numeric()
ntrs <- c(5, 10, 20, 50, 100, 150, 200)
for (ntr in ntrs) {
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = ntr, max_depth = 10, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	auc.xval <- c(auc.xval, h2o.auc(h2o.performance(model_drf, xval = T)))
	predicted <- as.data.frame(h2o.predict(model_drf, train.set.h2o))
	auc.train <- c(auc.train, as.numeric(roc(response = training.set$regulatory, predictor = predicted$p1)$auc))
}

plot(1:length(ntrs), auc.xval, ylim = range(c(auc.train, auc.xval)), pch = 16, col = "darkorange", xaxt = "n", xlab = "Number of trees", ylab = "AUC")
points(1:length(ntrs), auc.train, pch = 16, col = "forestgreen")
axis(1, 1:length(ntrs), ntrs)
legend("bottomright", legend = c("Cross-validation", "Training"), fill = c("darkorange", "forestgreen"), bty = "n")

## Optimise  tree depth

auc.xval <- auc.train <- numeric()
ntrs <- c(5, 10, 20, 50, 100)
for (ntr in ntrs) {
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 50, max_depth = ntr, fold_assignment = "Stratified", keep_cross_validation_predictions = TRUE, seed = 16052017)
	auc.xval <- c(auc.xval, h2o.auc(h2o.performance(model_drf, xval = T)))
	predicted <- as.data.frame(h2o.predict(model_drf, train.set.h2o))
	auc.train <- c(auc.train, as.numeric(roc(response = training.set$regulatory, predictor = predicted$p1)$auc))
}

plot(1:length(ntrs), auc.xval, ylim = range(c(auc.train, auc.xval)), pch = 16, col = "darkorange", xaxt = "n", xlab = "Number of trees", ylab = "AUC")
points(1:length(ntrs), auc.train, pch = 16, col = "forestgreen")
axis(1, 1:length(ntrs), ntrs)
legend("bottomright", legend = c("Cross-validation", "Training"), fill = c("darkorange", "forestgreen"), bty = "n")

summary(model_drf)
?h2o.randomForest
