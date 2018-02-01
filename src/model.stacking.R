require(pROC)
require(RColorBrewer)
rm(list = ls())
clrs <- brewer.pal(5, "Dark2")
nfolds <- 5

source("./lib/prepare.training.set.R")
set.seed(15052017)
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)

## prepare test set
test.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")


## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}


## initiate H2O
require(h2o)
localH2O = h2o.init(nthreads=4)

## Convert dataset
train.set.h2o <- as.h2o(training.set)
test.set.h2o <- as.h2o(test.set)
# train.set.h2o <- as.h2o(training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

### Build a random forest model ###
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)

### Build gradient boosting machine model ###
model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, distribution = "bernoulli", ntrees = 100, max_depth = 5, min_rows = 2, learn_rate = 0.2, nfolds = nfolds, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_gbm)
h2o.performance(model_gbm, xval = T)
h2o.varimp_plot(model_gbm, num_of_features = 50)


## Train ensemble model
model_ensemble <- h2o.stackedEnsemble(x = features, y = target,
								training_frame = train.set.h2o,
								validation_frame = test.set.h2o,
								model_id = "h2o_ensemble",
								base_models = list(model_gbm@model_id, model_drf@model_id))
summary(model_ensemble)
h2o.performance(model_ensemble, newdata = test.set.h2o, valid = T)
h2o.performance(model_drf, newdata = test.set.h2o, valid = T)
h2o.performance(model_gbm, newdata = test.set.h2o, valid = T)

h2o.performance(model_ensemble, xval = T)
h2o.performance(model_drf, xval = T)
h2o.performance(model_gbm, xval = T)

h2o.varimp_plot(model_drf, num_of_features = 50)
h2o.varimp_plot(model_gbm, num_of_features = 50)


plot(model_dr)
h2o.auc(model_drf, xval = T, train = T)
h2o.auc(model_gbm, xval = T, train = T)
h2o.auc(model_ensemble, train = T, valid = T)
h2o.varimp_plot(model_drf)
h2o.cross_validation_predictions(model_drf)
h2o.specificity(model_drf, thresholds = 0.1)
perf <- h2o.performance(model_drf, newdata = test.set.h2o)
h2o.specificity(perf, 0.1)

test.set.h2o$regulatory

table(predicted$predict)
table(as.vector(test.set.h2o$regulatory))

plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Model Ensembling")
predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])
predicted <- as.data.frame(h2o.predict(model_gbm, test.set.h2o))
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[2])
predicted <- as.data.frame(h2o.predict(model_ensemble, test.set.h2o))
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[3])

ntrees <- c(10, 50, 100, 500, 1000)
aucs <- numeric(5)

# i <- 1
for (i in 1:5){
	print(ntrees[i])
	rf2 <- randomForest(as.factor(regulatory) ~ ., data = no.dnase.training.set, importance = T, ntree = ntrees[i])
	roc1 <- roc(response = rf2$y, predictor = rf2$votes[,2], plot = T, smooth = T, add = T, col = clrs[i])
	aucs[i] <- round(roc1$auc, 2)
}
legend("bottomright", legend = paste0(ntrees, " trees; AUC=", aucs), fill = clrs[1:5], bty = "n")
dev.off()

### Plot ROC curves for individual models and predictor
as.data.frame(h2o.cross_validation_holdout_predictions(model_drf)
h2o.cross_validation_holdout_predictions
