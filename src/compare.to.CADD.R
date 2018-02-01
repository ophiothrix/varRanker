#### Compare predictions by the classifier and CADD score

### Take two similar, ideally the same, test sets. Run the predictive model to generate one ranking. Then also rank them by CADD score and see which one performs better.


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

### Prepare the input data ###
source("./lib/prepare.training.set.R")
# training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")
training.set <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
# training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "../toy.RF.run/cache/E126.annotated.fSkin_fibro.negative.set.rds")

head(training.set)
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

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
pdf("./graphs/feature.importance.pdf", 8, 6)
h2o.varimp_plot(model_drf, num_of_features = 30)
dev.off()

h2o.auc(model_drf, xval = T) # 0.7752843 with cadd; 0.7729757 without

## prepare test set
test.set.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.negative.set.rds")
test.set.semi.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.fMuscle.negative.set.rds")
test.set.unmatched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E119.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E119.annotated.fMuscle.negative.set.rds")


pdf("./graphs/comparison.with.cadd.score.figure.pdf")
plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Training (Cross-validation) set")
### Plot Cross validation AUC
predicted <- as.data.frame(h2o.cross_validation_holdout_predictions(model_drf))
roc1 <- roc(response = training.set$regulatory, predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])
roc2 <- roc(response = as.vector(train.set.h2o$regulatory), predictor = as.vector(train.set.h2o$cadd.score), plot = T, smooth = T, add = T, col = clrs[2])
abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Classifier",  "CADD score"), round(c(roc1$auc, roc2$auc), 3), sep = "; AUC="), fill = clrs[1:2], bty = "n")

## Plot closely matched test set
test.set.h2o <- as.h2o(test.set.matched)
## Predict the response
predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))

plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Test set, closely matched")
## Generate ROC curve for classifier predictions
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])
## Generate ROC curve for CADD score predictions
roc2 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = as.vector(test.set.h2o$cadd.score), plot = T, smooth = T, add = T, col = clrs[2])
abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Classifier",  "CADD score"), round(c(roc1$auc, roc2$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")




## Plot less closely matched test set
test.set.h2o <- as.h2o(test.set.semi.matched)
## Predict the response
predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))

plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Test set, loosely matched")
## Generate ROC curve for classifier predictions
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])
## Generate ROC curve for CADD score predictions
roc2 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = as.vector(test.set.h2o$cadd.score), plot = T, smooth = T, add = T, col = clrs[2])
abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Classifier",  "CADD score"), round(c(roc1$auc, roc2$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")


## Plot unmatched test set
test.set.h2o <- as.h2o(test.set.unmatched)

## Predict the response
predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))

plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Test set, mis-matched")
## Generate ROC curve for classifier predictions
roc1 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])
## Generate ROC curve for CADD score predictions
roc2 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = as.vector(test.set.h2o$cadd.score), plot = T, smooth = T, add = T, col = clrs[2])
abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Classifier",  "CADD score"), round(c(roc1$auc, roc2$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")

dev.off()


