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
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
# training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "../toy.RF.run/cache/E126.annotated.fSkin_fibro.negative.set.rds")

head(training.set)
dim(training.set)

summary(training.set$regulatory)

dim(training.set)
## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}


## Try h2o random forest with DNase
train.set.h2o <- as.h2o(training.set)
# train.set.h2o <- as.h2o(no.dnase.training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.varimp_plot(model_drf, num_of_features = 50)

## Using full damage scores - AUC - 0.7990888
## Using absolute damage scores - AUC - 


## vs old set full value = 0.837168
## vs new set full value = 0.8079437
## vs new set absolute values = 0.7990063
## vs old set absolute values = 0.8014541




## prepare test set
test.set.matched <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")
if( length(grep("IDEAS.state", colnames(training.set))) == 1){
	test.set.matched$IDEAS.state <- as.factor(test.set.matched$IDEAS.state)
	# table(training.set$IDEAS.state, training.set$regulatory)
}
## convert damage scores to absolute values
test.set.matched$hocomoco.damage.score.abs <- abs(test.set.matched$hocomoco.damage.score)
test.set.matched$jaspar.damage.score.abs <- abs(test.set.matched$jaspar.damage.score)

test.set.semi.matched <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/E120.annotated.fMuscle.negative.set.rds")
test.set.unmatched <- prepare.training.set(path.to.positive.set = "./cache/E119.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/E119.annotated.fMuscle.negative.set.rds")


pdf("./graphs/PDS.comparison.figure.pdf")
plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity", main = "Cross-validation and test AUC")
### Plot Cross validation AUC
predicted <- as.data.frame(h2o.cross_validation_holdout_predictions(model_drf))
roc1 <- roc(response = training.set$regulatory, predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[1])


## Plot closely matched test set
test.set.h2o <- as.h2o(test.set.matched)
h2o.performance(model_drf, newdata = test.set.h2o)

predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))
roc2 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[2])

hist(test.set.matched$hocomoco.damage.score, breaks = 100)



## Plot less closely matched test set
test.set.h2o <- as.h2o(test.set.semi.matched)
h2o.performance(model_drf, newdata = test.set.h2o)

predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))
roc3 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[3])

## Plot unmatched test set
test.set.h2o <- as.h2o(test.set.unmatched)
h2o.performance(model_drf, newdata = test.set.h2o)

predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))
roc4 <- roc(response = as.vector(test.set.h2o$regulatory), predictor = predicted$p1, plot = T, smooth = T, add = T, col = clrs[4])

lines(c(1,0), c(0,1), lwd = 2, lty = 2, col = "darkgrey")

legend("bottomright", legend = paste(c("Cross-validation data",  "Closely matched test set", "Loosely matched test set", "Unmatched test set"), round(c(roc1$auc, roc2$auc, roc3$auc, roc4$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")
dev.off()
