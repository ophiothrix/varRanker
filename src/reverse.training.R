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
training.set <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
test.set.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
test.set.h2o <- as.h2o(test.set.matched)

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

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Stratified", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
plot(h2o.performance(model_drf, xval = T))
h2o.performance(model_drf, train = T)
h2o.performance(model_drf, newdata = test.set.h2o)
pdf("./graphs/feature.importance.pdf", 8, 6)
h2o.varimp_plot(model_drf, num_of_features = 30)
dev.off()

## prepare test set
test.set.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.HSMM.negative.set.rds")
