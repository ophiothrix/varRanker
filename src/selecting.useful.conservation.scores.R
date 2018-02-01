rm(list = ls())
gc()
require(pROC)
require(h2o)
## Initiate H2O instance
localH2O = h2o.init(nthreads=-1)
nfolds <- 5

require(RColorBrewer)
clrs <- brewer.pal(5, "Dark2")

#### Check variant importance on loosely matched set

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds")
head(training.set)
dim(training.set)
summary(training.set$regulatory)
colnames(training.set)
## Remove qval features, as they are redundant with pval
training.set <- training.set[-grep("qval", colnames(training.set))]

## Run h2o random forest
train.set.h2o <- as.h2o(training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.varimp_plot(model_drf, num_of_features = 30)

varimp.loose.match <- h2o.varimp(model_drf)
write.csv(varimp.loose.match, "./reports/variant.importance.RF.loose.match.csv")

#### Check variant importance on loosely matched set

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")
head(training.set)
dim(training.set)
summary(training.set$regulatory)
colnames(training.set)
## Remove qval features, as they are redundant with pval
training.set <- training.set[-grep("qval", colnames(training.set))]

## Run h2o random forest
train.set.h2o <- as.h2o(training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.varimp_plot(model_drf, num_of_features = 30)

varimp.close.match <- h2o.varimp(model_drf)
write.csv(varimp.close.match, "./reports/variant.importance.RF.close.match.csv")
