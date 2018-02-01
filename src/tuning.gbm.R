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

## Train a gbm model
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
# model_gbm <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
## Specifying a validation frame
# model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 200, max_depth = 40, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3)

model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 100, max_depth = 5, learn_rate = 0.05, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3)

summary(model_gbm)

## Define parameter space. Hold the rest default as one is optimised
max.depths <- c(2, 5, 7, 10, 12, 15, 20, 25, 30)
learn.rates <- seq(0.13, 0.18, 0.01)
ntreess <- c(5, 10, 20, 50, 60, 80, 100)
max.depths <- seq(20, 30, 2)


max.depth.auc <- numeric()
for (max.depth in max.depths) {
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 40, max_depth = 10, learn_rate = 0.17, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3, stopping_metric = "AUC", stopping_tolerance = 0.001)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	max.depth.auc <- c(max.depth.auc, h2o.auc(model_gbm, valid = T))
	
}
names(max.depth.auc) <- max.depths
max.depth.auc
## At learning rate of 0.1 !! not the same ntrees !!
# 2         5        10        15        20 
# 0.8587870 0.8728655 0.8921220 0.8885545 0.8908820 


learn.rate.auc <- numeric()
for (learn.rate in learn.rates) {
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 50, max_depth = 10, learn_rate = learn.rate, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)#, stopping_rounds = 5, stopping_tolerance = 0.005)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	learn.rate.auc <- c(learn.rate.auc, h2o.auc(model_gbm, valid = T))
}
names(learn.rate.auc) <- learn.rates
plot(learn.rate.auc)

ntrees.auc <- numeric()
for (ntrees in ntreess) {
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = ntrees, max_depth = 10, learn_rate = 0.2,  keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)#, stopping_rounds = 3, stopping_tolerance = 0.01)
	print(h2o.auc(model_gbm, valid = T))
	ntrees.auc <- c(ntrees.auc, h2o.auc(model_gbm, valid = T))
	
}
names(ntrees.auc) <- ntreess


model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 100, max_depth = 10, learn_rate = 0.05, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5, stopping_tolerance = 0.005, sample_rate = 0.75, col_sample_rate = 0.75)
print(h2o.auc(model_gbm, valid = T))
print(nrow(h2o.scoreHistory(model_gbm)))
