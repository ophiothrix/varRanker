## Prepare the sets
rm(list = ls())
gc()
source("./src/prep.sets.R")



## Train a random forest model
require(pROC)
require(h2o)
require(RColorBrewer)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
nfolds <- 5
clrs <- brewer.pal(5, "Dark2")



train.set.h2o <- as.h2o(train.set)
valid.set.h2o <- as.h2o(valid.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)
## Remove qval features, as they are redundant with pval
features <- features[-grep("qval", features)]
## remove auxiliary columns from features
features <- features[!features %in% c("varID", "source")]

## With cross-validation
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, max_depth = 20, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)

## Specifying a validation frame
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = 100, max_depth = 20, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)

h2o.auc(model_drf, xval = T)
h2o.auc(model_drf, valid = T)
h2o.auc(model_drf, train = T)

####################################
### Optimise for number of trees ###
####################################
md <- 20
ntrees <- c(5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150)
ntrees.Niter <- ntrees.auc <- rep(NA, length(ntrees))
names(ntrees.auc) <- names(ntrees.Niter) <- ntrees
for (nt in ntrees) {
	print(nt)
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)
	ntrees.auc[as.character(nt)] <- h2o.auc(model_drf, valid = T)
	ntrees.Niter[as.character(nt)] <- nrow(h2o.scoreHistory(model_drf))
}
rbind(ntrees.auc, ntrees.Niter)
plot(ntrees.auc~names(ntrees.auc))
plot(ntrees.auc~ntrees.Niter)
## Optimisation stops around 70 trees. Might be worth exploring the space between 60 and 80

##################################
### Optimise for maximum depth ###
##################################
nt <- 100
mdepth <- seq(5, 50, 5)
mdepth.Niter <- mdepth.auc <- rep(NA, length(mdepth))
names(mdepth.auc) <- names(mdepth.Niter) <- mdepth
for (md in mdepth) {
	print(md)
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5, stopping_tolerance = 0.0001)
	mdepth.auc[as.character(md)] <- h2o.auc(model_drf, valid = T)
	mdepth.Niter[as.character(md)] <- nrow(h2o.scoreHistory(model_drf))
}
rbind(mdepth.auc, mdepth.Niter)
plot(mdepth.auc~names(mdepth.auc))
plot(mdepth.auc~mdepth.Niter)
## AUC maxes out at max depth between 10 and 20

## Search grid
sg <- matrix(NA, 3, 7)
colnames(sg) <- seq(67, 73, 1)
rownames(sg) <- seq(9, 11, 1)

for (nt in as.numeric(colnames(sg))) {
	print(nt)
	for (md in as.numeric(rownames(sg))) {
		print(md)
		model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)
		sg[as.character(md), as.character(nt)] <- h2o.auc(model_drf, valid = T)
	}
}

sg
heatmap(sg, scale = "none", Rowv = NA, Colv = NA)

###################
### Final model ###
###################
## ntrees = 71
nt = 71
## max_depth = 10
md = 10

## With cross-validation
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = nt, max_depth = md, fold_assignment = "Modulo", keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)


h2o.auc(model_drf, xval = T)
h2o.auc(model_drf, valid = T)
h2o.auc(model_drf, train = T)
