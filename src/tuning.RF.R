## Prepare the sets
rm(list = ls())
gc()
source("./src/prep.sets.R")
sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", path.to.test.sets = "./cache/", test.set.tissues = c("HSMM"))
names(sets.list)

train.set <- sets.list$training
valid.set <- sets.list$validation

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

table(train.set$source)

## With cross-validation
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 25, max_depth = 8, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)

## Specifying a validation frame
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = 100, max_depth = 20, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5)

h2o.auc(model_drf, xval = T)
h2o.auc(model_drf, valid = T)
h2o.auc(model_drf, train = T)

####################################
### Optimise for number of trees ###
####################################
md <- 10
ntrees <- c(30, 40, 50, 60, 70, 80, 100)#, 150)
# ntrees <- seq(21, 39, 2)
ntrees.Niter <- ntrees.auc <- ntrees.auc.train <- rep(NA, length(ntrees))
names(ntrees.auc) <- names(ntrees.auc.train) <- names(ntrees.Niter) <- ntrees
for (nt in ntrees) {
	print(nt)
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)#, stopping_rounds = 5, stopping_metric = "AUC", stopping_tolerance = 0.0001)
	ntrees.auc[as.character(nt)] <- h2o.auc(model_drf, valid = T)
	ntrees.auc.train[as.character(nt)] <- h2o.auc(model_drf, train = T)
	ntrees.Niter[as.character(nt)] <- nrow(h2o.scoreHistory(model_drf))
}
rbind(ntrees.auc.train, ntrees.auc, ntrees.Niter)
plot(ntrees.auc~names(ntrees.auc), ylim = range(c(ntrees.auc.train, ntrees.auc)))
points(ntrees.auc.train~names(ntrees.auc.train), col = "red")
plot(ntrees.auc~ntrees.Niter)
## Optimisation stops around 70 trees. Might be worth exploring the space between 60 and 80

##################################
### Optimise for maximum depth ###
##################################
nt <- 100
mdepth <- seq(4, 12, 1)
mdepth.Niter <- mdepth.auc <- mdepth.auc.train <- rep(NA, length(mdepth))
names(mdepth.auc) <- names(mdepth.auc.train) <- names(mdepth.Niter) <- mdepth
for (md in mdepth) {
	print(md)
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5, stopping_metric = "AUC", stopping_tolerance = 0.00001)
	mdepth.auc.train[as.character(md)] <- h2o.auc(model_drf, train = T)
	mdepth.auc[as.character(md)] <- h2o.auc(model_drf, valid = T)
	mdepth.Niter[as.character(md)] <- nrow(h2o.scoreHistory(model_drf))
	print(mdepth.auc[as.character(md)])
	print(mdepth.auc.train[as.character(md)])
	print(mdepth.Niter[as.character(md)])
}
rbind(mdepth.auc.train, mdepth.auc, mdepth.Niter)
plot(mdepth.auc~names(mdepth.auc), ylim = range(c(mdepth.auc.train, mdepth.auc)))
points(mdepth.auc.train~names(mdepth.auc.train), col = "red")
plot(mdepth.auc~mdepth.Niter)
## AUC maxes out at max depth between 10 and 20

## Search grid
sg <- matrix(NA, nrow = 6, ncol = 15)
colnames(sg) <- c(seq(36, 50, 1))#, 100, 120, 150)
rownames(sg) <- seq(5, 10, 1)

for (nt in as.numeric(colnames(sg))) {
	print(nt)
	for (md in as.numeric(rownames(sg))) {
		print(md)
		model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 0, ntrees = nt, max_depth = md, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)
		sg[as.character(md), as.character(nt)] <- h2o.auc(model_drf, valid = T)
		print(sg[as.character(md), as.character(nt)])
	}
}

sg
heatmap(sg, scale = "none", Rowv = NA, Colv = NA)
which.max(sg)
max(sg)
# top is at md=7 nt=100. But it's almost as good at md=6 nt=20.

###################
### Final model ###
###################
## ntrees = 71
	nt = 40
	## max_depth = 10
	md = 8
	
	## With cross-validation
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = nt, max_depth = md, fold_assignment = "Modulo", keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)
	
	
	h2o.auc(model_drf, xval = T)
	h2o.auc(model_drf, valid = T)
	h2o.auc(model_drf, train = T)
