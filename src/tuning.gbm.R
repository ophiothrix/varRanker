## Prepare the sets
rm(list = ls())
gc()
source("./src/prep.sets.R")
sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", path.to.test.sets = "./cache", test.set.tissues = c("fMuscle", "HSMM"))
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


## With cross-validation
# model_gbm <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
## Specifying a validation frame
# model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 200, max_depth = 40, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3)

model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 30, max_depth = 5, learn_rate = 0.05, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3)

summary(model_gbm)
h2o.auc(model_gbm, valid = T)

## Define parameter space. Hold the rest default as one is optimised
ntreess <- c(5, 10, 20, 50, 60, 80, 100)
max.depths <- seq(20, 30, 2)

##########################
### Optimise max depth ###
##########################
nt <- 100
lr <- 0.1
sr <- 0.9 # variant sample rate
cr <- 0.9 # feature sample rate
max.depths <- c(2, 5, 7, 10, 12, 15, 20, 25, 30)
max.depths <- seq(2, 12, 1)
mdepth.Niter <- mdepth.auc <- rep(NA, length(max.depths))
names(mdepth.Niter) <- names(mdepth.auc) <- max.depths
for (md in max.depths) {
	print(md)
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = nt, max_depth = md, learn_rate = lr, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr, stopping_metric = "AUC", stopping_rounds = 3, stopping_tolerance = 0.001)
	plot(h2o.scoreHistory(model_gbm)$validation_auc, ylim = c(0.7, 0.77), main = md)
	abline(h = 0.765)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	mdepth.auc[as.character(md)] <- h2o.auc(model_gbm, valid = T)
	mdepth.Niter[as.character(md)] <- nrow(h2o.scoreHistory(model_gbm))
}

md <- 5
rbind(mdepth.auc, mdepth.Niter)
plot(mdepth.auc~names(mdepth.auc))
plot(mdepth.auc~mdepth.Niter)

### AUC maxes out between depth 5 and 10. Note that the number or trees is not the same for all values of md because of early stopping. But it holds if forcing n trees to 50.
## Check the space between max depth 5 and 9
## At learning rate of 0.1 !! not the same ntrees !!
# 2         5        10        15        20 
# 0.8587870 0.8728655 0.8921220 0.8885545 0.8908820 



##############################
### Optimise learning rate ###
##############################
## Set a stopping criteria to get an idea of a good ntrees
nt <- 100
md <- 4
sr <- 0.9 # variant sample rate
cr <- 0.9 # feature sample rate
learn.rates <- seq(0.1, 0.18, 0.01)
lrate.Niter <- lrate.auc <- rep(NA, length(learn.rates))
names(lrate.Niter) <- names(lrate.auc) <- learn.rates
for (lr in learn.rates) {
	print(lr)
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = nt, max_depth = md, learn_rate = lr, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr, stopping_metric = "AUC", stopping_rounds = 3, stopping_tolerance = 0.001)
	plot(h2o.scoreHistory(model_gbm)$validation_auc, ylim = c(0.7, 0.77), main = lr)
	abline(h = 0.765)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	lrate.auc[as.character(lr)] <- h2o.auc(model_gbm, valid = T)
	lrate.Niter[as.character(lr)] <- nrow(h2o.scoreHistory(model_gbm))
}

rbind(lrate.auc, lrate.Niter)
plot(lrate.auc~names(lrate.auc))
plot(lrate.auc~lrate.Niter)

## The AUC maxes out at around 0.89 across a wide range of learning rates. What it seems to affect most (predictably) is the number of trees required. LR from 0.05 to 0.1 reaches AUC of about 0.89 at between 40 & 50 trees. Stick to LR of 0.1 and explore the ntrees between 40 and 50.



##############################
### Optimise ntrees ###
##############################
## Do not set a stopping criteria
nt <- 100
md <- 4
lr <- 0.16
sr <- 0.9 # variant sample rate
cr <- 0.9 # feature sample rate
ntrees <- seq(34, 75, 3)
ntree.Niter <- ntree.auc <- rep(NA, length(ntrees))
names(ntree.Niter) <- names(ntree.auc) <- ntrees
for (nt in ntrees) {
	print(nt)
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = nt, max_depth = md, learn_rate = lr, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr)
	plot(h2o.scoreHistory(model_gbm)$validation_auc, ylim = c(0.7, 0.77), main = lr)
	abline(h = 0.765)
	abline(v = 13)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	ntree.auc[as.character(nt)] <- h2o.auc(model_gbm, valid = T)
	ntree.Niter[as.character(nt)] <- nrow(h2o.scoreHistory(model_gbm))
}

rbind(ntree.auc, ntree.Niter)
plot(ntree.auc~names(ntree.auc))
plot(ntree.auc~ntree.Niter)

## AUC fluctuates around 0.89 and maxes out at ntrees = 43

## Let's see if we can improve it by subsampling
nt <- 70
md <- 4
lr <- 0.1

srates <- seq(0.7, 1, 0.05)
srate.Niter <- srate.auc <- rep(NA, length(srates))
names(srate.Niter) <- names(srate.auc) <- srates
for (sr in srates) {
	print(sr)
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = nt, max_depth = md, learn_rate = lr, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = 0.9, col_sample_rate = sr, stopping_metric = "AUC", stopping_rounds = 5, stopping_tolerance = 0.0001)
	print(h2o.auc(model_gbm, valid = T))
	print(nrow(h2o.scoreHistory(model_gbm)))
	srate.auc[as.character(sr)] <- h2o.auc(model_gbm, valid = T)
	srate.Niter[as.character(sr)] <- nrow(h2o.scoreHistory(model_gbm))
}

rbind(srate.auc, srate.Niter)
plot(srate.auc~names(srate.auc))
plot(srate.auc~srate.Niter)
points(srate.auc~names(srate.auc), col = "red")

## Optimise column and row sample rates simultaneously
## Search grid
sg <- matrix(NA, nrow = 3, ncol = 3)
colnames(sg) <- rownames(sg) <- seq(0.89, 0.91, 0.01)

for (sr in as.numeric(colnames(sg))) {
	print(sr)
	for (cr in as.numeric(rownames(sg))) {
		print(cr)
		model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = nt, max_depth = md, learn_rate = lr, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr, stopping_metric = "AUC", stopping_rounds = 5, stopping_tolerance = 0.0001)
		sg[as.character(cr), as.character(sr)] <- h2o.auc(model_gbm, valid = T)
		print(sg[as.character(cr), as.character(sr)])
	}
}

round(sg, 4)
heatmap(sg, scale = "none", Rowv = NA, Colv = NA)
which.max(sg)
max(sg)


###################
### Final model ###
###################
## ntrees = 43
nt = 13
## max_depth = 7
md = 4
## learn_rate = 0.1
lr = 0.16
# variant sample rate
sr <- 0.9
# feature sample rate
cr <- 0.9

model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", ntrees = nt, max_depth = md, learn_rate = lr, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr, nfolds = nfolds, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE)

h2o.auc(model_gbm, train = T)
h2o.auc(model_gbm, xval = T)
h2o.auc(model_gbm, valid = T)


for (l in 1:(length(names(sets.list)) - 2)) {
	print(names(sets.list)[l])
	test.set <- as.h2o(sets.list[[l]])
	print(h2o.auc(h2o.performance(model_gbm, newdata = test.set, valid = T)))
}

## Starting model (last)
# > h2o.auc(model_gbm, train = T)
# [1] 0.8843931
# > h2o.auc(model_gbm, xval = T)
# [1] 0.7740912
# > h2o.auc(model_gbm, valid = T)
# [1] 0.7575425

# model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 0, ntrees = 100, max_depth = 10, learn_rate = 0.05, keep_cross_validation_predictions = F, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 5, stopping_tolerance = 0.005, sample_rate = 0.75, col_sample_rate = 0.75)
# print(h2o.auc(model_gbm, valid = T))
# print(nrow(h2o.scoreHistory(model_gbm)))
