## Prepare the sets
rm(list = ls())
gc()
source("./src/prep.sets.R")
# sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", test.set.tissues = c("fMuscle", "HSMM"))
sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.k27.annotation.rds", path.to.test.sets = "./cache", test.set.tissues = "HSMM")
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

##########################
### Random fores model ###
##########################
## ntrees = 71
nt = 40
## max_depth = 10
md = 8

## With cross-validation
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = nt, max_depth = md, fold_assignment = "Modulo", keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)


h2o.auc(model_drf, train = T)
h2o.auc(model_drf, xval = T)
h2o.auc(model_drf, valid = T)

#################
### GBM model ###
#################
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


model_ensemble <- h2o.stackedEnsemble(x = features, y = target,
					training_frame = train.set.h2o,
					validation_frame = valid.set.h2o, 
					model_id = "h2o_ensemble4", 
					base_models = list(model_gbm@model_id, model_drf@model_id))


h2o.auc(h2o.performance(model_drf, newdata = valid.set.h2o, valid = T))
h2o.auc(h2o.performance(model_gbm, newdata = valid.set.h2o, valid = T))
h2o.auc(h2o.performance(model_ensemble, newdata = valid.set.h2o, valid = T))

h2o.auc(h2o.performance(model_drf, xval = T))
h2o.auc(h2o.performance(model_gbm, xval = T))
h2o.auc(h2o.performance(model_ensemble, xval = T))



for (l in 1:(length(names(sets.list)) - 2)) {
	print(names(sets.list)[l])
	test.set <- as.h2o(sets.list[[l]])
	print(h2o.auc(h2o.performance(model_drf, newdata = test.set, valid = T)))
	print(h2o.auc(h2o.performance(model_gbm, newdata = test.set, valid = T)))
	print(h2o.auc(h2o.performance(model_ensemble, newdata = test.set, valid = T)))
}


test.perfect <- as.h2o(readRDS("./cache/test.perfect3.rds"))
h2o.auc(h2o.performance(model_drf, newdata = test.perfect, valid = T))
h2o.auc(h2o.performance(model_gbm, newdata = test.perfect, valid = T))
h2o.auc(h2o.performance(model_ensemble, newdata = test.perfect, valid = T))


test.imperfect <- as.h2o(readRDS("./cache/test.imperfect4.rds"))
h2o.auc(h2o.performance(model_drf, newdata = test.imperfect, valid = T))
h2o.auc(h2o.performance(model_gbm, newdata = test.imperfect, valid = T))
h2o.auc(h2o.performance(model_ensemble, newdata = test.imperfect, valid = T))
#1 - 0.805871
#2 - 0.8128859
#3 - 0.8215855
####4 - 0.8222107
#5 - 0.7833855

test.mismatched <- as.h2o(readRDS("./cache/test.mismatched4.rds"))
h2o.auc(h2o.performance(model_drf, newdata = test.mismatched, valid = T))
h2o.auc(h2o.performance(model_gbm, newdata = test.mismatched, valid = T))
h2o.auc(h2o.performance(model_ensemble, newdata = test.mismatched, valid = T))
#1 - 0.6695983
#2 - 0.8105716
#3 - 0.6045804
#4 - 0.6073204
####5 - 0.5723777


test.cadd <- as.h2o(readRDS("./cache/test.cadd3.rds"))
h2o.auc(h2o.performance(model_drf, newdata = test.cadd, valid = T))
h2o.auc(h2o.performance(model_gbm, newdata = test.cadd, valid = T))
h2o.auc(h2o.performance(model_ensemble, newdata = test.cadd, valid = T))
#1 - 0.5723777
#2 - 0.8128859
#3 - 0.8222107
#4 - 
####5 - 


plot(h2o.performance(model_drf, newdata = test.perfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
plot(h2o.performance(model_gbm, newdata = test.perfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
plot(h2o.performance(model_ensemble, newdata = test.perfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
abline(v = c(0, 0.05), col = "#666666", lwd = 2, lty = 2)

plot(h2o.performance(model_drf, newdata = test.imperfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
plot(h2o.performance(model_gbm, newdata = test.imperfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
plot(h2o.performance(model_ensemble, newdata = test.imperfect, valid = T), col = "maroon", pch = 16, cex = 0.5)
abline(v = c(0, 0.05), col = "#666666", lwd = 2, lty = 2)

plot(h2o.performance(model_lreg, newdata = test.imperfect, valid = T), col = "maroon", pch = 16, cex = 0.5)


## Questions:
# Should I use a balanced training set?
# Should I play with ratios of negative variants
# 





## Predict the response
pdf("./graphs/test.vs.cadd.score.comparison.pdf")
plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity")#, main = "Test set, closely matched")
## Generate ROC curve for validation set predictions
predicted <- as.data.frame(h2o.predict(model_ensemble, valid.set.h2o))
roc1 <- roc(response = as.vector(valid.set.h2o$regulatory), predictor = as.vector(predicted$p1), plot = T, smooth = F, add = T, col = clrs[1])

## Generate ROC curve for correctly annotated test set
test.cadd <- as.h2o(readRDS("./cache/test.cadd3.rds"))
predicted <- as.data.frame(h2o.predict(model_ensemble, test.cadd))
roc2 <- roc(response = as.vector(test.cadd$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[2])

## Generate ROC curve for mismatched test set
test.cadd <- as.h2o(readRDS("./cache/test.cadd1.rds"))
predicted <- as.data.frame(h2o.predict(model_ensemble, test.cadd))
roc3 <- roc(response = as.vector(test.cadd$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[3])

## Generate ROC curve for CADD score predictions
roc4 <- roc(response = as.vector(test.cadd$regulatory), predictor = as.vector(test.cadd$cadd.score), plot = T, smooth = F, add = T, col = clrs[4])

abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Validation set",  "Matched test set", "Mismatched test set", "CADD score"), round(c(roc1$auc, roc2$auc, roc3$auc, roc4$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")
dev.off()

roc4$auc
roc2$auc
par(mfrow = c(1, 2))



plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity")#, main = "Test set, closely matched")
## Generate ROC curve for validation set predictions
predicted <- as.data.frame(h2o.predict(model_ensemble, valid.set.h2o))
roc1 <- roc(response = as.vector(valid.set.h2o$regulatory), predictor = as.vector(predicted$p1), plot = T, smooth = F, add = T, col = clrs[1])

## Generate ROC curve for correctly annotated test set
predicted <- as.data.frame(h2o.predict(model_ensemble, test.perfect))
roc2 <- roc(response = as.vector(test.perfect$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[2])

## Generate ROC curve for loosely matched test set
predicted <- as.data.frame(h2o.predict(model_ensemble, test.imperfect))
roc3 <- roc(response = as.vector(test.imperfect$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[3])

## Generate ROC curve for mismatched test set
predicted <- as.data.frame(h2o.predict(model_ensemble, test.mismatched))
roc4 <- roc(response = as.vector(test.mismatched$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[4])

abline(a = 1, b = -1, col = "darkgrey", lty = 2, lwd = 2)
legend("bottomright", legend = paste(c("Validation set",  "Well matched test set", "Matched test set", "Mismatched test set"), round(c(roc1$auc, roc2$auc, roc3$auc, roc4$auc), 3), sep = "; AUC="), fill = clrs[1:5], bty = "n")



#
#
#
#
#
#

source("./src/prep.sets.R")
# sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", test.set.tissues = c("fMuscle", "HSMM"))
sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", test.set.tissues = "fMuscle")
sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", path.to.test.sets = "./cache/funseq", test.set.tissues = "HSMM")

names(sets.list)


## Generate ROC curve for CADD score predictions
test.cato <- as.h2o(sets.list[[1]])
test.cato <- as.h2o(sets.list[[4]])
test.cato <- as.h2o(valid.set.cato)


plot(1, 1, xlim = c(1, 0), ylim = c(0,1), type = "n", ylab = "Sensitivity", xlab = "Specificity")#, main = "Test set, closely matched")

predicted <- as.data.frame(h2o.predict(model_ensemble, test.cato))
roc1 <- roc(response = as.vector(test.cato$regulatory), predictor = predicted$p1, plot = T, smooth = F, add = T, col = clrs[1])
roc1$auc

roc2 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$cadd.score), plot = T, smooth = F, add = T, col = clrs[2])
roc2$auc

## funseq
roc3 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$funseq.score), plot = T, smooth = F, add = T, col = clrs[3])
roc3$auc

## eigen
roc3 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$eigen.score), plot = T, smooth = F, add = T, col = clrs[3])
roc3$auc

## gwava.random
roc3 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$gwava.random), plot = T, smooth = F, add = T, col = clrs[3])
roc3$auc

## gwava.tss
roc4 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$gwava.tss), plot = T, smooth = F, add = T, col = clrs[4])
roc4$auc

## gwava.variant
roc5 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$gwava.variant), plot = T, smooth = F, add = T, col = clrs[5])
roc5$auc

table(as.vector(test.cato$regulatory))

plot(roc1)


roc2 <- roc(response = as.vector(test.cato$regulatory), predictor = as.vector(test.cato$cato.score), plot = T, smooth = F, add = T, col = clrs[2])
roc2$auc

roc2 <- roc(response = as.vector(test.cato$regulatory), predictor = runif(nrow(test.cato)), plot = T, smooth = F, add = T, col = "darkgrey")
roc2$auc


predicted <- h2o.predict(model_ensemble, valid.set.h2o)
predicted$regulatory <- valid.set.h2o$regulatory
predicted$source <- valid.set.h2o$source
predicted <- as.data.frame(predicted)
predicted$predict[predicted$p1 < 0.7] <- 0
head(predicted)
predicted %>%
	group_by(source) %>%
	summarise(wrong = sum(predict != regulatory), right = sum(predict == regulatory), total = sum(predict != regulatory) + sum(predict == regulatory)) %>%
	mutate(error.rate = wrong / total)
min(predicted$p1[predicted$predict == 1])

test.set <- as.h2o(sets.list$E120.annotated.HSMM.variants.rds)
plot(h2o.performance(model_ensemble, test.set))

table(test.set$source)
plot(h2o.performance(model_ensemble, valid.set.h2o))
