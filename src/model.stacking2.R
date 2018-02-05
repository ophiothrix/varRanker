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



model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", nfolds = 5, ntrees = 100, max_depth = 10, learn_rate = 0.1, keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3, stopping_tolerance = 0.005, sample_rate = 0.75, col_sample_rate = 0.75)
print(h2o.auc(model_gbm, valid = T))
print(nrow(h2o.scoreHistory(model_gbm)))

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = 5, ntrees = 200, max_depth = 30, keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, stopping_rounds = 3)
print(h2o.auc(model_drf, valid = T))
print(nrow(h2o.scoreHistory(model_drf)))


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
