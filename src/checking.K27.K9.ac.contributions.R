### Background: we are missing H3K27ac and H3K9ac marks for many tissues in encode. Either both or either of them. One way to deal with them would be to use imputed peaks for those marks. To make the beahaviour compatible for tissues, where the marks are available and not, would be to use only binary information from the available marks - exclude all signal information. The trouble with that approach is that the broad peak files seem to be the most useful, at least according to variant impact. But imputed data is only available as narrow or gapped peaks.

### Objective: Using either closely or loosely matched training sets, check the model performance using:
# all available metrics
# no K27ac or K9ac marks
# broad peak binary
# gapped peak binary
# narrow peak binary
# imputed gapped peak
# imputed narrow peak

### Note: A completely different approach would be to re-train the model each time based on the marks available for the set of interest. While somewhat applicable, the problem with that approach is the size of the training set. I.e. if we want to train a model with the maximum number of marks available, we will have relatively small number of training examples. Relaxing marks requirements will increase the size of the training set.

### To check: We can probably achieve quite a good result including only H3K4me1 and me3 marks on top of DNase data. Let's check if we would include more tissues by restricting the requirement to these marks only.

rm(list = ls())
gc()
require(pROC)
require(h2o)
## Initiate H2O instance
localH2O = h2o.init(nthreads=-1)
nfolds <- 5

require(RColorBrewer)
clrs <- brewer.pal(5, "Dark2")

#### Check feature importance on loosely matched set
## make a matrix to store the metric values
modes <- c("all.marks", "no.marks", "broad.bin", "gapped.bin", "narrow.bin", "gapped.imp", "narrow.imp")
metrics <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
sum.table <- matrix(NA, ncol = 7, nrow = length(modes), dimnames = list(modes, metrics))

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds")
head(training.set)
dim(training.set)
summary(training.set$regulatory)
colnames(training.set)
## Remove qval features, as they are redundant with pval
training.set <- training.set[-grep("qval", colnames(training.set))]

train.set.h2o <- as.h2o(training.set)
target <- "regulatory"

### all available metrics
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
## Save a table of variant importance
varimp.loose.match <- h2o.varimp(model_drf)
write.csv(varimp.loose.match, "./reports/variant.importance.RF.loose.match.core.marks.csv")

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["all.marks",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)


# no K27ac or K9ac marks
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["no.marks",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### broad peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("broadPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["broad.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

# gapped peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("gappedPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["gapped.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)


### narrow peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("narrowPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["narrow.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### imputed gapped peak
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("imputed.gappedPeak.bed.gPk.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["gapped.imp",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### imputed narrow peak
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("imputed.narrowPeak.bed.nPk.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["narrow.imp",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

## Save the summary table
write.csv(sum.table, "./reports/effect.of.acetylation.features.loose.match.csv")

##########################################################################################

#### Check feature importance on closely matched set ####

## make a matrix to store the metric values
modes <- c("all.marks", "no.marks", "broad.bin", "gapped.bin", "narrow.bin", "gapped.imp", "narrow.imp")
metrics <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
sum.table <- matrix(NA, ncol = 7, nrow = length(modes), dimnames = list(modes, metrics))

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds")

head(training.set)
dim(training.set)
summary(training.set$regulatory)
colnames(training.set)
## Remove qval features, as they are redundant with pval
training.set <- training.set[-grep("qval", colnames(training.set))]

train.set.h2o <- as.h2o(training.set)
target <- "regulatory"

### all available metrics
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
## Save a table of variant importance
varimp.loose.match <- h2o.varimp(model_drf)
write.csv(varimp.loose.match, "./reports/variant.importance.RF.close.match.core.marks.csv")

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["all.marks",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)


# no K27ac or K9ac marks
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["no.marks",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### broad peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("broadPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["broad.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

# gapped peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("gappedPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["gapped.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)


### narrow peak binary
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("narrowPeak.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["narrow.bin",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### imputed gapped peak
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("imputed.gappedPeak.bed.gPk.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["gapped.imp",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

### imputed narrow peak
features <- setdiff(colnames(train.set.h2o), target)
excl <- c(grep("K27ac", features), grep("K9ac", features))#
excl <- excl[-grep("imputed.narrowPeak.bed.nPk.bin", features[excl])]
features[excl]
features <- features[-excl]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["narrow.imp",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

## Save the summary table
write.csv(sum.table, "./reports/effect.of.acetylation.features.close.match.csv")
