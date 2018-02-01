#### Train a model with the bare minimum of marks:
# - DNase
# - H3K4me1
# - H3K4me3
par(mfrow = c(1, 2))

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
modes <- c("all.marks", "bare.minimum")
metrics <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
sum.table <- matrix(NA, ncol = 7, nrow = length(modes), dimnames = list(modes, metrics))

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds", prop.dnase = 3)
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

plot(h2o.performance(model_drf, xval = T))

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


# bare minimum marks
features <- setdiff(colnames(train.set.h2o), target)
## Find all epigenetic marks
epimarks <- c(grep("broad", features), grep("narrow", features), grep("gapped", features))
## Find the marks to retain
keepers <- unique(c(grep("DNase", features[epimarks]),
			 grep("H3K4me3", features[epimarks]),
			 grep("H3K4me1", features[epimarks])))
features[epimarks][keepers]
features <- features[-epimarks[-keepers]]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

plot(h2o.performance(model_drf, xval = T))

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["bare.minimum",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

## Save the summary table
write.csv(sum.table, "./reports/bare.core.marks.loose.match.csv")


#### Check feature importance on closely matched set
## make a matrix to store the metric values
modes <- c("all.marks", "bare.minimum")
metrics <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
sum.table <- matrix(NA, ncol = 7, nrow = length(modes), dimnames = list(modes, metrics))

### Prepare the input data ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds", prop.dnase = 3)
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

plot(h2o.performance(model_drf, xval = T))

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


# bare minimum marks
features <- setdiff(colnames(train.set.h2o), target)
## Find all epigenetic marks
epimarks <- c(grep("broad", features), grep("narrow", features), grep("gapped", features))
## Find the marks to retain
keepers <- unique(c(grep("DNase", features[epimarks]),
					grep("H3K4me3", features[epimarks]),
					grep("H3K4me1", features[epimarks])))
features[epimarks][keepers]
features <- features[-epimarks[-keepers]]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)

plot(h2o.performance(model_drf, xval = T))

## Extract perfomance metrics
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
sum.table["bare.minimum",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

## Save the summary table
write.csv(sum.table, "./reports/bare.core.marks.close.match.csv")



plot(h2o.performance(model_drf, xval = T))
