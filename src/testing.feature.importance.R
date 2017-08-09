#### Objective: look at the importance of different features in the prediction accuracy. ####

## E.g. for H3K9ac and H3K27ac, how much accuracy would we lose if we used only FC, or only pvalue or only binary compared to all of them, compared to each other. The idea being that we can substitute real experimental data with imputed.
## How much accuracy would we lose by excluding some of the marks? E.g. H3k4me2 or H2Z

rm(list = ls())
gc()
## Check the accuracy on a separate test set

require(pROC)
require(RColorBrewer)
## initiate H2O
require(h2o)
localH2O = h2o.init(nthreads=4)

clrs <- brewer.pal(8, "Dark2")
set.seed(15052017)
# training.set <- readRDS("./cache/prepared.training.set.rds")
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/E120.annotated.HSMM.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.50, prop.random = 0.05)
training.set$regulatory <- as.factor(training.set$regulatory)
table(training.set$regulatory)
## Remove qval features, as they are redundant with pval
training.set <- training.set[-grep("qval", colnames(training.set))]
dim(training.set)


## Convert dataset
train.set.h2o <- as.h2o(training.set)
target <- "regulatory"
all.features <- setdiff(colnames(train.set.h2o), target)
nfolds <- 5

metrics <- matrix(NA, nrow = 5, ncol = 7)
colnames(metrics) <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
rownames(metrics) <- c("all", "FC", "pval", "qval", "bin")
## Build the model using ALL features
model_drf <- h2o.randomForest(x = all.features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics["all",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

for (current.feature in c("FC", "pval", "qval", "bin")) {
	spare.features <- setdiff(c("FC", "pval", "qval", "bin"), current.feature)
	print(current.feature)
	features <- all.features[-c(grep(spare.features[1], all.features), grep(spare.features[2], all.features), grep(spare.features[3], all.features))]
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	model.perf <- h2o.performance(model_drf, xval = T)
	f1 <- h2o.F1(model.perf)
	metrics[current.feature,] <-c(
	AUC = h2o.auc(model.perf),
	Gini = h2o.giniCoef(model.perf),
	max.F1 = max(h2o.F1(model.perf)$f1),
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
	)
}
## also try a combination of the best two - FC and q value
features <- all.features[-c(grep("pval", all.features), grep("bin", all.features))]
model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics <- rbind(metrics, "FC+qval"=c(
	AUC = h2o.auc(model.perf),
	Gini = h2o.giniCoef(model.perf),
	max.F1 = max(h2o.F1(model.perf)$f1),
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
))
metrics <- metrics[order(metrics[,1], decreasing = T),]

pdf("./graphs/effect.of.epigenetic.metric.on.model.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()

## Compare building models with broad peak, narrow peak and gapped peak
metrics <- matrix(NA, nrow = 4, ncol = 7)
colnames(metrics) <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
rownames(metrics) <- c("all", "broad", "narrow", "gapped")
## Build the model using ALL features
model_drf <- h2o.randomForest(x = all.features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics["all",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

for (current.feature in c("broad", "narrow", "gapped")) {
	spare.features <- setdiff(c("broad", "narrow", "gapped"), current.feature)
	print(current.feature)
	features <- all.features[-c(grep(spare.features[1], all.features), grep(spare.features[2], all.features))]
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	model.perf <- h2o.performance(model_drf, xval = T)
	metrics[current.feature,] <-c(
		AUC = h2o.auc(model.perf),
		Gini = h2o.giniCoef(model.perf),
		max.F1 = max(h2o.F1(model.perf)$f1),
		sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
		specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
		recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
		accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
	)
	
}

pdf("./graphs/effect.of.mark.definition.on.model.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()

## When excluding all narrow marks, model performs substantially worse. Most likely because we exclude DNase, which only has narrow marks. include it back, after removing all other narrow marks. 
for (current.feature in c("broad", "narrow", "gapped")) {
	spare.features <- setdiff(c("broad", "narrow", "gapped"), current.feature)
	print(current.feature)
	features <- all.features[-c(grep(spare.features[1], all.features), grep(spare.features[2], all.features))]
	## Add DNase features back
	features <- unique(c(features, all.features[grep("DNase", all.features)]))
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	model.perf <- h2o.performance(model_drf, xval = T)
	metrics[current.feature,] <-c(
		AUC = h2o.auc(model.perf),
		Gini = h2o.giniCoef(model.perf),
		max.F1 = max(h2o.F1(model.perf)$f1),
		sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
		specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
		recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
		accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
	)
	
}

pdf("./graphs/effect.of.mark.definition.on.model.keep.dnase.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()

## Excluding narrow marks also excludes summits. Retain the summits, while discarding all other narrowPeak features
for (current.feature in c("broad", "narrow", "gapped")) {
	spare.features <- setdiff(c("broad", "narrow", "gapped"), current.feature)
	print(current.feature)
	features <- all.features[-c(grep(spare.features[1], all.features), grep(spare.features[2], all.features))]
	features <- unique(c(features, all.features[grep("summit", all.features)]))
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	model.perf <- h2o.performance(model_drf, xval = T)
	metrics[current.feature,] <-c(
		AUC = h2o.auc(model.perf),
		Gini = h2o.giniCoef(model.perf),
		max.F1 = max(h2o.F1(model.perf)$f1),
		sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
		specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
		recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
		accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
	)
	
}

pdf("./graphs/effect.of.mark.definition.on.model.keep.summits.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()

## Are summits useful?
metrics <- matrix(NA, nrow = 2, ncol = 7)
colnames(metrics) <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
rownames(metrics) <- c("w.summits", "wo.summits")
## Build the model using ALL features
model_drf <- h2o.randomForest(x = all.features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics["w.summits",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

## Remove summits
features <- all.features[-grep("summit", all.features)]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics["wo.summits",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)


## It looks like distance to summit has a very limited, although a positive effect on the model. 
pdf("./graphs/effect.of.summit.inclusion.on.model.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()

### Now concentrate on specific marks. Use all available metrics and peak types, just to simplify. Keep DNase in all cases
marks <- unique(gsub(".narrow.*", "", all.features[grep("narrow", all.features)]))
marks <- marks[-grep("DNase", marks)]
metrics <- matrix(NA, nrow = length(marks)+1, ncol = 7)
colnames(metrics) <- c("AUC", "Gini", "max.F1", "sensitivity", "specificity", "recall", "accuracy")
rownames(metrics) <- c("all", marks)
## Build the model using ALL features
model_drf <- h2o.randomForest(x = all.features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
model.perf <- h2o.performance(model_drf, xval = T)
f1 <- h2o.F1(model.perf)
metrics["all",] <- c(
	AUC=h2o.auc(model.perf),
	Gini=h2o.giniCoef(model.perf),
	max.F1=f1[which.max(f1$f1),2],
	sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),])[[1]],
	specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
	recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
	accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
)

for (current.feature in marks) {
	spare.features <- setdiff(marks, current.feature)
	print(current.feature)
	spare.ids <- numeric()
	for (i in 1:length(spare.features)) {
		spare.ids <- c(spare.ids, grep(spare.features[i], all.features))
	}
	features <- all.features[-spare.ids]
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
	model.perf <- h2o.performance(model_drf, xval = T)
	metrics[current.feature,] <-c(
		AUC = h2o.auc(model.perf),
		Gini = h2o.giniCoef(model.perf),
		max.F1 = max(h2o.F1(model.perf)$f1),
		sensitivity = h2o.sensitivity(model.perf, f1[which.max(f1$f1),1])[[1]],
		specificity = h2o.specificity(model.perf, f1[which.max(f1$f1),1])[[1]],
		recall = h2o.recall(model.perf, f1[which.max(f1$f1),1])[[1]],
		accuracy = h2o.accuracy(model.perf, f1[which.max(f1$f1),1])[[1]]
	)
	
}

pdf("./graphs/effect.of.mark.on.model.pdf", 8, 6)
barplot(metrics, beside = T, col = clrs[1:nrow(metrics)])
legend("topright", rownames(metrics), fill = clrs[1:nrow(metrics)], bty = "n")
dev.off()
