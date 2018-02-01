require(pROC)
require(RColorBrewer)
rm(list = ls())
gc()
require(h2o)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
clrs <- brewer.pal(5, "Dark2")
nfolds <- 5

set.seed(15052017)
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.7, prop.random = 0.05)

## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}

## Remove epigenetic data for all but 5 core marks:
core.marks <- c("H3K4me1", "H3K4me3", "H3K36me3", "DNase.macs2", "H3K9me3", "H3K27me3", "H3K27ac", "H3K9ac")
gsub("^(.*)\\..*\\..*", "\\1", colnames(training.set)[grep("FC", colnames(training.set))])

## get column ids of the marks we want to keep
marks.to.keep <- numeric()
for (mark in core.marks) {
	print(mark)
	marks.to.keep <- c(marks.to.keep, grep(mark, colnames(training.set)))
}
colnames(training.set)[marks.to.keep]
colnames(training.set)[-marks.to.keep]
## get id columns of all roadmap marks
all.marks <- numeric()
for (metric in c("bin", "FC", "pval", "qval", "dist.to.summit")) {
	print(metric)
	all.marks <- c(all.marks, grep(metric, colnames(training.set)))
}
length(all.marks)
length(marks.to.keep)
training.set.core <- training.set[, -all.marks[!(all.marks %in% marks.to.keep)]]
dim(training.set.core)
head(training.set.core)

## Train a model with all core marks
train.set.h2o <- as.h2o(training.set.core)
# train.set.h2o <- as.h2o(no.dnase.training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.varimp_plot(model_drf, num_of_features = 50)

## For the two intermittently available marks (H3K27ac and H3K9ac) leave only binary metric, which can be supposedly obtained from the imputed data.
training.set.core.imp <- training.set.core
training.set.core.imp <- training.set.core.imp[,-setdiff(grep("H3K27ac", colnames(training.set.core.imp)), setdiff(grep("H3K27ac.*bin", colnames(training.set.core.imp)), grep("H3K27ac.*broad.*bin", colnames(training.set.core.imp))))]
training.set.core.imp <- training.set.core.imp[,-setdiff(grep("H3K9ac", colnames(training.set.core.imp)), setdiff(grep("H3K9ac.*bin", colnames(training.set.core.imp)), grep("H3K9ac.*broad.*bin", colnames(training.set.core.imp))))]

colnames(training.set.core.imp)

## Train a model with all core marks
train.set.h2o <- as.h2o(training.set.core.imp)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Stratified", keep_cross_validation_predictions = TRUE, seed = 16052017)
summary(model_drf)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.varimp_plot(model_drf, num_of_features = 50)

### Prepare test set ###
test.set.semi.matched <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E120.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E120.annotated.fMuscle.negative.set.rds", prop.matched = 0.45, prop.dnase = 0.5, prop.random = 0.05)
test.set.h2o <- as.h2o(test.set.semi.matched)
table(test.set.semi.matched$regulatory)
h2o.performance(model_drf, newdata = test.set.h2o)
plot(h2o.performance(model_drf, newdata = test.set.h2o))

# hist(training.set$jaspar.abs.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088")
# hist(training.set$jaspar.abs.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
# 
# hist(training.set$jaspar.loss.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088")
# hist(training.set$jaspar.loss.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
# 
# hist(training.set$jaspar.gain.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088")
# hist(training.set$jaspar.gain.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
# 
# 
# hist(training.set$hocomoco.abs.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088", ylim = c(0,2))
# hist(training.set$hocomoco.abs.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
# 
# hist(training.set$hocomoco.gain.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088", ylim = c(0,2.2))
# hist(training.set$hocomoco.gain.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
# 
# hist(training.set$hocomoco.loss.score[training.set$regulatory == 0], breaks = 100, freq = F, col = "#00FF0088", ylim = c(0,2.2))
# hist(training.set$hocomoco.loss.score[training.set$regulatory == 1], breaks = 100, add = T, col = "#FF000088", freq = F)
