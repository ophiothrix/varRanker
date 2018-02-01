dir.create("./cache/phastCons")
source("./lib/phastCons.R")
variants <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds")
variants <- annotate.phastCons(variants)
saveRDS(variants, "./cache/phastCons/E126.annotated.fSkin_fibro.variants.rds")

variants <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds")
variants <- annotate.phastCons(variants)
saveRDS(variants, "./cache/phastCons/E126.annotated.fSkin_fibro.negative.set.rds")

variants <- readRDS("./cache/cadd.annotated.variants/E120.annotated.fMuscle.variants.rds")
variants <- annotate.phastCons(variants)
saveRDS(variants, "./cache/phastCons/E120.annotated.fMuscle.variants.rds")

variants <- readRDS("./cache/cadd.annotated.variants/E120.annotated.fMuscle.negative.set.rds")
variants <- annotate.phastCons(variants)
saveRDS(variants, "./cache/phastCons/E120.annotated.fMuscle.negative.set.rds")



## Setting up
rm(list = ls())
gc()
require(pROC)
require(h2o)
require(RColorBrewer)
## Initiate H2O instance
localH2O <- h2o.init(nthreads=-1)
nfolds <- 5
clrs <- brewer.pal(5, "Dark2")

### Prepare the training set ###
source("./lib/prepare.training.set.R")
training.set <- prepare.training.set(path.to.positive.set = "./cache/phastCons/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/phastCons/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.3, prop.dnase = 0.6, prop.random = 0.05)
table(training.set$regulatory)

### Prepare test set ###
test.set.semi.matched <- prepare.training.set(path.to.positive.set = "./cache/phastCons/E120.annotated.fMuscle.variants.rds", path.to.negative.set = "./cache/phastCons/E120.annotated.fMuscle.negative.set.rds", prop.matched = 0.6, prop.dnase = 0.6, prop.random = 0.05)
test.set.h2o <- as.h2o(test.set.semi.matched)

## Remove qval features, as they are redundant with pval
if( length(grep("qval", colnames(training.set))) != 0){
	training.set <- training.set[-grep("qval", colnames(training.set))]
}

## Train a random forest model
train.set.h2o <- as.h2o(training.set)
# train.set.h2o <- as.h2o(no.dnase.training.set)
target <- "regulatory"
features <- setdiff(colnames(train.set.h2o), target)
## Remove cadd score from training features
features <- features[!(features %in% c("cadd.score"))]
features <- features[!(features %in% c("cadd.score", "phastCons100", "phastCons46.mammals", "phastCons46.primates"))]
features <- features[!(features %in% c("cadd.score", "phastCons100", "phastCons46.mammals"))]
features <- features[!(features %in% c("cadd.score", "phastCons100"))]

model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = 100, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE, seed = 16052017)
h2o.performance(model_drf, xval = T)
h2o.performance(model_drf, train = T)
h2o.performance(model_drf, newdata = test.set.h2o)
h2o.varimp_plot(model_drf, num_of_features = 30)
as.data.frame(h2o.varimp(model_drf))

plot(h2o.performance(model_drf, xval = T))


# 0.7805424 - with cadd
# 0.7762511 - no cadd
# 0.7811876 - no phastCons100
# 0.7820895 - no phastConsMammals
# 0.7746062 - no phastConsPrimates
# 0.7417185 no phastCons
# 0.7451066 phastCons Primates
# 0.7426197 phastCons Mammals
