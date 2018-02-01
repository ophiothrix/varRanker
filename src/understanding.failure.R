## I'm mostly worried about the specificity. I.e. I'd rather have false negatives than false positives. Interestingly, they fluctuate quite a bit. False negative error rate - between 20 and 30%. False negatie is quite stable around 30. Needs confirming.

## Since I'm mostly worried about false positives, let's take a closer look at what we're predicting as false positive

h2o.performance(model_drf, newdata = test.set.h2o)
predicts <- as.data.frame(h2o.predict(object = model_drf, newdata = test.set.h2o))
head(predicts)
FPs <- test.set.semi.matched[which(predicts$predict == 1 & test.set.semi.matched$regulatory == 0),]
dim(FPs)
head(FPs)

## Load the original test set
src.test <- as.data.frame(mcols(readRDS("./cache/phastCons/E120.annotated.fMuscle.negative.set.rds")))
head(src.test)
src.test <- src.test[match(FPs$id, src.test$id),]
table(src.test$source)
table(src.test$DNase.macs2.narrowPeak.bin[src.test$source == "dnase"])
table(src.test$DNase.macs2.narrowPeak.bin[src.test$source == "matched"])

table(src.test$chromHMM.state[src.test$source == "dnase"])
table(src.test$chromHMM.state[src.test$source == "matched"])

table(FPs$DNase.macs2.narrowPeak.bin)

table(test.set.semi.matched$DNase.macs2.narrowPeak.bin[test.set.semi.matched$regulatory == 1])

colnames(FPs[,-ncol(FPs)]) == colnames(src.test[,-c(1,2,3,5)])


## 
