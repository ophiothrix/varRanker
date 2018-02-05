## Load training set
# rm(list = ls())
train.set <- readRDS("./cache/all.variants.partial.annotation.rds")
dim(train.set)
table(train.set$regulatory)

## We want to avoid using the same variants as in the test sets, so remove all the variants with matching variants IDs from the training set
test.set <- readRDS("./cache/test.perfect2.rds")
train.set <- train.set[!(train.set$varID %in% test.set$varID),]
table(train.set$regulatory)

# test.set <- readRDS("./cache/test.imperfect4.rds")
# train.set <- train.set[!(train.set$varID %in% test.set$varID),]
# table(train.set$regulatory)
# 
# test.set <- readRDS("./cache/test.mismatched5.rds")
# train.set <- train.set[!(train.set$varID %in% test.set$varID),]
# table(train.set$regulatory)


## Get validation set
set.seed(17012018)
pos.ids <- sample(which(train.set$regulatory == 1), 1000, replace = F)
neg.ids <- sample(which(train.set$regulatory == 0), 1000, replace = F)

valid.set <- train.set[c(pos.ids, neg.ids),]
train.set <- train.set[-c(pos.ids, neg.ids),]

## Remove duplicated variants from the validation set
valid.set <- valid.set[!duplicated(valid.set$varID),]
table(valid.set$regulatory)
table(valid.set$source)

## Remove variants from the training set that match variants IDs in the validation set
train.set <- train.set[!(train.set$varID %in% valid.set$varID),]
table(train.set$source)


## We've deliberately over-sampled the negative examples, so we can have some freedom to remove duplicated. Let's use that freedom
duplicated(train.set$varID[train.set$regulatory == 0])
train.set <- train.set[-which(duplicated(train.set$varID) & train.set$regulatory == 0),]

## Finally, let's subsample the training set to go back to the original ratios of positive and negative variants
dnase.ids <- sample(x = which(train.set$source == "dnase"), size = round(sum(train.set$regulatory == 1) * 0.5, 0), replace = F)

matched.ids <- sample(x = which(train.set$source == "matched"), size = round(sum(train.set$regulatory == 1) * 0.45, 0), replace = F)

random.ids <- sample(x = which(train.set$source == "random"), size = round(sum(train.set$regulatory == 1) * 0.05, 0), replace = F)

pos.ids <- which(train.set$regulatory == 1)

## Combine the negative variants into a new training set
train.set <- train.set[c(pos.ids, dnase.ids, matched.ids, random.ids),]
table(train.set$source)
table(train.set$regulatory)
