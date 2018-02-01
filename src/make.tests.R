## Make test objects for model validation
source("./lib/prepare.training.set.R")
set.seed(19012018)

test.imperfect1 <- prepare.training.set(path.to.positive.set = "./cache/E055.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E055.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.imperfect1$regulatory == 0), -diff(table(test.imperfect1$regulatory)), replace = F)
test.imperfect1 <- test.imperfect1[-ids,]
saveRDS(test.imperfect1, "./cache/test.imperfect1.rds")

test.imperfect2 <- prepare.training.set(path.to.positive.set = "./cache/E056.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E056.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.imperfect2$regulatory == 0), -diff(table(test.imperfect2$regulatory)), replace = F)
test.imperfect2 <- test.imperfect2[-ids,]
saveRDS(test.imperfect2, "./cache/test.imperfect2.rds")

test.imperfect3 <- prepare.training.set(path.to.positive.set = "./cache/E128.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E128.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.imperfect3$regulatory == 0), -diff(table(test.imperfect3$regulatory)), replace = F)
test.imperfect3 <- test.imperfect3[-ids,]
saveRDS(test.imperfect3, "./cache/test.imperfect3.rds")


test.imperfect4 <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.imperfect4$regulatory == 0), -diff(table(test.imperfect4$regulatory)), replace = F)
test.imperfect4 <- test.imperfect4[-ids,]
saveRDS(test.imperfect4, "./cache/test.imperfect4.rds")

test.imperfect5 <- prepare.training.set(path.to.positive.set = "./cache/E017.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E017.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.imperfect5$regulatory == 0), -diff(table(test.imperfect5$regulatory)), replace = F)
test.imperfect5 <- test.imperfect5[-ids,]
saveRDS(test.imperfect5, "./cache/test.imperfect5.rds")

test.mismatched1 <- prepare.training.set(path.to.positive.set = "./cache/E119.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E119.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.mismatched1$regulatory == 0), -diff(table(test.mismatched1$regulatory)), replace = F)
test.mismatched1 <- test.mismatched1[-ids,]
saveRDS(test.mismatched1, "./cache/test.mismatched1.rds")

test.mismatched2 <- prepare.training.set(path.to.positive.set = "./cache/E120.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E120.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.mismatched2$regulatory == 0), -diff(table(test.mismatched2$regulatory)), replace = F)
test.mismatched2 <- test.mismatched2[-ids,]
saveRDS(test.mismatched2, "./cache/test.mismatched2.rds")

test.mismatched3 <- prepare.training.set(path.to.positive.set = "./cache/E085.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E085.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.mismatched3$regulatory == 0), -diff(table(test.mismatched3$regulatory)), replace = F)
test.mismatched3 <- test.mismatched3[-ids,]
saveRDS(test.mismatched3, "./cache/test.mismatched3.rds")

test.mismatched4 <- prepare.training.set(path.to.positive.set = "./cache/E021.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E021.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.mismatched4$regulatory == 0), -diff(table(test.mismatched4$regulatory)), replace = F)
test.mismatched4 <- test.mismatched4[-ids,]
saveRDS(test.mismatched4, "./cache/test.mismatched4.rds")


test.mismatched5 <- prepare.training.set(path.to.positive.set = "./cache/E034.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E034.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.mismatched5$regulatory == 0), -diff(table(test.mismatched5$regulatory)), replace = F)
test.mismatched5 <- test.mismatched5[-ids,]
saveRDS(test.mismatched5, "./cache/test.mismatched5.rds")


## Make test objects with cadd annotations
test.cadd1 <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E034.cadd.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E034.cadd.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.cadd1$regulatory == 0), -diff(table(test.cadd1$regulatory)), replace = F)
test.cadd1 <- test.cadd1[-ids,]
saveRDS(test.cadd1, "./cache/test.cadd1.rds")

## Make test objects with cadd annotations
test.cadd2 <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E056.cadd.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E056.cadd.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.cadd2$regulatory == 0), -diff(table(test.cadd2$regulatory)), replace = F)
test.cadd2 <- test.cadd2[-ids,]
saveRDS(test.cadd2, "./cache/test.cadd2.rds")

## Make test objects with cadd annotations
test.cadd3 <- prepare.training.set(path.to.positive.set = "./cache/cadd.annotated.variants/E126.cadd.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/cadd.annotated.variants/E126.cadd.annotated.fSkin_fibro.negative.set.rds", prop.matched = 0.5, prop.dnase = 0.5, prop.random = 0.05)
## Subsample negative set
ids <- sample(which(test.cadd3$regulatory == 0), -diff(table(test.cadd3$regulatory)), replace = F)
test.cadd3 <- test.cadd3[-ids,]
saveRDS(test.cadd3, "./cache/test.cadd3.rds")


