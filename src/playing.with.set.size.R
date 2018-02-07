source("./lib/prepare.training.set.R")
tester.set <- prepare.training.set(path.to.positive.set = "./cache/funseq/E120.cadd.annotated.HSMM.variants.rds", path.to.negative.set = "./cache/funseq/E120.cadd.annotated.HSMM.negative.set.rds", prop.matched = 1, prop.dnase = 1, prop.random = 1, negative.positive.ratio = 6)
table(tester.set$regulatory)
test.cato <- as.h2o(tester.set)
