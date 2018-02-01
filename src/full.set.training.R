## To start with, stick to the partially annotated variants - there's more of them

## We will need to generate 2 test sets: one closely matched and another loosely matched.
set.seed(17072018)

train.perfect <- readRDS("./cache/perfect.variants.partial.annotation.rds")
pos.ids <- sample(which(train.perfect$regulatory == 1), 500, replace = F)
neg.ids <- sample(which(train.perfect$regulatory == 0), 500, replace = F)
test.perfect <- train.perfect[c(pos.ids, neg.ids),]
train.perfect <- train.perfect[-c(pos.ids, neg.ids),]
saveRDS(train.perfect, "./cache/train.perfect.rds")
saveRDS(test.perfect, "./cache/test.perfect.rds")

train.imperfect <- readRDS("./cache/imperfect.variants.partial.annotation.rds")
pos.ids <- sample(which(train.imperfect$regulatory == 1), 500, replace = F)
neg.ids <- sample(which(train.imperfect$regulatory == 0), 500, replace = F)
test.imperfect <- train.imperfect[c(pos.ids, neg.ids),]
train.imperfect <- train.imperfect[-c(pos.ids, neg.ids),]
saveRDS(train.imperfect, "./cache/train.imperfect.rds")
saveRDS(test.imperfect, "./cache/test.imperfect.rds")

train.full <- rbind(train.perfect, train.imperfect)
saveRDS(train.full, "./cache/train.full.rds")

dim(train.full)
dim(test.imperfect)


