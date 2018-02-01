## The rationale:

## We have some variant sets that are annotated with perfectly matched epigenetic features and others that are imperfectly matched. These are likely to behave differently in both training and predicting variant impact. These will be refered to as perfect and imperfect annotations.

## Additionally, we have some tissues where H3K9ac or H3K27ac or both, or none are available. These will be referred to as full (both marks available) and partial (none of the marks available) annotations. To keep things simple, we will not consider separately the cases where only one of the marks available.

## So we can stratify the datasets in several ways:

# We can split them by matching fidelity, i.e. perfectly or imperfectly matched
# Each of the above two classes (as well as all the variants) we can subdivide into variant sets with full and partial annotation.

rm(list = ls())
gc()

## Get the list of matched tissues
mappings <- read.csv("./reports/ROADMAP.to.AS_DHS.map.csv", header = T, stringsAsFactors = F, comment.char = "#")
head(mappings)
## Remove unmatched tissues
mappings <- mappings[mappings$AS.DHS.tissue != "",]
## We're going to use fetal skin fibroblasts as a test set. So remove them from the training part
mappings <- mappings[mappings$AS.DHS.tissue != "",]
mappings <- mappings[-grep("fSkin_fibro", mappings$AS.DHS.tissue),]
sum(mappings$AS.variants)
sum(mappings$AS.variants[mappings$perfect.match])
sum(mappings$AS.variants[!mappings$perfect.match])

# positive.sets <- positive.sets[-grep("cato", positive.sets)]
# positive.sets <- positive.sets[-grep("E119.annotated.fMuscle.variants.rds", positive.sets)]
# length(positive.sets)

variant.summary <- as.data.frame(matrix(NA, 6, 2))
colnames(variant.summary) <- c("negative", "positive")
rownames(variant.summary) <- c("all.variants.partial.annotation", 
							   "imperfect.variants.partial.annotation",
							   "perfect.variants.partial.annotation",
							   "all.variants.full.annotation", 
							   "imperfect.variants.full.annotation",
							   "perfect.variants.full.annotation")

# all.variants.partial.annotation ------------------------------------
## Get the list of all matched variant-tissue annotations
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)

p.sets <- p.sets[basename(p.sets) %in% paste0(mappings$tissue.id, ".annotated.", mappings$AS.DHS.tissue, ".variants.rds")]

combined.set <- combine.training.sets(p.sets)
dim(combined.set)
variant.summary["all.variants.partial.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/all.variants.partial.annotation.rds")


# imperfect.variants.partial.annotation ------------------------------
## Get the list of imperfectly matched variant-tissue annotations
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)
imperfect.variants <- mappings[!mappings$perfect.match,]

p.sets <- p.sets[basename(p.sets) %in% paste0(imperfect.variants$tissue.id, ".annotated.", imperfect.variants$AS.DHS.tissue, ".variants.rds")]

combined.set <- combine.training.sets(p.sets)
dim(combined.set)
variant.summary["imperfect.variants.partial.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/imperfect.variants.partial.annotation.rds")



# perfect.variants.partial.annotation --------------------------------
## Get the list of perfectly matched variant-tissue annotations
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)
perfect.variants <- mappings[mappings$perfect.match,]

p.sets <- p.sets[basename(p.sets) %in% paste0(perfect.variants$tissue.id, ".annotated.", perfect.variants$AS.DHS.tissue, ".variants.rds")]

combined.set <- combine.training.sets(p.sets)
dim(combined.set)
variant.summary["perfect.variants.partial.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/perfect.variants.partial.annotation.rds")



### all.variants.full.annotation
## Get the list of all variant-tissue annotations
## Limit the set to the tissues with full set of annotations available
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)

p.sets <- p.sets[basename(p.sets) %in% paste0(mappings$tissue.id, ".annotated.", mappings$AS.DHS.tissue, ".variants.rds")]

## Check the list of annotations
features <- list()
for (p.set.name in p.sets) {
	print(p.set.name)
	features[[p.set.name]] <- colnames(mcols(readRDS(p.set.name)))
}
full.set <- unique(unlist(features))
full.anno <- unlist(lapply(features, function(x) all(full.set %in% x)))

combined.set <- combine.training.sets(p.sets[full.anno])
dim(combined.set)
variant.summary["all.variants.full.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/all.variants.full.annotation.rds")



# imperfect.variants.full.annotation ---------------------------------
## Get the list of perfectly matched variant-tissue annotations
## Limit the set to the tissues with full set of annotations available
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)
imperfect.variants <- mappings[!mappings$perfect.match,]

p.sets <- p.sets[basename(p.sets) %in% paste0(imperfect.variants$tissue.id, ".annotated.", imperfect.variants$AS.DHS.tissue, ".variants.rds")]

## Check the list of annotations
features <- list()
for (p.set.name in p.sets) {
	print(p.set.name)
	features[[p.set.name]] <- colnames(mcols(readRDS(p.set.name)))
}
full.set <- unique(unlist(features))
full.anno <- unlist(lapply(features, function(x) all(full.set %in% x)))

combined.set <- combine.training.sets(p.sets[full.anno])
dim(combined.set)
variant.summary["imperfect.variants.full.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/imperfect.variants.full.annotation.rds")


# perfect.variants.full.annotation -----------------------------------
## Get the list of perfectly matched variant-tissue annotations
## Limit the set to the tissues with full set of annotations available
source("./lib/combine.training.sets.R")
p.sets <- list.files("./cache", "variants.rds", full.names = T)
perfect.variants <- mappings[mappings$perfect.match,]

p.sets <- p.sets[basename(p.sets) %in% paste0(perfect.variants$tissue.id, ".annotated.", perfect.variants$AS.DHS.tissue, ".variants.rds")]

## Check the list of annotations
features <- list()
for (p.set.name in p.sets) {
	print(p.set.name)
	features[[p.set.name]] <- colnames(mcols(readRDS(p.set.name)))
}
full.set <- unique(unlist(features))
full.anno <- unlist(lapply(features, function(x) all(full.set %in% x)))

combined.set <- combine.training.sets(p.sets[full.anno])
dim(combined.set)
variant.summary["perfect.variants.full.annotation",] <- table(combined.set$regulatory)
saveRDS(combined.set, "./cache/perfect.variants.full.annotation.rds")


write.csv(variant.summary, "./reports/annotated.variant.summary.csv")