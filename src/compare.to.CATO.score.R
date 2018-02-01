## Read in dbSNP annotated with CATO
rm(list = ls())
gc()
require(GenomicRanges)
require(data.table)
cato.scores <- readRDS("./cache/cato.scores.in.training.rds")
cato.ids <- paste0(as.character(seqnames(cato.scores)), "-", start(cato.scores), "-", cato.scores$ref, "-", cato.scores$alt)

train.pos <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds")

train.neg <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds")

pos.ids <- paste0(as.character(seqnames(train.pos)), "-", start(train.pos), "-", train.pos$REF, "-", train.pos$ALT)
train.pos$cato.score <- NA
train.pos$cato.score <- cato.scores$cato.score[match(pos.ids, cato.ids)]
## For all variants that don't overlap a DNAse site (DNAse bin is 0), set CATO score to 0
train.pos$cato.score[train.pos$DNase.macs2.narrowPeak.bin == 0] <- 0
summary(!is.na(train.pos$cato.score))
summary(train.pos %in% cato.scores)

neg.ids <- paste0(as.character(seqnames(train.neg)), "-", start(train.neg), "-", train.neg$REF, "-", train.neg$ALT)
train.neg$cato.score <- NA
train.neg$cato.score <- cato.scores$cato.score[match(neg.ids, cato.ids)]
## For all variants that don't overlap a DNAse site (DNAse bin is 0), set CATO score to 0
train.pos$cato.score[train.pos$DNase.macs2.narrowPeak.bin == 0] <- 0
summary(is.na(train.neg$cato.score))
summary(train.neg %in% cato.scores)
table(is.na(train.neg$cato.score), train.neg$source)


## Make a list of scores with CATO annotation. Use the loosely annotated scores
comb.pos <- GRanges()
train.pos <- readRDS("./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds")
comb.pos <- c(comb.pos, train.pos)
train.pos <- readRDS("./cache/cadd.annotated.variants/E120.annotated.fMuscle.variants.rds")
comb.pos <- c(comb.pos, train.pos)
train.pos <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds")
comb.pos <- c(comb.pos, train.pos)

comb.pos <- comb.pos[!duplicated(comb.pos)]
pos.ids <- paste0(as.character(seqnames(comb.pos)), "-", start(comb.pos), "-", comb.pos$REF, "-", comb.pos$ALT)
comb.pos$cato.score <- NA
comb.pos$cato.score <- cato.scores$cato.score[match(pos.ids, cato.ids)]
summary(!is.na(comb.pos$cato.score))
summary(comb.pos %in% cato.scores)
comb.pos <- comb.pos[!is.na(comb.pos$cato.score)]
saveRDS(comb.pos, "./cache/cato.annotated.positive.set.rds")

comb.neg <- GRanges()
train.neg <- readRDS("./cache/cadd.annotated.variants/E120.annotated.HSMM.negative.set.rds")
comb.neg <- c(comb.neg, train.neg)
train.neg <- readRDS("./cache/cadd.annotated.variants/E120.annotated.fMuscle.negative.set.rds")
comb.neg <- c(comb.neg, train.neg)
train.neg <- readRDS("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.negative.set.rds")
comb.neg <- c(comb.neg, train.neg)

comb.neg <- comb.neg[!duplicated(comb.neg)]
neg.ids <- paste0(as.character(seqnames(comb.neg)), "-", start(comb.neg), "-", comb.neg$REF, "-", comb.neg$ALT)
comb.neg$cato.score <- NA
comb.neg$cato.score <- cato.scores$cato.score[match(neg.ids, cato.ids)]
summary(is.na(comb.neg$cato.score))
summary(comb.neg %in% cato.scores)
comb.neg <- comb.neg[!is.na(comb.neg$cato.score)]
saveRDS(comb.neg, "./cache/cato.annotated.negative.set.rds")


## Add cato scores to an object
path.to.positive.training.set <- "./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds"
add.cato.score <- function(path.to.positive.training.set) {
	cato.scores <- readRDS("./cache/cato.scores.in.training.rds")
	cato.ids <- paste0(as.character(seqnames(cato.scores)), "-", start(cato.scores), "-", cato.scores$ref, "-", cato.scores$alt)
	
	## Load, annotate and save positive training set
	path.to.negative <- gsub(".variants.rds", ".negative.set.rds", path.to.positive.training.set)
	comb.pos <- readRDS(path.to.positive.training.set)
	comb.pos <- comb.pos[!duplicated(comb.pos)]
	pos.ids <- paste0(as.character(seqnames(comb.pos)), "-", start(comb.pos), "-", comb.pos$REF, "-", comb.pos$ALT)
	comb.pos$cato.score <- NA
	comb.pos$cato.score <- cato.scores$cato.score[match(pos.ids, cato.ids)]
	## For all variants that don't overlap a DNAse site (DNAse bin is 0), set CATO score to 0
	comb.pos$cato.score[comb.pos$DNase.macs2.narrowPeak.bin == 0] <- 0
	summary(!is.na(comb.pos$cato.score))
	summary(comb.pos %in% cato.scores)
	comb.pos <- comb.pos[!is.na(comb.pos$cato.score)]
	saveRDS(comb.pos, paste0("./cache/", gsub(".annotated", ".cato.annotated", basename(path.to.positive.training.set))))
	
	comb.neg <- readRDS(path.to.negative)
	comb.neg <- comb.neg[!duplicated(comb.neg)]
	neg.ids <- paste0(as.character(seqnames(comb.neg)), "-", start(comb.neg), "-", comb.neg$REF, "-", comb.neg$ALT)
	comb.neg$cato.score <- NA
	comb.neg$cato.score <- cato.scores$cato.score[match(neg.ids, cato.ids)]
	## For all variants that don't overlap a DNAse site (DNAse bin is 0), set CATO score to 0
	comb.neg$cato.score[comb.neg$DNase.macs2.narrowPeak.bin == 0] <- 0
	summary(!is.na(comb.neg$cato.score))
	summary(comb.neg %in% cato.scores)
	comb.neg <- comb.neg[!is.na(comb.neg$cato.score)]
	saveRDS(comb.neg, paste0("./cache/", gsub(".annotated", ".cato.annotated", basename(path.to.negative))))
}

add.cato.score("./cache/cadd.annotated.variants/E120.annotated.HSMM.variants.rds")
add.cato.score("./cache/cadd.annotated.variants/E120.annotated.fMuscle.variants.rds")
add.cato.score("./cache/cadd.annotated.variants/E126.annotated.fSkin_fibro.variants.rds")
add.cato.score("./cache/cadd.annotated.variants/E119.annotated.fMuscle.variants.rds")
