## Look at how the well-described regulatory mutations would fare
require(GenomicRanges)
terts <- GRanges("chr5", IRanges(c(1295228, 1295250), width = 1), REF = rep("G", 2), ALT = rep("A", 2)) 

source("./lib/annotate.GRanges.R")
terts.anno <- annotate.GRanges(terts, tissue.id = "E126")
terts.anno.df <- as.data.frame(mcols(terts.anno))
training.set <- terts.anno
test.set.h2o <- as.h2o(test.set.matched)
test.set.h2o <- as.h2o(terts.anno.df)
merge(x = test.set.matched, y = terts.anno.df)
colnames(terts.anno.df)[!(colnames(terts.anno.df) %in% colnames(test.set.matched))]
colnames(test.set.matched)[!(colnames(test.set.matched) %in% colnames(terts.anno.df))]
new.test <- rbind(
	test.set.matched[,intersect(colnames(test.set.matched), colnames(terts.anno.df))],
	terts.anno.df[,intersect(colnames(test.set.matched), colnames(terts.anno.df))])
test.set.h2o <- as.h2o(new.test)

apply(terts.anno.df, 2, print)
num.terts <- apply(terts.anno.df, 2, as.numeric)
num.terts[,which(is.na(num.terts[1,]))] <- apply(terts.anno.df[,which(is.na(num.terts[1,]))], 2, as.factor)

num.terts[,1]
saveRDS(terts.anno, "./cache/terts.rds")

terts.test <- prepare.training.set("./cache/terts.rds", "./cache/E126.annotated.fSkin_fibro.negative.set.rds")
ids <- intersect(colnames(terts.test), colnames(test.set.matched))
new.test <- rbind(terts.test[,ids], test.set.matched[test.set.matched$regulatory == 0,ids])
test.set.h2o <- as.h2o(new.test)

predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))
head(predicted$p1)
hist(predicted$p1)
summary(colnames(test.set.h2o) %in% colnames(train.set.h2o))
colnames(test.set.h2o)[!(colnames(test.set.h2o) %in% colnames(train.set.h2o))]
colnames(train.set.h2o)[!(colnames(train.set.h2o) %in% colnames(test.set.h2o))]

## Check the motifs
alt.seq <- paste0(substr(ref.seq[,2], 1, flnk), terts$ALT, substr(ref.seq[,2], flnk+2, flnk*2+1))
alt.seq <- cbind(coords, alt.seq)
ref.seq <- rbind(ref.seq, alt.seq)
getSeq(BSgenome.Hsapiens.UCSC.hg19, terts)
