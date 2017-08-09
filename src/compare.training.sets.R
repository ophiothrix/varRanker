## compare training sets

new.training.set <- prepare.training.set(path.to.positive.set = "./cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "./cache/E126.annotated.fSkin_fibro.negative.set.rds")
old.training.set <- prepare.training.set(path.to.positive.set = "../toy.RF.run/cache/E126.annotated.fSkin_fibro.variants.rds", path.to.negative.set = "../toy.RF.run/cache/E126.annotated.fSkin_fibro.negative.set.rds")

write.csv(cbind(colnames(old.training.set), colnames(new.training.set)), "./reports/compare.features.csv")


new.positives <- readRDS("./cache/E126.annotated.fSkin_fibro.variants.rds")
old.positives <- readRDS("../toy.RF.run/cache/E126.annotated.fSkin_fibro.variants.rds")
new.positives <- new.positives[match(old.positives, new.positives)]

new.positives <- as.data.frame(mcols(new.positives))[,-(1:3)]
old.positives <- as.data.frame(mcols(old.positives))[,-(1:3)]

new.positives$feature <- as.factor(new.positives$feature)
old.positives$feature <- as.factor(old.positives$feature)

new.positives$chromHMM.state <- as.factor(new.positives$chromHMM.state)
old.positives$chromHMM.state <- as.factor(old.positives$chromHMM.state)

ids <- intersect(colnames(new.positives), colnames(old.positives))

pdf("./graphs/compare.features.pdf")
for (id in ids) {
	print(id)
	print(all(new.positives[,id] == old.positives[,id]))
	plot(new.positives[,id], old.positives[,id], pch = 16, cex = 0.7, col = "#44444444", main = id)
}				 
dev.off()
dim(new.positives)
dim(old.positives)


new.positives <- readRDS("./cache/E126.annotated.fSkin_fibro.variants.rds")
old.positives <- readRDS("../toy.RF.run/cache/E126.annotated.fSkin_fibro.variants.rds")

mcols(new.positives) <- mcols(new.positives)[colnames(mcols(new.positives)) == "GERP.score"]
mcols(old.positives) <- mcols(old.positives)[colnames(mcols(old.positives)) == "GERP.score"]
new.positives <- new.positives[match(old.positives, new.positives)]

new.positives$old.GERP.score <- old.positives$GERP.score

plot(new.positives$GERP.score, old.positives$GERP.score, pch = 16, cex = 0.7, col = "#44444444", main = id)

plot(new.positives$hocomoco.damage.score, old.positives$hocomoco.damage.score, pch = 16, cex = 0.7, col = "#44444444", main = id)


neg.set <- readRDS("./cache/E126.annotated.fSkin_fibro.negative.set.2kb.rds")
table(neg.set$source)
	  

new.positives <- readRDS("./cache/E126.annotated.fSkin_fibro.variants.rds")
old.positives <- readRDS("../toy.RF.run/cache/E126.annotated.fSkin_fibro.variants.rds")
new.positives <- new.positives[match(old.positives, new.positives)]

new.positives <- as.data.frame(mcols(new.positives))
old.positives <- as.data.frame(mcols(old.positives))

ids <- intersect(colnames(new.positives), colnames(old.positives))
sames <- unlist(lapply(ids, function(x) (all(new.positives[,x] == old.positives[,x]))))
names(sames) <- ids
sames
all(sames)
which(sames == F)


new.positives$hocomoco.damage.score.old <- old.positives$hocomoco.damage.score

mcols(new.positives) <- mcols(new.positives)[grep("damage.score", colnames(mcols(new.positives)))]
mcols(old.positives) <- mcols(old.positives)[grep("damage.score", colnames(mcols(old.positives)))]

new.positives$hocomoco.damage.score != old.positives$hocomoco.damage.score
hocomoco.damage.score <- cbind(new.positives$hocomoco.damage.score, old.positives$hocomoco.damage.score)
hocomoco.damage.score[new.positives$hocomoco.damage.score != old.positives$hocomoco.damage.score,]
plot(hocomoco.damage.score[,1], hocomoco.damage.score[,2])
