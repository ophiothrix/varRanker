old.negative <- readRDS("../toy.RF.run/cache/E126.annotated.fSkin_fibro.negative.set.rds")
new.negative <- readRDS("./cache/E126.annotated.fSkin_fibro.negative.set.2kb.rds")

old.positive <- readRDS("../toy.RF.run/cache/E126.annotated.fSkin_fibro.variants.rds")
new.positive <- readRDS("./cache/E126.annotated.fSkin_fibro.variants.rds")

par(mfrow = c(1,2))
hist(old.negative$hocomoco.damage.score, breaks = 100, freq = F, col = "#00FF0088")
hist(new.negative$hocomoco.damage.score, breaks = 100, col = "#FF000088", add = T, freq = F)
legend("topright", legend = c("old set", "new set"), bty = "n", fill = c("#00FF0088", "#FF000088"))
hist(old.positive$hocomoco.damage.score, breaks = 100, freq = F, col = "#00FF0088")
hist(new.positive$hocomoco.damage.score, breaks = 100, col = "#FF000088", add = T, freq = F)
legend("topright", legend = c("old set", "new set"), bty = "n", fill = c("#00FF0088", "#FF000088"))


par(mfrow = c(1,3))
for (src in unique(old.negative$source)) {
	# src <- "dnase"
	ids.new <- new.negative$source == src
	ids.old <- old.negative$source == src
	
	hist(old.negative$hocomoco.damage.score[ids.old], breaks = 100, col = "#00FF0088", main = src, freq = F, ylim = c(0, 0.75))
	hist(new.negative$hocomoco.damage.score[ids.new], breaks = 100, col = "#FF000088", add = T, freq = F)
legend("topright", legend = c("old set", "new set"), bty = "n", fill = c("#00FF0088", "#FF000088"))
}
## The difference, as expected, is mostly in DNAse variants

par(mfrow = c(1,2))
hist(new.negative$hocomoco.damage.score, breaks = 100, col = "#00FF0088", main = "Positive vs new negative", freq = F, ylim = c(0, 0.75))
hist(new.positive$hocomoco.damage.score, breaks = 100, col = "#FF000088", add = T, freq = F)

hist(old.negative$hocomoco.damage.score, breaks = 100, col = "#00FF0088", main = "Positive vs old negative", freq = F, ylim = c(0, 0.75))
hist(new.positive$hocomoco.damage.score, breaks = 100, col = "#FF000088", add = T, freq = F)

hist(abs(old.negative$hocomoco.damage.score)[ids.old], breaks = 100, col = "#00FF0088", main = src, freq = F)
hist(abs(new.negative$hocomoco.damage.score)[ids.new], breaks = 100, col = "#FF000088", add = T, freq = F)



hist(abs(new.negative$hocomoco.damage.score), breaks = 100, col = "#00FF0088", main = "Positive vs new negative", freq = F, ylim = c(0, 1.7))
hist(abs(new.positive$hocomoco.damage.score), breaks = 100, col = "#FF000088", add = T, freq = F)
legend("topright", legend = c("New negative set", "Positive set"), bty = "n", fill = c("#00FF0088", "#FF000088"))

hist(abs(old.negative$hocomoco.damage.score), breaks = 100, col = "#00FF0088", main = "Positive vs old negative", freq = F, ylim = c(0, 1.7))
hist(abs(new.positive$hocomoco.damage.score), breaks = 100, col = "#FF000088", add = T, freq = F)
legend("topright", legend = c("Old negative set", "Positive set"), bty = "n", fill = c("#00FF0088", "#FF000088"))



## Check if the reference allele in the old and new negative sets are actually reference alleles
require(BSgenome.Hsapiens.UCSC.hg19)
new.ref <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, new.negative))
summary(new.negative$REF == new.ref)

old.ref <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, old.negative))
summary(old.negative$REF == old.ref)



#### Conclusion ####

# In the original analysis we computed the damage scores for positive variants from the alleles as they are in the AS DHS table, i.e. a lot of the time the reference (or rather allele1) allele is not the same as human genome reference allele. As a rule, the first allele in the table is the one that has more reads in the DHS. Therefore it appears that for these alleles there is a shift towards motif gain (negative damage score). Due to an oversight, a bias has been introduced in the original workflow, when generating the negative training set. In that case the reference and alternate alleles have been re-set so that the reference allele matches human reference genome. As a result of this re-setting almost all variants that have an strong effect on motif were assigned a loss of motif (positive damage score).

## In the second iteration, the negative scores from AS DHS sites were left unchanged. As a result, the distribution of damage scores became more similar to the positive training set with the consequent small drop in the model performance and down-grading of motif damage scores as predictive factors.

## Since we don't care a lot about the direction of the motif change (loss or gain), but only the potential of the variant to impact motif recognition, we'll keep the current practive of not re-setting reference and alternate alleles. But we will use absolute motif damage score instead of the signed one. The current practice still has a disadvantage of potentially missing de-novo motifs created by the variant, when the original sequence is so different from the motif, that it does not get picked up by the motif scanner. This caveat will potentially be addressed in the future.

## Notably, keeping the alleles as is (without resetting the reference) results in a slight drop in model performance and a substantial drop in the predictive power of motif damage score. This is most likely driven by the substantially larger number of gain of motif variants in the training set. Note that about half of all the variants in the negative set are common variants, where the reference would be matched to human reference genome. Therefore, we may inadvertently introduce an artificial bias by NOT resetting AS DHS variants. Effectively we would be saying "positive variants tend to have negative damage score". So when we take an absolute (which is equivalent to matching positive variants to the reference genome) we lose this bias. Note that the variant sets that we are likely to encounter as discovery sets WILL have their reference allele derived from the human genome reference. Therefore I propose to reset all AS DHS variant alleles to match human genome reference.
