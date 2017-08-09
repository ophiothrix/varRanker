#### Check the effect of REF-ALT directionality on motif damage score
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)

dhs.variants <- readRDS("./cache/allele.biased.DHS.variants.all.rds")
head(dhs.variants)

imbalanced <- dhs.variants[dhs.variants$significance.level != "not_imbalanced"]

## Annotate the imbalanced variants with the damage scores
source("./lib/motif.damage.annotation.R")

## Add JASPAR damage score annotation to SNVs
imbalanced <- get.damage.scores.direct(database.path = "~/utils/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme", variants = imbalanced)
imbalanced$jaspar.damage.score <- imbalanced$log2.damage.score
imbalanced$jaspar.motif.hit <- imbalanced$motif.hit

## Add HOCOMOCO damage score annotation SNVs
imbalanced <- get.damage.scores.direct(database.path = "~/utils/motif_databases/HUMAN/HOCOMOCOv9.meme", variants = imbalanced)
imbalanced$hocomoco.damage.score <- imbalanced$log2.damage.score
imbalanced$hocomoco.motif.hit <- imbalanced$motif.hit

## Plot the distribution of damage score for the AS DHS variants as they are
par(mfrow = c(1,2))
hist(imbalanced$hocomoco.damage.score, freq = F, breaks = 100, main = "HOCOMOCO damage score", ylim = c(0, 0.9))
hist(imbalanced$jaspar.damage.score, freq = F, breaks = 100, main = "JASPAR damage score")
## There are more low score variants in JASPAR database, likely to its lower specificity to human motifs.
## HOCOMOCO damage score has approximately equal number of positive and negative damage score. With a slight shift towards positive (loss of function). JASPAR damage score has a slightly stronger shift towards positive damage scores and somewhat smoother distribution of low effect scores.

## Get reference allele
imbalanced.reset <- imbalanced
refs <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, imbalanced))

head(refs)
summary(imbalanced$REF == refs) ## For about 50% of all variants allele designated as reference does not match human reference genome.
## Check that for the mismatched variants, the alternate allele matches genome reference
summary(imbalanced$ALT[imbalanced$REF != refs] == refs[imbalanced$REF != refs]) ## All true

alts <- NA
## For variants where reference allele matches genome reference, assign the alternative allele as is
alts[imbalanced$REF == refs] <- imbalanced$ALT[imbalanced$REF == refs]

## For mismatched variants, take the designated reference allele as alternate.
alts[imbalanced$REF != refs] <- imbalanced$REF[imbalanced$REF != refs]

summary(is.na(alts))

## Assign re-set reference and alternate alleles to the GRanges object
imbalanced.reset$REF <- refs
imbalanced.reset$ALT <- alts

## Calculate motif damage scores for the re-set object

## Add JASPAR damage score annotation to SNVs
imbalanced.reset <- get.damage.scores.direct(database.path = "~/utils/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme", variants = imbalanced.reset)
imbalanced.reset$jaspar.damage.score <- imbalanced.reset$log2.damage.score
imbalanced.reset$jaspar.motif.hit <- imbalanced.reset$motif.hit

## Add HOCOMOCO damage score annotation SNVs
imbalanced.reset <- get.damage.scores.direct(database.path = "~/utils/motif_databases/HUMAN/HOCOMOCOv9.meme", variants = imbalanced.reset)
imbalanced.reset$hocomoco.damage.score <- imbalanced.reset$log2.damage.score
imbalanced.reset$hocomoco.motif.hit <- imbalanced.reset$motif.hit


## Plot the distribution of damage score for the AS DHS variants as they are and the reset variants on the same plot
par(mfrow = c(1,2))
hist(imbalanced$hocomoco.damage.score, freq = F, breaks = 100, main = "HOCOMOCO damage score", ylim = c(0, 0.9), col = "#00FF0088")
hist(imbalanced.reset$hocomoco.damage.score, freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Variants as is", "Reset variants"), fill = c("#00FF0088", "#FF000088"))

hist(imbalanced$jaspar.damage.score, freq = F, breaks = 100, main = "JASPAR damage score", ylim = c(0, 0.9), col = "#00FF0088")
hist(imbalanced.reset$jaspar.damage.score, freq = F, breaks = 100, col = "#FF000088", add = T)
legend("topright", legend = c("Variants as is", "Reset variants"), fill = c("#00FF0088", "#FF000088"))

## As expected, resetting the alleles so that the reference is matching human reference genome almost completely shifts high gain of morif scores to the loss of motif scores (from negative to positive values). A similar pattern is observed at the lower values that cluster around 0.

## As discussed earlier, let's stick to the current practice of NOT re-setting the alleles. However, not doing so, results in both positive 
