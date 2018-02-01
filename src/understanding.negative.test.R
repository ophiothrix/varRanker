## Understanding negative test on the model performance

## Negative test elements impacting the model performance:
# - Window size for collecting matched variants
# - Allowed allele frequency
# - Ratio between negative classes

## Test scenarios:
## MAF > 5% vs all MAFs
## 1 kb window vs 2 kb window

rm(list = ls())
gc()
source("./lib/annotate.negative.set.R")

## 1kb window MAFs > 0.05
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 1000, maf.cutoff = 0.05)
system("mv ./cache/E126.annotated.fSkin_fibro.negative.set.rds ./cache/E126.annotated.fSkin_fibro.negative.set.1kb.maf5.rds")
## 2 kb window MAF > 0.05
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 2000, maf.cutoff = 0.05)
system("mv ./cache/E126.annotated.fSkin_fibro.negative.set.rds ./cache/E126.annotated.fSkin_fibro.negative.set.2kb.maf5.rds")
## 1 kb window all MAF
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 1000, maf.cutoff = 0)
system("mv ./cache/E126.annotated.fSkin_fibro.negative.set.rds ./cache/E126.annotated.fSkin_fibro.negative.set.1kb.all.maf.rds")
## 2 kb window all MAF
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 2000, maf.cutoff = 0)
system("mv ./cache/E126.annotated.fSkin_fibro.negative.set.rds ./cache/E126.annotated.fSkin_fibro.negative.set.2kb.all.maf.rds")

## Another factor is the source of DNAse variants. Currently I'm taking pan-tissue DNAse variants. This might not be so optimal, as the DNase may not actually be matched. So even though the idea is to take DNase non biased variants, the algorithm may see them as non-DNase variants, if the variant doesn't actually overlap DNase site in the given tissue.
