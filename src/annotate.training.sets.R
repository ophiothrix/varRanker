### Annotate a given set of variants with features from a given tissue ### 
# .libPaths("~/tools/R.packages/")
# rm(list = ls())
# gc()

### !!! Have you specified the path to fimo binary? !!! ###
## Load helper scripts
source("./lib/annotate.AS.DHS.R")
source("./lib/annotate.eQTLs.R")
source("./lib/annotate.negative.set.R")

## Foetal muscle with HSMM annotation
variants <- annotate.AS.DHS("fMuscle", "E120")
negative.variant.set <- annotate.negative.set(variant.file.id = "fMuscle", tissue.id = "E120")
rm(list = c("variants", "negative.variant.set"))

## HSMM with HSMM annotation
variants <- annotate.AS.DHS("HSMM", "E120")
negative.variant.set <- annotate.negative.set("HSMM", tissue.id = "E120")
rm(list = c("variants", "negative.variant.set"))

## Foetal Muscle with HMEC annotation
variants <- annotate.AS.DHS("fMuscle", "E119")
negative.variant.set <- annotate.negative.set("fMuscle", tissue.id = "E119")
rm(list = c("variants", "negative.variant.set"))

## Foetal skin with Dermal Fibroblast annotation
variants <- annotate.AS.DHS(variant.file.id = "fSkin_fibro", tissue.id = "E126")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 3000)
rm(list = c("variants", "negative.variant.set"))

## eQTL-based annotation
variants <- annotate.eQTLs("Breast_Mammary_Tissue", "E119")
negative.set <- annotate.negative.set("Breast_Mammary_Tissue", "E119")
rm(list = c("variants", "negative.variant.set"))


