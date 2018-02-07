### Annotate a given set of variants with features from a given tissue ### 
# .libPaths("~/tools/R.packages/")
rm(list = ls())
gc()

### !!! Have you specified the path to fimo binary? !!! ###
## Load helper scripts
source("./lib/annotate.AS.DHS.R")
source("./lib/annotate.eQTLs.R")
source("./lib/annotate.negative.set.R")

### Annotate HSMM variants with:
## foetal leg muscle
variants <- annotate.AS.DHS("HSMM", "E089")
negative.variant.set <- annotate.negative.set(variant.file.id = "HSMM", tissue.id = "E089")
rm(list = c("variants", "negative.variant.set"))
gc()

## foetal trunk muscle
variants <- annotate.AS.DHS("HSMM", "E090")
negative.variant.set <- annotate.negative.set(variant.file.id = "HSMM", tissue.id = "E090")
rm(list = c("variants", "negative.variant.set"))
gc()

## Primary T cells
variants <- annotate.AS.DHS("HSMM", "E034")
negative.variant.set <- annotate.negative.set(variant.file.id = "HSMM", tissue.id = "E034")
rm(list = c("variants", "negative.variant.set"))
gc()



## Skin fibroblasts with muscle and mammary epithelium annotations
variants <- annotate.AS.DHS("fSkin_fibro", "E120")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E120")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E119")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E119")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E017")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E017")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E085")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E085")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E128")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E128")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E021")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E021")
rm(list = c("variants", "negative.variant.set"))
gc()

variants <- annotate.AS.DHS("fSkin_fibro", "E034")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E034")
rm(list = c("variants", "negative.variant.set"))
gc()


## Foetal muscle with HSMM annotation
variants <- annotate.AS.DHS("fMuscle", "E120")
negative.variant.set <- annotate.negative.set(variant.file.id = "fMuscle", tissue.id = "E120")
rm(list = c("variants", "negative.variant.set"))
gc()

## HSMM with HSMM annotation
variants <- annotate.AS.DHS("HSMM", "E120")
negative.variant.set <- annotate.negative.set("HSMM", tissue.id = "E120")
rm(list = c("variants", "negative.variant.set"))
gc()


# ## Foetal Muscle with HMEC annotation
# variants <- annotate.AS.DHS("fMuscle", "E119")
# negative.variant.set <- annotate.negative.set(variant.file.id = "fMuscle", tissue.id = "E119")
# rm(list = c("variants", "negative.variant.set"))

# Foetal skin with Dermal Fibroblast annotation
variants <- annotate.AS.DHS(variant.file.id = "fSkin_fibro", tissue.id = "E126")
negative.variant.set <- annotate.negative.set(variant.file.id = "fSkin_fibro", tissue.id = "E126", window.to.match = 2000, maf.cutoff = 0)
rm(list = c("variants", "negative.variant.set"))
gc()

# ## eQTL-based annotation
# variants <- annotate.eQTLs("Breast_Mammary_Tissue", "E119")
# negative.set <- annotate.negative.set("Breast_Mammary_Tissue", "E119")
# rm(list = c("variants", "negative.variant.set"))
# 
# 

## Load helper scripts
source("./lib/annotate.AS.DHS.R")
source("./lib/annotate.eQTLs.R")
source("./lib/annotate.negative.set.R")

### Annotate all tissues with sufficient ROADMAP data ###
mappings <- read.csv("./reports/ROADMAP.to.AS_DHS.map.csv", header = T, stringsAsFactors = F, comment.char = "#")
dim(mappings)
## Remove unmatched tissues
mappings <- mappings[mappings$AS.DHS.tissue != "",]
sum(mappings$AS.variants[mappings$perfect.match])
sum(mappings$AS.variants[!mappings$perfect.match])
sum(mappings$AS.variants)
sum(mappings$AS.variants[!duplicated(mappings$AS.DHS.tissue)])
sum(!duplicated(mappings$AS.DHS.tissue))
head(mappings)

## ChromHMM calls are not available for all tissues. Remove the tissues for which ChromHMM calls are not available
chromhmm.ids <- gsub("_.*", "", list.files("./cache/ENCODE/chromHMM.calls/", "mnemonics.bed.gz$"))
mappings <- mappings[mappings$tissue.id %in% chromhmm.ids,]
dim(mappings)

mappings <- mappings[mappings$AS.DHS.tissue == "vHMEC",]
nrow(mappings) - nrow(mappings[mappings$AS.DHS.tissue == "fSkin_fibro",])


## Annotate each of the matches
for (i in 1:nrow(mappings)) {
# for (i in 1:37) {
	print(paste0("Annotating ", mappings$AS.DHS.tissue[i], " with ", mappings$tissue.name[i], " features"))
	positive.variant.set <- annotate.AS.DHS(variant.file.id = mappings$AS.DHS.tissue[i], tissue.id = mappings$tissue.id[i])
	negative.variant.set <- annotate.negative.set(variant.file.id = mappings$AS.DHS.tissue[i], tissue.id = mappings$tissue.id[i], window.to.match = 2000, maf.cutoff = 0)
	rm(list = c("positive.variant.set", "negative.variant.set"))
	gc()
}


head(mappings[order(mappings$AS.variants, decreasing = T),], 30)
