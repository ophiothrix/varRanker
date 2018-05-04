### Reassign damage scores to eQTLs.

### Inside the eQTL file every given variant can have multiple associated genes. This seems to be throwing the script, as it now assigns multiple scores to each variant.

### In reality for a given chr:pos:REF:ALT combination, there may only be one score. This is something to take into account during the run. Notably, in this case there are 113094 duplicated variants (about a quarter of all). Which is a waste of computational resources
summary(duplicated(original.variants$variant_id))

variants <- readRDS("./cache/E119.annotated.allSigVars.variants.rds")
head(variants)
summary(duplicated(variants$varID))
variants <- variants[!duplicated(variants$varID)]
length(variants)

variants$varID <- paste0(variants$varID, ":", variants$REF, ":", variants$ALT)
variants$varID <- paste(as.character(seqnames(variants)), start(variants), variants$REF, variants$ALT, sep = ":")


original.variants <- readRDS("./cache/original.variants.rds")
dim(original.variants)
original.variants$varID <- 	paste(original.variants$chr, original.variants$start, original.variants$REF, original.variants$ALT, sep = ":")
original.variants$varID


## Reassign damage scores to variants
summary(duplicated(match(original.variants$varID, variants$varID)))
original.variants$p.regulatory <- variants$p.regulatory[match(original.variants$varID, variants$varID)]
head(original.variants)

## There need to be some checks. I.e. if the variant (in vcf) is multi-alleling, we can either collapse the scores for all alleles, or keep them expanded. Need to rethink.
