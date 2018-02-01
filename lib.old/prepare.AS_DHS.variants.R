#### Take the list of allele-specific DNase hypersensitive site variants from Maurano et al paper and converts it to GRanges object.
### This object will be used to extract DHS portion of the negative training set as well as to assign reference and alternate alleles to the positive training set.
### Since other negative training set variants (most notably location matched variants) have their reference allele synchronised with human reference assembly, we need to do the same for the AS DHS variant alleles.
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)

## Load the set and convert to a GRanges object
all.as.dhs <- GRanges(read.table("./data/ng.3432-S5.txt", header = T, stringsAsFactors = F))
## Change to 1-based coordinates
start(all.as.dhs) <- start(all.as.dhs) + 1

## Add reference and alternative alleles
all.as.dhs$REF <- all.as.dhs$allele.1
all.as.dhs$ALT <- all.as.dhs$allele.2

## Check if allele.1 matches reference assembly
refs <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, all.as.dhs))
summary(refs == all.as.dhs$REF) ## It's about 50/50 split with only half of the variants REF allele matching reference assembly.

## For the mismatched variants, check if the alternate allele matched reference assembly
summary(refs[refs != all.as.dhs$REF] == all.as.dhs$ALT[refs != all.as.dhs$REF])

## Turns out they do. For mismatched variants swap reference and alternate alleles
alts <- NA
alts[refs == all.as.dhs$REF] <- all.as.dhs$ALT[refs == all.as.dhs$REF]
alts[refs != all.as.dhs$REF] <- all.as.dhs$REF[refs != all.as.dhs$REF]

new.refs <- NA
new.refs[refs == all.as.dhs$REF] <- all.as.dhs$REF[refs == all.as.dhs$REF]
new.refs[refs != all.as.dhs$REF] <- all.as.dhs$ALT[refs != all.as.dhs$REF]

## Final check. Do new reference alleles match the reference assembly?
summary(new.refs == refs) ## ALL true! Add them back to the variants.
all.as.dhs$REF <- new.refs
all.as.dhs$ALT <- alts

## Save the intermediate file containing all the variants
saveRDS(all.as.dhs, file = "./cache/allele.biased.DHS.variants.all.rds")

