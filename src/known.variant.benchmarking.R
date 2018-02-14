rm(list = ls())
gc()
source("./lib/annotate.and.predict.vcf.R")
annotate.vcf(path.to.vcf = "./cache/terts.vcf", tissue.id = "E007")

annotate.vcf(path.to.vcf = "./cache/fto.proxies.vcf", tissue.id = "E118")


variants
getSeq(BSgenome.Hsapiens.UCSC.hg19, GRanges(variants))


## Variants in LD with  obtained from https://analysistools.nci.nih.gov/LDlink/?var=rs1421085&pop=CEU%2BTSI%2BFIN%2BGBR%2BIBS&r2_d=r2&tab=ldproxy
ftos <- read.table("./cache/fto.proxy.txt", header = T, stringsAsFactors = F) %>%
	filter(R2 >= 0.8) %>%
	mutate(chr = gsub(":.*", "", Coord)) %>%
	mutate(pos = gsub(".*:", "", Coord)) %>%
	mutate(REF = gsub("\\((.*)/.*", "\\1", Alleles)) %>%
	mutate(ALT = gsub(".*/(.*)\\)", "\\1", Alleles)) %>%
	filter(ALT != "-" & REF != "-") %>%
	select(c(chr, pos, RS_Number, REF, ALT, MAF))
head(ftos)
dim(ftos)
write.table(ftos, "./cache/fto.proxies.vcf", sep = "\t", quote = F, row.names = F)
