require(data.table)
## Download the dsQTL and eQTL files
dir.create("./cache/degner.benchmark", showWarnings = F)
download.file(url = "http://eqtl.uchicago.edu/dsQTL_data/QTLs/GSE31388_dsQtlTable.txt.gz", destfile = "./cache/degner.benchmark/GSE31388_dsQtlTable.txt.gz")
download.file(url = "http://eqtl.uchicago.edu/dsQTL_data/QTLs/GSE31388_dsQtlTableLong.txt.gz", destfile = "./cache/degner.benchmark/GSE31388_dsQtlTableLong.txt.gz")
download.file(url = "http://eqtl.uchicago.edu/dsQTL_data/QTLs/GSE31388_eQtlTable.txt.gz", destfile = "./cache/degner.benchmark/GSE31388_eQtlTable.txt.gz")


## Read in the list of dsQTLs
dsqtls <- read.table("./cache/degner.benchmark/GSE31388_dsQtlTableLong.txt.gz", sep = "\t", header = T, stringsAsFactors = F)
dsqtls$varID <- paste(dsqtls$Chr, dsqtls$SNP, sep = ".")
dsqtls$ALT <- dsqtls$REF <- as.character(NA)
head(dsqtls)
tail(dsqtls)

for (chr in unique(dsqtls$Chr)) {
	print(chr)
	## Download SNP for the chromosome, unless exists
	if (!file.exists(paste0("./cache/degner.benchmark/", chr, ".YRI.snpdata.txt.gz"))) {
		download.file(url = paste0("http://eqtl.uchicago.edu/dsQTL_data/GENOTYPES/", chr, ".YRI.snpdata.txt.gz"), destfile = paste0("./cache/degner.benchmark/", chr, ".YRI.snpdata.txt.gz"))
	}
	## Read in SNP date
	snpdt <- fread(paste0("gunzip -c ./cache/degner.benchmark/", chr, ".YRI.snpdata.txt.gz"), sep = "\t", header = T, stringsAsFactors = F, skip = 1)
	class(snpdt$chr) <- "character"
	snpdt$chr <- chr
	snpdt$varID <- paste(snpdt$chr, snpdt$pos, sep = ".")
	chr.snps <- dsqtls[dsqtls$Chr == chr,]
	chr.snps$REF <- snpdt$A[match(chr.snps$varID, snpdt$varID)]
	chr.snps$ALT <- snpdt$B[match(chr.snps$varID, snpdt$varID)]
	dsqtls$REF[dsqtls$Chr == chr] <- chr.snps$REF
	dsqtls$ALT[dsqtls$Chr == chr] <- chr.snps$ALT
}

## Check if the dsQTL is internal to the region it affects


## The snp data is in hg18 lift it over to hg19
coords.hg18 <- dsqtls[,c(1, 4, 4, 9)]
coords.hg18$SNP <- coords.hg18$SNP-1
write.table(coords.hg18, "dsqtls.hg18.bed", quote = F, sep = "\t", row.names = F, col.names = F)
source("~/tools/liftOver/general.liftOver.script.R")
liftOver(fname = "dsqtls.hg18.bed", from = "hg18", to = "hg19")

coords.hg19 <- read.table("dsqtls.hg19.bed", col.names = c("chr", "start", "end", "varID"))
head(coords.hg19)
head(dsqtls)
## There are some variants that don't transfer, remove them
dsqtls <- dsqtls[dsqtls$varID %in% coords.hg19$varID,]
dsqtls$pos <- coords.hg19$end[match(dsqtls$varID, coords.hg19$varID)]

## Remove temporary files
file.remove(c("dsqtls.hg18.bed", "dsqtls.hg19.bed", "unmapped.bed"))

## Check that either reference or alternate allele matches actual 
require(BSgenome.Hsapiens.UCSC.hg19)
variants <- GRanges(dsqtls$Chr, IRanges(dsqtls$pos, dsqtls$pos), REF = dsqtls$REF, ALT = dsqtls$ALT)
variants$real.ref <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, variants))
dsqtls$real.ref <- variants$real.ref
summary(dsqtls$real.ref == dsqtls$REF | dsqtls$real.ref == dsqtls$ALT)
## There are 364 variants where none of the two alleles matches reference. Remove them.
dsqtls <- dsqtls[(dsqtls$real.ref == dsqtls$REF | dsqtls$real.ref == dsqtls$ALT),]
## For the variants where ALT allele matches the reference, swap REF and ALT alleles
dsqtls$ALT[dsqtls$real.ref == dsqtls$ALT] <- dsqtls$REF[dsqtls$real.ref == dsqtls$ALT]
dsqtls$REF <- dsqtls$real.ref
dsqtls <- dsqtls[,-(which(colnames(dsqtls) == "real.ref"))]

## Write a vcf file with the dsQTLs
write.table(dsqtls[,c(1, 12, 10, 11, 2, 3, 5, 6, 7, 8, 9)], "./cache/degner.benchmark/dsQTLs.vcf", sep = "\t", quote = F, row.names = F)


source("./lib/annotate.and.predict.vcf.R")
annotate.vcf(path.to.vcf = "./cache/degner.benchmark/dsQTLs.vcf", tissue.id = "E116")


anno.dsqtls <- read.csv("./cache/degner.benchmark/dsQTLs.E116.annotated.csv", stringsAsFactors = F)
head(anno.dsqtls)
summary(anno.dsqtls$jaspar.abs.score == 0)
hist(as.numeric(anno.dsqtls$p.regulatory), breaks = 100)
table(anno.dsqtls$DNase.macs2.narrowPeak.bin)
summary(as.numeric(anno.dsqtls$p.regulatory) > 0.65)

anno.vcf <- read.table("./cache/degner.benchmark/dsQTLs.annotated.vcf", col.names = c("chr", "pos", "REF", "ALT", "start", "end", "slope", "intersect", "t.stat", "pVal", "p.regulatory"))
head(anno.vcf)
smoothScatter(-log10(anno.vcf$pVal), as.numeric(anno.vcf$p.regulatory))
summary(anno.vcf$pos >= anno.vcf$start & anno.vcf$pos <= anno.vcf$end)
