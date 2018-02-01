## Download GERP score from UCSC browser and store them as a GRanges object
# rm(list = ls())
# gc()
dir.create("./cache/conservation.scores/GERP", showWarnings = F, recursive = T)
require(curl)
require(GenomicRanges)
require(BSgenome.Hsapiens.UCSC.hg19)
require(rtracklayer)

si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
si <- si[paste0("chr", c(1:22, "X", "Y", "M"))]

## Download the full GERP file
if (!file.exists("./cache/conservation.scores/GERP/hg19.GERP_scores.tar.gz")) {
	download.file(url = "http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz", destfile = "./cache/conservation.scores/GERP/hg19.GERP_scores.tar.gz", quiet = F, method = "curl")
}

## Unpack the archive into individual files
system("tar -xvf ./cache/conservation.scores/GERP/hg19.GERP_scores.tar.gz -C ./cache/conservation.scores/GERP/")
file.remove("./cache/conservation.scores/GERP/hg19.GERP_scores.tar.gz")
gc()
## Extracting the score is far more efficient if we first convert it to a BigWig format
for (fname in list.files("./cache/conservation.scores/GERP", "maf.rates$", full.name = T)) {
	print(fname)
	chname <- gsub(".maf.rates", "", basename(fname))
	wigname <- gsub(".maf.rates", ".wigFix.gz", fname)
	## Convert to fixed step Wig
	system(paste0('printf "\tfixedStep chrom=', chname,' start=1 step=1\n" > ', fname, ".tmp"))
	system(paste0("cut -f2 ", fname, ".tmp ", fname, " | gzip -c > ", wigname))
	gc()
	## Convert to bigWig
	wigToBigWig(wigname, si)
	## Remove intermediate files
	file.remove(fname)
	file.remove(paste0(fname, ".tmp"))
	file.remove(wigname)
	gc()
}
