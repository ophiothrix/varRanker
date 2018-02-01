#### phyloP scores ####

annotate.phyloP <- function(variants) {
	require(curl)
	require(BSgenome.Hsapiens.UCSC.hg19)
	require(rtracklayer)
	
	si <- Seqinfo(seqnames=seqnames(Hsapiens), seqlengths=seqlengths(Hsapiens))
	si <- si[paste0("chr", c(1:22, "X", "Y"))]
	
	dname <- "./cache/conservation.scores/phyloP/"
	dir.create(dname, showWarnings = F)
	
	
	### 100 way alignment
	urlname <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP100way/hg19.100way.phyloP100way.bw"
	fname <- basename(urlname)
	if (!file.exists(paste0(dname, fname))) {
		## Download the file
		download.file(urlname, paste0(dname, fname), method = "curl")
	}
	variants$phyloP100 <- as.numeric(NaN)
	
	print("Getting 100 way phyloP scores for all variants")
	
	file.remove(c("positions.bed", "scores.txt"))
	write.table(as.data.frame(variants)[,1:3], "positions.bed", col.names = F, row.names = F, quote = F, sep = "\t")
	system(paste0("./lib/ExtractBigWigScore.py ", dname, fname))
	scores <- read.table("scores.txt", stringsAsFactors = F)[,1]
	scores <- gsub("\\[", "", scores)
	scores <- as.numeric(gsub("\\]", "", scores))
	variants$phyloP100 <- scores
	file.remove(c("positions.bed", "scores.txt"))
	
	variants$phyloP100[is.nan(variants$phyloP100)] <- 0
	gc()
	
	
	#### 46 way alignment Vertebrates
	urlname <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/vertebrate/"
	variants$phyloP46.vertebrate <- as.numeric(NaN)
	
	for (chr in seqlevels(variants)) {
		print(paste0("Getting vertebrate phyloP scores for ", chr))
		fname <- paste0(chr, ".phyloP46way.wigFix.gz")
		
		if (!file.exists(paste0(dname, gsub("wigFix.gz", "bw", fname)))) {
			## Download the file
			download.file(paste0(urlname, fname), paste0(dname, fname), method = "curl")
			
			## Convert to bigWig
			wigToBigWig(paste0(dname, fname), si)
			
			## Remove the fixed step Wig file
			file.remove(paste0(dname, fname))
			gc()
		}
		
		
		file.remove(c("positions.bed", "scores.txt"))
		write.table(as.data.frame(variants[seqnames(variants) == chr])[,1:3], "positions.bed", col.names = F, row.names = F, quote = F, sep = "\t")
		system(paste0("./lib/ExtractBigWigScore.py ", dname, gsub("wigFix.gz", "bw", fname)))
		scores <- read.table("scores.txt", stringsAsFactors = F)[,1]
		scores <- gsub("\\[", "", scores)
		scores <- as.numeric(gsub("\\]", "", scores))
		variants[seqnames(variants) == chr]$phyloP46.vertebrate <- scores
		file.remove(c("positions.bed", "scores.txt"))
	}
	variants$phyloP46.vertebrate[is.nan(variants$phyloP46.vertebrate)] <- 0
	gc()
	
	
	#### 46 way alignment Mammals
	urlname <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/placentalMammals/"
	variants$phyloP46.mammals <- as.numeric(NaN)
	
	for (chr in seqlevels(variants)) {
		print(paste0("Getting placental mammals phyloP scores for ", chr))
		fname <- paste0(chr, ".phyloP46way.placental.wigFix.gz")
		
		if (!file.exists(paste0(dname, gsub("wigFix.gz", "bw", fname)))) {
			## Download the file
			download.file(paste0(urlname, fname), paste0(dname, fname), method = "curl")
			
			## Convert to bigWig
			wigToBigWig(paste0(dname, fname), si)
			
			## Remove the fixed step Wig file
			file.remove(paste0(dname, fname))
			gc()
		}
		
		
		file.remove(c("positions.bed", "scores.txt"))
		write.table(as.data.frame(variants[seqnames(variants) == chr])[,1:3], "positions.bed", col.names = F, row.names = F, quote = F, sep = "\t")
		system(paste0("./lib/ExtractBigWigScore.py ", dname, gsub("wigFix.gz", "bw", fname)))
		scores <- read.table("scores.txt", stringsAsFactors = F)[,1]
		scores <- gsub("\\[", "", scores)
		scores <- as.numeric(gsub("\\]", "", scores))
		variants[seqnames(variants) == chr]$phyloP46.mammals <- scores
		file.remove(c("positions.bed", "scores.txt"))
	}
	variants$phyloP46.mammals[is.nan(variants$phyloP46.mammals)] <- 0
	gc()
	
	#### 46 way alignment Primates
	urlname <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/phyloP46way/primates/"
	variants$phyloP46.primates <- as.numeric(NaN)
	
	for (chr in seqlevels(variants)) {
		print(paste0("Getting primates phyloP scores for ", chr))
		fname <- paste0(chr, ".phyloP46way.primate.wigFix.gz")
		
		if (!file.exists(paste0(dname, gsub("wigFix.gz", "bw", fname)))) {
			## Download the file
			download.file(paste0(urlname, fname), paste0(dname, fname), method = "curl")
			
			## Convert to bigWig
			wigToBigWig(paste0(dname, fname), si)
			
			## Remove the fixed step Wig file
			file.remove(paste0(dname, fname))
			gc()
		}
		
		
		file.remove(c("positions.bed", "scores.txt"))
		write.table(as.data.frame(variants[seqnames(variants) == chr])[,1:3], "positions.bed", col.names = F, row.names = F, quote = F, sep = "\t")
		system(paste0("./lib/ExtractBigWigScore.py ", dname, gsub("wigFix.gz", "bw", fname)))
		scores <- read.table("scores.txt", stringsAsFactors = F)[,1]
		scores <- gsub("\\[", "", scores)
		scores <- as.numeric(gsub("\\]", "", scores))
		variants[seqnames(variants) == chr]$phyloP46.primates <- scores
		file.remove(c("positions.bed", "scores.txt"))
	}
	variants$phyloP46.primates[is.nan(variants$phyloP46.primates)] <- 0
	gc()
	
	return(variants)
}


variants <- annotate.phyloP(variants)
