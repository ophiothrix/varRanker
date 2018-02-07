## Script to extract available cato scores given a GRanges object

variants <- readRDS("./cache/E120.annotated.HSMM.variants.rds")


for (fls in list.files("./cache", "*fMuscle.*rds", full.names = T)) {
	print(fls)
	vrnts <- readRDS(fls)
	flsout <- gsub("annotated", "cato.annotated", fls)
	flsout <- gsub("cache", "cache/cato.annotated.variants", flsout)
	vrnts <- extract.cato.score(variants = vrnts)
	saveRDS(vrnts, flsout)
}


require(BSgenome.Hsapiens.UCSC.hg19)
valid.GR <- GRanges(gsub("(.*):.*:.*", "\\1", valid.set$varID), IRanges(as.numeric(gsub(".*:(.*):.*", "\\1", valid.set$varID)), width = 1), ALT = gsub(".*:.*:(.*)", "\\1", valid.set$varID))
valid.GR$REF <- as.character(getSeq(x = BSgenome.Hsapiens.UCSC.hg19, valid.GR))
valid.GR$varID <- valid.set$varID
valid.GR$DNase.macs2.narrowPeak.bin <- valid.set$DNase.macs2.narrowPeak.bin
valid.GR <- extract.cato.score(valid.GR)
summary
valid.set.cato <- valid.set
valid.set.cato$cato.score <- NA
valid.set.cato$cato.score[match(valid.GR$varID, valid.set.cato$varID)] <- valid.GR$cato.score
valid.set.cato <- valid.set.cato[!is.na(valid.set.cato$cato.score),]
head(valid.set.cato)

extract.cato.score <- function(variants) {
	require(GenomicRanges)
	system("rm variant.scores.tsv variants.file.tsv")
	variant.table <- as.data.frame(variants)
	variant.table <- variant.table[order(variant.table$seqnames, variant.table$start),]
	variant.table$coords <- paste(variant.table$start, variant.table$end, sep="-")
	# head(variant.table)
	# summary(duplicated(variant.table))
	## Remove duplicated records (multiple alternative alleles) - no need to look them up twice
	variants.file <- variant.table[,c("seqnames", "coords")]
	variants.file <- variants.file[!duplicated(variants.file),]
	write.table(variants.file, file="variants.file.tsv", row.names=F, col.names=F, quote=F, sep=":")
	
	# Get cato scores
	print("Extracting CATO scores")
	system("cat variants.file.tsv | xargs -I {} /Users/aholik/tools/miniconda2/bin/tabix ~/Downloads/dbSNP142.CATO.V1.1.txt.gz {} > variant.scores.tsv")
	
	variant.scores <- read.table("variant.scores.tsv", header = F, stringsAsFactors = F, col.names =  c("chr", "start", "end", "rsid", "score", "strand", "motif", "position", "ref", "alt", "cell.lines"))
	variant.scores$start <- variant.scores$end
	head(variant.scores)
	variant.scores <- GRanges(variant.scores$chr, IRanges(variant.scores$start, width = 1), 
							  REF = variant.scores$ref, ALT = variant.scores$alt, score = variant.scores$score)
	
	olaps <- findOverlaps(variants, variant.scores)
	matched.SNVs <- variant.scores[subjectHits(olaps)]
	matched.SNVs$variant.id <- queryHits(olaps)
	matched.SNVs$original.ALT <- variants$ALT[queryHits(olaps)]
	summary(matched.SNVs$ALT == matched.SNVs$original.ALT)
	matched.SNVs <- matched.SNVs[matched.SNVs$ALT == matched.SNVs$original.ALT]
	variants$cato.score <- as.numeric(NA)
	variants$cato.score[matched.SNVs$variant.id] <- matched.SNVs$score
	summary(is.na(variants$cato.score))
	
	## Variants that fall outside of DNase sites should be assigned CATO score of 0.
	variants$cato.score[variants$DNase.macs2.narrowPeak.bin == 0] <- 0
	
	## Remove variants without score annotation
	variants <- variants[!is.na(variants$cato.score)]
	
	system("rm variant.scores.tsv")
	return(variants)
}
