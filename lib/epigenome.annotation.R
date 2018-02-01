#### A function that takes a list of variants (as a GRanges object) and a ROADMAP tissue ID and annotates the above variants with epigenetic features from the given tissue ####
## tissue.id <- "E119"

epigenome.annotation <- function(variants, tissue.id) {
	### Download the ROADMAP data ###
	source("./lib/download.ROADMAP.data.R")
	download.ROADMAP.data(tissue.id)
	
	### Annotate the variants with ENCODE data
	annos <- list.files("./cache/ENCODE", pattern = tissue.id, full.names = T)
	imputed.annos <- annos[grep("imputed", annos)]
	annos <- annos[-grep("imputed", annos)]
	anno.matrix <- matrix(0, nrow = length(variants), ncol = length(annos) * 4)
	colnames(anno.matrix) <- paste(rep(gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", annos), each = 4), c("bin", "FC", "pval", "qval"), sep = ".")
	imputed.anno.matrix <- matrix(0, nrow = length(variants), ncol = length(imputed.annos))
	colnames(imputed.anno.matrix) <- paste(gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", imputed.annos), "bin", sep = ".")
	
	## Add "local distance to summit" columns for narrow peak marks
	dist.to.summit.local <- matrix(1000, nrow = length(variants), ncol = length(grep("narrow", annos)))
	colnames(dist.to.summit.local) <- paste(gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", annos[grep("narrow", annos)]), "dist.to.summit.local", sep = ".")
	
	## Add "global distance to summit" columns for narrow peak marks
	dist.to.summit.global <- matrix(NA, nrow = length(variants), ncol = length(grep("narrow", annos)))
	colnames(dist.to.summit.global) <- paste(gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", annos[grep("narrow", annos)]), "dist.to.summit.global", sep = ".")
	
	
	## Combine with the rest of the annotation table
	anno.matrix <- cbind(anno.matrix, dist.to.summit.local, dist.to.summit.global, imputed.anno.matrix)
	
	## Set default values for FC to 1
	anno.matrix[, grep("FC", colnames(anno.matrix))] <- 1
	
	# dim(anno.matrix)
	# head(anno.matrix, 23)
	
	for (anno in annos) {
		print(anno)
		anno.name <- gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", anno)
		print(anno.name)
		## Import a peaks file
		if (length(grep("broad", anno)) == 1) { 
			## Process broad peak file
			peaks <- GRanges(read.table(anno, sep = "\t", col.names = c("chr", "start", "end", "rank", "score", "strand", "FC", "pval", "qval"), stringsAsFactors = F))
		} else { 
			## Process narrow peak file
			if (length(grep("narrow", anno)) == 1) { ## assign correct column names for narrow peak file
				peaks <- GRanges(read.table(anno, sep = "\t", col.names = c("chr", "start", "end", "rank", "score", "strand", "FC", "pval", "qval", "summit"), stringsAsFactors = F))
				peaks$summit.coord <- start(peaks) + peaks$summit
			} else { 
				## Process gapped peak file
				gapped.peaks <- read.table(anno, sep = "\t", col.names = c("chr", "start", "end", "rank", "x10.qval", "id", "narrow.start", "narrow.end", "colour", "n.peaks", "block.length", "block.start", "FC", "pval", "qval"), stringsAsFactors = F)
				
				## make a table of narrow peak coordinates from the gapped peaks
				narrow.peaks <- as.data.frame(matrix(NA, nrow = sum(gapped.peaks$n.peaks), ncol = 6))
				colnames(narrow.peaks) <- c("chr", "start", "end", "FC", "pval", "qval")
				head(narrow.peaks, 20)
				narrow.peaks$chr <- rep(gapped.peaks$chr, gapped.peaks$n.peaks)
				narrow.peaks$FC <- rep(gapped.peaks$FC, gapped.peaks$n.peaks)
				narrow.peaks$pval <- rep(gapped.peaks$pval, gapped.peaks$n.peaks)
				narrow.peaks$qval <- rep(gapped.peaks$qval, gapped.peaks$n.peaks)
				narrow.peaks$start <- rep(gapped.peaks$start, gapped.peaks$n.peaks)
				narrow.peaks$end <- rep(gapped.peaks$end, gapped.peaks$n.peaks)
				starts <- as.numeric(unlist(strsplit(gapped.peaks$block.start, ",")))
				widths <- as.numeric(unlist(strsplit(gapped.peaks$block.length, ",")))
				narrow.peaks$end <- narrow.peaks$start + starts + widths
				narrow.peaks$start <- narrow.peaks$start + starts
				
				## Convert to GRanges
				peaks <- GRanges(narrow.peaks)
			}
		}
		peaks.meta <- mcols(peaks)
		olaps <- findOverlaps(variants, peaks)
		olaps
		for (metric in c("FC", "pval", "qval")) {
			print(metric)
			anno.matrix[queryHits(olaps), paste0(anno.name, ".", metric)] <- peaks.meta[subjectHits(olaps), metric]
		}
		## If the track file is for narrow peak, add distance to summit annotation
		if (length(grep("narrow", anno)) == 1) {
			print("Distance to summit")
			## Add local (within a peak) distance
			anno.matrix[queryHits(olaps), paste0(anno.name, ".", "dist.to.summit.local")] <- (abs(start(variants[queryHits(olaps)]) - peaks.meta$summit.coord[subjectHits(olaps)])) - 1
			## Add global (nearest summit) distance
			summit.anno <- GRanges(seqnames(peaks), IRanges(peaks$summit.coord, peaks$summit.coord))
			dist.to.sum <- as.data.frame(distanceToNearest(variants, summit.anno))
			anno.matrix[dist.to.sum$queryHits, paste0(anno.name, ".", "dist.to.summit.global")] <- dist.to.sum$distance
		# head(anno.matrix[,grep("DNase", colnames(anno.matrix))])
		}
		
		## Add binary yes / no annotation for each mark
		anno.matrix[queryHits(olaps), paste0(anno.name, ".bin")] <- 1
	}
	
	### Add imputed marks
	for (anno in imputed.annos) {
		print(anno)
		anno.name <- gsub(paste0(".*", tissue.id, "-(.*)\\..*"), "\\1", anno)
		print(anno.name)
		## Load peak coordinates
		peaks  <- GRanges(read.table(anno, sep = "\t", col.names = c("chr", "start", "end", "peakType"), stringsAsFactors = F))
		## Find overlaps between variants and annotations
		olaps <- findOverlaps(variants, peaks)
		## Mark variants overlapping the annotations
		anno.matrix[queryHits(olaps), paste0(anno.name, ".bin")] <- 1
	}
	
	### Add H3K4me1/me3 ratio for FC and p val
	anno.matrix <- as.data.frame(anno.matrix)
	anno.matrix$H3K4me1.me3.ratio.broad.FC <- anno.matrix$H3K4me1.broadPeak.FC / anno.matrix$H3K4me3.broadPeak.FC
	anno.matrix$H3K4me1.me3.ratio.broad.qval <- (anno.matrix$H3K4me1.broadPeak.qval + 0.1) / (anno.matrix$H3K4me3.broadPeak.qval + 0.1)
	anno.matrix$H3K4me1.me3.ratio.broad.pval <- (anno.matrix$H3K4me1.broadPeak.pval + 0.1) / (anno.matrix$H3K4me3.broadPeak.pval + 0.1)
	anno.matrix$H3K4me1.me3.ratio.narrow.FC <- anno.matrix$H3K4me1.narrowPeak.FC / anno.matrix$H3K4me3.narrowPeak.FC
	anno.matrix$H3K4me1.me3.ratio.narrow.qval <- (anno.matrix$H3K4me1.narrowPeak.qval + 0.1) / (anno.matrix$H3K4me3.narrowPeak.qval + 0.1)
	anno.matrix$H3K4me1.me3.ratio.narrow.pval <- (anno.matrix$H3K4me1.narrowPeak.pval + 0.1) / (anno.matrix$H3K4me3.narrowPeak.pval + 0.1)
	# head(anno.matrix)
	
	## Assign the annotations back to the GRanges object
	mcols(variants) <- cbind(mcols(variants), anno.matrix)
	
	
	### Add tissue-specific ChromHMM and GEP states ###
	# Add promoter and enhancer annotations -----------------------------------
	# source("~/tools/R.functions.library/annotate.variants.from.GRanges.R")
	
	## Add annotation for different chromHMM classes for a given tissue
	print("Adding tissue-specific ChromHMM annotation")
	state <- list.files("./cache/ENCODE/chromHMM.calls", paste0(tissue.id, "_18_core"), full.names = T)
	if (length(state) == 0) {
		stop("There is no ChromHMM annotation available for this tissue ID")
	}
	
	print(state)
	state <- GRanges(read.table(state, sep = "\t", col.names = c("chr", "start", "end", "state"), stringsAsFactors = F))
	# head(state)
	olaps <- findOverlaps(variants, state)
	## Remove duplicated query hits, just in case (this will only leave the first hit)
	olaps <- olaps[!duplicated(queryHits(olaps))]
	variants$chromHMM.state <- state[subjectHits(olaps)]$state
	rm(list=c("state"))
	gc()
	
	### Add annotation for IDEAS chromatin state calls
	## Add annotation for different chromHMM classes for a given tissue
	state <- list.files("./cache/ENCODE/IDEAS", tissue.id, full.names = T)
	print(state)
	state <- GRanges(read.table(state, sep = "\t", col.names = c("chr", "start", "end", "state", "width", "strand", "start2", "end2", "colour"), stringsAsFactors = F))
	
	olaps <- findOverlaps(variants, state)
	## Remove duplicated query hits, just in case (this will only leave the first hit)
	olaps <- olaps[!duplicated(queryHits(olaps))]
	## It seems that occasionally some variants don't get IDEAS annotation - there's no obvious pattern other than all of them come from "random" sampling
	## Check if all variants got IDEAS annotation. If not, annotate those that have a match and remove those that don't
	if (all(1:length(variants) %in% queryHits(olaps))) {
		variants$IDEAS.state <- state[subjectHits(olaps)]$state
	} else {
		variants$IDEAS.state <- as.character(NA)
		variants[queryHits(olaps)]$IDEAS.state <- state[subjectHits(olaps)]$state
		## Remove variants without IDEAS annotation
		variants <- variants[!is.na(variants$IDEAS.state)]
	}
	rm(list=c("state"))
	gc()
	
	## Return annotated variants
	return(variants)
}
