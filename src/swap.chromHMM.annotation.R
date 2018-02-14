### Replace 18-state ChromHMM annotation with 15-state one, that does not rely on H3K27ac data availability

fnames <- list.files("./cache", "^E[0-9]+", full.names = T)
for (i in 1:length(fnames)) {
	print(i)
	fname <- fnames[i]
	print(fname)
	variants <- readRDS(fname)
	print("Adding tissue-specific ChromHMM annotation")
	state <- list.files("./cache/ENCODE/chromHMM.core.calls", paste0(substr(basename(fname), 1, 4), "_15_core"), full.names = T)
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
	saveRDS(variants, fname)
	rm(list=c("state", "variants"))
	gc()
}
