#### GERP score ####

conservation.score.annotation <- function(variants) {
	## Adding GERP conservation score
	## If there are no conservation score files, download GERP scores
	if (length(list.files("./data/conservation.scores/GERP", "maf.rates.gz$")) != 25) {
		print("GERP score file is missing. Downloading it...")
		source("./lib/download.GERP.score.R")
	}
	score.dir <- "./data/conservation.scores/GERP/"
	
	# variants <- readRDS("./data/play.training.set.rds")
	# variants <- variants[sample(1:length(variants), 100, replace = F)]
	
	## Check that the variants object is sorted
	if (any(order(variants) != 1:length(variants))) {
		stop("Variant object is not sorted.")
	}
	
	scores <- character()
	for (chr in seqlevels(variants)) {
		print(paste0("Extracting GERP score from ", chr))
		if (length(variants[seqnames(variants) == chr]) == 0) {
			next
		}
		## write the file with lines to be extracted
		write.table(paste0(start(variants[seqnames(variants) == chr]), "p"), file = "linesfile", quote = F, col.names = F, row.names = F)
		write.table(paste0(start(variants[seqnames(variants) == chr])[sum(seqnames(variants) == chr)]+1, "q"), file = "linesfile", quote = F, col.names = F, row.names = F, append = T)
		x <- system(paste0("gunzip -c ", score.dir, chr, ".maf.rates.gz | sed -n -f linesfile"), intern = T)
		
		if (length(x) != sum(seqnames(variants) == chr)) {
			stop("The number of extracted scores is different from the number of input lines")
		} else {
			scores <- c(scores, x)
		}
	}
	
	score <- do.call(rbind, strsplit(scores, "\t"))[,2]
	return(score)
	file.remove("linesfile")
}
