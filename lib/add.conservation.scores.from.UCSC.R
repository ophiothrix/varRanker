#### GERP score ####

conservation.score.annotation <- function(variants) {
	source("./lib/annotate.variants.from.GRanges.R")
	if (!file.exists("./data/conservation.scores/dbNsfpGerpRs.rds")) {
		source("./lib/download.GERP.score.R")
		cons.scores <- readRDS("./data/conservation.scores/dbNsfpGerpRs.rds")
	} else {
		cons.scores <- readRDS("./data/conservation.scores/dbNsfpGerpRs.rds")
	}
	gerp.score <- annotate.variants.from.GRanges(variants = variants, regions = cons.scores)
	
	variants$GERP.score <- do.call(rbind, strsplit(scores, "\t"))[,2]
	return(variants)
}
