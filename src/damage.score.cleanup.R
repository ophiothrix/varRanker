## Take the annotated variants containing all damage scores and clean them up, so that for each of the scores, all variants missing that score are removed.

all.scores <- list.files("./cache/all.score.annotated.variants", full.names = T)
scores <- c("cadd", "eigen", "funseq", "gwava")

for (score in scores) {
	print(score)
	dirname <- paste0("./cache/", score)
	dir.create(dirname, showWarnings = F)
	for (fname in all.scores) {
		print(fname)
		variant.set <- readRDS(fname)
		df <- as.data.frame(mcols(variant.set))
		col.id <- grep(score, colnames(df))[1]
		na.ids <- is.na(df[,col.id])
		variant.set <- variant.set[!na.ids]
		outname <- gsub("all.score.annotated.variants", score, fname)
		saveRDS(variant.set, outname)
	}
}
