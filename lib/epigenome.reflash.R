### Given an annotated positive variant set, reannotate it with the epigenetic context from a specified tissue id.
# path.to.positive.set <- "./cache/E120.annotated.HSMM.variants.rds"
# tissue.id = "E034"
epigenome.reflash <- function(path.to.positive.set, tissue.id) {
	source("./lib/epigenome.annotation.R")
	variants <- readRDS(path.to.positive.set)
	variants <- epigenome.annotation(variants = variants, tissue.id = tissue.id)
	
	## Now we have a lot of double columns, which we need to remove (old annotations)
	## New annotations get assigned ".1" at the end
	anno.mat <- as.data.frame(mcols(variants))
	new.annos <- grep("\\.1$", colnames(anno.mat), value = T)
	old.annos <- which(colnames(anno.mat) %in% gsub("\\.1$", "", new.annos))
	anno.mat <- anno.mat[,-old.annos]
	colnames(anno.mat) <- gsub("\\.1$", "", colnames(anno.mat))
	mcols(variants) <- anno.mat
	
	## Save the reflashed object
	saveRDS(variants, gsub("E[0-9][0-9][0-9]", tissue.id, path.to.positive.set))
	
	## Repeat for the negative set
	path.to.negative.set <- gsub(".variants.rds", ".negative.set.rds", path.to.positive.set)
	variants <- readRDS(path.to.negative.set)
	variants <- epigenome.annotation(variants = variants, tissue.id = tissue.id)
	
	## Now we have a lot of double columns, which we need to remove (old annotations)
	## New annotations get assigned ".1" at the end
	anno.mat <- as.data.frame(mcols(variants))
	new.annos <- grep("\\.1$", colnames(anno.mat), value = T)
	old.annos <- which(colnames(anno.mat) %in% gsub("\\.1$", "", new.annos))
	anno.mat <- anno.mat[,-old.annos]
	colnames(anno.mat) <- gsub("\\.1$", "", colnames(anno.mat))
	mcols(variants) <- anno.mat
	
	## Save the reflashed object
	saveRDS(variants, gsub("E[0-9][0-9][0-9]", tissue.id, path.to.negative.set))
}
