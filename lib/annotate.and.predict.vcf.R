#### Annotate variants from a VCF file ####

### Function to annotate variants from a vcf file with a given tissue-specific annotation and convert them into a dataframe compatible with prediction.
# path.to.vcf <- "./cache/BRCA-US.vsmall.vcf.gz"
# tissue.id <- "E119"
# path.to.vcf <- "./cache/starrseq.AS.variants.vcf"
# tissue.id <- "E199"
annotate.vcf <- function(path.to.vcf, tissue.id) {
	require(GenomicRanges)
	require(data.table)
	if (length(grep("gz", path.to.vcf)) == 0) {
		variants <- fread(path.to.vcf, stringsAsFactors = F)
	} else {
		variants <- fread(paste0("gunzip -c ", path.to.vcf), stringsAsFactors = F)
	}
	# head(variants, 12)
	## Check if there is a SNP id column
	if (sum(length(grep("A", variants[,3][[1]])), length(grep("C", variants[,3][[1]])), length(grep("G", variants[,3][[1]])), length(grep("T", variants[,3][[1]]))) != 0) {
		colnames(variants)[1:4] <- c("chr", "start", "REF", "ALT")
	} else {
		colnames(variants)[1:5] <- c("chr", "start", "rsID", "REF", "ALT")
	}
	
	
	## Save the original copy
	variants$varID <- paste(variants$chr, variants$start, sep = ":")
	saveRDS(variants, "./cache/original.variants.rds")
	
	## Remove extra columns
	# variants[, grep("^V[0-9]+", colnames(variants)):=NULL]
	variants[, c("rsID"):=NULL]
	variants[,5:ncol(variants):=NULL]
	variants$varID <- paste(variants$chr, variants$start, sep = ":")
	
	## Convert to UCSC chromosome names, if not already
	if (length(grep("chr", variants$chr)) == 0) {
		variants[, chr:=paste0("chr", chr)]
	}
	
	## Check if it's a multi-allele vcf. If so, split into single alleles
	if (length(grep(",", variants$ALT)) != 0) {
		print("The vcf file is multi-allelic. Splitting into single alleles.")
		split.ALT <- strsplit(variants$ALT, ",")
		# split.af <- strsplit(variants$AF, ",")
		n.alleles <- unlist(lapply(split.ALT, length))
		variants <- variants[rep(1:nrow(variants), n.alleles)]
		variants$ALT <- unlist(split.ALT)
		# variants$AF <- unlist(split.af)
	}
	
	## Add end coordinate
	variants[,end:=start]
	head(variants, 12)
	
	## Convert to GRanges object
	variants <- GRanges(variants)
	gc()
	
	## Order variant object by genomic coordinates
	variants <- variants[order(variants)]
	
	## Remove indels (for now)
	if (any(nchar(variants$REF) > 1 | nchar(variants$ALT) > 1)) {
		print("Current implementation is unable to predict damage impact of indels. They will be removed from further analysis.")
		variants <- variants[nchar(variants$REF) == 1 & nchar(variants$ALT) == 1]
	}
	
	## Annotate variants with genomic features, poly-tissue GEP and ChromHMM sets, TF motif scores
	if (!file.exists(paste0(path.to.vcf, ".feature.motif.tmp.rds")) | fresh.run == T) {
		source("./lib/feature.motif.annotation.R")
		variants <- feature.motif.annotation(variants)
		## Save intermediate file
		saveRDS(variants, paste0(path.to.vcf, ".feature.motif.tmp.rds"))
	} else {
		variants <- readRDS(paste0(path.to.vcf, ".feature.motif.tmp.rds"))
	}
	
	## Annotate variants with conservation scores
	if (!file.exists(paste0(path.to.vcf, ".conservation.tmp.rds")) | fresh.run == T) {
		source("./lib/add.conservation.scores.R")
		variants <- annotate.GERP(variants)
		variants <- annotate.phastCons(variants)
		variants <- annotate.phyloP(variants)
		## Save intermediate file
		saveRDS(variants, paste0(path.to.vcf, ".conservation.tmp.rds"))
	} else {
		variants <- readRDS(paste0(path.to.vcf, ".conservation.tmp.rds"))
	}
	
	## Annotate variants with tissue-specific epigenetic marks and ChromHMM and GEP predicted states
	source("./lib/epigenome.annotation.R")
	variants <- epigenome.annotation(variants, tissue.id)
	
	## Save the variant file
	saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))
}

	###################################
	####### Predicting variants #######
	###################################
predict.vcf <- function(path.to.vcf, tissue.id) {
	variants <- readRDS(paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))
	source("./src/prep.sets.R")
	# sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", test.set.tissues = c("fMuscle", "HSMM"))
	sets.list <- prep.sets(path.to.full.set = "./cache/all.variants.partial.annotation.rds", path.to.test.sets = "./cache/", test.set.tissues = "HSMM")
	names(sets.list)
	
	train.set <- sets.list$training
	valid.set <- sets.list$validation
	
	## Train a random forest model
	require(pROC)
	require(h2o)
	require(RColorBrewer)
	## Initiate H2O instance
	localH2O <- h2o.init(nthreads=-1)
	nfolds <- 5
	clrs <- brewer.pal(5, "Dark2")
	
	train.set.h2o <- as.h2o(train.set)
	valid.set.h2o <- as.h2o(valid.set)
	target <- "regulatory"
	features <- setdiff(colnames(train.set.h2o), target)
	## Remove qval features, as they are redundant with pval
	features <- features[-grep("qval", features)]
	## remove auxiliary columns from features
	features <- features[!features %in% c("varID", "source")]
	
	####
	# features <- features[features %in% colnames(mcols(variants))]
	####
	
	##########################
	### Random fores model ###
	##########################
	## ntrees = 71
	nt = 40
	## max_depth = 10
	md = 8
	
	## With cross-validation
	model_drf <- h2o.randomForest(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_drf", nfolds = nfolds, ntrees = nt, max_depth = md, fold_assignment = "Modulo", keep_cross_validation_predictions = T, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T)
	
	
	# h2o.auc(model_drf, xval = T)
	# h2o.auc(model_drf, valid = T)
	# h2o.auc(model_drf, train = T)
	
	#################
	### GBM model ###
	#################
	## ntrees = 43
	nt = 13
	## max_depth = 7
	md = 4
	## learn_rate = 0.1
	lr = 0.16
	# variant sample rate
	sr <- 0.9
	# feature sample rate
	cr <- 0.9
	
	model_gbm <- h2o.gbm(x = features, y = target, training_frame = train.set.h2o, model_id = "h2o_gbm", ntrees = nt, max_depth = md, learn_rate = lr, seed = 16052017, validation_frame = valid.set.h2o, score_each_iteration = T, sample_rate = sr, col_sample_rate = cr, nfolds = nfolds, fold_assignment = "Modulo", keep_cross_validation_predictions = TRUE)
	
	# h2o.auc(model_gbm, train = T)
	# h2o.auc(model_gbm, xval = T)
	# h2o.auc(model_gbm, valid = T)
	
	
	model_ensemble <- h2o.stackedEnsemble(x = features, y = target,
										  training_frame = train.set.h2o,
										  validation_frame = valid.set.h2o, 
										  model_id = "h2o_ensemble", 
										  base_models = list(model_gbm@model_id, model_drf@model_id))
	
	predicted <- h2o.predict(model_ensemble, as.h2o(as.data.frame(variants)))
	predicted <- as.data.frame(predicted)
	dim(predicted)
	variants$p.regulatory <- as.character(round(predicted$p1, 4))
	## The score does not try to predict the functional impact of coding variants. To keep all the original variants, but avoid possible mis-intepreptation, assign NA score to coding variants
	variants$p.regulatory[variants$feature == "coding"] <- "NA_coding"
	hist(as.numeric(variants$p.regulatory), breaks = 20)
	summary(as.numeric(variants$p.regulatory) > 0.6)
	variants[which(as.numeric(variants$p.regulatory) > 0.6)]
	## The multi-allele variants have been expanded, compress them back to match original VCF file
	compressed.p.reg <- tapply(variants$p.regulatory, variants$varID, function(x) paste(x, collapse = ","))
	
	## Load original variant set
	original.variants <- readRDS("./cache/original.variants.rds")
	
	original.variants$p.regulatory <- "NA_indel"
	original.variants$p.regulatory[match(names(compressed.p.reg), original.variants$varID)] <- compressed.p.reg
	original.variants[1:20,]
	original.variants[, varID:=NULL]
	colnames(original.variants)[1] <- "#CHR"
	colnames(original.variants)[2] <- "POS"
	## Write out original vcf adding a column with damage annotation
	out.path.to.vcf <- gsub(".vcf.*", ".annotated.vcf", path.to.vcf)
	write.table(original.variants, out.path.to.vcf, sep = "\t", quote = F, row.names = F)
	
	## Additionally, write out a csv file containing all annotated variants, note that indels and coding variants are missing from this table
	# variants.dt <- as.data.frame(variants) %>%
	# 	select(c("seqnames", "start", "REF", "ALT", "p.regulatory")) %>%
	# 	filter(p.regulatory != "NA_coding") %>%
	# 	mutate(p.regulatory=as.numeric(p.regulatory))
	
	write.csv(variants[order(variants$p.regulatory, decreasing = T),], gsub(".vcf.*", paste0(".", tissue.id, ".annotated.csv"), path.to.vcf), quote = F, row.names = F)
	
	## Save the variant file
	saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))
	# return(variants)
}
