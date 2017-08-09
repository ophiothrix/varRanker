### This file specifies base functions aimed at annotating the variant list with a selection of features ###

#### Please specify correct path to FIMO installation here! ####
## You need to have MEME suite installed ##
## The version used for development is meme_4.12.0 ##
## The motif database version used is motif_databases.12.15 ##
fimo.path <- "~/utils/meme/bin/fimo"


#### Extract Position Frequency Matrices from a given MEME motif database ####

### The function takes the full path to the meme file and returns and saves a list with motif PFMs as elements, as well as a dataframe with motif information. 

extract.pfm <- function (database.path, force.run = F) {
	# Check if the file already exists
	if (file.exists(gsub("meme", "PFM.list.rds", paste0("./cache/", basename(database.path)))) & force.run == F) {
		print("Loading the existing motif PFM file")
		pfm.list <- readRDS(gsub("meme", "PFM.list.rds", paste0("./cache/", basename(database.path))))
		return(pfm.list)
	} else {
		print("Making a new motif PFM file")
		# Extract PFMs from the motif database -------------------------------------
		meme.file <- readLines(database.path)
		# Extract motif IDs and gene symbols
		motifs <- as.data.frame(matrix(unlist(strsplit(grep("MOTIF", meme.file, value = T), " ")), ncol = 3, byrow = T)[,2:3], stringsAsFactors = F)
		colnames(motifs) <- c("MotifID", "Symbol")
		# Get the motif lengths
		motifs$length <- matrix(unlist(strsplit(grep("alength", meme.file, value=4), " ")), ncol = 10, byrow = T)[,6]
		
		# Extract PFM for each motif. Get the line of each motif and the line of the subsequent motif and get the lines between them minus padding. Need some extra care with the last motif. Essentially, add a "pretend" motif at the end of the object). Also need to add an extra line to the object, or it cuts the last PFM short.
		pfms <- c(grep("MOTIF", meme.file), length(meme.file)+1)
		pfm.list <- list()
		
		## Unfortunately the motif DB files are not uniform, e.g. in gaps between motifs. Need to hard code it in.
		if (length(grep("URL", meme.file)) == 0) {
			gap <- 2
		} else {
			gap <- 4
		}
		for (i in 1:(length(pfms)-1)){
			pfm.list[[motifs$MotifID[i]]] <- matrix(as.numeric(unlist(strsplit(meme.file[(pfms[i]+3):(pfms[i+1]-gap)], "\t"))), ncol = 4, byrow = T)
		}
		# Add column names
		for (i in 1:length(pfm.list)) {
			colnames(pfm.list[[i]]) <- c("A", "C", "G", "T")
		}
		# Check that the motif length corresponds to the PFM matrices
		if (!all(unlist(lapply(pfm.list, nrow)) == motifs$length)) {
			stop("Motif length is inconsistent between the declared motif annotation and extracted PFM") 
		}
		
		saveRDS(motifs, file=gsub("meme", "motif.table.rds", paste0("./cache/", basename(database.path))))
		saveRDS(pfm.list, file=gsub("meme", "PFM.list.rds", paste0("./cache/", basename(database.path))))
		return(pfm.list)
	}
}

#### Generate motif damage score ####
### Function to generate motif damage score for each variant
## Takes the list of variants as a GRanges object and path to motif database as input.
## For each variant the function extracts the flanking genomic sequence based on the length of the longest motif
## Each sequence is then run against the selected motif database using FIMO
## Once motifs are mapped to the sequences, we map each variant to the motifs and take frequency ratio of Reference and Alternate alleles in the motif's PWM. The ratio is log2 transformed and acts as a damage score
## Additionally, we extract a binary feature reflecting whether a variant has been mapped to any motif.
get.damage.scores.direct <- function(database.path, variants) {
	## Specify the "default" path to fimo binary - need to think how best do it. Perhaps in the annotation script
	require(BSgenome.Hsapiens.UCSC.hg19)
	require(GenomicRanges)
	require(GenomicFeatures)
	require(OrganismDbi)
	compl.table <- c(A = "T", C = "G", G = "C", T = "A")
	variants$log2.damage.score <- 0
	variants$motif.hit <- 0
	
	### Get PFMs from the database file
	motif.PFMs <- extract.pfm(database.path, force.run = F)
	
	### Run FIMO on the target variants
	variants.tmp <- variants
	
	mcols(variants.tmp) <- as.data.frame(mcols(variants.tmp))[colnames(as.data.frame(mcols(variants.tmp))) %in% c("REF", "ALT")]
	head(variants.tmp)
	
	## Extend the region by flank bp either side of the variant, where flank should be the length of the longest motif + 1
	flnk <- max(unlist(lapply(motif.PFMs, nrow))) + 1
	start(variants.tmp) <- start(variants.tmp)-flnk
	end(variants.tmp) <- end(variants.tmp)+flnk
	
	## Get the reference sequence
	ref.seq <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, variants.tmp))
	## Check if the reference allele is the same in reference sequence and in the vcf
	ref.allele <- substr(ref.seq, flnk+1, nchar(ref.seq)-flnk)
	all(ref.allele == variants.tmp$REF)
	## To each sequence add corresponding genomic coordinates
	coords <- paste0(">", paste(as.character(seqnames(variants.tmp)), paste(start(variants.tmp), end(variants.tmp), sep = "-"), sep = ":"))
	ref.seq <- cbind(coords, ref.seq)
	
	## Eventually, it would be nice (especially for indels) to reconstruct the reference and alternate seqnences, run FIMO on those separately and compare the scores. But for now (SNP solution only) we can take the reference genome sequence. The difference will be calculated between the reference and alternate nuceotides
	
	## write sequence files
	write.table(ref.seq, "./ref.seq.file.txt", quote = F, col.names = F, row.names = F, sep = "\n")
	
	## Run FIMO on the saved sequences
	system(paste0(fimo.path, " --text --skip-matched-sequence --parse-genomic-coord ", database.path, " ref.seq.file.txt > motif.map.tmp.txt"))
	file.remove("ref.seq.file.txt")
	
	## Load the map of motifs
	mapped.motifs.chr <- GRanges(read.table("motif.map.tmp.txt", col.names = c("Motif", "Alt_motif_id", "Chr", "Start", "End", "Strand", "Score", "Pvalue"), stringsAsFactors = F))
	names(mapped.motifs.chr) <- mapped.motifs.chr$Motif
	file.remove("motif.map.tmp.txt")
	
	# Check that all motif IDs in the fimo_out file are present in meme database (i.e. don't care if there are some motifs in the database that didn't get called, as long as all that did get called are in the PFM database)
	all(unique(names(mapped.motifs.chr)) %in% names(motif.PFMs))
	
	## Map the variants to motif hits
	mapped.vars <- as.data.frame(mapToTranscripts(variants, mapped.motifs.chr))
	mapped.vars$REF <- variants$REF[mapped.vars$xHits]
	mapped.vars$ALT <- variants$ALT[mapped.vars$xHits]
	
	# Change the variants overlapping with reverse-strand motifs to the complementary
	if (any(mapped.vars$strand == "-")) {
		mapped.vars[mapped.vars$strand == "-", c("REF", "ALT")] <- names(compl.table)[match(unlist(mapped.vars[mapped.vars$strand == "-", c("REF", "ALT")]), compl.table)]
	}
	
	# Reset the variants as factors, so we can easily deal with them as numbers later
	mapped.vars$REF <- factor(mapped.vars$REF, levels=c("A", "C", "G", "T"))
	mapped.vars$ALT <- factor(mapped.vars$ALT, levels=c("A", "C", "G", "T"))
	mapped.vars$damage.score <- NA
	
	# For a given motif, extract frequencies for the reference allele and alternative allele and calculate the REF/ALT ratio
	for (i in as.character(unique(mapped.vars$seqnames))) {
		pfm <- motif.PFMs[[i]]
		mapped.vars$damage.score[ mapped.vars$seqnames == i] <- 
			(pfm[cbind(mapped.vars$start[ mapped.vars$seqnames == i], as.numeric(mapped.vars$REF[ mapped.vars$seqnames == i]))]+0.000001)/
			(pfm[cbind(mapped.vars$start[ mapped.vars$seqnames == i], as.numeric(mapped.vars$ALT[ mapped.vars$seqnames == i]))]+0.000001)
	}
	# The ratio can capture both the variants that potentially disrupt the binding (higher freq base -> lower freq base) and enhance the binding (lower freq base -> higher freq base). Take a log2 to make the damage score symmetrical
	mapped.vars$log2.damage.score <- log2(mapped.vars$damage.score)
	
	# Only some of the variants are mapped to the motifs. Additionally, each variant can map to more than one motif. So in order to assign a damage score to each of the original variables two things need to be done: 1) Get a single damage score per motif based on the highest absolute value. 2) Assign the motifs with damage scores back to the original list of variants. The variants without damage score get 0.
	max.log2.damage.score <- tapply(mapped.vars$log2.damage.score, mapped.vars$xHits, function(x) x[which.max(abs(x))])
	variants$log2.damage.score[as.numeric(names(max.log2.damage.score))] <- max.log2.damage.score
	motif.hit <- tapply(mapped.vars$log2.damage.score, mapped.vars$xHits, function(x) x[which.max(abs(x))])
	variants$motif.hit[unique(mapped.vars$xHits)] <- 1
	summary(variants$motif.hit == 1 & variants$log2.damage.score == 0)
	## Also make a separate feature to mark whether a variant hits a motif without specifying the damage
	gc()
	return(variants)
}
