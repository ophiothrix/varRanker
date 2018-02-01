#### Function to download ROADMAP data for a specific tissue. Takes ROADMAP tissue id ####
download.ROADMAP.data <- function(tissue.id) {
	require(GenomicRanges)
	require(data.table)
	require(curl)
	require(RCurl)
	### Download all available ChromHMM calls ###
	## Check if the pan-tissue matrix of ChromHMM calls exists, if not, download calls for all the tissues and build the matrix
	dir.create("./cache/ENCODE/chromHMM.calls", recursive = T, showWarnings = F)
	
	if (length(list.files("./cache/ENCODE/chromHMM.calls/", "mnemonics.bed.gz")) >= 98) {
		print("ChromHMM call files exist. Skipping download.")
	} else {
		print("Downloading ChromHMM calls for all available tissues")
		dir.create("./cache/ENCODE/chromHMM.calls", recursive = T)
		chromhmm.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz"
		curl_download(url = chromhmm.url, destfile = paste0("./cache/ENCODE/", basename(chromhmm.url)))
		## Extract the files
		system("tar -xvf ./cache/ENCODE/all.mnemonics.bedFiles.tgz -C ./cache/ENCODE/chromHMM.calls/")
	}
	### Build a pan-tissue matrix for each state
	### Take all available chromHMM called tissues from ROADMAP and for each state generate:
	# - union call - combine all tissues
	# - frequency matrix - for each region count in how many tissues it is present
	
	## Check if the files already exist and only generate them if they don't
	if (length(list.files("./cache/ENCODE/chromHMM.calls", "union.rds")) == 18) {
		print("Pan-tissie chromatin state files already exist. Moving on...")
	} else {
		print("Generating pan-tissue matrix of chromatin states...")
		
		## Combine all chromatin state files into a single file annotating each region with the tissue it came from
		
		system("bash ./lib/combine.chrom.states.sh")
		
		## Load the union file
		all.tissues <- fread("all.tissues.bed")
		head(all.tissues)
		colnames(all.tissues) <- c("chr", "start", "end", "state", "tissue.name")
		all.tissues.GR <- GRanges(all.tissues)
		rm(all.tissues)
		gc()
		saveRDS(all.tissues.GR, "./cache/ENCODE/all.tissues.rds")
		file.remove("all.tissues.bed")
		all.tissues.GR[all.tissues.GR$state == "12_ZNF/Rpts"]$state <- "12_ZNF_Rpts"
		gc()
		
		## Split into state-specific GRanges
		for (state in unique(all.tissues.GR$state)) {
			print(state)
			# Select all the regions for the given state and merge overlapping regions
			full.state.GR <- all.tissues.GR[all.tissues.GR$state ==  state]
			reduced.state.GR <- reduce(all.tissues.GR[all.tissues.GR$state ==  state])
			## Look in how many tissues each reduced region overlaps with the same state
			olaps <- as.data.frame(findOverlaps(reduced.state.GR, full.state.GR))
			olaps$tissue.name <- full.state.GR[olaps$subjectHits]$tissue.name
			reduced.state.GR$n.tissues <- tapply(olaps$tissue.name, olaps$queryHits, function(x) length(unique(x)))
			# Save the regions as a file
			saveRDS(reduced.state.GR, file = paste0("./cache/ENCODE/chromHMM.calls/", state, ".union.rds"))
		}
	}

	### Download ROADMAP peak calls data ###
	## Set a list of core histone marks
	core.marks <- c("DNase.macs2","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3")
	extra.marks <- c("H3K27ac","H3K9ac")
	
	## Download relevant ROADMAP files for a given tissue id
	narrow.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"
	narrow.files <- getURLContent(narrow.url)
	narrow.files <- strsplit(narrow.files, "<td>")[[1]]
	narrow.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, narrow.files, value = T))
	
	
	## get files for core marks
	narrow.files.core <- unlist(sapply(core.marks, function(x) grep(x, narrow.files)))
	## Check that all core marks are available for a given tissue
	if (length(narrow.files.core) != 6) {
		stop("One or more core marks are not available for this tissue ID. Please select a tissue ID that has at least the following marks available: DNase.macs2, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3.")
	}
	narrow.files.extra <- unlist(sapply(extra.marks, function(x) grep(x, narrow.files)))
	
	## Combine available marks
	narrow.files <- narrow.files[c(narrow.files.core, narrow.files.extra)]
	
	## Download files
	for (narrow.file in narrow.files) {
		if (file.exists(paste0("./cache/ENCODE/", narrow.file))) {
			print(paste0("./cache/ENCODE/", narrow.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", narrow.file, "..."))
			curl_download(url = paste0(narrow.url, narrow.file), destfile = paste0("./cache/ENCODE/", narrow.file))
		}
	}
	
	### Broad peaks ###
	broad.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/"
	broad.files <- getURLContent(broad.url)
	broad.files <- strsplit(broad.files, "<td>")[[1]]
	broad.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, broad.files, value = T))
	## get files for core marks
	broad.files.core <- unlist(sapply(core.marks, function(x) grep(x, broad.files)))
	## Check that all core marks are available for a given tissue
	if (length(broad.files.core) != 5) {
		stop("One or more core marks are not available for this tissue ID. Please select a tissue ID that has at least the following marks available: DNase.macs2, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3.")
	}
	broad.files.extra <- unlist(sapply(extra.marks, function(x) grep(x, broad.files)))
	## Combine available marks
	broad.files <- broad.files[c(broad.files.core, broad.files.extra)]
	
	for (broad.file in broad.files) {
		if (file.exists(paste0("./cache/ENCODE/", broad.file))) {
			print(paste0("./cache/ENCODE/", broad.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", broad.file, "..."))
			curl_download(url = paste0(broad.url, broad.file), destfile = paste0("./cache/ENCODE/", broad.file))
		}
	}
	
	gapped.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/gappedPeak/"
	gapped.files <- getURLContent(gapped.url)
	gapped.files <- strsplit(gapped.files, "<td>")[[1]]
	gapped.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, gapped.files, value = T))
	## get files for core marks
	gapped.files.core <- unlist(sapply(core.marks, function(x) grep(x, gapped.files)))
	## Check that all core marks are available for a given tissue
	if (length(gapped.files.core) != 5) {
		stop("One or more core marks are not available for this tissue ID. Please select a tissue ID that has at least the following marks available: DNase.macs2, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3.")
	}
	gapped.files.extra <- unlist(sapply(extra.marks, function(x) grep(x, gapped.files)))
	## Combine available marks
	gapped.files <- gapped.files[c(gapped.files.core, gapped.files.extra)]
	
	for (gapped.file in gapped.files) {
		if (file.exists(paste0("./cache/ENCODE/", gapped.file))) {
			print(paste0("./cache/ENCODE/", gapped.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", gapped.file, "..."))
			curl_download(url = paste0(gapped.url, gapped.file), destfile = paste0("./cache/ENCODE/", gapped.file))
		}
	}
	
	### Download imputed peaks for extra marks
	
	### Narrow imputed peaks
	imputed.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/"
	imputed.files <- getURLContent(imputed.url)
	imputed.files <- strsplit(imputed.files, "<td>")[[1]]
	imputed.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, imputed.files, value = T))
	imputed.files <- imputed.files[grep(".gz$", imputed.files)]
	imputed.files.extra <- sapply(extra.marks, function(x) grep(x, imputed.files))
	imputed.files <- imputed.files[imputed.files.extra]
	for (imputed.file in imputed.files) {
		if (file.exists(paste0("./cache/ENCODE/", imputed.file))) {
			print(paste0("./cache/ENCODE/", imputed.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", imputed.file, "..."))
			curl_download(url = paste0(imputed.url, imputed.file), destfile = paste0("./cache/ENCODE/", imputed.file))
		}
	}
	
	### gapped imputed peaks
	imputed.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/gappedPeak/"
	imputed.files <- getURLContent(imputed.url)
	imputed.files <- strsplit(imputed.files, "<td>")[[1]]
	imputed.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, imputed.files, value = T))
	imputed.files <- imputed.files[grep(".gz$", imputed.files)]
	imputed.files.extra <- sapply(extra.marks, function(x) grep(x, imputed.files))
	imputed.files <- imputed.files[imputed.files.extra]
	for (imputed.file in imputed.files) {
		if (file.exists(paste0("./cache/ENCODE/", imputed.file))) {
			print(paste0("./cache/ENCODE/", imputed.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", imputed.file, "..."))
			curl_download(url = paste0(imputed.url, imputed.file), destfile = paste0("./cache/ENCODE/", imputed.file))
		}
	}
	
	# Puff, imputed files only contain coordinates. So we can either forego the signal data from those files. Or stick to only using tissues with all marks available. BTW, we could use many more marks for the tissues with imputed data! Technically, we only need the signal data from H3K4me1 and me3 to take the ratio...
	
	### Download IDEAS chromatin states
	dir.create("./cache/ENCODE/IDEAS", showWarnings = F, recursive = T)
	
	## If doesn't already exist, generate tissue id to url mappings
	if (file.exists(paste0("./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed.gz"))) {
		print("IDEAS chromatin state file exists. Skipping download...")
	} else {
		print("Downloading and converting IDEAS chromatin states file")
		if (!file.exists("./cache/ENCODE/IDEAS/ideas.table.rds")) {	
		meta.file <- getURLContent("http://bx.psu.edu/~yuzhang/Roadmap_ideas/trackDb_test.txt")
		meta.file <- unlist(strsplit(meta.file, "\n"))
		urls <- gsub("bigDataUrl ", "", grep("bigDataUrl", meta.file, value = T))
		
		tissues <- gsub("shortLabel ", "", grep("shortLabel", meta.file, value = T))
		tissue.ids <- gsub("^(E[0-9]+)_.*", "\\1", tissues)
		
		ideas.table <- data.frame(tissue.id = tissue.ids, url = urls, stringsAsFactors = F)
		saveRDS(ideas.table, "./cache/ENCODE/IDEAS/ideas.table.rds")
		}
		### Download and convert IDEAS states for a given tissue id
		ideas.table <- readRDS("./cache/ENCODE/IDEAS/ideas.table.rds")
		curl_download(url = ideas.table$url[ideas.table$tissue.id == tissue.id], destfile = paste0("./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb"), quiet = F)
		
		## Convert bigBed to bed file using UCSC tool compiled for the appropriate tissue
		ucsc.path <- "~/tools/UCSC.tools/"
		system(paste0(ucsc.path, "/bigBedToBed ./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb ./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed"))
		system(paste0("gzip ./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed"))
		file.remove(paste0("./cache/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb"))
	} 
}

