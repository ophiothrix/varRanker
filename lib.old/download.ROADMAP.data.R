#### Function to download ROADMAP data for a specific tissue. Takes ROADMAP tissue id ####
download.ROADMAP.data <- function(tissue.id) {
	require(GenomicRanges)
	require(data.table)
	require(curl)
	require(RCurl)
	### Download all available ChromHMM calls ###
	## Check if the pan-tissue matrix of ChromHMM calls exists, if not, download calls for all the tissues and build the matrix
	dir.create("./data/ENCODE/chromHMM.calls", recursive = T, showWarnings = F)
	
	if (length(list.files("./data/ENCODE/chromHMM.calls/", "mnemonics.bed.gz")) >= 98) {
		print("ChromHMM call files exist. Skipping download.")
	} else {
		print("Downloading ChromHMM calls for all available tissues")
		dir.create("./data/ENCODE/chromHMM.calls", recursive = T)
		chromhmm.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz"
		curl_download(url = chromhmm.url, destfile = paste0("./data/ENCODE/", basename(chromhmm.url)))
		## Extract the files
		system("tar -xvf ./data/ENCODE/all.mnemonics.bedFiles.tgz -C ./data/ENCODE/chromHMM.calls/")
	}
	### Build a pan-tissue matrix for each state
	### Take all available chromHMM called tissues from ROADMAP and for each state generate:
	# - union call - combine all tissues
	# - frequency matrix - for each region count in how many tissues it is present
	
	## Check if the files already exist and only generate them if they don't
	if (length(list.files("./data/ENCODE/chromHMM.calls", "union.rds")) == 18) {
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
		saveRDS(all.tissues.GR, "./data/ENCODE/all.tissues.rds")
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
			saveRDS(reduced.state.GR, file = paste0("./data/ENCODE/chromHMM.calls/", state, ".union.rds"))
		}
	}

	### Download ROADMAP peak calls data ###
	## Set a list of core histone marks
	core.marks <- c("DNase","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9ac","H3K9me3")
	
	## Download relevant ROADMAP files for a given tissue id
	narrow.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"
	narrow.files <- getURLContent(narrow.url)
	narrow.files <- strsplit(narrow.files, "<td>")[[1]]
	narrow.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, narrow.files, value = T))
	for (narrow.file in narrow.files) {
		if (file.exists(paste0("./data/ENCODE/", narrow.file))) {
			print(paste0("./data/ENCODE/", narrow.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", narrow.file, "..."))
			curl_download(url = paste0(narrow.url, narrow.file), destfile = paste0("./data/ENCODE/", narrow.file))
		}
	}
	
	broad.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/broadPeak/"
	broad.files <- getURLContent(broad.url)
	broad.files <- strsplit(broad.files, "<td>")[[1]]
	broad.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, broad.files, value = T))
	for (broad.file in broad.files) {
		if (file.exists(paste0("./data/ENCODE/", broad.file))) {
			print(paste0("./data/ENCODE/", broad.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", broad.file, "..."))
			curl_download(url = paste0(broad.url, broad.file), destfile = paste0("./data/ENCODE/", broad.file))
		}
	}
	
	gapped.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/gappedPeak/"
	gapped.files <- getURLContent(gapped.url)
	gapped.files <- strsplit(gapped.files, "<td>")[[1]]
	gapped.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, gapped.files, value = T))
	for (gapped.file in gapped.files) {
		if (file.exists(paste0("./data/ENCODE/", gapped.file))) {
			print(paste0("./data/ENCODE/", gapped.file, " exists. Skipping download."))
		} else {
			print(paste0("Downloading ", gapped.file, "..."))
			curl_download(url = paste0(gapped.url, gapped.file), destfile = paste0("./data/ENCODE/", gapped.file))
		}
	}
	
	## Check that all core marks are available.
	missing.marks <- setdiff(core.marks, unique(gsub("\\..*", "", gsub(".*-", "", list.files("./data/ENCODE/", tissue.id)))))
	if (length(missing.marks) != 0) {
		imputed.url <- "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidatedImputed/narrowPeak/"
		imputed.files <- getURLContent(imputed.url)
		imputed.files <- strsplit(imputed.files, "<td>")[[1]]
		imputed.files <- gsub(".*>(E.*)<\\/a.*", "\\1", grep(tissue.id, imputed.files, value = T))
		imputed.files <- grep(missing.marks, imputed.files, value = T)
		for (imputed.file in imputed.files) {
			if (file.exists(paste0("./data/ENCODE/", imputed.file))) {
				print(paste0("./data/ENCODE/", imputed.file, " exists. Skipping download."))
			} else {
				print(paste0("Downloading ", imputed.file, "..."))
				curl_download(url = paste0(imputed.url, imputed.file), destfile = paste0("./data/ENCODE/", imputed.file))
			}
		}
		
	}
	
	# Puff, imputed files only contain coordinates. So we can either forego the signal data from those files. Or stick to only using tissues with all marks available. BTW, we could use many more marks for the tissues with imputed data! Technically, we only need the signal data from H3K4me1 and me3 to take the ratio...
	
	### Download IDEAS chromatin states
	dir.create("./data/ENCODE/IDEAS", showWarnings = F, recursive = T)
	
	## If doesn't already exist, generate tissue id to url mappings
	if (file.exists(paste0("./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed.gz"))) {
		print("IDEAS chromatin state file exists. Skipping download...")
	} else {
		print("Downloading and converting IDEAS chromatin states file")
		if (!file.exists("./data/ENCODE/IDEAS/ideas.table.rds")) {	
		meta.file <- getURLContent("http://bx.psu.edu/~yuzhang/Roadmap-25state/trackDb.txt")
		meta.file <- unlist(strsplit(meta.file, "\n"))
		urls <- gsub("bigDataUrl ", "", grep("bigDataUrl", meta.file, value = T))
		
		tissues <- gsub("shortLabel ", "", grep("shortLabel", meta.file, value = T))
		tissue.ids <- gsub("^(E[0-9]+) .*", "\\1", tissues)
		
		ideas.table <- data.frame(tissue.id = tissue.ids, url = urls, stringsAsFactors = F)
		saveRDS(ideas.table, "./data/ENCODE/IDEAS/ideas.table.rds")
		}
		### Download and convert IDEAS states for a given tissue id
		ideas.table <- readRDS("./data/ENCODE/IDEAS/ideas.table.rds")
		curl_download(url = ideas.table$url[ideas.table$tissue.id == tissue.id], destfile = paste0("./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb"), quiet = F)
		
		## Convert bigBed to bed file using UCSC tool compiled for the appropriate tissue
		ucsc.path <- "~/utils/UCSC.tools/"
		system(paste0(ucsc.path, "/bigBedToBed ./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb ./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed"))
		system(paste0("gzip ./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bed"))
		file.remove(paste0("./data/ENCODE/IDEAS/", tissue.id, ".ideas.states.bb"))
	} 
}

