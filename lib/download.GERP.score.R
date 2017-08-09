## Download GERP score from UCSC browser and store them as a GRanges object
rm(list = ls())
gc()
dir.create("./data/conservation.scores/GERP", showWarnings = F, recursive = T)
require(curl)
require(GenomicRanges)
require(data.table)


## Download the full GERP file
if (!file.exists("./data/conservation.scores/GERP/hg19.GERP_scores.tar.gz")) {
	curl_download(url = "http://mendel.stanford.edu/SidowLab/downloads/gerp/hg19.GERP_scores.tar.gz", destfile = "./data/conservation.scores/GERP/hg19.GERP_scores.tar.gz", quiet = F)
}

system("tar -xvf ./data/conservation.scores/GERP/hg19.GERP_scores.tar.gz -C ./data/conservation.scores/GERP/")

## Unpack the archive into individual files and gzip them
for (fname in list.files("./data/conservation.scores/GERP", "maf.rates$", full.name = T)) {
	print(fname)
	system(paste0("gzip ", fname))
}

file.remove("./data/conservation.scores/GERP/hg19.GERP_scores.tar.gz")
