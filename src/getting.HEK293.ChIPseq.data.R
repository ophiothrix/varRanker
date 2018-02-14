library(data.table)
library(tidyverse)
library(curl)
tracks <- fread("gunzip -c ENCFF135NOY.bed.gz", stringsAsFactors = F)
head(tracks)
table(tracks$V6)
table(tracks$V7)
table(tracks$V8)


epitable <- read.table("metadata (6).tsv", stringsAsFactors = F, sep = "\t", header = T)
head(epitable)
table(epitable$Experiment.target)
download.table <- epitable %>%
    filter(Library.made.from == "DNA" &
               Experiment.target != "Control-human" &
               File.format %in% c("bed broadPeak", "bed narrowPeak") &
               Assay == "ChIP-seq" &
               Output.type == "replicated peaks" &
               Assembly == "GRCh38") %>%
    mutate(Experiment.target = gsub("-human", "", Experiment.target)) %>%
    select(c(File.format, Output.type, Experiment.target, Assembly, File.download.URL))

download.table
for (i in 1:nrow(download.table)) {
    download.file(url = download.table$File.download.URL[i], destfile = paste0("E199-", download.table$Experiment.target[i], ".GRCh38.narrowPeak.gz"))
}

for (fname in list.files(".", "GRCh38.narrowPeak.gz")) {
    print(fname)
    outname <- gsub(".GRCh38", "", fname)
    hg38 <- read.table(fname)
    write.table(hg38[,1:4], "tmp.hg38.bed", sep = "\t", quote = F, row.names = F, col.names = F)
    system("~/tools/liftOver/liftOver tmp.hg38.bed ~/tools/liftOver/hg38ToHg19.over.chain.gz tmp.hg19.bed unmapped.bed")
    hg19 <- read.table("tmp.hg19.bed")
    hg19 <- cbind(hg19, hg38[match(hg19$V4, hg38$V4),5:ncol(hg38)])
    gz1 <- gzfile(outname, "w")
    write.table(hg19, gz1, sep = "\t", quote = F, row.names = F, col.names = F)
    close(gz1)
    system("rm tmp.hg19.bed tmp.hg38.bed unmapped.bed")
}
