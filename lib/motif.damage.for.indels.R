saveRDS(variants, "./cache/variants.from.vcf.w.indels.rds")
variants <- readRDS("./cache/variants.from.vcf.w.indels.rds")
head(variants, 12)
## Subset to indels only
variants <- variants[nchar(variants$REF) != 1 | nchar(variants$ALT) != 1]

### Adjust variant coordinates ###

## For deletions (REF is longer than ALT) the affected locus is equal to the number of bases being deleted

base.dif <- nchar(variants$REF) - nchar(variants$ALT)
base.dif[base.dif < 0] <- 0
head(base.dif)
variants[base.dif == 0]
width(variants) <- base.dif + 1

variants$check.REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, variants))


dels <- GRanges("chr8", IRanges(rep(2000000, 3), width = 4))
dels$REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, dels))
width(dels) <- 1
dels$ALT <- c("C", "CA", "CAT")

dels <- variants

base.dif <- nchar(dels$REF) - nchar(dels$ALT)

# dels <- dels[base.dif > 0]
# base.dif <- base.dif[base.dif > 0]

end.aff <- start(dels) + nchar(dels$REF) - 1
start.aff <- end.aff - base.dif + 1

flnk <- 22

new.dels <- dels
end(new.dels) <- end.aff
start(new.dels) <- start.aff
new.dels
## Extend by flnk
start(new.dels) <- start(new.dels) - flnk
end(new.dels) <- end(new.dels) + flnk
real.ref <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, new.dels))
new.dels
coords <- paste0(">", paste(as.character(seqnames(new.dels)), paste(start(new.dels), end(new.dels), sep = "-"), sep = ":"))
ref.seq <- cbind(coords, real.ref)

real.alt <- paste0(substr(ref.seq[,2], 1, flnk), substr(real.ref, nchar(real.ref) - flnk + 1, nchar(real.ref)))
alt.seq <- cbind(coords, real.alt)

write.table(ref.seq, "./ref.seq.file.txt", quote = F, col.names = F, row.names = F, sep = "\n")
write.table(alt.seq, "./alt.seq.file.txt", quote = F, col.names = F, row.names = F, sep = "\n")


hist(nchar(variants$REF), breaks = 100)
hist(nchar(variants$ALT), breaks = 100)
## True reference needs to be repositioned so that it:
## starts at "start.aff - flnk"
## ends at "end.aff + flnk


## Run FIMO on the saved sequences
system(paste0(fimo.path, " --text --skip-matched-sequence --parse-genomic-coord ", database.path, " alt.seq.file.txt > motif.map.tmp.txt"), ignore.stderr = T)

## Load the map of motifs
mapped.motifs.alt <- GRanges(read.table("motif.map.tmp.txt", col.names = c("Motif", "Alt_motif_id", "Chr", "Start", "End", "Strand", "Score", "Pvalue"), stringsAsFactors = F))
names(mapped.motifs.alt) <- mapped.motifs.alt$Motif

## Map the variants to motif hits
mapped.vars <- as.data.frame(mapToTranscripts(new.dels, mapped.motifs.chr))
mapped.vars$REF <- variants$REF[mapped.vars$xHits]
mapped.vars$ALT <- variants$ALT[mapped.vars$xHits]
head(mapped.vars)

mapped.vars.alt <- as.data.frame(mapToTranscripts(new.dels, mapped.motifs.alt))
mapped.vars.alt$REF <- variants$REF[mapped.vars.alt$xHits]
mapped.vars.alt$ALT <- variants$ALT[mapped.vars.alt$xHits]
head(mapped.vars.alt)

head(mapped.vars)
mapped.vars[mapped.vars$xHits == 2,]
mapped.vars.alt[mapped.vars.alt$xHits == 2,]



olaps <- as.data.frame(findOverlaps(new.dels, mapped.motifs.chr))
olaps[olaps$queryHits == 1,]
new.dels[1]
mapped.motifs.chr[olaps$subjectHits[olaps$queryHits == 1]]

olaps.alt <- as.data.frame(findOverlaps(new.dels, mapped.motifs.alt))
olaps.alt[olaps.alt$queryHits == 1,]
new.dels[1]
mapped.motifs.alt[olaps.alt$subjectHits[olaps.alt$queryHits == 1]]

olaps.ref <- findOverlaps(new.dels, mapped.motifs.chr)
full.map.ref <- mapped.motifs.chr[subjectHits(olaps.ref)]
full.map.ref$xHits <- queryHits(olaps.ref)
names(full.map.ref) <- NULL
tst.ref <- as.data.frame(full.map.ref)

tst.ref %>%
	group_by(.dots=c("xHits", "Motif")) %>%
	summarise(max_score = max(Score))

olaps.alt <- findOverlaps(new.dels, mapped.motifs.chr)
full.map.alt <- mapped.motifs.chr[subjectHits(olaps.alt)]
full.map.alt$xHits <- queryHits(olaps.alt)
names(full.map.alt) <- NULL
tst.alt <- as.data.frame(full.map.alt)

tst.alt %>%
	group_by(.dots=c("xHits", "Motif")) %>%
	summarise(max_score = max(Score))



## Get ref seq
new.dels$new.REF <- getSeq(BSgenome.Hsapiens.UCSC.hg19, new.dels)

downstream <- upstream <- new.dels
start(upstream) <- start(new.dels) - 1 - flnk
end(upstream) <- start(new.dels) - 1
upstream$new.REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, upstream))

end(downstream) <- end(new.dels) + 1 + flnk
start(downstream) <- end(new.dels) + 1
downstream$new.REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, downstream))

full.stream <- c(upstream, downstream)
full.stream <- GRanges(seqlevels(dels), IRanges(min(start(full.stream)), max(end(full.stream))))
full.stream$REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, full.stream))

paste0(upstream$new.REF, downstream$new.REF)

real.ref <- new.dels$real.ref

substr(real.ref, 1, flnk)
paste0(substr(real.ref, 1, flnk),
substr(real.ref, nchar(real.ref) - flnk + 1, nchar(real.ref)))


### GCATGTG ###
### GC---TG ###
### GCA--TG ###
### GCAT-TG ###
max(nchar(variants$REF))
max(nchar(variants$ALT))

## You have 3 bases that are to be deleted
### 1 2 3 4 5 6 7 ### 
### G C - - - T G ###

## Reference - no problem, the numbers are as they are
## Alternate - puffssssst


library(dplyr)
tst <- as.data.frame(mapped.motifs.chr)
tst.alt <- as.data.frame(mapped.motifs.alt)

tst %>%
	group_by(Motif) %>%
	summarise(max_score = max(Score))

tst.alt %>%
	group_by(Motif) %>%
	summarise(max_score = max(Score))
