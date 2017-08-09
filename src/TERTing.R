## Look at how the well-described regulatory mutations would fare
require(GenomicRanges)
terts <- GRanges("chr5", IRanges(c(1295228, 1295250), width = 1), REF = rep("G", 2), ALT = rep("A", 2)) 

source("./lib/annotate.GRanges.R")
terts.anno <- annotate.GRanges(terts, tissue.id = "E126")
training.set <- terts.anno

test.set.h2o <- as.h2o(training.set)

predicted <- as.data.frame(h2o.predict(model_drf, test.set.h2o))


## Check the motifs
alt.seq <- paste0(substr(ref.seq[,2], 1, flnk), terts$ALT, substr(ref.seq[,2], flnk+2, flnk*2+1))
alt.seq <- cbind(coords, alt.seq)
ref.seq <- rbind(ref.seq, alt.seq)
getSeq(BSgenome.Hsapiens.UCSC.hg19, terts)
