require(tidyverse)
## Convert eQTL file into a suitable format
eqtls <- read.table("./cache/allSigVars.tsv", header = T, stringsAsFactors = F, comment.char = "")
head(eqtls)

hist(table(eqtls$X.gene_id))

chrs <- gsub("^(.*)_(.*)_(.*)_(.*)_.*", "\\1", eqtls$variant_id) %>%
	paste0("chr", .)
head(chrs)
tail(chrs)

pos <- as.numeric(gsub("^(.*)_(.*)_(.*)_(.*)_.*", "\\2", eqtls$variant_id))
head(pos)
tail(pos)

refs <- gsub("^(.*)_(.*)_(.*)_(.*)_.*", "\\3", eqtls$variant_id)
alts <- gsub("^(.*)_(.*)_(.*)_(.*)_.*", "\\4", eqtls$variant_id)
head(refs)
head(alts)

eqtls.vcf <- cbind(CHR = chrs, POS = pos, REF = refs, ALT = alts, eqtls)
head(eqtls.vcf)
write.table(eqtls.vcf, "./cache/allSigVars.vcf", quote = F, sep = "\t", row.names = F)

source("./lib/annotate.and.predict.vcf.R")
annotate.vcf(path.to.vcf = "./cache/allSigVars.vcf", tissue.id = "E119")
