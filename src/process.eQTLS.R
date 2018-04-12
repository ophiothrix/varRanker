library(tidyverse)
eqtls <- read.table("./cache/leadsnps_100kb_header_fdr.tsv.gz")
head(eqtls)
eqtls.vcf <- eqtls %>%
	mutate(chr = paste0("chr", gsub("^(.*)_.*_.*_.*_.*", "\\1", V2))) %>%
	mutate(pos = gsub("^.*_(.*)_.*_.*_.*", "\\1", V2)) %>%
	mutate(REF = gsub("^.*_.*_(.*)_.*_.*", "\\1", V2)) %>%
	mutate(ALT = gsub("^.*_.*_.*_(.*)_.*", "\\1", V2)) %>%
	select(c(chr, pos, REF, ALT))

write.table(eqtls.vcf, "./cache/pancanrisk.eqtls.vcf", sep = "\t", quote = F, row.names = F, col.names = F)

source("./lib/annotate.and.predict.vcf.R")
annotate.vcf(path.to.vcf = "./cache/pancanrisk.eqtls.vcf", tissue.id = "E029")
variants <- readRDS("./cache/E034.annotated.pancanrisk.eqtls.variants.rds")
source("./lib/epigenome.annotation.R")
variants <- epigenome.annotation(variants = variants, tissue.id = "E056")
tissue.id = "E056"
anno.mat <- as.data.frame(mcols(variants))
anno.mat <- anno.mat[,-(43:142)]
colnames(anno.mat) <- gsub("\\.1$", "", colnames(anno.mat))

head(anno.mat)
mcols(variants) <- anno.mat

breast.eqtls <- fread("./cache/pancanrisk.eqtls.E119.annotated.csv", stringsAsFactors = F)
breast.eqtls$p.regulatory <- as.numeric(breast.eqtls$p.regulatory)
breast.eqtls <- breast.eqtls[order(breast.eqtls$varID),]

liver.eqtls <- fread("./cache/pancanrisk.eqtls.E118.annotated.csv", stringsAsFactors = F)
liver.eqtls$p.regulatory <- as.numeric(liver.eqtls$p.regulatory)
liver.eqtls <- liver.eqtls[order(liver.eqtls$varID),]

tcell.eqtls <- fread("./cache/pancanrisk.eqtls.E034.annotated.csv", stringsAsFactors = F)
tcell.eqtls$p.regulatory <- as.numeric(tcell.eqtls$p.regulatory)
tcell.eqtls <- tcell.eqtls[order(tcell.eqtls$varID),]

mono.eqtls <- fread("./cache/pancanrisk.eqtls.E029.annotated.csv", stringsAsFactors = F)
mono.eqtls$p.regulatory <- as.numeric(mono.eqtls$p.regulatory)
mono.eqtls <- mono.eqtls[order(mono.eqtls$varID),]

lung.fibro <- fread("./cache/pancanrisk.eqtls.E056.annotated.csv", stringsAsFactors = F)
lung.fibro$p.regulatory <- as.numeric(lung.fibro$p.regulatory)
lung.fibro <- lung.fibro[order(lung.fibro$varID),]

skin.fibro <- fread("./cache/pancanrisk.eqtls.E126.annotated.csv", stringsAsFactors = F)
skin.fibro$p.regulatory <- as.numeric(skin.fibro$p.regulatory)
skin.fibro <- skin.fibro[order(skin.fibro$varID),]


predicted.eqtls$p.regulatory <- as.numeric(predicted.eqtls$p.regulatory)
summary((predicted.eqtls$p.regulatory))
hist(predicted.eqtls$p.regulatory, breaks = 100)
summary(predicted.eqtls$p.regulatory > .73)
length(unique(predicted.eqtls$varID))
dim(predicted.eqtls)

pdf("./graphs/regulatory.score.scatters.pdf", 10, 6)
par(mfrow = c(2,2))
par(mar = c(4, 4, 0.5, 1))
smoothScatter(tcell.eqtls$p.regulatory, breast.eqtls$p.regulatory, xlab = "T cell regulatory score", ylab = "Breast regulatory score")
# abline(0, 1, lty = 2, lwd = 2, col = "#666666")
smoothScatter(breast.eqtls$p.regulatory, liver.eqtls$p.regulatory, xlab = "Breast regulatory score", ylab = "Liver regulatory score")
# abline(0, 1, lty = 2, lwd = 2, col = "#666666")
smoothScatter(tcell.eqtls$p.regulatory, liver.eqtls$p.regulatory, xlab = "T cell regulatory score", ylab = "Liver regulatory score")
# abline(0, 1, lty = 2, lwd = 2, col = "#666666")
smoothScatter(tcell.eqtls$p.regulatory, mono.eqtls$p.regulatory, xlab = "T cell regulatory score", ylab = "Monocyte regulatory score")
dev.off()

plot(liver.eqtls$p.regulatory, breast.eqtls$p.regulatory, xlab = "T cell regulatory score", ylab = "Breast regulatory score", cex = 0.5, pch = 16, col = "#88888888")

smoothScatter(lung.fibro$p.regulatory, skin.fibro$p.regulatory, xlab = "T cell regulatory score", ylab = "Monocyte regulatory score")

pdf("./graphs/feature.importance.pdf", 10)
h2o.varimp_plot(model_drf, 30)
dev.off()
