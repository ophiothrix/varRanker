## Get cut off for 5% FDR
predict.valid <- as.data.frame(h2o.predict(model_ensemble, valid.set.h2o))
predict.valid$truth <- as.numeric(as.vector(valid.set.h2o$regulatory))
predict.valid <- predict.valid[order(predict.valid$p1, decreasing = T),]
predict.valid$fdr <- cumsum(predict.valid$truth == 0) / sum(predict.valid$truth == 0)
head(predict.valid)
cut.off <- max(which(predict.valid$fdr < 0.1))
cut.off <- predict.valid$p1[cut.off]
predict.valid[(cut.off - 5):(cut.off + 5),]



path.to.vcf <- "./cache/BRCA-US.small.vcf.gz"
tissue.id <- "E119"

if (file.exists(paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))) {
	variants <- readRDS(paste0("./cache/", tissue.id, ".annotated.", gsub(".vcf.*", "", basename(path.to.vcf)), ".variants.rds"))
}

variants <- as.data.frame(variants)
variants <- variants[!variants$feature == "coding",]
variants$AF <- as.numeric(variants$AF)
variants$AF[variants$AF > 0.5] <- 1 - variants$AF[variants$AF > 0.5]
head(variants)
target.set <- as.h2o(as.data.frame(variants))
table(variants$feature)
table(train.set$feature)

## Predict functionality
predicted <- as.data.frame(h2o.predict(model_ensemble, target.set))
variants$p1 <- predicted$p1
variants$predicted <- predicted$predict
table(variants$predicted)

hist(variants$p1, breaks = 100)
hist(variants$AF, breaks = 100, freq = F)

hist(variants$AF[variants$p1 >= cut.off], breaks = 500, freq = F, xlim = c(0, 0.1), col = "#ff000044")
hist(variants$AF[variants$p1 < cut.off], breaks = 500, freq = F, add = T, col = "#00ffff44")

hist(variants$AF[variants$p1 >= cut.off], breaks = 100, freq = F, xlim = c(0, 0.5), col = "#ff000044")
hist(variants$AF[variants$p1 < cut.off], breaks = 100, freq = F, add = T, col = "#00ffff44")
summary(variants$p1 >= cut.off)

hist(variants$AF[variants$predicted == 1], breaks = 500, freq = F, xlim = c(0, 0.1), col = "#ff000044")
hist(variants$AF[variants$predicted == 0], breaks = 500, freq = F, add = T, col = "#00ffff44")


hist(variants$AF[variants$predicted == 1], breaks = 100, freq = F, xlim = c(0, 0.5), col = "#ff000044")
hist(variants$AF[variants$predicted == 0], breaks = 100, freq = F, add = T, col = "#00ffff44")
