library(dplyr)
starrseq.vars <- fread("~/Downloads/starr-seq.variants.csv", header = T, stringsAsFactors = F)
starrseq.vars[, OR_fdr := as.numeric(OR_fdr)]

AS.variants <- starrseq.vars %>%
	mutate(OR_fdr = as.numeric(OR_fdr)) %>%
	filter(label != "Inactive") %>%
	filter(OR_fdr < 0.1)
table(AS.variants$label)

write.csv(AS.variants, "./reports/starrseq.AS.variants.csv", row.names = F, quote = F)

AS.variants <- fread("./reports/starrseq.AS.variants.updated.csv", header = T, stringsAsFactors = F)


## Check that they got pooled correctly
dbsnp.vars <- read.table("./cache/extracted.snps.vcf.recode.vcf", col.names = c("CHROM", "POS", "SNP", "REF", "ALT", "QUAL", "FILTER", "INFO"), stringsAsFactors = F)

head(dbsnp.vars, 2)

combo.table <- as.data.table(left_join(AS.variants, dbsnp.vars, by = "SNP"))
summary(combo.table$Chr == combo.table$CHROM)
summary(combo.table$Coordinate == combo.table$POS)

combo.table[is.na(POS) ,c("SNP", "Chr", "CHROM", "Coordinate", "POS")]
combo.table[,c("SNP", "Chr", "CHROM", "Coordinate", "POS")]

## Write out a vcf
write.table(combo.table[,c("CHROM", "POS", "SNP", "REF", "ALT")], "./cache/starrseq.AS.variants.vcf", row.names = F, col.names = F, sep = "\t", quote = F)
