### Requires GTEx_Analysis_v6p_eQTL.tar file from GTEx portal. Can't be downloaded automatically, at least not easily.
require(data.table)
require(GenomicRanges)
dir.create("./data/GTEx", showWarnings = F, recursive = T)
## Place the file in "./data/GTEx"

system("tar xvf ./data/GTEx/GTEx_Analysis_v6p_eQTL.tar -C ./data/GTEx")
file.remove("./data/GTEx/GTEx_Analysis_v6p_eQTL.tar")
## Combine variant ids from all tissues into a single file
system("gunzip -c ./data/GTEx/GTEx_Analysis_v6p_eQTL/*.signif_snpgene_pairs.txt.gz | cut -f1 | gzip -c > ./data/GTEx/all.variant.ids.txt.gz")

## Read in variant ids
variant.ids <- fread(input = "gunzip -c ./data/GTEx/all.variant.ids.txt.gz", header = T, stringsAsFactors = F)
variant.ids <- variant.ids[variant_id != "variant_id"]
variant.ids <- unique(variant.ids$variant_id)

## Split variant ids into the constituent components
variant.table <- do.call(rbind, strsplit(variant.ids, split = "_"))
gc()
colnames(variant.table) <- c("chr", "pos", "REF", "ALT", "refgen")
variant.table <- as.data.table(variant.table)
if (nrow(variant.table) != length(variant.ids)) {
	stop("The number of rows in the variant table is not the same as the length of variant ids vector")
}
## Convert the table into GRanges object
variant.table$chr <- paste0("chr", variant.table$chr)
variant.table$end <- variant.table$start <- as.numeric(variant.table$pos)
head(variant.table)
variant.table <- GRanges(variant.table[,-c(2, 5)])

## For deletions, any of the underlying variants could be an eQTL, extend the ranges to cover the size of deletion
ref.len <- nchar(variant.table$REF)
alt.len <- nchar(variant.table$ALT)
summary(ref.len > 1 & alt.len > 1) ## There are no variants with both reference and alternate alleles longer than 1
summary(ref.len > 1) # 135000 deletions
width(variant.table) <- ref.len

## Save the full list of eQTLs
saveRDS(variant.table, "./data/GTEx/all.eQTLs.rds")
