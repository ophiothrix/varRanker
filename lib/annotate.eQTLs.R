#### Function to annotate eQTL variants from GTEx ####
## Takes as arguments eQTL file id (tissue name only) and ROADMAP tissue id to be used for variant annotation
## Returns and saves GRanges object containing annotated variants.

annotate.eQTLs <- function(variant.file.id, tissue.id) {
  require(GenomicRanges)
  variants <- read.table(paste0("/no_backup/so/aholik/RegVar/data/GTEx/", variant.file.id, "_Analysis.snpgenes"), stringsAsFactors = F, sep = "\t", header = T)
  ## Subset to chosen SNPs
  table(variants$is_chosen_snp)
  variants <- variants[variants$is_chosen_snp == 1,]
  head(variants)
  variants <- GRanges(variants$snp_chrom, IRanges(variants$snp_pos, variants$snp_pos), REF = variants$ref, ALT = variants$alt)
  seqlevels(variants) <- paste0("chr", seqlevels(variants))
  ## At the moment, can't deal with indels. Remove them
  require(stringr)
  variants <- variants[nchar(variants$REF) == 1 & nchar(variants$ALT) == 1]
  
  ## Order the varians
  variants <- variants[order(variants)]
  
  ## Annotate variants with genomic features, poly-tissue GEP and ChromHMM sets, TF motif scores, pre-computed damage scores
  source("./lib/feature.motif.annotation.R")
  variants <- feature.motif.annotation(variants, tissue.id)

  ## Annotate variants with conservation scores
  source("./lib/add.conservation.scores.R")
  variants <- annotate.GERP(variants)
  variants <- annotate.phastCons(variants)
  variants <- annotate.phyloP(variants)
  
  ## Annotate variants with tissue-specific epigenetic marks and ChromHMM and GEP predicted states
  source("./lib/epigenome.annotation.R")
  variants <- epigenome.annotation(variants, tissue.id)

  ## Save the variant file
  saveRDS(variants, file = paste0("./cache/", tissue.id, ".annotated.", variant.file.id, ".variants.rds"))
  return(variants)
}
