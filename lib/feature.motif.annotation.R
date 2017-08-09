feature.motif.annotation <- function(variants, tissue_id) {
  ## Add feature annotation
  require(GenomicRanges)
  source("./lib/add.features.to.variant.sets.R")
  variants <- assign.features(variants)

  variants$id <- 1:length(variants)
  
  # Adding motif damage scores ----------------------------------------------
  
  
  ### Add motif damage score for ALL the variants. It can't be explicitly calculated for indels, so we'll make an assumption that if indel hits a motif, it destroys it and hence is highly damaging. We could make a more refined assumption based on how central the indel is. But let's keep it simple for now.
  
  ## Start with SNPs and generate explicit motif damage.
  print("Calculating motif damage scores...")
  # Make a new object only containing SNV - Both REF and ALT are single nucleotide
  snp.variants <- variants[nchar(variants$REF) == 1 & nchar(variants$ALT) == 1]
  source("./lib/motif.damage.annotation.R")
  
  ## Add JASPAR damage score annotation to SNVs
  snp.variants <- get.damage.scores.direct(database.path = "~/utils/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme", variants = snp.variants)
  snp.variants$jaspar.damage.score <- snp.variants$log2.damage.score
  snp.variants$jaspar.motif.hit <- snp.variants$motif.hit
  
  ## Add HOCOMOCO damage score annotation SNVs
  snp.variants <- get.damage.scores.direct(database.path = "~/utils/motif_databases/HUMAN/HOCOMOCOv9.meme", variants = snp.variants)
  snp.variants$hocomoco.damage.score <- snp.variants$log2.damage.score
  snp.variants$hocomoco.motif.hit <- snp.variants$motif.hit
  
  ## Add the damage scores back to the original variant set
  variants$hocomoco.damage.score <- variants$jaspar.damage.score <- NA
  test1 <- (paste0(snp.variants$id, snp.variants$ALT))
  test2 <- (paste0(variants$id, variants$ALT))
  variants$jaspar.damage.score <- snp.variants$jaspar.damage.score[match(test2, test1)]
  variants$hocomoco.damage.score <- snp.variants$hocomoco.damage.score[match(test2, test1)]
  
  # Clean up
  rm(snp.variants)
  gc()
  
  # Add promoter and enhancer annotations -----------------------------------
  source("./lib/annotate.variants.from.GRanges.R")
  
  ## Add annotation for different chromHMM classes
  print("Adding pan-tissue chromatin state annotation...")
  all.states <- data.frame(matrix(nrow = length(variants)))
  head(all.states)
  for (state in list.files("./data/ENCODE/chromHMM.calls", "rds", full.names = T)) {
    print(state)
    state.name <- gsub(".union.rds", "", gsub(".*_", "", state))
    state <- readRDS(state)
    anno <- annotate.variants.from.GRanges(variants, state)
    all.states[,state.name] <- anno$n.tissues
    rm(list=c("anno", "state"))
    gc()
  }
  all.states <- all.states[,-1]
  mcols(variants) <- cbind(values(variants), all.states)
  rm(all.states)
  gc()
  
  ## Add annotation for GEP predicted enhancers
  print("Adding GEP enhancer predictions...")
  gep.enhancers <- readRDS("./data/GEP.enhancer.union.rds")
  anno <- annotate.variants.from.GRanges(variants, gep.enhancers)
  variants$GEP <- anno$n.tissues
  rm(anno)
  gc()
  
  ## Add annotation for Core promoters
  print("Adding core promoter annotation...")
  core.proms <- read.table("./data/gc19_pc.promCore.expanded.bed")
  colnames(core.proms) <- c("chr", "start", "end", "gene", "score", "strand")
  core.proms <- GRanges(core.proms)
  seqlevels(core.proms) <- paste0("chr", seqlevels(core.proms))
  core.proms
  anno <- annotate.variants.from.GRanges(variants, core.proms)
  variants$prom.Core <- anno
  rm(list = c("anno", "core.proms"))
  gc()
  
  ## Add annotation for Domain Promoters
  print("Adding domain promoter annotation...")
  domain.proms <- read.table("./data/gc19_pc.promDomain.expanded.bed")
  colnames(domain.proms) <- c("chr", "start", "end", "gene", "score", "strand")
  domain.proms <- GRanges(domain.proms)
  seqlevels(domain.proms) <- paste0("chr", seqlevels(domain.proms))
  anno <- annotate.variants.from.GRanges(variants, domain.proms)
  variants$prom.Domain <- anno
  rm(list = c("anno", "domain.proms"))
  gc()
  
  ## The problem with most of the pre-computed damage scores is that they are not available for ALL possible variants, with the exception of CADD. So don't include them in the model. We could include CADD and it actually performs quite well, but the files are massive and it's very resource intensive extracting the scores. So not so good for production.
  # ## Add Funseq Scores
  # variant.ids <- paste(seqnames(variants), start(variants), variants$REF, variants$ALT, sep = ".")
  # 
  # funseq.scores <- read.table("/no_backup/so/aholik/marker.paper/data/funseq.scores.simple.bed", sep = "\t", header = F, stringsAsFactors = F)
  # colnames(funseq.scores) <- c("chr", "start", "end", "ref", "alt", "score")
  # head(funseq.scores)
  # funseq.scores$chr <- paste0("chr", funseq.scores$chr)
  # funseq.ids <- paste(funseq.scores$chr, funseq.scores$start, funseq.scores$ref, funseq.scores$alt, sep=".")
  # 
  # variants$funseq.PCAGW.score <- funseq.scores$score[match(variant.ids, funseq.ids)]
  # summary(is.na(variants$funseq.PCAGW.score))
  # 
  # saveRDS(variants, "./cache/variants.tmp.rds")
  # 
  # 
  # ## Extract pre-computed scores, e.g. CADD, GWAVA, FunSeq, etc. Note that some may not be available.
  # source("~/RegVar/src/variant.scoring.prototype.R")
  # 
  # ## The problem with AS DHS is that first allele is not necessarily a reference. So we can't always get a CADD score.
  # require(BSgenome.Hsapiens.UCSC.hg19)
  # refs <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, variants))
  # summary(refs == variants$REF)
  # ## Make a temporary variants object and swap reference and ALT alleles
  # ref.variants <- variants
  # ref.variants[refs != variants$REF]$REF <- variants[refs != variants$REF]$ALT
  # ref.variants[refs != variants$REF]$ALT <- variants[refs != variants$REF]$REF
  # 
  # ref.variants <- extract.scores(ref.variants)
  # variants <- ref.variants

# summary(is.na(ref.variants$funseq.score))
# summary(is.na(ref.variants$cadd.score))
# summary(is.na(ref.variants$gwava.tss))
# summary(is.na(ref.variants$gwava.random))
  
  # test.variants <- variants[1:100]
  # test.variants <- extract.scores(test.variants)
  # 
  # ### Save the annotated object
  # saveRDS(variants, file = "./data/chrom.anno.play.training.set.rds")
  return(variants)
}
