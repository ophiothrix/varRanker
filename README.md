#### IMPORTANT ####
The package relies on having a current installation of MEME suite and motif databases for calculating motif damage score. Make sure to specify the correct path to FIMO function at the top of *./lib/base.functions.R* script


# varRank

List of data files
- Training set files from Maurano paper:
	* ng.3432-S5.txt downloadable as paper supplementary - this file lists SNPs with reference and alternate allele with evidence for allele imbalance across all tissues - used to extract REF and ALT alleles for tissue-specific variants
	* snps.multicell.bySample downloadable from CATO score page: http://www.mauranolab.org/CATO/snps.multicell.bySample.tar.gz - one file per tissue indicating whether a given allele showed evidence for allele imbalance


Dependencies
UCSC's bigBedtoBed
MEME suite
